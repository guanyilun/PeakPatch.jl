# Large-scale RNG null test for the production white noise (Threefry-2x64-20).
#
# Question (J.R. Bond): can the RNG imprint spurious correlations on large
# (angular) scales when ~2.3e11 deviates are drawn? The Fortran random.f90 LCG
# guards against this by partitioning its 2^46 cycle into disjoint streams.
#
# Two tests, both against independent Xoshiro white noise pushed through
# IDENTICAL statistics and geometry (so box replication, pixelization and
# mode counting cancel and only the generator is tested):
#
#  A. ENSEMBLE (generator defect test): K Threefry seeds vs K Xoshiro seeds of
#     the coarse block-averaged noise (`_downsample_noise`, block 12^3 as in
#     production) at N=1536 → two-sample z and KS on every statistic. Plus K
#     fine-resolution 256^3 subcubes taken deep in the N=6144 counter space
#     (incl. the Box-Muller pairing along x).
#  B. PRODUCTION REALIZATION: `_downsample_noise(6144, 512, 12345)` — the exact
#     coarse noise feeding every production octant — ranked within a Xoshiro
#     512^3 ensemble. One realization: tells whether our seed is a typical
#     draw, not whether the generator is sound (that is A).
#
# Statistics per field: moments, low-k and axis-aligned 3D power, shell P(k),
# lag correlations along each axis, max |ξ| over ALL 3D separations, and
# angular C_ℓ of radial-shell HEALPix maps from an observer at the box origin
# (== the 8-octant full-sky production geometry, with periodic replication).
#
# Usage: julia --project=validation -t 32 validation/rng/rng_large_scale_test.jl OUTDIR [K]
#        RNG_SMOKE=1 → tiny grids (functional check only)

using PeakPatch, FFTW, Healpix, Random, Statistics, Printf, DelimitedFiles

const OUT   = ARGS[1]
const K     = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 40
const SMOKE = get(ENV, "RNG_SMOKE", "0") == "1"
const BLOCK = 12                                     # production 6144/512
const NP, MP, SEEDP = SMOKE ? (768, 64, 12345) : (6144, 512, 12345)   # production
const NE = SMOKE ? 384 : 1536; const ME = NE ÷ BLOCK                   # ensemble
const NSIDE_P, NSIDE_E = SMOKE ? (16, 8) : (64, 32)
const NSUB = SMOKE ? 32 : 256
const SHELLS = [(0.30, 0.50), (0.50, 0.75), (0.75, 1.00), (0.30, 1.00)]  # r / L
const LBINS = [(2, 5), (6, 10), (11, 20), (21, 40), (41, 80)]
const LAGS = [1, 2, 3, 4, 6, 8, 16, 32, 64, 128, 256]
const PART = get(ENV, "RNG_PART", "AB")            # any of "A", "B", "C"
mkpath(OUT)
# Threaded FFTW segfaulted (spawn_apply) on 512^3 under Julia 1.12, job 5635142;
# RNG_FFTW_THREADS=1 avoids it at a few s per transform.
FFTW.set_num_threads(parse(Int, get(ENV, "RNG_FFTW_THREADS", string(Threads.nthreads()))))

# ------------------------------------------------------------------ statistics
"""3D statistics of a unit-variance white field → ordered name => value pairs."""
function stats3d(A::Array{Float32,3})
    n = size(A, 1); n3 = n^3; x = Float64.(A)
    mu = mean(x); v = var(x)
    out = Pair{String,Float64}["mean*sqrt(n3)" => mu * sqrt(n3), "variance" => v,
        "skewness" => mean((x .- mu) .^ 3) / v^1.5, "exkurtosis" => mean((x .- mu) .^ 4) / v^2 - 3]
    F = rfft(x); P = abs2.(F) ./ n3           # white → <P> = 1 per mode
    kf(i) = i - 1 <= n ÷ 2 ? i - 1 : i - 1 - n
    edges = [0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 32.5, 64.5, 128.5, 256.5]
    nb = count(<=(n ÷ 2), edges) - 1
    psum = zeros(nb); pcnt = zeros(nb); pmax = 0.0
    for k in 1:n, j in 1:n, i in 1:size(P, 1)
        kk = sqrt((i - 1)^2 + kf(j)^2 + kf(k)^2)
        0 < kk <= 32 && (pmax = max(pmax, P[i, j, k]))
        b = searchsortedlast(edges, kk)
        (1 <= b <= nb) || continue
        psum[b] += P[i, j, k]; pcnt[b] += 1
    end
    for b in 1:nb
        push!(out, @sprintf("P(|k| %g-%g)", edges[b] + 0.5, edges[b+1] - 0.5) => psum[b] / pcnt[b])
    end
    push!(out, "max P(0<|k|<=32)" => pmax)
    na = min(32, n ÷ 2)
    push!(out, "n axis modes P>6 (|k|<=$na)" => count(>(6), vcat(P[2:na+1, 1, 1], P[1, 2:na+1, 1], P[1, 1, 2:na+1])))
    push!(out, "P kx-axis 1-$na" => mean(P[2:na+1, 1, 1]), "P ky-axis 1-$na" => mean(P[1, 2:na+1, 1]),
               "P kz-axis 1-$na" => mean(P[1, 1, 2:na+1]))
    xi = irfft(abs2.(F), n) ./ n3; xi ./= xi[1, 1, 1]
    for l in filter(<=(n ÷ 2), LAGS)
        push!(out, "xi_x($l)*sqrt(n3)" => xi[l+1, 1, 1] * sqrt(n3),
                   "xi_y($l)*sqrt(n3)" => xi[1, l+1, 1] * sqrt(n3),
                   "xi_z($l)*sqrt(n3)" => xi[1, 1, l+1] * sqrt(n3))
    end
    xi[1, 1, 1] = 0.0
    push!(out, "max|xi(r!=0)|*sqrt(n3)" => maximum(abs, xi) * sqrt(n3))
    return out
end

"""Mean-noise HEALPix maps in radial shells from an observer at the origin of
the periodic box (== 8-octant production geometry) → C_ℓ averaged in ℓ bins."""
function shell_stats(A::Array{Float32,3}, nside::Int)
    n = size(A, 1); L = Float64(n)
    res = Healpix.Resolution(nside); npix = 12 * nside^2; ns = length(SHELLS)
    nt = Threads.maxthreadid()                 # 1.12: threadid() includes interactive pool
    sums = [zeros(npix, ns) for _ in 1:nt]; cnts = [zeros(npix, ns) for _ in 1:nt]
    Threads.@threads :static for k in 1:n
        acc_s = sums[Threads.threadid()]; acc_c = cnts[Threads.threadid()]
        @inbounds for j in 1:n, i in 1:n
            val = Float64(A[i, j, k])
            for sx in (0, -1), sy in (0, -1), sz in (0, -1)
                px = i - 0.5 + sx * L; py = j - 0.5 + sy * L; pz = k - 0.5 + sz * L
                r = sqrt(px^2 + py^2 + pz^2) / L
                (r < SHELLS[1][1] || r >= 1.0) && continue
                pix = Healpix.vec2pixRing(res, px, py, pz)
                for (s, (r0, r1)) in enumerate(SHELLS)
                    if r0 <= r < r1
                        acc_s[pix, s] += val; acc_c[pix, s] += 1
                    end
                end
            end
        end
    end
    S = sum(sums); C = sum(cnts)
    lmax = min(3 * nside - 1, LBINS[end][2])
    out = Pair{String,Float64}[]
    for (s, (r0, r1)) in enumerate(SHELLS)
        @assert all(C[:, s] .> 0) "empty pixels in shell $s (nside too high)"
        m = Healpix.HealpixMap{Float64, Healpix.RingOrder}(S[:, s] ./ C[:, s])
        cl = Healpix.alm2cl(Healpix.map2alm(m; lmax = lmax, niter = 0))  # same treatment for both ensembles
        for (l0, l1) in LBINS
            l1 <= lmax || continue
            push!(out, @sprintf("Cl shell %.2f-%.2f l%d-%d", r0, r1, l0, l1) => mean(cl[l0+1:l1+1]))
        end
    end
    return out
end

analyze(A, nside) = vcat(stats3d(A), shell_stats(A, nside))
vals(st) = last.(st)

# ------------------------------------------------------------------ comparison
"""Two-sample KS statistic and asymptotic p-value."""
function ks2(a, b)
    a = sort(a); b = sort(b); na, nb = length(a), length(b); d = 0.0
    for x in vcat(a, b)
        d = max(d, abs(searchsortedlast(a, x) / na - searchsortedlast(b, x) / nb))
    end
    ne = na * nb / (na + nb); λ = (sqrt(ne) + 0.12 + 0.11 / sqrt(ne)) * d
    p = 2 * sum((-1)^(j - 1) * exp(-2 * j^2 * λ^2) for j in 1:100)
    return d, clamp(p, 0.0, 1.0)
end

io = open(joinpath(OUT, PART == "AB" ? "REPORT.md" : "REPORT_$(PART).md"), "w")
pr(a...) = (println(io, a...); println(a...); flush(io))
tsec(t) = @sprintf("%.0f s", t)

pr("# RNG large-scale null test — Threefry-2x64-20 vs Xoshiro (Julia default RNG)\n")
pr("Block $(BLOCK)³ coarse noise as in production (`_downsample_noise`). Shells r/L = ",
   "$(SHELLS); ℓ bins $(LBINS). Observer at the periodic-box origin (8-octant geometry).\n")

function ensemble_table(names, tf, xo, label)
    T = reduce(hcat, tf); X = reduce(hcat, xo); nstat = length(names)
    zmax = 0.0; pmin = 1.0; rows = Any[]
    pr("| statistic | Threefry mean | Xoshiro mean | z | KS p |\n|---|---|---|---|---|")
    for i in 1:nstat
        a = T[i, :]; b = X[i, :]
        z = (mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))
        _, p = ks2(a, b); zmax = max(zmax, abs(z)); pmin = min(pmin, p)
        pr(@sprintf("| %s | %.5g | %.5g | %+.2f | %.3f |", names[i], mean(a), mean(b), z, p))
        push!(rows, (names[i], mean(a), std(a), mean(b), std(b), z, p))
    end
    pr(@sprintf("\n**%s summary**: %d statistics; max |z| = %.2f (typical max of %d Gaussians ≈ %.1f); ",
                label, nstat, zmax, nstat, sqrt(2 * log(nstat))),
       @sprintf("min KS p = %.4f → Bonferroni p = %.3f (> 0.05 ⇒ no detectable difference).\n", pmin,
                min(1.0, pmin * nstat)))
    writedlm(joinpath(OUT, "ensemble_$(label).csv"),
             vcat(permutedims(["stat", "tf_mean", "tf_std", "xo_mean", "xo_std", "z", "ks_p"]),
                  reduce(vcat, [permutedims(collect(r)) for r in rows])), ',')
end

# ---- A. ensemble ----------------------------------------------------------
if occursin("A", PART)
pr("## A. Ensemble test (generator soundness): K = $K Threefry seeds vs K Xoshiro seeds\n")
pr("### Coarse fields: N = $NE → M = $ME (block $(BLOCK)³), Nside $NSIDE_E\n")
tA = @elapsed begin
    stA = analyze(PeakPatch.MultiResolution._downsample_noise(NE, ME, 1), NSIDE_E)
    tf = [s == 1 ? vals(stA) : vals(analyze(PeakPatch.MultiResolution._downsample_noise(NE, ME, s), NSIDE_E)) for s in 1:K]
    xo = [vals(analyze(randn(Xoshiro(10_000 + s), Float32, ME, ME, ME), NSIDE_E)) for s in 1:K]
end
ensemble_table(first.(stA), tf, xo, "coarse")
pr("(coarse ensemble: $(tsec(tA)))\n")

# Fine subcubes at random offsets deep in the production counter space (N = NP)
subpos(s) = let r = Xoshiro(30_000 + s); ntuple(_ -> rand(r, 1:NP-NSUB+1), 3) end
function fine_threefry(s)
    A = Array{Float32,3}(undef, NSUB, NSUB, NSUB); o = subpos(s)
    PeakPatch.RandomField.fill_noise_threefry_region!(A, ntuple(d -> o[d]:o[d]+NSUB-1, 3), NP, s)
    return A
end
pr(@sprintf("### Fine-resolution subcubes: %d³ at random offsets in the N = %d grid (linear counter up to %.2e)\n",
            NSUB, NP, Float64(NP)^3))
tF = @elapsed begin
    stF = stats3d(fine_threefry(1))
    ftf = [s == 1 ? vals(stF) : vals(stats3d(fine_threefry(s))) for s in 1:K]
    fxo = [vals(stats3d(randn(Xoshiro(20_000 + s), Float32, NSUB, NSUB, NSUB))) for s in 1:K]
end
ensemble_table(first.(stF), ftf, fxo, "fine")
pr("(fine ensemble: $(tsec(tF)))\n")
end  # part A

# ---- B. production realization ---------------------------------------------
if occursin("B", PART)
pr("## B. Production realization: `_downsample_noise($NP, $MP, $SEEDP)` ranked in K = $K Xoshiro $(MP)³ fields (Nside $NSIDE_P)\n")
coarse_path = joinpath(OUT, "coarse_threefry_N$(NP)_M$(MP)_seed$(SEEDP).f32")
tP = @elapsed begin
    global TP
    if isfile(coarse_path)
        TP = Array{Float32,3}(undef, MP, MP, MP); read!(coarse_path, TP)
    else
        TP = PeakPatch.MultiResolution._downsample_noise(NP, MP, SEEDP); write(coarse_path, TP)
    end
end
pr(@sprintf("(production coarse noise, %.2e deviates: %s)\n", Float64(NP)^3, tsec(tP)))
tB = @elapsed begin
    stP = analyze(TP, NSIDE_P)
    XP = reduce(hcat, [vals(analyze(randn(Xoshiro(40_000 + s), Float32, MP, MP, MP), NSIDE_P)) for s in 1:K])
end
nflag = 0
pr("| statistic | production | Xoshiro mean ± std | two-sided rank p |\n|---|---|---|---|")
for (i, (nm, v)) in enumerate(stP)
    b = XP[i, :]; lo = count(<(v), b); hi = count(>(v), b)
    p = min(1.0, 2 * (min(lo, hi) + 1) / (length(b) + 1))
    global nflag += p < 0.05
    pr(@sprintf("| %s | %.5g | %.5g ± %.2g | %.3f |", nm, v, mean(b), std(b), p))
end
pr(@sprintf("\n**Production summary**: %d / %d statistics with rank p < 0.05 (≈ %.1f expected by chance; ",
            nflag, length(stP), 0.05 * length(stP)),
   "rank p floor with K = $K is $(round(2 / (K + 1), digits = 3))). (analysis: $(tsec(tB)))")
end  # part B

# ---- C. production-scale ensemble -------------------------------------------
# B found seed 12345 ~1-in-40 on two statistics (ky-axis power, max mode). A ran
# its coarse ensemble at N=1536, so a defect appearing only when all 2.3e11
# production counters feed the large-scale modes was untested: KC seeds at the
# full production N=6144 → M=512 vs the same Xoshiro 512³ ensemble as B.
if occursin("C", PART)
    KC = parse(Int, get(ENV, "RNG_KC", "12"))
    pr("## C. Production-scale ensemble: $KC Threefry seeds at N = $NP → M = $MP vs K = $K Xoshiro $(MP)³ (Nside $NSIDE_P)\n")
    function coarse_cached(seed)
        f = joinpath(OUT, "coarse_threefry_N$(NP)_M$(MP)_seed$(seed).f32")
        isfile(f) && (A = Array{Float32,3}(undef, MP, MP, MP); read!(f, A); return A)
        A = PeakPatch.MultiResolution._downsample_noise(NP, MP, seed); write(f, A); return A
    end
    tC = @elapsed begin
        stC = analyze(coarse_cached(1), NSIDE_P)
        ctf = [vals(stC)]
        for s in 2:KC
            push!(ctf, vals(analyze(coarse_cached(s), NSIDE_P))); @info "C: Threefry seed $s/$KC"
        end
        cxo = [vals(analyze(randn(Xoshiro(40_000 + s), Float32, MP, MP, MP), NSIDE_P)) for s in 1:K]
    end
    ensemble_table(first.(stC), ctf, cxo, "production_scale")
    pr("(production-scale ensemble: $(tsec(tC)))\n")
end
close(io)
@info "done" OUT
