# Field-level test of the multi-resolution (MUSIC-style) synthesis, plan item A5(b)
# (docs/paper_comparison_plan_2026-09.md). CPU only.
#
# For one global grid (N, cellsize = Websky 0.85221 Mpc/h, seed) and a FIXED coarse grid M
# (block = N/M = 12, as in production cf=32), synthesizes δ and ψ₁ₓ with several tilings
# (ntile = 1, 2, 4, 8; nbuff = 16 fixed) and compares each against the exact global periodic
# FFT field built from the same Threefry noise and the same √P kernel convention:
#   1. noise / residual identity across tilings (exact, bitwise)
#   2. real-space rms error vs global, and its profile vs distance to the nearest tile face
#   3. tiling-vs-tiling difference
#   4. P(k) ratio and cross-correlation r(k) of the stitched core field vs global (same
#      region, same taper) — looks for features at the coarse Nyquist and tile scales
#   5. a 2-D cross-section through tile boundaries (coarse / self / total / global) for F2
#
#   julia --project=validation -t 16 validation/tiling/field_split_test.jl [outdir]
using PeakPatch
using FFTW, Printf, Statistics
import PeakPatch.MultiResolution as MR
import PeakPatch.RandomField: fill_noise_threefry!
import PeakPatch.PowerSpectrum: load_pk

FFTW.set_num_threads(parse(Int, get(ENV, "FIELD_FFTW_THREADS", "1")))
const OUT = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "results")
mkpath(OUT)
const DATA = joinpath(@__DIR__, "..", "websky_6144", "data")

const N = parse(Int, get(ENV, "FIELD_N", "480"))
const NBUFF = 16
const BLOCK = parse(Int, get(ENV, "FIELD_BLOCK", "12"))
const M = N ÷ BLOCK
const ALATT = 0.85221
const BOX = N * ALATT
const SEED = 12345
const NTILES = parse.(Int, split(get(ENV, "FIELD_NTILES", "1,2,4,8"), r"[,:]"))
const NC = N - 2NBUFF                      # region covered by tile cores
@assert N % BLOCK == 0
for nt in NTILES; @assert NC % nt == 0 "NC=$NC not divisible by ntile=$nt"; end
pk = load_pk(joinpath(DATA, "pk_websky.dat"))

@printf("N=%d  box=%.2f Mpc/h  cell=%.5f  M=%d (block %d, coarse cell %.2f Mpc/h)  nbuff=%d  NC=%d\n",
        N, BOX, ALATT, M, BLOCK, BLOCK * ALATT, NBUFF, NC)
@printf("coarse Nyquist k = %.4f h/Mpc\n", π / (BLOCK * ALATT))
flush(stdout)

# ---------------- exact global reference (periodic N³ FFT, same kernels) ----------------
t0 = time()
noise = Array{Float32,3}(undef, N, N, N)
fill_noise_threefry!(noise, N, SEED)
noise_k = rfft(noise)
dk_ = copy(noise_k); MR._periodic_convolve!(dk_, pk, N, BOX)
dglob = irfft(dk_, N)
pk_ = copy(noise_k); MR._periodic_convolve!(pk_, pk, N, BOX; kernel_fn=MR._kernel_1lpt(1))
pglob = irfft(pk_, N)
dk_ = nothing; pk_ = nothing; noise_k = nothing
coreR = 1:NC                               # global indices of the tiled region (1..NC)
dglob_c = dglob[coreR, coreR, coreR]; pglob_c = pglob[coreR, coreR, coreR]
dglob = nothing; pglob = nothing; GC.gc()
@printf("global reference: σ_δ=%.4f  σ_ψx=%.4f Mpc/h  (%.0f s)\n", std(dglob_c), std(pglob_c), time() - t0)

# ---------------- coarse fields (shared by all tilings: depend on M only) ----------------
coarse_noise = MR._downsample_noise(N, M, SEED)
ck = rfft(coarse_noise)
# FIELD_COMP=1: coarse splice compensation D/T (MultiResolution._splice_compensation)
const COMP = get(ENV, "FIELD_COMP", "0") == "1" ? MR._splice_compensation(M, BLOCK) : nothing
@info "splice compensation" on=(COMP !== nothing)
tmp = copy(ck); MR._periodic_convolve!(tmp, pk, M, BOX; comp=COMP); dcoarse = irfft(tmp, M)
tmp = copy(ck); MR._periodic_convolve!(tmp, pk, M, BOX; kernel_fn=MR._kernel_1lpt(1), comp=COMP); pcoarse = irfft(tmp, M)
ck = nothing; tmp = nothing

# ---------------- 1. noise / residual identity across tilings ----------------
# Global cell g is covered by the tile it contains in every tiling; compare the tile-local
# arrays at a set of global cells (including buffer cells shared between neighbouring tiles).
function tile_of(g, nsub); it = clamp((g - 1) ÷ nsub + 1, 1, NC ÷ nsub); it; end
let nbad_noise = 0, nbad_res = 0, nchk = 0
    ref = Dict{NTuple{3,Int},Tuple{Float32,Float32}}()
    gcells = [(g, 7, NC ÷ 2 - 3) for g in unique((1, NBUFF + 1, NC ÷ 8, NC ÷ 8 + 1, NC ÷ 4, NC ÷ 4 + 1, NC ÷ 2, NC ÷ 2 + 1, 2NC ÷ 3, NC))]
    for nt in NTILES
        nsub = NC ÷ nt; nmesh = nsub + 2NBUFF
        cache = Dict{NTuple{3,Int},Tuple{Array{Float32,3},Array{Float32,3}}}()
        for g in gcells
            t = (tile_of(g[1], nsub), tile_of(g[2], nsub), tile_of(g[3], nsub))
            haskey(cache, t) || (cache[t] = (MR._generate_tile_noise(t..., nsub, nmesh, N, SEED),
                                             MR._generate_extended_residual(t..., nsub, nmesh, N, SEED, coarse_noise, M, 0)))
            nz, rs = cache[t]
            l = ntuple(d -> g[d] - ((t[d] - 1) * nsub + 1 - NBUFF) + 1, 3)
            v = (nz[l...], rs[l...]); nchk += 1
            if haskey(ref, g)
                v[1] === ref[g][1] || (nbad_noise += 1)
                v[2] === ref[g][2] || (nbad_res += 1)
            else
                ref[g] = v
            end
            @assert v[1] === noise[g...] "tile noise ≠ global noise at $g (ntile=$nt)"
        end
    end
    @printf("\n[1] noise identity: %d cell checks, noise mismatches %d, residual mismatches %d (tile noise also == global noise bitwise)\n",
            nchk, nbad_noise, nbad_res)
end
noise = nothing; GC.gc()

# ---------------- 2–5. synthesize each tiling and compare ----------------
# Stitched core fields; for ntile=4 also keep the coarse and self components for the F2 slice.
const SLICE_NT = 4 in NTILES ? 4 : NTILES[end]
const SLICE_K = min(50, NC)                           # z-plane (global index) of the cross-section
dmax = maximum(NC ÷ nt for nt in NTILES) ÷ 2
results = Dict{Int,Any}()
stitched = Dict{Int,Tuple{Array{Float32,3},Array{Float32,3}}}()

for nt in NTILES
    local t0 = time()
    nsub = NC ÷ nt; nmesh = nsub + 2NBUFF; boxl = nmesh * ALATT
    dst = Array{Float32,3}(undef, NC, NC, NC); pst = similar(dst)
    slc = nt == SLICE_NT ? zeros(Float32, NC, NC, 2) : nothing     # (coarse, self) at SLICE_K
    for kt in 1:nt, jt in 1:nt, it in 1:nt
        res = MR._generate_extended_residual(it, jt, kt, nsub, nmesh, N, SEED, coarse_noise, M, 0)
        dself = MR._isolated_convolve(res, pk, boxl, nmesh)
        dlong = MR._interpolate_to_tile(dcoarse, it, jt, kt, nsub, nmesh, N, M)
        pself = MR._isolated_convolve(res, pk, boxl, nmesh; kernel_fn=MR._kernel_1lpt(1))
        plong = MR._interpolate_to_tile(pcoarse, it, jt, kt, nsub, nmesh, N, M)
        c = NBUFF+1:NBUFF+nsub
        gi = (it-1)*nsub+1:it*nsub; gj = (jt-1)*nsub+1:jt*nsub; gk = (kt-1)*nsub+1:kt*nsub
        dst[gi, gj, gk] .= dself[c, c, c] .+ dlong[c, c, c]
        pst[gi, gj, gk] .= pself[c, c, c] .+ plong[c, c, c]
        if slc !== nothing && SLICE_K in gk
            lk = NBUFF + SLICE_K - gk[1] + 1
            slc[gi, gj, 1] .= dlong[c, c, lk]; slc[gi, gj, 2] .= dself[c, c, lk]
        end
    end
    stitched[nt] = (dst, pst)

    # rms error vs global, overall and vs distance to nearest tile face (core cells)
    ed = dst .- dglob_c; ep = pst .- pglob_c
    rd = sqrt(mean(abs2, ed)) / std(dglob_c); rp = sqrt(mean(abs2, ep)) / std(pglob_c)
    cd = cor(vec(dst), vec(dglob_c)); cp = cor(vec(pst), vec(pglob_c))
    prof_d = zeros(nsub ÷ 2); prof_p = zeros(nsub ÷ 2); cnt = zeros(Int, nsub ÷ 2)
    @inbounds for k in 1:NC, j in 1:NC, i in 1:NC
        d = minimum(min(mod(x - 1, nsub), nsub - 1 - mod(x - 1, nsub)) for x in (i, j, k)) + 1
        d > nsub ÷ 2 && continue
        prof_d[d] += ed[i, j, k]^2; prof_p[d] += ep[i, j, k]^2; cnt[d] += 1
    end
    prof_d = sqrt.(prof_d ./ max.(cnt, 1)) ./ std(dglob_c)
    prof_p = sqrt.(prof_p ./ max.(cnt, 1)) ./ std(pglob_c)
    results[nt] = (; nsub, nmesh, rd, rp, cd, cp, prof_d, prof_p, maxd=maximum(abs, ed) / std(dglob_c))
    @printf("\n[2] ntile=%d (nsub=%d, nmesh=%d): δ rms err/σ = %.4f%%  corr = %.7f  max|err|/σ = %.3f;  ψx rms err/σ = %.4f%%  corr = %.7f  (%.0f s)\n",
            nt, nsub, nmesh, 100rd, cd, results[nt].maxd, 100rp, cp, time() - t0)
    @printf("    δ err/σ vs distance to tile face (cells 1,2,4,8,16,%d): %s\n", nsub ÷ 2,
            join([@sprintf("%.4f%%", 100prof_d[min(d, end)]) for d in (1, 2, 4, 8, 16, nsub ÷ 2)], "  "))
    open(joinpath(OUT, "field_edge_profile_ntile$(nt).csv"), "w") do io
        println(io, "dist_cells,derr_over_sigma,psierr_over_sigma,ncell")
        for d in eachindex(cnt); cnt[d] > 0 && @printf(io, "%d,%.6e,%.6e,%d\n", d, prof_d[d], prof_p[d], cnt[d]); end
    end
    if slc !== nothing
        open(joinpath(OUT, "field_slice_ntile$(nt)_k$(SLICE_K).f32"), "w") do io
            write(io, slc[:, :, 1], slc[:, :, 2], dst[:, :, SLICE_K], dglob_c[:, :, SLICE_K])
        end
        open(joinpath(OUT, "field_slice_README.txt"), "w") do io
            println(io, "field_slice_ntile$(nt)_k$(SLICE_K).f32: 4 consecutive Float32 $(NC)x$(NC) column-major planes at global z-index $(SLICE_K):")
            println(io, "  1 coarse (interpolated δ_coarse), 2 self (isolated residual convolution), 3 total = 1+2, 4 exact global FFT δ")
            println(io, "  cellsize $(ALATT) Mpc/h; tile faces at multiples of nsub=$(NC ÷ nt) cells; coarse cell = $(BLOCK) cells")
        end
    end
    GC.gc()
end

# ---------------- 3. tiling vs tiling ----------------
println()
for a in NTILES, b in NTILES
    a < b || continue
    dd = stitched[a][1] .- stitched[b][1]; dp = stitched[a][2] .- stitched[b][2]
    @printf("[3] ntile %d vs %d: δ rms diff/σ = %.4f%%   ψx rms diff/σ = %.4f%%\n", a, b,
            100sqrt(mean(abs2, dd)) / std(dglob_c), 100sqrt(mean(abs2, dp)) / std(pglob_c))
end

# ---------------- 4. P(k) ratio and r(k) on the tiled region ----------------
# Same region (NC³, non-periodic), same cosine taper (width = nbuff cells) for all fields,
# so the window/leakage cancels in ratios. Physical P = V/N⁶ |δ_k|² (= (2π)³ × pk file).
function taper(n, w)
    t = ones(Float32, n)
    for i in 1:w; t[i] = t[n+1-i] = Float32(0.5 - 0.5cos(π * (i - 0.5) / w)); end
    t
end
const TP = taper(NC, NBUFF)
tapered(f) = f .* reshape(TP, :, 1, 1) .* reshape(TP, 1, :, 1) .* reshape(TP, 1, 1, :)
kf = 2π / (NC * ALATT); nbin = 48
kedges = exp.(range(log(kf), log(π / ALATT), length=nbin + 1))
function binned(fa, fb)
    A = rfft(tapered(fa)); B = rfft(tapered(fb))
    kx = FFTW.rfftfreq(NC, NC * kf); ky = FFTW.fftfreq(NC, NC * kf)
    Paa = zeros(nbin); Pbb = zeros(nbin); Pab = zeros(nbin); nk = zeros(Int, nbin); ks = zeros(nbin)
    @inbounds for iz in 1:NC, iy in 1:NC, ix in 1:size(A, 1)
        k = sqrt(kx[ix]^2 + ky[iy]^2 + ky[iz]^2)
        (k < kedges[1] || k >= kedges[end]) && continue
        b = searchsortedlast(kedges, k); w = (ix == 1 || ix == size(A, 1)) ? 1 : 2
        Paa[b] += w * abs2(A[ix, iy, iz]); Pbb[b] += w * abs2(B[ix, iy, iz])
        Pab[b] += w * real(A[ix, iy, iz] * conj(B[ix, iy, iz])); nk[b] += w; ks[b] += w * k
    end
    norm = (NC * ALATT)^3 / Float64(NC)^6 / mean(TP .^ 2)^3
    (; k=ks ./ max.(nk, 1), Pa=norm .* Paa ./ max.(nk, 1), Pb=norm .* Pbb ./ max.(nk, 1),
       r=Pab ./ sqrt.(Paa .* Pbb), nk)
end
println("\n[4] spectra of the stitched fields vs the exact global fields on the same tapered region")
kcn = π / (BLOCK * ALATT)
open(joinpath(OUT, "field_pk_ratio.csv"), "w") do io
    println(io, "field,ntile,k_hMpc,P_split,P_global,P_input_phys,ratio,transfer,r,nmodes")
    for (fname, fi, gref) in (("delta", 1, dglob_c), ("psix", 2, pglob_c)), nt in NTILES
        s = binned(stitched[nt][fi], gref)
        T = [s.r[b] * sqrt(s.Pa[b] / s.Pb[b]) for b in 1:nbin]      # Re⟨a b*⟩/⟨|b|²⟩
        for b in 1:nbin
            s.nk[b] > 0 || continue
            @printf(io, "%s,%d,%.6e,%.6e,%.6e,%.6e,%.8f,%.8f,%.8f,%d\n", fname, nt, s.k[b], s.Pa[b], s.Pb[b],
                    (2π)^3 * pk(s.k[b]), s.Pa[b] / s.Pb[b], T[b], s.r[b], s.nk[b])
        end
        @printf("  %-5s ntile=%d  k/k_Nc: ", fname, nt)
        for b in 1:nbin
            (s.nk[b] > 0 && (b % 4 == 0)) || continue
            @printf("%.2f:%.3f/%.4f ", s.k[b] / kcn, s.Pa[b] / s.Pb[b], s.r[b])
        end
        println("   (ratio/r)")
    end
    s = binned(dglob_c, dglob_c)
    @printf("  global/input P (0.05<k<1): %s\n",
            join([@sprintf("%.3f", s.Pa[b] / ((2π)^3 * pk(s.k[b]))) for b in 1:nbin if s.nk[b] > 0 && 0.05 < s.k[b] < 1][1:4:end], " "))
end

# σ(R) of the top-hat-smoothed δ for the filter-bank radii: the quantity that sets peak counts.
# Periodic FFT of the tapered NC³ region (identical treatment for both fields → ratio meaningful).
println("\n[4b] σ(R) top-hat, stitched / global")
filt = [parse(Float64, split(l)[3]) for l in readlines(joinpath(DATA, "filters_websky_finecell.dat"))[2:end] if !isempty(strip(l))]
Rs = filt[1:3:end]
function sigmaR(f, Rs)
    A = rfft(tapered(f)); kx = FFTW.rfftfreq(NC, NC * kf); ky = FFTW.fftfreq(NC, NC * kf)
    acc = zeros(length(Rs))
    @inbounds for iz in 1:NC, iy in 1:NC, ix in 1:size(A, 1)
        k = sqrt(kx[ix]^2 + ky[iy]^2 + ky[iz]^2); k == 0 && continue
        w = (ix == 1 || ix == size(A, 1)) ? 1 : 2; p = abs2(A[ix, iy, iz])
        for (q, R) in enumerate(Rs)
            x = k * R; W = 3 * (sin(x) - x * cos(x)) / x^3
            acc[q] += w * p * W^2
        end
    end
    sqrt.(acc)
end
sg = sigmaR(dglob_c, Rs)
open(joinpath(OUT, "field_sigmaR_ratio.csv"), "w") do io
    println(io, "ntile,R_Mpch,sigma_ratio")
    for nt in NTILES
        ss = sigmaR(stitched[nt][1], Rs)
        @printf("  ntile=%d: %s\n", nt, join([@sprintf("R=%.1f:%.4f", Rs[q], ss[q] / sg[q]) for q in eachindex(Rs)], " "))
        for q in eachindex(Rs); @printf(io, "%d,%.4f,%.6f\n", nt, Rs[q], ss[q] / sg[q]); end
    end
end

open(joinpath(OUT, "field_summary.csv"), "w") do io
    println(io, "ntile,nsub,nmesh,delta_rms_err_over_sigma,delta_corr,delta_max_err_over_sigma,psix_rms_err_over_sigma,psix_corr")
    for nt in NTILES
        r = results[nt]
        @printf(io, "%d,%d,%d,%.6e,%.9f,%.4e,%.6e,%.9f\n", nt, r.nsub, r.nmesh, r.rd, r.cd, r.maxd, r.rp, r.cp)
    end
end
println("\nwrote results to $OUT")
