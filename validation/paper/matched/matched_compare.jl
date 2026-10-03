#!/usr/bin/env julia
# Same-field comparison (MATCHED_FORTRAN_2026-10.md): Fortran hpkvd + merge_pkvd vs Julia, all on ONE linear
# field (fields/Fvec_jmatched), so cosmic variance cancels in every ratio.
#   F    : Fortran raw → Fortran merge_pkvd            (the reference chain)
#   FJ   : Fortran raw → Julia merge (exclusion + volume reduction)   (isolates the merge)
#   JE   : Julia exact global field (CPU run_multitile) → Julia merge   (isolates the finder vs F/FJ)
#   JS   : Julia production multires split (GPU), same noise → Julia merge   (isolates the split vs JE)
# Statistics per catalog, after rank-matching each to the Tinker08 abundance of the analysis region:
#   * raw-peak cross-match F_raw ↔ JE_raw (Lagrangian cell + filter), R_TH agreement
#   * merged halo-by-halo match (Lagrangian position within 1 cell) by mass
#   * cross bias with the shared linear field, Lagrangian and Eulerian positions, per mass bin (k < 0.1 h/Mpc),
#     8-subvolume jackknife errors on the ratios
#   * ξ(3–15 Mpc/h) for M > 5e12 (natural estimator against randoms in the same region)
#   usage: julia --project=validation -t 32 matched_compare.jl <rundir>
using PeakPatch, FFTW, Printf, Statistics, Random, DelimitedFiles
import PeakPatch.Cosmology: CosmologyParams, growth_factor
import PeakPatch.MassFunction: precompute_sigma, tinker_dndlnM

const RUN = ARGS[1]
# geometry: defaults = the 1056³ matched box; MC_N / MC_NB / MC_A override (e.g. the 1536³ production-tile box)
const N = parse(Int, get(ENV, "MC_N", "1056")); const NB = parse(Int, get(ENV, "MC_NB", "26"))
const A = parse(Float64, get(ENV, "MC_A", string(258.2207 / 303))); const L = N * A
const RHO_M = 2.775e11 * 0.31
mass(R) = 4 / 3 * π * RHO_M * Float64(R)^3
gidx(x) = x / A + (N + 1) / 2                         # fine-grid (1-based, cell-centred) coordinate of position x
# analysis region: coarse (4-cell) cells C0.. inside the non-periodic cores NB+1..N−NB (N=1056: cells 8..257 =
# fine 29..1028 inside 27..1030)
const CB = 4; const C0 = 8; const NC = (N - 2 * (CB * (C0 - 1))) ÷ CB
const XLO = (CB * (C0 - 1) + 0.5 - (N + 1) / 2) * A; const XHI = XLO + NC * CB * A
const VREG = (XHI - XLO)^3
inreg(x, y, z) = XLO <= x < XHI && XLO <= y < XHI && XLO <= z < XHI
const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
const ZS = 0.7; const DZ = growth_factor(ZS, COSMO)

struct Cat
    xL::Vector{Float64}; yL::Vector{Float64}; zL::Vector{Float64}
    xE::Vector{Float64}; yE::Vector{Float64}; zE::Vector{Float64}
    R::Vector{Float64}; Rf::Vector{Float64}
end
# pksc in the Julia/hpkvd raw layout (Lagrangian + 1LPT/2LPT displacements)
function cat_lagr(halos)
    h = filter(q -> q.RTHL > 0, halos)
    Cat([Float64(q.x) for q in h], [Float64(q.y) for q in h], [Float64(q.z) for q in h],
        [Float64(q.x + q.vx + q.vx2) for q in h], [Float64(q.y + q.vy + q.vy2) for q in h], [Float64(q.z + q.vz + q.vz2) for q in h],
        [Float64(q.RTHL) for q in h], [Float64(q.Rf) for q in h])
end
# merge_pkvd output: Eulerian x (1–3), velocity (4–6), R (7), Lagrangian xL (8–10), harray (Rf at col 28)
cat_mpkvd(h) = Cat([Float64(q.vx2) for q in h], [Float64(q.vy2) for q in h], [Float64(q.vz2) for q in h],
                   [Float64(q.x) for q in h], [Float64(q.y) for q in h], [Float64(q.z) for q in h],
                   [Float64(q.RTHL) for q in h], [Float64(q.Rf) for q in h])
sub(c::Cat, s) = Cat(c.xL[s], c.yL[s], c.zL[s], c.xE[s], c.yE[s], c.zE[s], c.R[s], c.Rf[s])

# ---- Tinker08 abundance in the region, rank matching ----
function tinker_counts()
    pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "..", "..", "websky_6144", "data", "pk_websky_RAW_unnormalized.dat"))
    Mg = 10 .^ range(11.5, 16; length=1801); lnMg = log.(Mg); sg = precompute_sigma(Mg, pk, 0.31)
    dn = [tinker_dndlnM(sqrt(Mg[i] * Mg[i+1]), DZ * sqrt(sg[i] * sg[i+1]), (log(sg[i+1]) - log(sg[i])) / (lnMg[i+1] - lnMg[i]), ZS, 0.31) *
          (lnMg[i+1] - lnMg[i]) for i in 1:length(Mg)-1]
    Mg, reverse(cumsum(reverse(dn))) .* VREG
end
const MG, NGT = tinker_counts()
Mof(k) = begin i = findlast(>=(k), NGT); i === nothing ? MG[1] : MG[i] end
function am_masses(c::Cat)                            # rank-matched masses of the halos inside the region (Lagrangian)
    s = findall(i -> inreg(c.xL[i], c.yL[i], c.zL[i]), eachindex(c.R))
    o = sortperm(c.R[s]; rev=true); M = zeros(length(s))
    for (k, i) in enumerate(o); M[i] = Mof(k); end
    sub(c, s), M
end

# ---- shared linear field on the 4-cell analysis grid (block average), and its FFT ----
function linear_field()
    fn = joinpath(RUN, "fields", "Fvec_jmatched"); @assert filesize(fn) == 4N^3
    d = Array{Float32}(undef, N, N, N); open(io -> read!(io, d), fn)
    g = zeros(Float64, NC, NC, NC); g0 = CB * (C0 - 1)
    Threads.@threads for K in 1:NC
        for J in 1:NC, I in 1:NC
            s = 0.0
            for dk in 1:CB, dj in 1:CB, di in 1:CB
                s += d[g0+CB*(I-1)+di, g0+CB*(J-1)+dj, g0+CB*(K-1)+dk]
            end
            g[I, J, K] = s / CB^3
        end
    end
    g .- mean(g)
end
function deposit(x, y, z)                              # NGP counts on the analysis grid → overdensity
    g = zeros(Float64, NC, NC, NC); h = CB * A; n = 0
    @inbounds for i in eachindex(x)
        inreg(x[i], y[i], z[i]) || continue
        I = clamp(floor(Int, (x[i] - XLO) / h) + 1, 1, NC); J = clamp(floor(Int, (y[i] - XLO) / h) + 1, 1, NC)
        K = clamp(floor(Int, (z[i] - XLO) / h) + 1, 1, NC); g[I, J, K] += 1; n += 1
    end
    g ./ (n / NC^3) .- 1
end
# cross bias Σ Re(H D*) / Σ |D|² over 0 < k < KMAX, total and per jackknife-deleted octant of the region
const KMAX = 0.1
function kmask(n, Lb)
    kf = 2π / Lb; m = falses(n ÷ 2 + 1, n, n)
    for k3 in 1:n, k2 in 1:n, k1 in 1:n÷2+1
        q1 = k1 - 1; q2 = k2 <= n ÷ 2 + 1 ? k2 - 1 : k2 - 1 - n; q3 = k3 <= n ÷ 2 + 1 ? k3 - 1 : k3 - 1 - n
        kk = kf * sqrt(q1^2 + q2^2 + q3^2); m[k1, k2, k3] = 0 < kk < KMAX
    end
    m
end
const KM = kmask(NC, NC * CB * A)
function bx(H, D)
    num = 0.0; den = 0.0
    @inbounds for i in eachindex(KM)
        KM[i] || continue; num += real(H[i] * conj(D[i])); den += abs2(D[i])
    end
    num / den
end
# jackknife: delete one of 8 octant sub-cubes (zero both fields there), recompute
const HALF = NC ÷ 2
function bx_jk(h, d)
    vals = Float64[]
    for o in 0:7
        hh = copy(h); dd = copy(d)
        r1 = (o & 1 == 0 ? (1:HALF) : (HALF+1:NC)); r2 = (o & 2 == 0 ? (1:HALF) : (HALF+1:NC)); r3 = (o & 4 == 0 ? (1:HALF) : (HALF+1:NC))
        hh[r1, r2, r3] .= 0; dd[r1, r2, r3] .= 0
        push!(vals, bx(rfft(hh), rfft(dd)))
    end
    vals
end
jkerr(v) = sqrt((length(v) - 1) / length(v) * sum(abs2, v .- mean(v)))

# ---- ξ by pair counts in the region (chaining mesh, non-periodic), natural estimator with randoms ----
function dd_counts(x, y, z, edges)
    nb = length(edges) - 1; rmax = edges[end]; nc = max(3, floor(Int, (XHI - XLO) / rmax)); cs = (XHI - XLO) / nc
    cell(v) = clamp(floor(Int, (v - XLO) / cs) + 1, 1, nc)
    heads = zeros(Int, nc, nc, nc); nxt = zeros(Int, length(x))
    for i in eachindex(x); c = (cell(x[i]), cell(y[i]), cell(z[i])); nxt[i] = heads[c...]; heads[c...] = i; end
    le = log10(edges[1]); dl = (log10(rmax) - le) / nb
    function worker(i0, i1)
        h = zeros(nb)
        @inbounds for i in i0:i1
            ci = (cell(x[i]), cell(y[i]), cell(z[i]))
            for a in -1:1, b in -1:1, c in -1:1
                (1 <= ci[1] + a <= nc && 1 <= ci[2] + b <= nc && 1 <= ci[3] + c <= nc) || continue
                j = heads[ci[1]+a, ci[2]+b, ci[3]+c]
                while j > 0
                    if j > i
                        d2 = (x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2
                        if 0 < d2 < rmax^2
                            k = floor(Int, (0.5log10(d2) - le) / dl) + 1
                            1 <= k <= nb && (h[k] += 1)
                        end
                    end
                    j = nxt[j]
                end
            end
        end
        h
    end
    n = length(x); nw = Threads.nthreads(); bn = [1 + div(n * (w - 1), nw) for w in 1:nw+1]; bn[end] = n + 1
    reduce(+, fetch.([Threads.@spawn worker(bn[w], bn[w+1] - 1) for w in 1:nw]))
end
const EDGES = collect(10.0 .^ range(0, log10(50); length=11)); const RC = [sqrt(EDGES[i] * EDGES[i+1]) for i in 1:10]
const WX = [i for i in 1:10 if 3 <= RC[i] <= 15]
const RR = begin
    rng = MersenneTwister(7); nr = 2_000_000
    rx = XLO .+ (XHI - XLO) .* rand(rng, nr); ry = XLO .+ (XHI - XLO) .* rand(rng, nr); rz = XLO .+ (XHI - XLO) .* rand(rng, nr)
    dd_counts(rx, ry, rz, EDGES) ./ (nr * (nr - 1) / 2)
end
function xi(x, y, z)
    s = findall(i -> inreg(x[i], y[i], z[i]), eachindex(x)); n = length(s)
    dd_counts(x[s], y[s], z[s], EDGES) ./ (n * (n - 1) / 2) ./ RR .- 1
end

# ---- halo-by-halo match on Lagrangian position (within `tol` Mpc/h), via a cell hash ----
function match(a::Cat, b::Cat; tol=A)
    key(x, y, z) = (floor(Int, x / tol), floor(Int, y / tol), floor(Int, z / tol))
    H = Dict{NTuple{3,Int},Vector{Int}}()
    for j in eachindex(b.R); push!(get!(H, key(b.xL[j], b.yL[j], b.zL[j]), Int[]), j); end
    m = zeros(Int, length(a.R))
    for i in eachindex(a.R)
        k = key(a.xL[i], a.yL[i], a.zL[i]); best = 0; bd = tol^2
        for p in -1:1, q in -1:1, r in -1:1, j in get(H, (k[1] + p, k[2] + q, k[3] + r), Int[])
            d2 = (a.xL[i] - b.xL[j])^2 + (a.yL[i] - b.yL[j])^2 + (a.zL[i] - b.zL[j])^2
            d2 <= bd && (bd = d2; best = j)
        end
        m[i] = best
    end
    m
end

function main()
    rd(f) = read_pksc(joinpath(RUN, f))[1]
    Fraw = rd("output/fortran_raw.pksc.12345")
    fr = filter(q -> q.RTHL > 0, Fraw)
    @info "Fortran raw" total = length(Fraw) kept = length(fr)
    cats = Dict{String,Cat}()
    cats["F"] = cat_mpkvd(rd("output/fortran_merge.pksc.12345"))
    cats["FJ"] = cat_lagr(merge_catalog(fr; volume_reduction=true))
    cats["JE"] = cat_lagr(rd("julia/julia_exact_merged.pksc"))
    isfile(joinpath(RUN, "julia/julia_split_merged.pksc")) && (cats["JS"] = cat_lagr(rd("julia/julia_split_merged.pksc")))
    labs = filter(l -> haskey(cats, l), ["F", "FJ", "JE", "JS"])
    out = open(joinpath(@__DIR__, "..", "results", "matched_fortran.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    say(@sprintf("same-field box: N=%d, cell %.5f Mpc/h, z=%.2f; analysis region %.1f Mpc/h cube (core interior)", N, A, ZS, XHI - XLO))

    # -- 0. grid convention check: Lagrangian positions sit on cell centres in both codes
    for (lab, c) in (("F_raw", cat_lagr(fr)), ("JE_raw", cat_lagr(rd("julia/julia_exact_raw.pksc"))))
        fx = [abs(gidx(x) - round(gidx(x))) for x in c.xL[1:min(end, 100000)]]
        say(@sprintf("%-7s n=%d  max |offset from cell centre| = %.3g cells", lab, length(c.R), maximum(fx)))
    end

    # -- 1. raw peaks: F_raw vs JE_raw, matched by Lagrangian cell (+ same filter)
    Fr = cat_lagr(fr); Jr = cat_lagr(rd("julia/julia_exact_raw.pksc"))
    m = match(Fr, Jr; tol=0.5A)
    ok = findall(>(0), m); sameRf = count(i -> abs(Fr.Rf[i] - Jr.Rf[m[i]]) < 1e-3, ok)
    say("\n== 1. raw peaks (R_TH>0): Fortran ", length(Fr.R), "  Julia exact ", length(Jr.R))
    say(@sprintf("  Fortran peaks with a Julia peak in the same cell: %.4f  (same filter: %.4f)", length(ok) / length(Fr.R), sameRf / length(Fr.R)))
    rr = [Jr.R[m[i]] / Fr.R[i] for i in ok]
    say(@sprintf("  R_TH Julia/Fortran over matched: median %.5f  mean %.5f  |Δ|>1%%: %.4f  |Δ|>10%%: %.4f",
                 median(rr), mean(rr), count(r -> abs(r - 1) > 0.01, rr) / length(rr), count(r -> abs(r - 1) > 0.1, rr) / length(rr)))
    dE = [sqrt((Jr.xE[m[i]] - Fr.xE[i])^2 + (Jr.yE[m[i]] - Fr.yE[i])^2 + (Jr.zE[m[i]] - Fr.zE[i])^2) for i in ok]
    dsp = [sqrt((Fr.xE[i] - Fr.xL[i])^2 + (Fr.yE[i] - Fr.yL[i])^2 + (Fr.zE[i] - Fr.zL[i])^2) for i in ok]
    say(@sprintf("  Eulerian position difference: median %.4f Mpc/h (displacement median %.3f)", median(dE), median(dsp)))
    for (lo, hi) in ((5e12, 1.3e13), (1.3e13, 5e13), (5e13, 1e16))
        s = filter(i -> lo <= mass(Fr.R[i]) < hi, ok); isempty(s) && continue
        r2 = [Jr.R[m[i]] / Fr.R[i] for i in s]
        say(@sprintf("  M_F %.1e–%.1e: n=%d  R ratio median %.4f  p05 %.4f  p95 %.4f  frac(R_J > 1.02 R_F) %.4f  frac(R_J < 0.98 R_F) %.4f",
                     lo, hi, length(s), median(r2), quantile(r2, 0.05), quantile(r2, 0.95), count(>(1.02), r2) / length(r2), count(<(0.98), r2) / length(r2)))
    end
    cap = count(i -> Fr.R[i] <= 1.75 * Fr.Rf[i] + 0.01 && Jr.R[m[i]] > 1.75 * Fr.Rf[i] + A, ok)
    say(@sprintf("  matched peaks where Fortran R_TH sits at the 1.75 R_f hunt cap and Julia is larger: %.4f", cap / length(ok)))

    # -- 2. merged catalogs, AM, halo-by-halo
    say("\n== 2. merged catalogs (rank-matched to Tinker08 in the region)")
    AM = Dict(l => am_masses(cats[l]) for l in labs)
    for l in labs
        c, M = AM[l]
        say(@sprintf("  %-3s merged %d, in region %d, M>5e12: %d, raw R_TH of the 5e12 rank: %.3f Mpc/h", l, length(cats[l].R), length(c.R),
                     count(>(5e12), M), minimum(c.R[M.>5e12])))
    end
    MB = [5e12, 8e12, 1.3e13, 2e13, 3.2e13, 7.9e13, 2.5e14]
    for (a, b) in (("F", "FJ"), ("F", "JE"), ("JE", "JS"))
        (haskey(AM, a) && haskey(AM, b)) || continue
        ca, Ma = AM[a]; cb, Mb = AM[b]; mm = match(ca, cb)
        say("  match $a → $b (Lagrangian position within 1 cell), fraction of $a halos with a $b counterpart:")
        for i in 1:length(MB)-1
            s = findall(m -> MB[i] <= m < MB[i+1], Ma); isempty(s) && continue
            f = count(j -> mm[j] > 0, s) / length(s)
            fb = count(j -> mm[j] > 0 && MB[i] <= Mb[mm[j]] < MB[i+1], s) / length(s)
            say(@sprintf("    M %.1e–%.1e: n=%d  matched %.4f  same bin %.4f", MB[i], MB[i+1], length(s), f, fb))
        end
    end

    # -- 3. cross bias with the shared linear field
    say("\n== 3. cross bias with the linear field, k < $(KMAX) h/Mpc (b_L: Lagrangian positions, b_E: Eulerian); δ_lin scaled to z=$(ZS)")
    dlin = linear_field() .* DZ; Dk = rfft(dlin)
    B = Dict{Tuple{String,Int,Symbol},Tuple{Float64,Vector{Float64}}}()
    for l in labs
        c, M = AM[l]
        for i in 0:length(MB)-1
            s = i == 0 ? findall(>(5e12), M) : findall(m -> MB[i] <= m < MB[i+1], M)
            for (sym, X, Y, Z) in ((:L, c.xL, c.yL, c.zL), (:E, c.xE, c.yE, c.zE))
                h = deposit(X[s], Y[s], Z[s])
                B[(l, i, sym)] = (bx(rfft(h), Dk), bx_jk(h, dlin))
            end
        end
    end
    binlab(i) = i == 0 ? "M>5e12      " : @sprintf("%.1e–%.1e", MB[i], MB[i+1])
    for sym in (:L, :E)
        say("  b_$(sym):            " * join([@sprintf("%-8s", l) for l in labs]) * "   ratios (jackknife σ)")
        for i in 0:length(MB)-1
            vals = join([@sprintf("%-8.4f", B[(l, i, sym)][1]) for l in labs])
            rats = String[]
            for (a, b) in (("F", "FJ"), ("F", "JE"), ("JE", "JS"))
                (haskey(cats, a) && haskey(cats, b)) || continue
                va, ja = B[(a, i, sym)]; vb, jb = B[(b, i, sym)]
                push!(rats, @sprintf("%s/%s %.4f±%.4f", a, b, va / vb, jkerr(ja ./ jb)))
            end
            say("  " * binlab(i) * "  " * vals * "   " * join(rats, "  "))
        end
    end

    # -- 4. ξ(3–15) for M > 5e12 (Eulerian)
    say("\n== 4. ξ(3–15 Mpc/h), M>5e12, Eulerian, region randoms")
    X = Dict{String,Vector{Float64}}()
    for l in labs
        c, M = AM[l]; s = findall(>(5e12), M); X[l] = xi(c.xE[s], c.yE[s], c.zE[s])
        say(@sprintf("  %-3s ξ(3–15) %.4f   ξ(1–3) %.4f", l, mean(X[l][WX]), mean(X[l][[i for i in 1:10 if 1 <= RC[i] <= 3]])))
    end
    for (a, b) in (("F", "FJ"), ("F", "JE"), ("JE", "JS"))
        (haskey(X, a) && haskey(X, b)) || continue
        say(@sprintf("  %s/%s ξ(3–15) %.4f   by r: %s", a, b, mean(X[a][WX]) / mean(X[b][WX]),
                     join([@sprintf("%.3f", X[a][i] / X[b][i]) for i in 1:10], " ")))
    end
    say("  [full-sky Websky/ours v4fs: ξ(3–15) 0.948 ± 0.002; b ratio ≈ 0.95–0.96 below 8e13]")
    close(out)
end
if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
