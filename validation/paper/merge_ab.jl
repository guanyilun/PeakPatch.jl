#!/usr/bin/env julia
# A/B test of the merge step: Lagrangian exclusion only (every catalog so far) vs exclusion +
# Fortran 'shared' volume reduction (merge_pkvd.f90:139-140; missing in our merge). Question: does
# the missing reduction explain why v3fs is ~7% more biased than Websky at fixed number density,
# with ~17% more close pairs (TIERA_V3FS_2026-09-28.md, tierA_diag.jl)?
#
# One periodic snapshot box at z = 0.7 (the Tier-A shell redshift), v3 physics: N = 1536
# (ntile 4, n 434, nbuff 25, periodic cores), cell 0.852213 Mpc/h, cf 32. The halo finder runs
# ONCE; the same raw catalog is merged both ways, finalized, and rank-matched to the Tinker08
# abundance of the box. Periodic natural-estimator ξ (no randoms).
#   usage: julia --project=validation -t 32 validation/paper/merge_ab.jl <config.toml>
using TOML, Printf, Statistics, DelimitedFiles
using PeakPatch
using CUDA
import PeakPatch.Cosmology: CosmologyParams, growth_factor
import PeakPatch.MassFunction: precompute_sigma, tinker_dndlnM

const RHO_M = 2.775e11 * 0.31
mass(R) = 4 / 3 * π * RHO_M * Float64(R)^3
Rof(M) = cbrt(3M / (4π * RHO_M))

# ---- periodic pair counts (chaining mesh) ----
function pairs_periodic(x, y, z, L, edges; wrapfn=nothing)
    nb = length(edges) - 1; rmax = edges[end]; nc = max(3, floor(Int, L / rmax)); cs = L / nc
    cell(v) = mod(floor(Int, v / cs), nc) + 1
    heads = zeros(Int, nc, nc, nc); nxt = zeros(Int, length(x))
    for i in eachindex(x); c = (cell(x[i]), cell(y[i]), cell(z[i])); nxt[i] = heads[c...]; heads[c...] = i; end
    le = log10(edges[1]); dl = (log10(rmax) - le) / nb
    function worker(i0, i1)
        h = zeros(nb)
        @inbounds for i in i0:i1
            ci = (cell(x[i]), cell(y[i]), cell(z[i]))
            for a in -1:1, b in -1:1, c in -1:1
                j = heads[mod1(ci[1] + a, nc), mod1(ci[2] + b, nc), mod1(ci[3] + c, nc)]
                while j > 0
                    if j > i
                        dx = x[i] - x[j]; dx -= L * round(dx / L)
                        dy = y[i] - y[j]; dy -= L * round(dy / L)
                        dz = z[i] - z[j]; dz -= L * round(dz / L)
                        d2 = dx^2 + dy^2 + dz^2
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
function xi_periodic(x, y, z, L, edges)
    DD = pairs_periodic(x, y, z, L, edges); n = length(x)
    [DD[k] / (n * (n - 1) / 2 * (4π / 3) * (edges[k+1]^3 - edges[k]^3) / L^3) - 1 for k in eachindex(DD)]
end
function upairs_periodic(x, y, z, R, L, uedges)
    h = zeros(length(uedges) - 1); Rm = maximum(R); nc = max(3, floor(Int, L / (2Rm * uedges[end]))); cs = L / nc
    cell(v) = mod(floor(Int, v / cs), nc) + 1
    heads = zeros(Int, nc, nc, nc); nxt = zeros(Int, length(x))
    for i in eachindex(x); c = (cell(x[i]), cell(y[i]), cell(z[i])); nxt[i] = heads[c...]; heads[c...] = i; end
    @inbounds for i in eachindex(x)
        ci = (cell(x[i]), cell(y[i]), cell(z[i]))
        for a in -1:1, b in -1:1, c in -1:1
            j = heads[mod1(ci[1] + a, nc), mod1(ci[2] + b, nc), mod1(ci[3] + c, nc)]
            while j > 0
                if j > i
                    dx = x[i] - x[j]; dx -= L * round(dx / L); dy = y[i] - y[j]; dy -= L * round(dy / L)
                    dz = z[i] - z[j]; dz -= L * round(dz / L)
                    u = sqrt(dx^2 + dy^2 + dz^2) / (R[i] + R[j]); k = searchsortedlast(uedges, u)
                    1 <= k < length(uedges) && (h[k] += 1)
                end
                j = nxt[j]
            end
        end
    end
    h ./ length(x)
end

function main()
    cfgd = TOML.parsefile(ARGS[1]); cfg = PipelineConfig(cfgd); rc = cfgd["run"]
    ntile = rc["ntile"]; zsn = Float64(cfg.z_out)
    nsub, N = grid_layout(cfg, ntile); L = N * cfg.boxsize / cfg.n
    @info "merge A/B snapshot" N L zsn
    raw = run_multitile_split(cfg; ntile=ntile, seed=rc["seed"], coarse_factor=rc["coarse_factor"], use_gpu=true,
                              devices=collect(0:length(CUDA.devices())-1), verbose=true)
    @info "raw halos" length(raw)
    cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
    cats = Dict{String,Any}()
    VARIANTS = (("A_excl", false, false), ("B_excl+red", true, false), ("C_red+fortran_ties", true, true),
                ("D_red+fortran_cap", true, false))
    # D: emulate the local Fortran get_homel outward-search bug (peakvoidsubs.f90:447 reuses jp, so the
    # search at :487-489 never runs): R_TH capped at the hunt radius ≈ 1.75 R_f before the merge.
    capR(h) = typeof(h)(ntuple(i -> i == 7 ? min(h.RTHL, Float32(1.75) * h.Rf) : getfield(h, i), fieldcount(typeof(h)))...)
    rawcap = capR.(raw)
    @info "Fortran-cap emulation" changed = count(i -> rawcap[i].RTHL < raw[i].RTHL, eachindex(raw)) of = length(raw)
    for (lab, red, ties) in VARIANTS
        src = startswith(lab, "D_") ? rawcap : raw
        m = merge_catalog(src; verbose=true, volume_reduction=red, fortran_ties=ties)
        cats[lab] = finalize_eulerian(m, cosmo, (0.0, 0.0, 0.0); ievol=0, z_out=zsn)
    end
    # Tinker08 N(>M) in the box at zsn, for rank matching
    pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky_RAW_unnormalized.dat"))
    Mg = 10 .^ range(11.5, 16; length=1801); lnMg = log.(Mg); sg = precompute_sigma(Mg, pk, 0.31); Dg = growth_factor(zsn, cosmo)
    dn = [tinker_dndlnM(sqrt(Mg[i] * Mg[i+1]), Dg * sqrt(sg[i] * sg[i+1]), (log(sg[i+1]) - log(sg[i])) / (lnMg[i+1] - lnMg[i]), zsn, 0.31) * (lnMg[i+1] - lnMg[i]) for i in 1:length(Mg)-1]
    Ngt = reverse(cumsum(reverse(dn))) .* L^3                 # N(>Mg[i]) for i = 1..end-1
    Mof(k) = begin i = findlast(>=(k), Ngt); i === nothing ? Mg[1] : Mg[i] end
    # linear ξ_mm and Tinker10 bias
    pkd = readdlm(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky.dat"); comments=true, comment_char='#')
    KK = Float64.(pkd[:, 1]); PK = Float64.(pkd[:, 2]) .* (2π)^3
    spec(f) = sum(0.5 * (f(i) + f(i + 1)) * (KK[i+1] - KK[i]) for i in 1:length(KK)-1)
    Wth(x) = x < 1e-3 ? 1.0 - x^2 / 10 : 3 * (sin(x) - x * cos(x)) / x^3
    sigma2(R) = spec(i -> PK[i] * Wth(KK[i] * R)^2 * KK[i]^2) / (2π^2)
    xi_mm(r) = Dg^2 * spec(i -> PK[i] * KK[i]^2 * sin(KK[i] * r) / (KK[i] * r)) / (2π^2)
    yT = log10(200.0); AT = 1 + 0.24yT * exp(-(4 / yT)^4); aT = 0.44yT - 0.88; CT = 0.019 + 0.107yT + 0.19exp(-(4 / yT)^4)
    tb(M) = begin ν = 1.686 / (Dg * sqrt(sigma2(Rof(M)))); 1 - AT * ν^aT / (ν^aT + 1.686^aT) + 0.183ν^1.5 + CT * ν^2.4 end
    edges = collect(10.0 .^ range(0, log10(50); length=11)); rcs = [sqrt(edges[i] * edges[i+1]) for i in 1:10]
    wb = [i for i in 1:10 if 6 <= rcs[i] <= 18]; wx = [i for i in 1:10 if 3 <= rcs[i] <= 15]; w13 = [i for i in 1:10 if 1 <= rcs[i] <= 3]
    ximm = mean(xi_mm(rcs[i]) for i in wb)
    uedges = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    MB = [5e12, 8e12, 1.3e13, 2e13, 3.2e13, 7.9e13, 2.5e14]
    out = open(joinpath(@__DIR__, "results", get(ENV, "MAB_OUT", "merge_ab.txt")), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    say(@sprintf("snapshot z=%.2f, N=%d, L=%.1f Mpc/h, raw halos %d; <ξ_mm,lin>(6–18)=%.4f", zsn, N, L, length(raw), ximm))
    res = Dict{String,Any}()
    for (lab, _, _) in VARIANTS
        h = cats[lab]; Mr = [mass(q.RTHL) for q in h]; o = sortperm(Mr; rev=true)
        Mam = similar(Mr); for (k, i) in enumerate(o); Mam[i] = Mof(k); end
        x = [mod(Float64(q.x), L) for q in h]; y = [mod(Float64(q.y), L) for q in h]; z = [mod(Float64(q.z), L) for q in h]
        s5 = findall(>(5e12), Mam)
        xi5 = xi_periodic(x[s5], y[s5], z[s5], L, edges)
        up = upairs_periodic(x[s5], y[s5], z[s5], Rof.(Mam[s5]), L, uedges)
        bs = Float64[]
        for i in 1:length(MB)-1
            s = findall(m -> MB[i] <= m < MB[i+1], Mam)
            push!(bs, sqrt(max(mean(xi_periodic(x[s], y[s], z[s], L, edges)[wb]), 0) / ximm))
        end
        res[lab] = (n=length(h), n5=length(s5), xi=xi5, up=up, b=bs, Mr=Mr, o=o)
        say(@sprintf("\n[%s] halos %d, AM M>5e12: %d; ξ(3–15) %.4f  ξ(1–3) %.4f", lab, length(h), length(s5), mean(xi5[wx]), mean(xi5[w13])))
        for i in 1:length(MB)-1
            say(@sprintf("  b(M %.1e-%.1e) = %.3f   Tinker10 %.3f   b/Tk %.3f", MB[i], MB[i+1], bs[i], tb(sqrt(MB[i] * MB[i+1])), bs[i] / tb(sqrt(MB[i] * MB[i+1]))))
        end
    end
    A = res["A_excl"]; B = res["B_excl+red"]; Cv = res["C_red+fortran_ties"]; Dv = res["D_red+fortran_cap"]
    say("\n== D/B (Fortran R_TH cap at 1.75 R_f, on top of reduction) — residual to explain: Websky/ours full sky ξ(3–15) 0.948 ± 0.002")
    say(@sprintf("ξ(3–15) D/B %.4f   ξ(1–3) D/B %.4f", mean(Dv.xi[wx]) / mean(B.xi[wx]), mean(Dv.xi[w13]) / mean(B.xi[w13])))
    for i in 1:length(MB)-1; say(@sprintf("b(M %.1e-%.1e) D/B %.4f", MB[i], MB[i+1], Dv.b[i] / B.b[i])); end
    say("\n== C/B (Fortran double reduction of equal-radius pairs, on top of reduction) — residual to explain: Websky/ours full sky ξ(3–15) 0.948 ± 0.002")
    say(@sprintf("ξ(3–15) C/B %.4f   ξ(1–3) C/B %.4f", mean(Cv.xi[wx]) / mean(B.xi[wx]), mean(Cv.xi[w13]) / mean(B.xi[w13])))
    for i in 1:length(MB)-1; say(@sprintf("b(M %.1e-%.1e) C/B %.4f", MB[i], MB[i+1], Cv.b[i] / B.b[i])); end
    for k in 1:length(uedges)-1; say(@sprintf("close pairs/halo u %.2f-%.2f: C/B %.3f", uedges[k], uedges[k+1], Cv.up[k] / B.up[k])); end
    say("\n== B/A (reduction on / off) — compare the Websky/ours Tier-A ratios in brackets")
    say(@sprintf("ξ(3–15) B/A %.3f   [W/ours 0.877 ± 0.029]", mean(B.xi[wx]) / mean(A.xi[wx])))
    say(@sprintf("ξ(1–3)  B/A %.3f   [W/ours 0.79]", mean(B.xi[w13]) / mean(A.xi[w13])))
    for k in 1:length(uedges)-1
        say(@sprintf("close pairs/halo u %.2f-%.2f: A %.4f B %.4f  B/A %.3f", uedges[k], uedges[k+1], A.up[k], B.up[k], B.up[k] / A.up[k]))
    end
    say("  [W/ours in the same u bins: 0.827 0.838 0.894 0.928 0.933 0.954 0.975]")
    for i in 1:length(MB)-1
        say(@sprintf("b(M %.1e-%.1e) B/A %.3f", MB[i], MB[i+1], B.b[i] / A.b[i]))
    end
    say("  [W/ours b(M) in the tierA_diag bins: 0.923 0.989 0.930 1.017 0.950 1.080]")
    close(out)
end
main()
