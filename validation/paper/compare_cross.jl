#!/usr/bin/env julia
# Cross-spectrum comparisons against the released Websky maps (paper_comparison_plan B1,
# Websky Figs 10–11 and Table 1): the SAME pairs measured on our full-sky maps and on
# theirs; ratio of cross spectra, and the correlation coefficient r_ℓ of each sky.
#
#   CAMPAIGN=prod|v2 julia --project=validation -t 32 compare_cross.jl
#
# Pairs: κ×CIB_ν (Fig 10; κ×φ has the same ours/ref ratio), y×CIB_ν and κ×y (Fig 11 —
# the y×(CIB+y) combination follows from these and the autos with the Table 2 factors),
# CIB_ν×CIB_ν' decoherence averaged over 150<ℓ<1000 (Table 1). y pairs at Nside 2048
# (the released y map), κ/CIB pairs at Nside 4096 (lmax 4096).
include(joinpath(@__DIR__, "spectra.jl"))
const C = get(ENV, "CAMPAIGN", "prod")
const F = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/fullsky_$(C)"
const W = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"
const NUS = ["0143", "0217", "0353", "0545", "0857"]

ours(n) = joinpath(F, "$(n)_$(C)_fullsky_nside4096.fits")
function alms_of(paths; nside=0, lmax)
    Dict(k => almof(loadmap(p; nside=nside), lmax) for (k, p) in paths)
end

function report(tag, AO, AR, pairs, edges; out)
    rows = []
    @printf("\n== %s\n%-8s", tag, "ell")
    for (a, b) in pairs; @printf(" %-22s", "$(a)x$(b): ratio  r_o  r_r"); end
    println()
    cl = Dict{Any,Any}()
    for s in (:o, :r), k in unique(vcat(first.(pairs), last.(pairs)))
        A = s == :o ? AO : AR
        cl[(s, k, k)] = binned(clof(A[k]), edges)[1]
    end
    for s in (:o, :r), (a, b) in pairs
        A = s == :o ? AO : AR
        cl[(s, a, b)] = binned(clof(A[a], A[b]), edges)[1]
    end
    le = binned(clof(first(values(AO))), edges)[2]
    cols = Any[le]; hdr = "ell_eff"
    for (a, b) in pairs
        ratio = cl[(:o, a, b)] ./ cl[(:r, a, b)]
        ro = cl[(:o, a, b)] ./ sqrt.(cl[(:o, a, a)] .* cl[(:o, b, b)])
        rr = cl[(:r, a, b)] ./ sqrt.(cl[(:r, a, a)] .* cl[(:r, b, b)])
        push!(cols, cl[(:o, a, b)], cl[(:r, a, b)], ratio, ro, rr)
        hdr *= " C_$(a)x$(b)_ours C_$(a)x$(b)_ref ratio r_ours r_ref"
    end
    for i in eachindex(le)
        @printf("%-8.0f", le[i])
        for j in 0:length(pairs)-1
            @printf(" %7.3f %6.3f %6.3f       ", cols[4+5j][i], cols[5+5j][i], cols[6+5j][i])
        end
        println()
    end
    write_table(out, hdr, cols...)
    cl
end

edges = lbins(4096; lmin=30)
# κ and CIB at Nside 4096
po = Dict("k" => ours("kappa_lt4.5")); pr = Dict("k" => joinpath(W, "kap_lt4.5.fits"))
for ν in NUS
    po["c$ν"] = ours("cib_nu$(ν)_wcut"); pr["c$ν"] = joinpath(W, "cib_nu$(ν).fits")
end
AO = alms_of(po; lmax=4096); AR = alms_of(pr; lmax=4096)
report("κ × CIB (Websky Fig 10)", AO, AR, [("k", "c$ν") for ν in NUS], edges;
       out=joinpath(@__DIR__, "results", "cross_$(C)_kappa_cib.txt"))
# decoherence (Table 1): mean r over 150<ℓ<1000 for all CIB pairs
cpairs = [("c$a", "c$b") for (i, a) in enumerate(NUS) for b in NUS[i+1:end]]
cl = report("CIB × CIB", AO, AR, cpairs, edges; out=joinpath(@__DIR__, "results", "cross_$(C)_cib_cib.txt"))
le = binned(clof(AO["k"]), edges)[2]; sel = findall(l -> 150 < l < 1000, le)
@printf("\nCIB decoherence <C^νν'/√(C^νν C^ν'ν')> over 150<ℓ<1000 (ours | Websky):\n")
for (a, b) in cpairs
    f(s) = mean((cl[(s, a, b)] ./ sqrt.(cl[(s, a, a)] .* cl[(s, b, b)]))[sel])
    @printf("  %s×%s: %.3f | %.3f\n", a[2:end], b[2:end], f(:o), f(:r))
end
AO = nothing; AR = nothing; GC.gc()

# y pairs at Nside 2048
edges2 = lbins(3000; lmin=30)
po = Dict("y" => ours("tsz_y"), "k" => ours("kappa_lt4.5"))
pr = Dict("y" => joinpath(W, "tsz_2048.fits"), "k" => joinpath(W, "kap_lt4.5.fits"))
for ν in NUS
    po["c$ν"] = ours("cib_nu$(ν)_wcut"); pr["c$ν"] = joinpath(W, "cib_nu$(ν).fits")
end
AO = alms_of(po; nside=2048, lmax=3000); AR = alms_of(pr; nside=2048, lmax=3000)
report("y × CIB and κ × y (Websky Fig 11 ingredients)", AO, AR,
       vcat([("y", "c$ν") for ν in NUS], [("k", "y")]), edges2;
       out=joinpath(@__DIR__, "results", "cross_$(C)_y.txt"))
