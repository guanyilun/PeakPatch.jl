#!/usr/bin/env julia
# ours-vs-reference full-sky auto spectrum with the A3 error model.
#
#   julia --project=validation -t 32 compare_auto.jl <label> <ours.fits> <ref.fits> <lmax> \
#         [scale_ours] [scale_ref] [nside_common] [--jk]
#
# --jk adds the octant jackknife for BOTH maps (8 masked transforms each). Output:
# validation/paper/results/auto_<label>.txt with ℓ_eff, C_ours, C_ref, ratio, σ_ratio.
include(joinpath(@__DIR__, "spectra.jl"))

label, pours, pref = ARGS[1], ARGS[2], ARGS[3]
lmax = parse(Int, ARGS[4])
so = length(ARGS) >= 5 && !startswith(ARGS[5], "--") ? parse(Float64, ARGS[5]) : 1.0
sr = length(ARGS) >= 6 && !startswith(ARGS[6], "--") ? parse(Float64, ARGS[6]) : 1.0
nsc = length(ARGS) >= 7 && !startswith(ARGS[7], "--") ? parse(Int, ARGS[7]) : 0
JK = "--jk" in ARGS

t0 = time()
mo = loadmap(pours; scale=so, nside=nsc); mr = loadmap(pref; scale=sr, nside=nsc)
mo.resolution.nside == mr.resolution.nside || error("Nside differs: pass nside_common")
lmax <= 2mo.resolution.nside || @warn "lmax > 2·Nside"
@printf("means: ours %.4e  ref %.4e  (ratio %.4f)\n", mean(mo.pixels), mean(mr.pixels),
        mean(mo.pixels) / mean(mr.pixels))
edges = lbins(lmax)
co, le, nm = binned(clof(almof(mo, lmax)), edges)
cr, _, _ = binned(clof(almof(mr, lmax)), edges)
r = co ./ cr
σg = r .* sqrt.(4 ./ max.(nm, 1))                   # two independent realizations
σ = copy(σg); σjo = zeros(length(r)); σjr = zeros(length(r))
if JK
    S = octant_spectra(Dict("o" => mo, "r" => mr), [("o", "o"), ("r", "r")], lmax, edges)
    σjo = jk_sigma(S[("o", "o")]); σjr = jk_sigma(S[("r", "r")])
    σj = r .* sqrt.((σjo ./ co) .^ 2 .+ (σjr ./ cr) .^ 2)
    σ = max.(σg, σj)
end
@printf("\n%-8s %-12s %-12s %-8s %-8s %-8s\n", "ell", "C_ours", "C_ref", "ratio", "sig", "sig_g")
for b in eachindex(r)
    @printf("%-8.1f %-12.4e %-12.4e %-8.4f %-8.4f %-8.4f\n", le[b], co[b], cr[b], r[b], σ[b], σg[b])
end
write_table(joinpath(@__DIR__, "results", "auto_$(label).txt"),
            "ell_eff C_ours C_ref ratio sigma_ratio sigma_gauss sigma_jk_ours sigma_jk_ref Nmodes ($(pours) vs $(pref), lmax=$(lmax))",
            le, co, cr, r, σ, σg, σjo, σjr, nm)
@printf("done in %.1f min\n", (time() - t0) / 60)
