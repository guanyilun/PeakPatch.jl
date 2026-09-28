#!/usr/bin/env julia
# Full-sky maps vs theory (paper_comparison_plan A4), independent of the Websky maps:
#   κ (ours v2 total + field, Websky kap_lt4.5) vs Limber linear / Halofit, z<4.5
#   field kSZ (ours v2, and the frozen campaign for before/after the splice fix) vs exact
#   linear LOS Doppler (z_max 4.5, ℓ≤1000) + linear Ostriker-Vishniac.
# Theory tables: validation/paper_theory/results/ (THEORY_ANCHORS_2026-09-26.md).
# Gaussian pixel-window approximation applied to the maps (Nside 4096).
include(joinpath(@__DIR__, "spectra.jl"))
using DelimitedFiles
const D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const W = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"
const T = joinpath(@__DIR__, "..", "paper_theory", "results")
const C = get(ENV, "CAMPAIGN", "v2")               # campaign under test; kSZ "before" column = frozen (v2) or v2 (v3)
const CB = C == "v2" ? "prod" : C == "v3" ? "v2" : "v3"

kl = readdlm(joinpath(T, "kappa_limber.txt"); comments=true)          # ell, lin..., halofit...
kl_ell = kl[:, 1]; kl_lin45 = kl[:, 3]; kl_hf45 = kl[:, 6]
dop = readdlm(joinpath(T, "ksz_doppler.txt"); comments=true)          # ell, eq4.3(4.5), eq4.3(4.6), fullLOS, ...
ov = readdlm(joinpath(T, "ksz_ov.txt"); comments=true)                # ell, OV
interp(x, y, x0) = (i = clamp(searchsortedlast(x, x0), 1, length(x) - 1);
                    t = (x0 - x[i]) / (x[i+1] - x[i]); y[i] * (1 - t) + y[i+1] * t)
bin_theory(f, edges) = [sum((2l + 1) * f(l) for l in edges[b]:edges[b+1]-1) /
                        sum(2l + 1 for l in edges[b]:edges[b+1]-1) for b in 1:length(edges)-1]

edges = lbins(3000; lmin=20)
function measured(path; scale=1.0)
    m = loadmap(path; scale=scale); ns = m.resolution.nside
    c, le, nm = binned(clof(almof(m, 3000)), edges)
    pw = bin_theory(l -> pixwin2_gauss(l, ns), edges)
    c ./ pw, le, nm
end

# ---- κ ----
co, le, nm = measured(joinpath(D, "fullsky_$(C)", "kappa_lt4.5_$(C)_fullsky_nside4096.fits"))
cf, _, _ = measured(joinpath(D, "fullsky_$(C)", "kappa_field_$(C)_fullsky_nside4096.fits"))
cw, _, _ = measured(joinpath(W, "kap_lt4.5.fits"))
tlin = bin_theory(l -> interp(kl_ell, kl_lin45, l), edges)
thf = bin_theory(l -> interp(kl_ell, kl_hf45, l), edges)
σ = sqrt.(2 ./ nm)
@printf("\n== κ (z<4.5) vs Limber theory; fractional Gaussian error per map in last column\n")
@printf("%-7s %-11s %-11s %-11s %-12s %-12s %-6s\n", "ell", "ours/HF", "websky/HF", "field/lin", "ours/websky", "HF/lin", "σ")
for b in eachindex(le)
    @printf("%-7.0f %-11.3f %-11.3f %-11.3f %-12.3f %-12.3f %-6.3f\n", le[b], co[b] / thf[b],
            cw[b] / thf[b], cf[b] / tlin[b], co[b] / cw[b], thf[b] / tlin[b], σ[b])
end
write_table(joinpath(@__DIR__, "results", "theory_$(C)_kappa.txt"),
            "ell_eff C_ours C_field_ours C_websky C_halofit_0-4.5 C_linear_0-4.5 (maps pixwin-deconvolved, Gaussian approx)",
            le, co, cf, cw, thf, tlin)

# ---- field kSZ ----
e2 = lbins(1000; lmin=20)
edges = e2
theo = bin_theory(l -> 2π / (l * (l + 1)) * (interp(dop[:, 1], dop[:, 4], l) + interp(ov[:, 1], ov[:, 2], l)), e2)
cv2, le2, nm2 = measured(joinpath(D, "fullsky_$(C)", "ksz_field_uK_$(C)_fullsky_nside4096.fits"))
cp, _, _ = measured(joinpath(D, "fullsky_$(CB)", "ksz_field_uK_$(CB)_fullsky_nside4096.fits"))
@printf("\n== field kSZ vs linear (exact-LOS Doppler + OV), z<4.5\n")
@printf("%-7s %-12s %-14s %-14s %-6s\n", "ell", "D_theory", "$(C)/theory", "$(CB)/theory", "σ")
for b in eachindex(le2)
    @printf("%-7.0f %-12.4f %-14.3f %-14.3f %-6.3f\n", le2[b], theo[b] * le2[b] * (le2[b] + 1) / 2π,
            cv2[b] / theo[b], cp[b] / theo[b], sqrt(2 / nm2[b]))
end
write_table(joinpath(@__DIR__, "results", "theory_$(C)_ksz_field.txt"),
            "ell_eff C_$(C) C_$(CB) C_theory(Doppler_fullLOS+OV_lin) [uK^2]", le2, cv2, cp, theo)
