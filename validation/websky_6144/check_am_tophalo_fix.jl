#!/usr/bin/env julia
# Quantify the AM top-halo fix (build_abundance_table, 2026-09-24) on a FROZEN production
# octant: rebuild the table with the fixed code from the raw catalog and compare, halo by
# halo, with the committed _AM catalog (made with the pre-fix code). Writes the fixed AM
# catalog to scratch only — the frozen /project catalogs are never touched.
# Usage: julia --project=validation -t 32 validation/websky_6144/check_am_tophalo_fix.jl [OCT]
using PeakPatch, Printf
OCT = length(ARGS) >= 1 ? ARGS[1] : "000"
D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
bits = OCT; obs = Tuple(b == '1' ? 2618.0 : -2618.0 for b in reverse(bits))   # octZYX → (x,y,z)
raw, _, zout = read_pksc(joinpath(D, "catalog_websky_6144_prod_oct$(OCT).pksc"))
old, _, _    = read_pksc(joinpath(D, "catalog_websky_6144_prod_oct$(OCT)_AM.pksc"))
cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "data", "pk_websky_RAW_unnormalized.dat"))
table = build_abundance_table(raw, cosmo, pk; hmf=:tinker, z_max=4.5, obs=obs, fsky=1/8)
new = abundance_match(raw, table, cosmo; obs=obs)
rho = 2.775e11 * 0.31; M(h) = (4π/3) * rho * Float64(h.RTHL)^3
@assert length(new) == length(old)
dM = [M(n) / M(o) - 1 for (n, o) in zip(new, old)]
chg = findall(x -> abs(x) > 1e-4, dM)
@printf("octant %s obs=%s: %d halos; %d changed (|ΔM/M|>1e-4)\n", OCT, obs, length(new), length(chg))
Mo = M.(old[chg]); Mn = M.(new[chg])
for i in sortperm(Mn; rev=true)[1:min(15, end)]
    @printf("  M_old=%.3e  M_new=%.3e  ΔM/M=%+.3f\n", Mo[i], Mn[i], dM[chg[i]])
end
@printf("changed: M_new range %.2e–%.2e; median ΔM/M %+.3f; max %+.3f\n",
        minimum(Mn), maximum(Mn), sort(dM[chg])[cld(end, 2)], maximum(dM[chg]))
for Mc in (1e14, 5e14, 1e15)
    @printf("N(>%.0e): old %d  new %d\n", Mc, count(h -> M(h) > Mc, old), count(h -> M(h) > Mc, new))
end
out = "/home/yguan/scratch/websky_6144/catalog_prod_oct$(OCT)_AMfix.pksc"
write_pksc(out, new, Float32(maximum(h.RTHL for h in new)), Float32(zout)); @info "wrote" out
