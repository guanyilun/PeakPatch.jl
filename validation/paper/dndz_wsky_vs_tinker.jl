#!/usr/bin/env julia
# Which side is off in the Tier-A dN/dz M>1e13 at z 2.25–3 (Websky/ours 0.90–0.93, unchanged
# v3fs→v4fs; V4_RESULTS_2026-10-01.md)? Our AM catalogs match Tinker08 there to ~1% (A4), so compare
# Websky's patch (inscribed cap, as tierA_v3fs.jl) and our 8-cap mean directly with the Tinker08
# expectation per steradian in the same bins.
using Printf, PeakPatch
import PeakPatch.Cosmology: CosmologyParams, chi, build_chi_to_z, chi_to_z, growth_factor
import PeakPatch.MassFunction: precompute_sigma, tinker_dndlnM
const C = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81); const c2z = build_chi_to_z(C; z_max=6.0)
pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky_RAW_unnormalized.dat"))
Mg = 10 .^ range(12, 16; length=801); lM = log.(Mg); sg = precompute_sigma(Mg, pk, 0.31)
function nsr(Mc, za, zb)                       # Tinker08 N(>Mc) per steradian in za ≤ z < zb
    ra, rb = chi(za, C), chi(zb, C); ns = 40; acc = 0.0
    for j in 1:ns
        r0 = ra + (j - 1) * (rb - ra) / ns; r1 = ra + j * (rb - ra) / ns
        z = chi_to_z(c2z, (r0 + r1) / 2); D = growth_factor(z, C); n = 0.0
        for i in 1:length(Mg)-1
            Mm = sqrt(Mg[i] * Mg[i+1]); Mm < Mc && continue
            n += tinker_dndlnM(Mm, D * sqrt(sg[i] * sg[i+1]), (log(sg[i+1]) - log(sg[i])) / (lM[i+1] - lM[i]), z, 0.31) * (lM[i+1] - lM[i])
        end
        acc += n * (r1^3 - r0^3) / 3
    end
    acc
end
# Websky / ours per sr from the Tier-A output (dN/dz M>1e13 rows)
rows = filter(l -> occursin("dN/dz M>1e+13", l), readlines(joinpath(@__DIR__, "results", "tierA_v4fs.txt")))
@printf("%-14s %10s %10s %10s %8s %8s\n", "z bin", "Websky", "ours", "Tinker", "W/Tk", "ours/Tk")
for l in rows
    m = match(r"z ([\d.]+)-([\d.]+)\s+Websky\s+(\S+) \| ours\s+(\S+)", l); m === nothing && continue
    za, zb = parse(Float64, m[1]), parse(Float64, m[2]); w, o = parse(Float64, m[3]), parse(Float64, m[4])
    t = nsr(1e13, za, zb)
    @printf("z %.2f-%.2f  %10.4g %10.4g %10.4g %8.3f %8.3f\n", za, zb, w, o, t, w / t, o / t)
end
