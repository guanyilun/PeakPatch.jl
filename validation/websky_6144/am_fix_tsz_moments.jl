#!/usr/bin/env julia
# tSZ-weighted impact of the AM top-halo fix on a frozen production octant:
# Σ M^(5/3) (∝ mean Compton-y) and Σ M^(10/3) (∝ 1-halo tSZ power), old vs fixed, per z-bin.
using PeakPatch, Printf
D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
old, _, _ = read_pksc(joinpath(D, "catalog_websky_6144_prod_oct000_AM.pksc"))
new, _, _ = read_pksc("/home/yguan/scratch/websky_6144/catalog_prod_oct000_AMfix.pksc")
cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81); c2z = PeakPatch.Cosmology.build_chi_to_z(cosmo; z_max=6.0)
rho = 2.775e11 * 0.31; M(h) = (4π/3) * rho * Float64(h.RTHL)^3; o = -2618.0
zed(h) = PeakPatch.Cosmology.chi_to_z(c2z, sqrt((h.x-o)^2 + (h.y-o)^2 + (h.z-o)^2))
edges = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.5]
S = zeros(length(edges)-1, 4)                  # old 5/3, new 5/3, old 10/3, new 10/3
for (a, b) in zip(old, new)
    ma = M(a); ma < 1e13 && M(b) < 1e13 && continue
    iz = clamp(searchsortedlast(edges, zed(a)), 1, length(edges)-1)
    S[iz, 1] += ma^(5/3); S[iz, 2] += M(b)^(5/3); S[iz, 3] += ma^(10/3); S[iz, 4] += M(b)^(10/3)
end
@printf("%-10s %14s %14s\n", "z bin", "ΣM^5/3 new/old", "ΣM^10/3 new/old")
for i in 1:size(S, 1)
    @printf("%4.1f–%-4.1f  %14.5f %14.5f\n", edges[i], edges[i+1], S[i,2]/S[i,1], S[i,4]/S[i,3])
end
T = sum(S; dims=1); @printf("all z      %14.5f %14.5f\n", T[2]/T[1], T[4]/T[3])
