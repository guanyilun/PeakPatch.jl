#!/usr/bin/env julia
# Production-scale closure: read the real 6144^3 lightcone octant, bin by redshift
# shell, and compare each shell's n(>1.69e12) to the VERIFIED z=0 finder density
# (1.975e-3 /Mpc³ from the 256^3 snapshot). The lowest-z shell MUST match that if the
# full production assembly (multi-res GPU lightcone) is correct. Faster z-drop than
# the verified ievol evolution would expose a production-scale loss.
using PeakPatch, Printf
import PeakPatch.Cosmology: CosmologyParams

D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
cat = joinpath(D, "catalog_websky_6144_oct000_pkfix.pksc")
@info "reading production catalog (7GB)..."
halos, RTHLmax, z_out = read_pksc(cat)
@info "read" n=length(halos)

cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)  # Om total=0.31
chi2z = PeakPatch.build_chi_to_z(cosmo; z_max=5.0)
obs = (-3850.0, -3850.0, -3850.0)
rho_m = 2.775e11 * 0.31
Mfloor = 1.69e12

# z=0 finder anchor (verified 256^3 snapshot): N(>1.69e12)=65223 over box 320.833^3
n0_finder = 65223 / 320.8333^3
@printf("verified z=0 finder n(>1.69e12) = %.3e /Mpc³\n\n", n0_finder)

# bin halos above floor by redshift; accumulate count per shell
zedges = collect(0.0:0.25:4.75)
nshell = length(zedges)-1
cnt = zeros(Int, nshell)
for h in halos
    M = 4/3*pi*Float64(h.RTHL)^3*rho_m
    M > Mfloor || continue
    r = sqrt((Float64(h.x)-obs[1])^2 + (Float64(h.y)-obs[2])^2 + (Float64(h.z)-obs[3])^2)
    z = chi2z(r)
    (z < zedges[1] || z >= zedges[end]) && continue
    b = clamp(floor(Int,(z-zedges[1])/0.25)+1, 1, nshell)
    cnt[b] += 1
end

# comoving distance for shell volumes (octant = 1/8 sky)
chi(z) = (2.998e5/100.0) * sum(1/sqrt(0.31*(1+zz)^3+0.69)*0.001 for zz in 0.0005:0.001:z)
@printf("%-12s %-10s %-12s %-12s %-8s\n","z-shell","N(>floor)","V_oct Mpc³","n(>floor)","n/n0")
for b in 1:nshell
    zlo,zhi = zedges[b],zedges[b+1]
    rlo,rhi = chi(zlo),chi(zhi)
    V = (1/8)*(4pi/3)*(rhi^3-rlo^3)
    n = cnt[b]/V
    @printf("%.2f-%.2f    %-10d %-12.3e %-12.3e %.3f\n", zlo,zhi,cnt[b],V,n,n/n0_finder)
end
@printf("\nLowest-z shell n/n0 ≈ 1.0 => production assembly matches verified finder at z~0.\n")
