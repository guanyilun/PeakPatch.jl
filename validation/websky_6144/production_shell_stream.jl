#!/usr/bin/env julia
# Streaming version (memory-light): bin the 6144^3 lightcone octant by z-shell,
# compare each shell's n(>1.69e12) to the verified z=0 finder density 1.975e-3/Mpc³.
using PeakPatch, Printf
import PeakPatch.Cosmology: CosmologyParams

cat = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_pkfix.pksc"
cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
chi2z = PeakPatch.build_chi_to_z(cosmo; z_max=5.0)
ox,oy,oz = -3850.0,-3850.0,-3850.0
rho_m = 2.775e11*0.31
Mfloor = 1.69e12
# RTHL floor (avoid M calc per halo): M>Mfloor <=> RTHL > Rfloor
Rfloor = (3*Mfloor/(4pi*rho_m))^(1/3)

zedges = collect(0.0:0.25:4.75); nshell=length(zedges)-1
cnt = zeros(Int, nshell)

function stream!(cnt, cat, chi2z, ox,oy,oz, Rfloor, zedges, nshell)
    io = open(cat,"r")
    N = read(io, Int32); RTHLmax = read(io, Float32); z_out = read(io, Float32)
    @info "streaming" N RTHLmax z_out Rfloor
    nf = 33; chunk = 1_000_000
    buf = Vector{Float32}(undef, chunk*nf)
    ztop = zedges[end]
    ndone = 0
    while ndone < N
        m = min(chunk, Int(N)-ndone)
        read!(io, view(buf, 1:m*nf))
        @inbounds for h in 1:m
            b0 = (h-1)*nf
            buf[b0+7] > Rfloor || continue
            x=Float64(buf[b0+1]); y=Float64(buf[b0+2]); z=Float64(buf[b0+3])
            r = sqrt((x-ox)^2+(y-oy)^2+(z-oz)^2)
            zz = PeakPatch.chi_to_z(chi2z, r)
            (zz < 0.0 || zz >= ztop) && continue
            bb = clamp(floor(Int, zz/0.25)+1, 1, nshell)
            cnt[bb] += 1
        end
        ndone += m
    end
    close(io)
    return Int(N)
end
Ntot = stream!(cnt, cat, chi2z, ox,oy,oz, Rfloor, zedges, nshell)

n0 = 65223/320.8333^3   # verified z=0 finder density at floor
chi(z) = z <= 0.001 ? 0.0 : (2.998e5/100.0)*sum(1/sqrt(0.31*(1+zz)^3+0.69)*0.002 for zz in 0.001:0.002:z)
@printf("\nverified z=0 finder n(>1.69e12)=%.3e ; total above floor=%d (=%.1fM)\n\n",
        n0, sum(cnt), sum(cnt)/1e6)
@printf("%-11s %-11s %-12s %-11s %-7s\n","z-shell","N(>floor)","V_oct","n(>floor)","n/n0")
for b in 1:nshell
    zlo,zhi=zedges[b],zedges[b+1]; rlo,rhi=chi(zlo),chi(zhi)
    V=(1/8)*(4pi/3)*(rhi^3-rlo^3); n=cnt[b]/V
    @printf("%.2f-%.2f   %-11d %-12.3e %-11.3e %.3f\n", zlo,zhi,cnt[b],V,n,n/n0)
end
