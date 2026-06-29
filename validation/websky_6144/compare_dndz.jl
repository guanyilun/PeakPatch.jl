#!/usr/bin/env julia
# Separate factor 1 (solid angle / normalization) from factor 2 (radial/z completeness).
# dN/dz per deg² for our octant vs the Websky patch, plus octant angular coverage & max r.
#  - ratio flat across z  -> uniform normalization (solid angle / global density)
#  - ratio falls at high z -> radial truncation / z-incompleteness
using PeakPatch, Printf
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams

D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
cosmo = CosmologyParams(0.31+0.049, 0.049, 0.69, 0.68, 0.965, 0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0)
rho_m = 2.775e11*0.31; hub=0.68
Mcut = parse(Float64, get(ARGS, 1, "3e12"))   # pass 0 for AM-invariant total-count dN/dz
zedges = collect(0.0:0.4:4.8)
zc = [(zedges[i]+zedges[i+1])/2 for i in 1:length(zedges)-1]

binz!(h,z) = (i=searchsortedfirst(zedges,z)-1; (1<=i<=length(zc)) && (h[i]+=1.0))

# our octant: stream, obs at corner, M from RTHL, also track angular signs & max r
function our_dndz(cat, obs, Mcut)
    io=open(cat); N=read(io,Int32); read(io,Float32); read(io,Float32)
    nf=33; chunk=2_000_000; buf=Vector{Float32}(undef,chunk*nf)
    h=zeros(length(zc)); ndone=0; rmax=0.0; nq=zeros(Int,8); ntot=0
    while ndone<N
        m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            R=Float64(buf[(k-1)*nf+7]); R>0 || continue
            M=4/3*pi*rho_m*R^3
            dx=Float64(buf[(k-1)*nf+1])-obs[1]; dy=Float64(buf[(k-1)*nf+2])-obs[2]; dz=Float64(buf[(k-1)*nf+3])-obs[3]
            r=sqrt(dx^2+dy^2+dz^2); r>rmax && (rmax=r)
            # octant sign histogram (which of 8 octants the direction falls in)
            q=1+(dx<0)+2*(dy<0)+4*(dz<0); nq[q]+=1; ntot+=1
            M>Mcut || continue
            binz!(h, r>0 ? chi_to_z(chi2z,r) : 0.0)
        end
        ndone+=m
    end
    close(io); h, rmax, nq, ntot
end

# websky patch: obs at origin, R(field7)*h
function wsky_dndz(wf, Mcut)
    io=open(wf); Nw=read(io,Int32); read(io,Int32); read(io,Int32)
    nf=10; buf=Vector{Float32}(undef,Int(Nw)*nf); read!(io,buf); close(io)
    h=zeros(length(zc))
    for i in 1:Nw
        R=Float64(buf[(i-1)*nf+7])*hub; M=4/3*pi*rho_m*R^3; M>Mcut || continue
        # Websky positions are in Mpc (comoving); our chi2z table is in Mpc/h -> convert (*hub).
        x=Float64(buf[(i-1)*nf+1]); y=Float64(buf[(i-1)*nf+2]); z=Float64(buf[(i-1)*nf+3])
        r=sqrt(x^2+y^2+z^2)*hub; binz!(h, r>0 ? chi_to_z(chi2z,r) : 0.0)
    end
    h
end

@info "our RAW octant dN/dz..."
ho, rmax, nq, ntot = our_dndz(joinpath(D,"catalog_websky_6144_oct000_pkfix.pksc"), (-3850.0,-3850.0,-3850.0), Mcut)
@info "our AM octant dN/dz..."
ha, _, _, _ = our_dndz(joinpath(D,"catalog_websky_6144_oct000_pkfix_AM.pksc"), (-3850.0,-3850.0,-3850.0), Mcut)
@info "websky patch dN/dz..."
hw = wsky_dndz("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc", Mcut)

oct_deg2 = 41253.0/8; patch_deg2 = 100.0
@printf("\nmax r in our octant = %.0f Mpc/h  (chi(z=4.6) ~ 5300)\n", rmax)
@printf("octant sign distribution: octant1(+++)=%.2f%%  (others should be ~0)\n", 100*nq[1]/max(1,ntot))
@printf("\n%-8s %-11s %-11s %-11s %-9s %-9s\n","z","RAW/deg²","AM/deg²","wsky/deg²","RAW/wsky","AM/wsky")
for i in 1:length(zc)
    nr=ho[i]/oct_deg2; na=ha[i]/oct_deg2; nw=hw[i]/patch_deg2
    @printf("%-8.2f %-11.3f %-11.3f %-11.3f %-9.3f %-9.3f\n", zc[i], nr, na, nw, nr/max(nw,1e-9), na/max(nw,1e-9))
end
@printf("\nM>%.1e cut. If AM/wsky still crashes at high z -> genuine lightcone count deficit (not mass calib).\n", Mcut)
