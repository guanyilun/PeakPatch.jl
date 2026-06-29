#!/usr/bin/env julia
# Nail the displacement->velocity conversion (Fortran merge_pkvd.f90:164-170):
#   v[km/s] = a·H(a)·f(a) · (disp1 + 2·disp2),  with H = 100·E(a) (the h cancels for disp in Mpc/h)
# Our catalog stores disp1 = vx,vy,vz (fields 4-6) and disp2 = vx2,vy2,vz2 (fields 8-10).
# Unknown: are disp in lattice units (×alatt to get Mpc/h) or already Mpc/h? Test both vs Websky σ_vr.
using Printf
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams
rho_m = 2.775e11*0.31; hub=0.68; alatt=0.85221
cosmo = CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo; z_max=6.0)
Efac(a)= sqrt(0.31*a^-3 + 0.69)               # E(a)=H/H0
fgr(a) = (0.31*a^-3/(0.31*a^-3+0.69))^0.55     # f≈Ωm(a)^0.55
D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
Medges=[10.0^l for l in 12.5:0.5:14.0]; Mc=[sqrt(Medges[i]*Medges[i+1]) for i in 1:length(Medges)-1]
binM(M)=(i=searchsortedfirst(Medges,M)-1; (1<=i<=length(Mc)) ? i : 0)

# Websky reference σ_vr(M)
function wsky()
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc"); Nw=read(io,Int32);read(io,Int32);read(io,Int32)
    nf=10; buf=Vector{Float32}(undef,Int(Nw)*nf); read!(io,buf); close(io)
    s=zeros(length(Mc)); s2=zeros(length(Mc)); n=zeros(Int,length(Mc))
    for i in 1:Nw
        b=(i-1)*nf; R=Float64(buf[b+7])*hub; ib=binM(4/3*pi*rho_m*R^3); ib==0 && continue
        x=Float64(buf[b+1]);y=Float64(buf[b+2]);z=Float64(buf[b+3]); r=sqrt(x^2+y^2+z^2)
        vr= r>0 ? (x*buf[b+4]+y*buf[b+5]+z*buf[b+6])/r : 0.0
        s[ib]+=vr; s2[ib]+=vr^2; n[ib]+=1
    end
    [n[i]>1 ? sqrt(s2[i]/n[i]-(s[i]/n[i])^2) : 0.0 for i in 1:length(Mc)]
end
# our converted σ_vr(M) under unit option `scale` (1.0 = disp already Mpc/h, alatt = lattice units)
function ours(cat, obs, scale)
    io=open(cat); N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33; chunk=1_000_000; buf=Vector{Float32}(undef,chunk*nf)
    s=zeros(length(Mc)); s2=zeros(length(Mc)); n=zeros(Int,length(Mc)); ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            b=(k-1)*nf; R=Float64(buf[b+7]); R>0||continue; ib=binM(4/3*pi*rho_m*R^3); ib==0&&continue
            # Eulerian position = q + (disp1+disp2)*scale
            d1=(Float64(buf[b+4]),Float64(buf[b+5]),Float64(buf[b+6]))
            d2=(Float64(buf[b+8]),Float64(buf[b+9]),Float64(buf[b+10]))
            ex=Float64(buf[b+1])+(d1[1]+d2[1])*scale; ey=Float64(buf[b+2])+(d1[2]+d2[2])*scale; ez=Float64(buf[b+3])+(d1[3]+d2[3])*scale
            dx=ex-obs; dy=ey-obs; dz=ez-obs; r=sqrt(dx^2+dy^2+dz^2)
            z = r>0 ? chi_to_z(chi2z,r) : 0.0; a=1.0/(1.0+z)
            fac = a*100.0*Efac(a)*fgr(a)*scale     # km/s per (disp unit)
            vx=fac*(d1[1]+2*d2[1]); vy=fac*(d1[2]+2*d2[2]); vz=fac*(d1[3]+2*d2[3])
            vr = r>0 ? (dx*vx+dy*vy+dz*vz)/r : 0.0
            s[ib]+=vr; s2[ib]+=vr^2; n[ib]+=1
        end
        ndone+=m
    end
    close(io)
    [n[i]>1 ? sqrt(s2[i]/n[i]-(s[i]/n[i])^2) : 0.0 for i in 1:length(Mc)]
end

@info "Websky..."; W=wsky()
@info "ours scale=alatt (disp in lattice units)..."; Oa=ours(joinpath(D,"catalog_websky_6144_oct000_finecell.pksc"),-2618.0, alatt)
@info "ours scale=1 (disp already Mpc/h)..."; O1=ours(joinpath(D,"catalog_websky_6144_oct000_finecell.pksc"),-2618.0, 1.0)
@printf("\n%-10s %-10s | %-10s %-7s | %-10s %-7s\n","M","σvr_wsky","σvr(×alatt)","ratio","σvr(×1)","ratio")
for i in 1:length(Mc)
    @printf("%-10.2e %-10.1f | %-10.1f %-7.3f | %-10.1f %-7.3f\n", Mc[i],W[i],Oa[i],Oa[i]/max(W[i],1e-9),O1[i],O1[i]/max(W[i],1e-9))
end
@printf("\nWhichever scale gives ratio~1 is the correct disp unit. Then v=a·100·E·f·scale·(disp1+2disp2).\n")
