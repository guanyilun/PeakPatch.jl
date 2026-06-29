#!/usr/bin/env julia
# Implement the missing merge_pkvd stage (Eulerian position + km/s velocity) with the CORRECT
# growth factor, as a post-process on the existing raw catalog. Verifies σ_vr vs Websky inline.
#
# Two fixes vs the stored catalog:
#  (1) D_pk bug: stored disp uses D/a (Dlinear_ab 3rd return); true displacement uses D (1st return).
#      Our stored disp1 = Sbar·(D/a) [lattice]; true 1LPT disp = Sbar·D = stored·a; 2LPT ∝ D² so ×a².
#  (2) merge_pkvd: x_E = q + (disp1+disp2)[Mpc/h];  v[km/s] = a·100·E(a)·f(a)·(disp1 + 2·disp2)[Mpc/h]
#      (H=100h·E, the h cancels for disp in Mpc/h).
# args: [1]=in catalog, [2]=out catalog (optional; if omitted, verify only), [3]=obs, [4]=alatt
using Printf, LinearAlgebra
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams, Dlinear_ab, Dlinear_tables
rho_m = 2.775e11*0.31; hub=0.68
cat   = get(ARGS,1, "/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_finecell.pksc")
obs   = parse(Float64, get(ARGS,3,"-2618.0"))
alatt = parse(Float64, get(ARGS,4,"0.85221"))
cosmo = CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0); gt = Dlinear_tables(cosmo)
Om, OL = cosmo.Om, cosmo.OL
Efac(a)=sqrt(Om*a^-3+OL)
Om_a(a)=Om*a^-3/(Om*a^-3+OL)

# Websky reference σ_vr(M)
Medges=[10.0^l for l in 12.5:0.5:14.0]; Mc=[sqrt(Medges[i]*Medges[i+1]) for i in 1:length(Medges)-1]
binM(M)=(i=searchsortedfirst(Medges,M)-1; (1<=i<=length(Mc)) ? i : 0)
function wsky()
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc");read(io,Int32);read(io,Int32);read(io,Int32)
    seekstart(io); Nw=read(io,Int32);read(io,Int32);read(io,Int32); nf=10
    buf=Vector{Float32}(undef,Int(Nw)*nf); read!(io,buf); close(io)
    s=zeros(length(Mc));s2=zeros(length(Mc));n=zeros(Int,length(Mc))
    for i in 1:Nw
        b=(i-1)*nf;R=Float64(buf[b+7])*hub;ib=binM(4/3*pi*rho_m*R^3);ib==0&&continue
        x=Float64(buf[b+1]);y=Float64(buf[b+2]);z=Float64(buf[b+3]);r=sqrt(x^2+y^2+z^2)
        vr=r>0 ? (x*buf[b+4]+y*buf[b+5]+z*buf[b+6])/r : 0.0; s[ib]+=vr;s2[ib]+=vr^2;n[ib]+=1
    end
    [n[i]>1 ? sqrt(s2[i]/n[i]-(s[i]/n[i])^2) : 0.0 for i in 1:length(Mc)]
end

@info "Websky..."; W=wsky()
# Process one contiguous halo range [k0,k1] of the buffer -> (s,s2,n). Task-local accumulators.
function _range(buf, k0, k1, nf, obs, alatt, chi2z, gt, Efac, Medges, rho_m, nb)
    s=zeros(nb); s2=zeros(nb); n=zeros(Int,nb)
    bM(M)=(i=searchsortedfirst(Medges,M)-1; (1<=i<=nb) ? i : 0)
    for k in k0:k1
        b=(k-1)*nf; R=Float64(buf[b+7]); R>0 || continue
        ib=bM(4/3*pi*rho_m*R^3); ib==0 && continue
        q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3])
        ds11=Float64(buf[b+4]);ds12=Float64(buf[b+5]);ds13=Float64(buf[b+6])   # stored = Sbar·(D/a) [lattice]
        ds21=Float64(buf[b+8]);ds22=Float64(buf[b+9]);ds23=Float64(buf[b+10])  # stored 2LPT ∝ (D/a)²
        rL=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2)
        z = rL>0 ? chi_to_z(chi2z,rL) : 0.0; a=1.0/(1.0+z)
        # fix (1): true disp (Mpc/h): d1=Sbar·D=stored·a; d2 ∝ D² -> ×a². Sbar is ALREADY Mpc/h
        # (kernel kᵢ/k²), NOT lattice cells -> do NOT multiply by alatt (RTHL is cell-based, Sbar is not).
        d11=ds11*a; d12=ds12*a; d13=ds13*a
        d21=ds21*a^2; d22=ds22*a^2; d23=ds23*a^2
        # fix (2): Eulerian position + velocity (exact f from Dlinear_ab[2])
        ex=q1+d11+d21; ey=q2+d12+d22; ez=q3+d13+d23
        rE=sqrt((ex-obs)^2+(ey-obs)^2+(ez-obs)^2); zE=rE>0 ? chi_to_z(chi2z,rE) : 0.0; aE=1.0/(1.0+zE)
        _, f_aE, _ = Dlinear_ab(aE, gt)
        vfac = aE*100.0*Efac(aE)*f_aE
        vx=vfac*(d11+2*d21); vy=vfac*(d12+2*d22); vz=vfac*(d13+2*d23)
        dx=ex-obs;dy=ey-obs;dz=ez-obs; vr=rE>0 ? (dx*vx+dy*vy+dz*vz)/rE : 0.0
        s[ib]+=vr; s2[ib]+=vr^2; n[ib]+=1
    end
    (s,s2,n)
end
function process(cat, obs, alatt, chi2z, gt, Efac, Mc, Medges, rho_m)
    nb=length(Mc); nw=max(Threads.nthreads(),1)
    s=zeros(nb); s2=zeros(nb); n=zeros(Int,nb)
    io=open(cat); N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33; chunk=4_000_000; buf=Vector{Float32}(undef,chunk*nf); ndone=0
    @info "streaming with $nw workers..."
    while ndone<N
        m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
        # split [1,m] into nw contiguous ranges, one @spawn each with its own accumulator
        bnds=[1 + div((m)*(w-1),nw) for w in 1:nw+1]; bnds[end]=m+1
        tasks=[Threads.@spawn _range(buf, bnds[w], bnds[w+1]-1, nf, obs, alatt, chi2z, gt, Efac, Medges, rho_m, nb) for w in 1:nw]
        for t in tasks
            ls,ls2,ln=fetch(t); s.+=ls; s2.+=ls2; n.+=ln
        end
        ndone+=m
    end
    close(io)
    (s,s2,n)
end
@info "applying Eulerian+velocity conversion to ours..." cat
s, s2, n = process(cat, obs, alatt, chi2z, gt, Efac, Mc, Medges, rho_m)
@printf("\n%-10s %-10s %-10s %-7s\n","M","σvr_wsky","σvr_ours","ratio")
for i in 1:length(Mc)
    so=n[i]>1 ? sqrt(s2[i]/n[i]-(s[i]/n[i])^2) : 0.0
    @printf("%-10.2e %-10.1f %-10.1f %-7.3f\n", Mc[i], W[i], so, so/max(W[i],1e-9))
end
@printf("\nratio~1 across mass -> corrected velocity (D-fix + merge_pkvd conversion) matches Websky.\n")
