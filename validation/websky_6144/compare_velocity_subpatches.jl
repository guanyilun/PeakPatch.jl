#!/usr/bin/env julia
# Does a 100 deg² σ_vr measurement have ~20% cosmic variance? Carve OUR octant (corrected velocities)
# into ~54 patches of ~100 deg² and report the patch-to-patch scatter of σ_vr (M>1e13).
#   scatter ~20% -> the single Websky 10x10 patch at 0.77 is consistent with cosmic variance.
#   scatter ~5%  -> the 0.77 is a real systematic, not variance.
# Uses the SAME corrected velocity as apply_eulerian_velocity.jl (D-fix + merge_pkvd conversion).
using Printf, Statistics
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams, Dlinear_ab, Dlinear_tables
const rho_m=2.775e11*0.31
cat = get(ARGS,1,"/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_finecell_AM.pksc")
obs = parse(Float64,get(ARGS,2,"-2618.0")); alatt=parse(Float64,get(ARGS,3,"0.85221")); Mcut=parse(Float64,get(ARGS,4,"1e13"))
const nphi=9; const ncos=6; const npatch=nphi*ncos          # ~54 patches of ~100 deg² over the octant
cosmo=CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo;z_max=6.0); gt=Dlinear_tables(cosmo)
Efac(a)=sqrt(0.31*a^-3+0.69)

function _range(buf,k0,k1,nf,obs,alatt,chi2z,gt,Mcut)
    s=zeros(npatch);s2=zeros(npatch);n=zeros(Int,npatch)
    for k in k0:k1
        b=(k-1)*nf; R=Float64(buf[b+7]); R>0||continue; (4/3*pi*rho_m*R^3)>Mcut||continue
        q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3])
        ds11=Float64(buf[b+4]);ds12=Float64(buf[b+5]);ds13=Float64(buf[b+6])
        ds21=Float64(buf[b+8]);ds22=Float64(buf[b+9]);ds23=Float64(buf[b+10])
        rL=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2); z=rL>0 ? chi_to_z(chi2z,rL) : 0.0; a=1.0/(1.0+z)
        d11=ds11*a*alatt;d12=ds12*a*alatt;d13=ds13*a*alatt
        d21=ds21*a^2*alatt;d22=ds22*a^2*alatt;d23=ds23*a^2*alatt
        ex=q1+d11+d21;ey=q2+d12+d22;ez=q3+d13+d23
        dx=ex-obs;dy=ey-obs;dz=ez-obs; rE=sqrt(dx^2+dy^2+dz^2); rE>0||continue
        zE=chi_to_z(chi2z,rE); aE=1.0/(1.0+zE); _,f,_=Dlinear_ab(aE,gt); vf=aE*100.0*Efac(aE)*f
        vx=vf*(d11+2*d21);vy=vf*(d12+2*d22);vz=vf*(d13+2*d23); vr=(dx*vx+dy*vy+dz*vz)/rE
        nx=dx/rE;ny=dy/rE;nz=dz/rE; (nx>0&&ny>0&&nz>0)||continue   # +++ octant
        iphi=clamp(floor(Int,atan(ny,nx)/(pi/2/nphi)),0,nphi-1); ic=clamp(floor(Int,nz/(1.0/ncos)),0,ncos-1)
        p=iphi*ncos+ic+1; s[p]+=vr;s2[p]+=vr^2;n[p]+=1
    end
    (s,s2,n)
end
function run(cat)
    nw=max(Threads.nthreads(),1); s=zeros(npatch);s2=zeros(npatch);n=zeros(Int,npatch)
    io=open(cat);N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33;chunk=4_000_000;buf=Vector{Float32}(undef,chunk*nf);ndone=0
    @info "streaming with $nw workers..."
    while ndone<N
        m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
        bnds=[1+div(m*(w-1),nw) for w in 1:nw+1]; bnds[end]=m+1
        ts=[Threads.@spawn _range(buf,bnds[w],bnds[w+1]-1,nf,obs,alatt,chi2z,gt,Mcut) for w in 1:nw]
        for t in ts; ls,ls2,ln=fetch(t); s.+=ls;s2.+=ls2;n.+=ln; end
        ndone+=m
    end
    close(io); (s,s2,n)
end
@info "subpatch σ_vr..." cat Mcut
s,s2,n=run(cat)
sig=[n[p]>100 ? sqrt(max(s2[p]/n[p]-(s[p]/n[p])^2,0.0)) : NaN for p in 1:npatch]
valid=filter(!isnan,sig)
@printf("\n%d patches of ~100 deg² with >100 halos (M>%.0e). per-patch σ_vr (km/s):\n", length(valid), Mcut)
for (i,v) in enumerate(valid); @printf("%7.1f%s", v, i%9==0 ? "\n" : " "); end
m=mean(valid); sd=std(valid)
@printf("\n\nmean σ_vr = %.1f km/s,  scatter = %.1f km/s = %.1f%% ,  min=%.1f max=%.1f\n", m, sd, 100*sd/m, minimum(valid), maximum(valid))
@printf("Websky 10x10 patch σ_vr ≈ 261; ours octant-mean ≈ 200 (ratio 0.77).\n")
@printf("If patch scatter >~20%%, the single Websky patch at 0.77 is within cosmic variance. If ~5%%, it's systematic.\n")
