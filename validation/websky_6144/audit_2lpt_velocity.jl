#!/usr/bin/env julia
# Audit the 2LPT velocity term. v = a·H·f·(disp1 + c·disp2). merge_pkvd uses c=2 (EdS f₂=2f₁).
# Compute σ_vr for c ∈ {0 (1LPT only), 1, 2, -2} and see which matches Websky (~261).
#   c=0 isolates the 1LPT velocity; the spread shows the 2LPT contribution + its sign.
using Printf
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams, Dlinear_ab, Dlinear_tables
const rho_m=2.775e11*0.31; const hub=0.68
const COEFS=(0.0, 1.0, 2.0, -2.0); const NC=length(COEFS)
cat=get(ARGS,1,"/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_finecell_AM.pksc")
obs=parse(Float64,get(ARGS,2,"-2618.0"))
cosmo=CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo;z_max=6.0); gt=Dlinear_tables(cosmo)
Efac(a)=sqrt(0.31*a^-3+0.69)
const Medges=[10.0^l for l in 12.5:0.5:14.0]; const Mc=[sqrt(Medges[i]*Medges[i+1]) for i in 1:length(Medges)-1]; const NB=length(Mc)

# Websky reference σ_vr(M)
function wsky()
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc");read(io,Int32);read(io,Int32);read(io,Int32)
    seekstart(io);Nw=read(io,Int32);read(io,Int32);read(io,Int32);nf=10
    buf=Vector{Float32}(undef,Int(Nw)*nf);read!(io,buf);close(io)
    s=zeros(NB);s2=zeros(NB);n=zeros(Int,NB)
    bM(M)=(i=searchsortedfirst(Medges,M)-1;(1<=i<=NB) ? i : 0)
    for i in 1:Nw
        b=(i-1)*nf;R=Float64(buf[b+7])*hub;ib=bM(4/3*pi*rho_m*R^3);ib==0&&continue
        x=Float64(buf[b+1]);y=Float64(buf[b+2]);z=Float64(buf[b+3]);r=sqrt(x^2+y^2+z^2)
        vr=r>0 ? (x*buf[b+4]+y*buf[b+5]+z*buf[b+6])/r : 0.0;s[ib]+=vr;s2[ib]+=vr^2;n[ib]+=1
    end
    [n[i]>1 ? sqrt(s2[i]/n[i]-(s[i]/n[i])^2) : 0.0 for i in 1:NB]
end

function _range(buf,k0,k1,nf,obs,chi2z,gt)
    s=zeros(NC,NB);s2=zeros(NC,NB);n=zeros(Int,NB)
    bM(M)=(i=searchsortedfirst(Medges,M)-1;(1<=i<=NB) ? i : 0)
    for k in k0:k1
        b=(k-1)*nf;R=Float64(buf[b+7]);R>0||continue;ib=bM(4/3*pi*rho_m*R^3);ib==0&&continue
        q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3])
        ds11=Float64(buf[b+4]);ds12=Float64(buf[b+5]);ds13=Float64(buf[b+6])
        ds21=Float64(buf[b+8]);ds22=Float64(buf[b+9]);ds23=Float64(buf[b+10])
        rL=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2);z=rL>0 ? chi_to_z(chi2z,rL) : 0.0;a=1.0/(1.0+z)
        d11=ds11*a;d12=ds12*a;d13=ds13*a; d21=ds21*a^2;d22=ds22*a^2;d23=ds23*a^2
        ex=q1+d11+d21;ey=q2+d12+d22;ez=q3+d13+d23
        dx=ex-obs;dy=ey-obs;dz=ez-obs;rE=sqrt(dx^2+dy^2+dz^2);rE>0||continue
        zE=chi_to_z(chi2z,rE);aE=1.0/(1.0+zE);_,f,_=Dlinear_ab(aE,gt);vf=aE*100.0*Efac(aE)*f
        A=(dx*d11+dy*d12+dz*d13)/rE; B=(dx*d21+dy*d22+dz*d23)/rE   # 1LPT & 2LPT radial displacement
        for c in 1:NC
            vr=vf*(A+COEFS[c]*B); s[c,ib]+=vr; s2[c,ib]+=vr^2
        end
        n[ib]+=1
    end
    (s,s2,n)
end
function run(cat)
    nw=max(Threads.nthreads(),1);s=zeros(NC,NB);s2=zeros(NC,NB);n=zeros(Int,NB)
    io=open(cat);N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33;chunk=4_000_000;buf=Vector{Float32}(undef,chunk*nf);ndone=0
    @info "streaming with $nw workers..."
    while ndone<N
        m=min(chunk,Int(N)-ndone);read!(io,view(buf,1:m*nf))
        bnds=[1+div(m*(w-1),nw) for w in 1:nw+1];bnds[end]=m+1
        ts=[Threads.@spawn _range(buf,bnds[w],bnds[w+1]-1,nf,obs,chi2z,gt) for w in 1:nw]
        for t in ts;ls,ls2,ln=fetch(t);s.+=ls;s2.+=ls2;n.+=ln;end
        ndone+=m
    end
    close(io);(s,s2,n)
end
@info "Websky..."; W=wsky()
@info "auditing 2LPT velocity coefficient..." cat
s,s2,n=run(cat)
@printf("\n%-10s %-9s", "M", "σvr_wsky"); for c in COEFS; @printf(" | c=%-4g ratio", c); end; println()
for i in 1:NB
    @printf("%-10.2e %-9.1f", Mc[i], W[i])
    for c in 1:NC
        so=n[i]>1 ? sqrt(max(s2[c,i]/n[i]-(s[c,i]/n[i])^2,0.0)) : 0.0
        @printf(" | %-6.1f %-5.3f", so, so/max(W[i],1e-9))
    end
    println()
end
@printf("\nc=0 is 1LPT-only. Whichever c gives ratio~1 is the right 2LPT velocity coefficient/sign.\n")
