#!/usr/bin/env julia
# Halo bias b(M): the 2-halo amplitude that sets the large-ℓ tSZ/CIB/lensing map spectra.
# In matched geometry (inscribed angular cap), per mass bin:
#   relative   b_ours(M)/b_wsky(M) = sqrt(<ξ_ours>/<ξ_wsky>)   [theory-free, robust]
#   absolute   b(M) = sqrt(<ξ_hh(M)>/<ξ_mm>)  using linear ξ_mm from the CAMB P(k), D(z)²-scaled
# overlaid with analytic Tinker(2010, Δ=200) b(M). Window 6-18 Mpc/h (linear-biased, geometry-safe).
using Printf, Random, Statistics, LinearAlgebra, DelimitedFiles
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams, Dlinear_ab, Dlinear_tables
const rho_m=2.775e11*0.31; const hub=0.68
const r1=1600.0; const r2=2000.0; const Mlow=5e12
const Medges=[10.0^l for l in (12.7,13.1,13.5,13.9,14.4)]; const Mc=[sqrt(Medges[i]*Medges[i+1]) for i in 1:length(Medges)-1]
const redges=collect(10.0 .^ range(log10(1.0),log10(50.0);length=11)); const RMAX=redges[end]
const WIN=(6.0,18.0)
cosmo=CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo;z_max=6.0); gt=Dlinear_tables(cosmo)
Dr="/home/yguan/projects/aip-aspuru-ab/yguan/websky"
seedrng()=MersenneTwister(12345)
ortho_basis(a)=(t=abs(a[1])<0.9 ? [1.0,0.0,0.0] : [0.0,1.0,0.0]; e1=normalize(cross(a,t)); e2=cross(a,e1); (e1,e2))
incap(nx,ny,nz, ax,ay,az, cosw)= (nx*ax+ny*ay+nz*az) >= cosw
binM(M)=(i=searchsortedfirst(Medges,M)-1; (1<=i<=length(Mc)) ? i : 0)

# ---------- linear theory from CAMB P(k) (pre-divided by (2π)³) ----------
pk=readdlm("/home/yguan/work/PeakPatch.jl/validation/websky_6144/data/pk_websky.dat",comments=true,comment_char='#')
kk=Float64.(pk[:,1]); Pk=Float64.(pk[:,2]).*(2pi)^3                  # undo /(2π)³ -> standard P(k)
ok=(kk.>1e-4).&(Pk.>0); kk=kk[ok]; Pk=Pk[ok]
Wth(x)= x<1e-3 ? 1.0-x^2/10 : 3*(sin(x)-x*cos(x))/x^3
function spec_int(f)                                                # ∫ f(k) dk via trapezoid in k
    s=0.0; for i in 1:length(kk)-1; s+=0.5*(f(i)+f(i+1))*(kk[i+1]-kk[i]); end; s
end
sigma2(R)=1/(2pi^2)*spec_int(i->Pk[i]*Wth(kk[i]*R)^2*kk[i]^2)
xi_mm0(r)=1/(2pi^2)*spec_int(i->(kr=kk[i]*r; Pk[i]*kk[i]^2*(kr<1e-3 ? 1.0 : sin(kr)/kr)))
sig8=sqrt(sigma2(8.0)); @info "linear-theory check" sigma8=round(sig8;digits=4)
# Tinker 2010 b(ν), Δ=200
δc=1.686; yT=log10(200.0)
A=1.0+0.24*yT*exp(-(4/yT)^4); aT=0.44*yT-0.88; B=0.183; bT=1.5; C=0.019+0.107*yT+0.19*exp(-(4/yT)^4); cT=2.4
M2R(M)=(3*M/(4pi*rho_m))^(1/3)
tinker_b(M,Dz)=(ν=δc/(Dz*sqrt(sigma2(M2R(M)))); 1-A*ν^aT/(ν^aT+δc^aT)+B*ν^bT+C*ν^cT)

# ---------- pair counts / ξ (chaining mesh, threaded) ----------
function pairhist(qx,qy,qz, tx,ty,tz, auto::Bool)
    nb=length(redges)-1; L=RMAX; xmin=minimum(tx);ymin=minimum(ty);zmin=minimum(tz)
    grid=Dict{NTuple{3,Int},Vector{Int}}()
    @inbounds for j in eachindex(tx); c=(floor(Int,(tx[j]-xmin)/L),floor(Int,(ty[j]-ymin)/L),floor(Int,(tz[j]-zmin)/L)); push!(get!(grid,c,Int[]),j); end
    nw=max(Threads.nthreads(),1); rmax2=L^2; le=log10(redges[1]); lr=log10(redges[end])
    function worker(i0,i1)
        h=zeros(Float64,nb)
        @inbounds for i in i0:i1
            cx=floor(Int,(qx[i]-xmin)/L);cy=floor(Int,(qy[i]-ymin)/L);cz=floor(Int,(qz[i]-zmin)/L)
            for dx=-1:1,dy=-1:1,dz=-1:1
                cc=(cx+dx,cy+dy,cz+dz); haskey(grid,cc)||continue
                for j in grid[cc]
                    auto && j<=i && continue
                    d2=(qx[i]-tx[j])^2+(qy[i]-ty[j])^2+(qz[i]-tz[j])^2; (d2<rmax2 && d2>0.0)||continue
                    bin=floor(Int,(0.5*log10(d2)-le)/((lr-le)/nb))+1; (1<=bin<=nb)&&(h[bin]+=1.0)
                end
            end
        end
        h
    end
    nq=length(qx); bnds=[1+div(nq*(w-1),nw) for w in 1:nw+1]; bnds[end]=nq+1
    reduce(+, fetch.([Threads.@spawn worker(bnds[w],bnds[w+1]-1) for w in 1:nw]))
end
function xi_LS(dx,dy,dz, rx,ry,rz)
    Nd=length(dx);Nr=length(rx); DD=pairhist(dx,dy,dz,dx,dy,dz,true); RR=pairhist(rx,ry,rz,rx,ry,rz,true); DR=pairhist(dx,dy,dz,rx,ry,rz,false)
    @. (DD/(Nd*(Nd-1)/2) - 2*DR/(Float64(Nd)*Nr) + RR/(Nr*(Nr-1)/2))/max(RR/(Nr*(Nr-1)/2),1e-30)
end
function rand_cap(nr, obs, ax,ay,az,cosw, rsamp, rng)
    e1,e2=ortho_basis([ax,ay,az]); nrs=length(rsamp); x=Float64[];y=Float64[];z=Float64[]
    for _ in 1:nr
        ct=cosw+rand(rng)*(1-cosw);st=sqrt(max(1-ct^2,0.0));ph=2pi*rand(rng)
        nx=ct*ax+st*cos(ph)*e1[1]+st*sin(ph)*e2[1];ny=ct*ay+st*cos(ph)*e1[2]+st*sin(ph)*e2[2];nz=ct*az+st*cos(ph)*e1[3]+st*sin(ph)*e2[3]
        r=rsamp[rand(rng,1:nrs)]; push!(x,obs+r*nx);push!(y,obs+r*ny);push!(z,obs+r*nz)
    end
    (x,y,z)
end

# ---------- loaders (cap), carrying mass ----------
function load_ours(cat,obs, ax,ay,az,cosw)
    x=Float64[];y=Float64[];z=Float64[];M=Float64[]
    io=open(cat);N=read(io,Int32);read(io,Float32);read(io,Float32); nf=33;chunk=4_000_000;buf=Vector{Float32}(undef,chunk*nf);ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone);read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            b=(k-1)*nf;R=Float64(buf[b+7]);R>0||continue; Mh=4/3*pi*rho_m*R^3; Mh>Mlow||continue
            q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3]); rq=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2);(r1<=rq<=r2)||continue
            a=1.0/(1.0+chi_to_z(chi2z,rq))
            d11=Float64(buf[b+4])*a;d12=Float64(buf[b+5])*a;d13=Float64(buf[b+6])*a; d21=Float64(buf[b+8])*a^2;d22=Float64(buf[b+9])*a^2;d23=Float64(buf[b+10])*a^2
            ex=q1+d11-d21;ey=q2+d12-d22;ez=q3+d13-d23; dx=ex-obs;dy=ey-obs;dz=ez-obs;rr=sqrt(dx^2+dy^2+dz^2)
            incap(dx/rr,dy/rr,dz/rr,ax,ay,az,cosw)||continue
            push!(x,ex);push!(y,ey);push!(z,ez);push!(M,Mh)
        end
        ndone+=m
    end
    close(io);(x,y,z,M)
end
function load_wsky()
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc");read(io,Int32);read(io,Int32);read(io,Int32)
    seekstart(io);Nw=read(io,Int32);read(io,Int32);read(io,Int32);nf=10; buf=Vector{Float32}(undef,Int(Nw)*nf);read!(io,buf);close(io)
    x=Float64[];y=Float64[];z=Float64[];M=Float64[]
    for i in 1:Nw
        b=(i-1)*nf;R=Float64(buf[b+7])*hub; Mh=4/3*pi*rho_m*R^3; Mh>Mlow||continue
        X=Float64(buf[b+1])*hub;Y=Float64(buf[b+2])*hub;Z=Float64(buf[b+3])*hub; r=sqrt(X^2+Y^2+Z^2);(r1<=r<=r2)||continue
        push!(x,X);push!(y,Y);push!(z,Z);push!(M,Mh)
    end
    (x,y,z,M)
end

# ---- websky -> inscribed cap ----
@info "loading websky..."; wx,wy,wz,wM=load_wsky()
wax=mean(wx);way=mean(wy);waz=mean(wz);wn=sqrt(wax^2+way^2+waz^2);wax/=wn;way/=wn;waz/=wn
cosmin=minimum((wx[i]*wax+wy[i]*way+wz[i]*waz)/sqrt(wx[i]^2+wy[i]^2+wz[i]^2) for i in eachindex(wx))
cosw=cos(acos(cosmin)/sqrt(2)); keep=[(wx[i]*wax+wy[i]*way+wz[i]*waz)/sqrt(wx[i]^2+wy[i]^2+wz[i]^2)>=cosw for i in eachindex(wx)]
wx=wx[keep];wy=wy[keep];wz=wz[keep];wM=wM[keep]
ax=1/sqrt(3);ay=1/sqrt(3);az=1/sqrt(3)
@info "loading ours..."; ox,oy,oz,oM=load_ours(joinpath(Dr,"catalog_websky_6144_oct000_finecell_AM.pksc"),-2618.0,ax,ay,az,cosw)
rng=seedrng()
# linear ξ_mm at the shell redshift
zsh=chi_to_z(chi2z,0.5*(r1+r2)); Dz,_,_=Dlinear_ab(1.0/(1.0+zsh),gt)
rc=[sqrt(redges[i]*redges[i+1]) for i in 1:length(redges)-1]; wsel=findall(r->WIN[1]<=r<=WIN[2],rc)
ximm=[Dz^2*xi_mm0(r) for r in rc]; ximm_w=mean(ximm[wsel])
@info "shell" z=round(zsh;digits=3) D=round(Dz;digits=3) xi_mm_win=round(ximm_w;digits=4)

# randoms once per catalog (footprint only)
orq=[sqrt((ox[i]+2618.0)^2+(oy[i]+2618.0)^2+(oz[i]+2618.0)^2) for i in eachindex(ox)]
wrq=[sqrt(wx[i]^2+wy[i]^2+wz[i]^2) for i in eachindex(wx)]
orx,ory,orz=rand_cap(3*length(ox),-2618.0,ax,ay,az,cosw,orq,rng)
wrx,wry,wrz=rand_cap(3*length(wx),0.0,wax,way,waz,cosw,wrq,rng)

@printf("\n%-10s %-7s %-7s %-9s %-9s %-9s %-9s %-9s\n","M_bin","N_our","N_wsk","b_ours","b_wsky","b_rel","b_Tink","ours/Tk")
for ib in 1:length(Mc)
    oi=findall(m->binM(m)==ib,oM); wi=findall(m->binM(m)==ib,wM)
    (length(oi)<50 || length(wi)<50) && (@printf("%-10.2e %-7d %-7d   (too few)\n",Mc[ib],length(oi),length(wi)); continue)
    xo=xi_LS(ox[oi],oy[oi],oz[oi],orx,ory,orz); xw=xi_LS(wx[wi],wy[wi],wz[wi],wrx,wry,wrz)
    bo=sqrt(max(mean(xo[wsel]),0)/ximm_w); bw=sqrt(max(mean(xw[wsel]),0)/ximm_w)
    bt=tinker_b(Mc[ib],Dz)
    @printf("%-10.2e %-7d %-7d %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f\n",Mc[ib],length(oi),length(wi),bo,bw,bo/max(bw,1e-9),bt,bo/max(bt,1e-9))
end
@printf("\nb_rel~1 -> our bias(M) matches Websky; b_ours~b_Tink -> matches analytic Tinker. Sets the 2-halo map amplitude.\n")
