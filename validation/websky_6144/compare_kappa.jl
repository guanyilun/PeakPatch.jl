#!/usr/bin/env julia
# CMB lensing convergence κ from halos: the first Tier-B (painting-level) test, and the
# cheapest — pure mass projection, no gas physics. XGPaint.jl (the extragalactic painting
# reference) has no lensing module, so this is self-contained.
#
# Method: Born-approximation halo-only κ maps painted IDENTICALLY from both catalogs
# (matched footprint = square inscribed in the shared angular cap, gnomonic flat-sky,
# CMB source plane z*=1100), then C_ℓ^κκ via 2D FFT. Because painting/footprint/estimator
# are identical, C_ours/C_wsky is a pure catalog test (missing-field-mass, χ*(no-radiation),
# and pixel-window effects all cancel in the ratio). A Limber 2-halo overlay
# (empirical f_coll(χ)·b_eff(χ) from the catalog × Tinker b × linear P) gives context
# for the low-ℓ absolute level; halo shot noise is reported per catalog.
# Ours: Eulerian positions reconstructed from the OLD-convention catalog (q+d1·a-d2·a²,
# stored 2LPT carries +3/7 so true term enters with minus). Websky: positions ×h (Mpc→Mpc/h).
using Printf, Statistics, LinearAlgebra, DelimitedFiles, FFTW
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams, Dlinear_ab, Dlinear_tables, chi
const rho_m=2.775e11*0.31; const hub=0.68
const Mcut=2e12                       # both catalogs complete & MF-matched above this
const rmin=50.0; const rmax=5150.0    # common radial range (ours reaches z≈4.57)
const NPIX=512
cosmo=CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo;z_max=6.0); gt=Dlinear_tables(cosmo)
const chistar=chi(1100.0,cosmo)       # Mpc/h; no radiation term (few % high) — cancels in the ratio
const pref=1.5*cosmo.Om/2997.92^2     # (3/2)Ωm(H0/c)², lengths in Mpc/h
Dr="/home/yguan/projects/aip-aspuru-ab/yguan/websky"
ortho_basis(a)=(t=abs(a[1])<0.9 ? [1.0,0.0,0.0] : [0.0,1.0,0.0]; e1=normalize(cross(a,t)); e2=cross(a,e1); (e1,e2))

# κ weight of one halo: W_L(χ)·(M/ρ̄)/(χ²ΔΩ), W_L = pref·χ(χ*-χ)/χ*·(1+z)
kweight(M,r,z,dom)=pref*r*(chistar-r)/chistar*(1+z)*M/(rho_m*r^2*dom)

# ---------- linear theory (CAMB P(k) pre-divided by (2π)³) + Tinker b ----------
pk=readdlm("/home/yguan/work/PeakPatch.jl/validation/websky_6144/data/pk_websky.dat",comments=true,comment_char='#')
kk=Float64.(pk[:,1]); Pk=Float64.(pk[:,2]).*(2pi)^3
ok=(kk.>1e-4).&(Pk.>0); kk=kk[ok]; Pk=Pk[ok]; lkk=log.(kk); lPk=log.(Pk)
Plin(k)= (k<=kk[1]||k>=kk[end]) ? 0.0 : exp(lPk[searchsortedfirst(kk,k)-1]+(lPk[searchsortedfirst(kk,k)]-lPk[searchsortedfirst(kk,k)-1])*(log(k)-lkk[searchsortedfirst(kk,k)-1])/(lkk[searchsortedfirst(kk,k)]-lkk[searchsortedfirst(kk,k)-1]))
Wth(x)= x<1e-3 ? 1.0-x^2/10 : 3*(sin(x)-x*cos(x))/x^3
spec_int(f)=(s=0.0; for i in 1:length(kk)-1; s+=0.5*(f(i)+f(i+1))*(kk[i+1]-kk[i]); end; s)
sigma2(R)=1/(2pi^2)*spec_int(i->Pk[i]*Wth(kk[i]*R)^2*kk[i]^2)
@info "theory check" sigma8=round(sqrt(sigma2(8.0));digits=4)
δc=1.686; yT=log10(200.0)
A=1.0+0.24*yT*exp(-(4/yT)^4); aT=0.44*yT-0.88; B=0.183; bT=1.5; C=0.019+0.107*yT+0.19*exp(-(4/yT)^4); cT=2.4
M2R(M)=(3*M/(4pi*rho_m))^(1/3)
# σ²(M) lookup so per-halo Tinker b is cheap
lMtab=collect(range(log10(5e11),16.0;length=160)); s2tab=[sigma2(M2R(10.0^l)) for l in lMtab]
sig2M(M)=(l=log10(M); l<=lMtab[1] ? s2tab[1] : l>=lMtab[end] ? s2tab[end] : (i=searchsortedfirst(lMtab,l); s2tab[i-1]+(s2tab[i]-s2tab[i-1])*(l-lMtab[i-1])/(lMtab[i]-lMtab[i-1])))
tinker_b(M,Dz)=(ν=δc/(Dz*sqrt(sig2M(M))); 1-A*ν^aT/(ν^aT+δc^aT)+B*ν^bT+C*ν^cT)

# ---------- flat-sky painting into the inscribed square ----------
# footprint: square of gnomonic half-width t inscribed in the cap of radius θc (t=tan(θc)/√2)
mutable struct Painter
    map::Matrix{Float64}; t::Float64; dom::Float64
    fM::Vector{Float64}; fMb::Vector{Float64}; redges::Vector{Float64}   # per-shell ΣM, ΣM·b for f_coll·b_eff
    shot::Float64; n::Int
end
Painter(t,nsh)=Painter(zeros(NPIX,NPIX),t,(2t/NPIX)^2,zeros(nsh),zeros(nsh),collect(range(rmin,rmax;length=nsh+1)),0.0,0)
function paint!(P::Painter, M,r, gx,gy)   # gx,gy: gnomonic coords; returns true if painted
    (abs(gx)<=P.t && abs(gy)<=P.t) || return false
    z=chi_to_z(chi2z,r); w=kweight(M,r,z,P.dom)
    ix=min(NPIX,1+floor(Int,(gx+P.t)/(2P.t)*NPIX)); iy=min(NPIX,1+floor(Int,(gy+P.t)/(2P.t)*NPIX))
    P.map[ix,iy]+=w; P.shot+=(w*P.dom)^2; P.n+=1
    ish=1+floor(Int,(r-rmin)/(P.redges[2]-P.redges[1])); Dz,_,_=Dlinear_ab(1.0/(1+z),gt)
    (1<=ish<=length(P.fM)) && (P.fM[ish]+=M; P.fMb[ish]+=M*tinker_b(M,Dz))
    true
end

function load_paint_ours(cat,obs, ax,ay,az, e1,e2, t, nsh)
    P=Painter(t,nsh)
    io=open(cat);N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33;chunk=4_000_000;buf=Vector{Float32}(undef,chunk*nf);ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone);read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            b=(k-1)*nf;R=Float64(buf[b+7]);R>0||continue; Mh=4/3*pi*rho_m*R^3; Mh>Mcut||continue
            q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3])
            rq=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2); (rmin-100<=rq<=rmax+100)||continue
            a=1.0/(1.0+chi_to_z(chi2z,rq))
            ex=q1+Float64(buf[b+4])*a-Float64(buf[b+8])*a^2; ey=q2+Float64(buf[b+5])*a-Float64(buf[b+9])*a^2; ez=q3+Float64(buf[b+6])*a-Float64(buf[b+10])*a^2
            dx=ex-obs;dy=ey-obs;dz=ez-obs; r=sqrt(dx^2+dy^2+dz^2); (rmin<=r<=rmax)||continue
            na=(dx*ax+dy*ay+dz*az); na>0||continue
            paint!(P,Mh,r,(dx*e1[1]+dy*e1[2]+dz*e1[3])/na,(dx*e2[1]+dy*e2[2]+dz*e2[3])/na)
        end
        ndone+=m
    end
    close(io); P
end
function load_paint_wsky(ax,ay,az, e1,e2, t, nsh)
    P=Painter(t,nsh)
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc")
    Nw=read(io,Int32);read(io,Int32);read(io,Int32);nf=10
    buf=Vector{Float32}(undef,Int(Nw)*nf);read!(io,buf);close(io)
    for i in 1:Nw
        b=(i-1)*nf;R=Float64(buf[b+7])*hub; Mh=4/3*pi*rho_m*R^3; Mh>Mcut||continue
        X=Float64(buf[b+1])*hub;Y=Float64(buf[b+2])*hub;Z=Float64(buf[b+3])*hub
        r=sqrt(X^2+Y^2+Z^2); (rmin<=r<=rmax)||continue
        na=X*ax+Y*ay+Z*az; na>0||continue
        paint!(P,Mh,r,(X*e1[1]+Y*e1[2]+Z*e1[3])/na,(X*e2[1]+Y*e2[2]+Z*e2[3])/na)
    end
    P
end

# ---------- C_ℓ of a map (flat sky) ----------
function cl_map(P::Painter, ledges)
    m=P.map.-mean(P.map)                             # weights already per-pixel: map holds κ directly
    L=2*P.t                                          # radians (flat-sky)
    F=fft(m); nb=length(ledges)-1; S=zeros(nb); Nm=zeros(nb)
    freq=fftfreq(NPIX,NPIX/L).*(2pi)                 # ℓ values per axis
    for j in 1:NPIX, i in 1:NPIX
        l=sqrt(freq[i]^2+freq[j]^2); (ledges[1]<=l<ledges[end])||continue
        b=searchsortedlast(ledges,l); S[b]+=abs2(F[i,j]); Nm[b]+=1
    end
    [Nm[b]>0 ? S[b]/Nm[b]*L^2/NPIX^4 : 0.0 for b in 1:nb]
end

# ---------- Limber 2-halo theory: C_ℓ = ∫dχ [W_L(χ)f_coll(χ)b_eff(χ)]² D²P(ℓ/χ)/χ² ----------
function cl_2halo(P::Painter, dom_foot, lcens)
    re=P.redges; nsh=length(P.fM); cl=zeros(length(lcens))
    for s in 1:nsh
        P.fM[s]>0 || continue
        r=0.5*(re[s]+re[s+1]); dchi=re[s+1]-re[s]; V=dom_foot*(re[s+1]^3-re[s]^3)/3
        f=P.fM[s]/(rho_m*V); beff=P.fMb[s]/P.fM[s]
        z=chi_to_z(chi2z,r); Dz,_,_=Dlinear_ab(1.0/(1+z),gt)
        W=pref*r*(chistar-r)/chistar*(1+z)*f*beff
        for (i,l) in enumerate(lcens); cl[i]+=W^2*Dz^2*Plin((l+0.5)/r)/r^2*dchi; end
    end
    cl
end

# ---- footprint from the Websky patch (same inscribed cap as previous scripts) ----
@info "scanning websky for footprint..."
io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc")
Nw=read(io,Int32);read(io,Int32);read(io,Int32);nf=10
buf=Vector{Float32}(undef,Int(Nw)*nf);read!(io,buf);close(io)
sx=0.0;sy=0.0;sz=0.0
for i in 1:Nw; b=(i-1)*nf; global sx+=buf[b+1]; global sy+=buf[b+2]; global sz+=buf[b+3]; end
wn=sqrt(sx^2+sy^2+sz^2); wax=sx/wn;way=sy/wn;waz=sz/wn
cosmin=minimum((Float64(buf[(i-1)*nf+1])*wax+Float64(buf[(i-1)*nf+2])*way+Float64(buf[(i-1)*nf+3])*waz)/sqrt(Float64(buf[(i-1)*nf+1])^2+Float64(buf[(i-1)*nf+2])^2+Float64(buf[(i-1)*nf+3])^2) for i in 1:Nw)
thc=acos(cosmin)/sqrt(2)                       # inscribed-cap radius
t=tan(thc)/sqrt(2)                             # half-width of square inscribed in that cap
@info "footprint" cap_deg=round(rad2deg(thc);digits=2) square_deg=round(2*rad2deg(atan(t));digits=2)
we1,we2=ortho_basis([wax,way,waz])
oax=1/sqrt(3);oay=1/sqrt(3);oaz=1/sqrt(3); oe1,oe2=ortho_basis([oax,oay,oaz])
NSH=25
@info "painting websky..."; Pw=load_paint_wsky(wax,way,waz,we1,we2,t,NSH)
@info "painting ours...";  Po=load_paint_ours(joinpath(Dr,"catalog_websky_6144_oct000_finecell_AM.pksc"),-2618.0,oax,oay,oaz,oe1,oe2,t,NSH)
@info "counts" n_ours=Po.n n_wsky=Pw.n

ledges=[10.0^l for l in range(log10(150.0),log10(10000.0);length=15)]
lc=[sqrt(ledges[i]*ledges[i+1]) for i in 1:length(ledges)-1]
co=cl_map(Po,ledges); cw=cl_map(Pw,ledges)
A=(2t)^2; sho=Po.shot/A; shw=Pw.shot/A
c2h=cl_2halo(Po,A,lc)
@printf("\nshot noise: ours=%.3e  websky=%.3e (flat C_shot; set by the few NEAREST massive halos,\n",sho,shw)
@printf("a realization-dependent quantity — our octant and Websky's patch are DIFFERENT sky)\n")
@printf("%-8s %-11s %-11s %-8s %-11s %-9s %-9s\n","ell","Cl_ours","Cl_wsky","raw_rat","Cl_2halo","ours/2h","wsky/2h")
for i in eachindex(lc)
    (co[i]>0 && cw[i]>0) || continue
    @printf("%-8.0f %-11.3e %-11.3e %-8.3f %-11.3e %-9.3f %-9.3f\n",
            lc[i],co[i],cw[i],co[i]/cw[i],c2h[i],(co[i]-sho)/c2h[i],(cw[i]-shw)/c2h[i])
end
sel=findall(l->300<=l<=1400,lc)   # 2-halo-dominated window; above ~1400 both maps are pure shot
@printf("\nshot-subtracted vs Limber 2-halo, mean over 300<=ell<=1400: ours/theory=%.3f  websky/theory=%.3f\n",
        mean((co[sel].-sho)./c2h[sel]), mean((cw[sel].-shw)./c2h[sel]))
@printf("PASS criterion: each catalog's shot-subtracted C_l tracks the 2-halo Limber level (~1 within\n")
@printf("single-cap sample variance). Raw ours/wsky != 1 at high ell is EXPECTED (shot of different sky).\n")
