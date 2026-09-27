#!/usr/bin/env julia
# Mean pairwise (streaming) velocity v12(r) = <(v_i - v_j)·r̂_ij> over pairs at separation r.
# This is the kSZ driver (Websky Fig 6): pairwise infall sources the pairwise-kSZ signal. σ_vr only
# tests the 1-point velocity variance; v12(r) tests the CORRELATED infall. Bulk flow cancels in the
# difference, so this isolates relative streaming (negative = convergent infall).
# Matched geometry (same inscribed angular cap for both). Ours: Eulerian pos + km/s velocity
# reconstructed from the OLD-convention catalog (q+d1-d2; v=a·100·E·f·(d1-2·d2), d2=stored·a² carries
# the +3/7 sign so true 2LPT term = -d2). Websky: km/s velocities read directly.
using Printf, Statistics, LinearAlgebra
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams, Dlinear_ab, Dlinear_tables
const rho_m=2.775e11*0.31; const hub=0.68
const r1=1600.0; const r2=2000.0; const Mcut=1e13
const redges=collect(range(2.0,60.0;length=16)); const RMAX=redges[end]   # linear r-bins, kSZ scales
cosmo=CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo;z_max=6.0); gt=Dlinear_tables(cosmo)
Efac(a)=sqrt(cosmo.Om*a^-3+cosmo.OL)
Dr="/home/yguan/projects/aip-aspuru-ab/yguan/websky"
ortho_basis(a)=(t=abs(a[1])<0.9 ? [1.0,0.0,0.0] : [0.0,1.0,0.0]; e1=normalize(cross(a,t)); e2=cross(a,e1); (e1,e2))

# ---------- pairwise-velocity histograms (chaining mesh, threaded) ----------
function v12hist(x,y,z, vx,vy,vz)
    nb=length(redges)-1; L=RMAX
    xmin=minimum(x);ymin=minimum(y);zmin=minimum(z)
    grid=Dict{NTuple{3,Int},Vector{Int}}()
    @inbounds for j in eachindex(x)
        c=(floor(Int,(x[j]-xmin)/L),floor(Int,(y[j]-ymin)/L),floor(Int,(z[j]-zmin)/L)); push!(get!(grid,c,Int[]),j)
    end
    nw=max(Threads.nthreads(),1); rmax2=L^2; r0=redges[1]; dr=(redges[end]-redges[1])/nb
    function worker(i0,i1)
        s=zeros(Float64,nb); n=zeros(Float64,nb)
        @inbounds for i in i0:i1
            cx=floor(Int,(x[i]-xmin)/L);cy=floor(Int,(y[i]-ymin)/L);cz=floor(Int,(z[i]-zmin)/L)
            for ax=-1:1,ay=-1:1,az=-1:1
                cc=(cx+ax,cy+ay,cz+az); haskey(grid,cc)||continue
                for j in grid[cc]
                    j<=i && continue
                    dx=x[i]-x[j];dy=y[i]-y[j];dz=z[i]-z[j]; d2=dx*dx+dy*dy+dz*dz
                    (d2<rmax2 && d2>0.0)||continue; r=sqrt(d2)
                    proj=((vx[i]-vx[j])*dx+(vy[i]-vy[j])*dy+(vz[i]-vz[j])*dz)/r    # Δv·r̂
                    bin=floor(Int,(r-r0)/dr)+1; (1<=bin<=nb)||continue
                    s[bin]+=proj; n[bin]+=1.0
                end
            end
        end
        (s,n)
    end
    nq=length(x); bnds=[1+div(nq*(w-1),nw) for w in 1:nw+1]; bnds[end]=nq+1
    ts=[Threads.@spawn worker(bnds[w],bnds[w+1]-1) for w in 1:nw]
    S=zeros(nb);N=zeros(nb); for t in ts; ls,ln=fetch(t); S.+=ls; N.+=ln; end
    [N[b]>0 ? S[b]/N[b] : 0.0 for b in 1:nb]
end

# ---------- loaders (cap-filtered), with km/s velocities ----------
function load_ours(cat,obs, ax,ay,az,cosw)
    x=Float64[];y=Float64[];z=Float64[];vx=Float64[];vy=Float64[];vz=Float64[]
    io=open(cat);N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33;chunk=4_000_000;buf=Vector{Float32}(undef,chunk*nf);ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone);read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            b=(k-1)*nf;R=Float64(buf[b+7]);R>0||continue;(4/3*pi*rho_m*R^3)>Mcut||continue
            q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3])
            rq=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2);(r1<=rq<=r2)||continue
            a=1.0/(1.0+chi_to_z(chi2z,rq))
            d11=Float64(buf[b+4])*a;d12=Float64(buf[b+5])*a;d13=Float64(buf[b+6])*a    # 1LPT disp (D-fix)
            d21=Float64(buf[b+8])*a^2;d22=Float64(buf[b+9])*a^2;d23=Float64(buf[b+10])*a^2  # 2LPT (stored +3/7 sign)
            ex=q1+d11-d21;ey=q2+d12-d22;ez=q3+d13-d23                                   # Eulerian (d2_true=-d21)
            dx=ex-obs;dy=ey-obs;dz=ez-obs;rr=sqrt(dx^2+dy^2+dz^2)
            (incap(dx/rr,dy/rr,dz/rr,ax,ay,az,cosw))||continue
            aE=1.0/(1.0+chi_to_z(chi2z,rr)); _,f,_=Dlinear_ab(aE,gt); vf=aE*100.0*Efac(aE)*f
            push!(x,ex);push!(y,ey);push!(z,ez)
            push!(vx,vf*(d11-2*d21));push!(vy,vf*(d12-2*d22));push!(vz,vf*(d13-2*d23))  # km/s
        end
        ndone+=m
    end
    close(io);(x,y,z,vx,vy,vz)
end
incap(nx,ny,nz, ax,ay,az, cosw)= (nx*ax+ny*ay+nz*az) >= cosw
function load_wsky_full()
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc");read(io,Int32);read(io,Int32);read(io,Int32)
    seekstart(io);Nw=read(io,Int32);read(io,Int32);read(io,Int32);nf=10
    buf=Vector{Float32}(undef,Int(Nw)*nf);read!(io,buf);close(io)
    x=Float64[];y=Float64[];z=Float64[];vx=Float64[];vy=Float64[];vz=Float64[]
    for i in 1:Nw
        b=(i-1)*nf;R=Float64(buf[b+7])*hub;(4/3*pi*rho_m*R^3)>Mcut||continue
        X=Float64(buf[b+1])*hub;Y=Float64(buf[b+2])*hub;Z=Float64(buf[b+3])*hub
        r=sqrt(X^2+Y^2+Z^2);(r1<=r<=r2)||continue
        push!(x,X);push!(y,Y);push!(z,Z);push!(vx,Float64(buf[b+4]));push!(vy,Float64(buf[b+5]));push!(vz,Float64(buf[b+6]))
    end
    (x,y,z,vx,vy,vz)
end
sigvr(x,y,z,vx,vy,vz,obs)=(vr=[((x[i]-obs)*vx[i]+(y[i]-obs)*vy[i]+(z[i]-obs)*vz[i])/sqrt((x[i]-obs)^2+(y[i]-obs)^2+(z[i]-obs)^2) for i in eachindex(x)]; std(vr))

# ---- websky first -> matched inscribed cap ----
@info "loading websky..."; wx,wy,wz,wvx,wvy,wvz=load_wsky_full()
wax=mean(wx);way=mean(wy);waz=mean(wz);wn=sqrt(wax^2+way^2+waz^2);wax/=wn;way/=wn;waz/=wn
cosmin=minimum((wx[i]*wax+wy[i]*way+wz[i]*waz)/sqrt(wx[i]^2+wy[i]^2+wz[i]^2) for i in eachindex(wx))
cosw=cos(acos(cosmin)/sqrt(2))
keep=[(wx[i]*wax+wy[i]*way+wz[i]*waz)/sqrt(wx[i]^2+wy[i]^2+wz[i]^2)>=cosw for i in eachindex(wx)]
wx=wx[keep];wy=wy[keep];wz=wz[keep];wvx=wvx[keep];wvy=wvy[keep];wvz=wvz[keep]
ax=1/sqrt(3);ay=1/sqrt(3);az=1/sqrt(3)
@info "loading ours..."; ox,oy,oz,ovx,ovy,ovz=load_ours(joinpath(Dr,"catalog_websky_6144_oct000_finecell_AM.pksc"),-2618.0,ax,ay,az,cosw)
@info "counts" n_ours=length(ox) n_wsky=length(wx) cap_deg=round(acosd(cosw);digits=2)
@printf("sanity σ_vr [km/s]: ours=%.1f  websky=%.1f  ratio=%.3f\n",
        sigvr(ox,oy,oz,ovx,ovy,ovz,-2618.0), sigvr(wx,wy,wz,wvx,wvy,wvz,0.0),
        sigvr(ox,oy,oz,ovx,ovy,ovz,-2618.0)/sigvr(wx,wy,wz,wvx,wvy,wvz,0.0))
@info "computing v12(r)..."
vo=v12hist(ox,oy,oz,ovx,ovy,ovz); vw=v12hist(wx,wy,wz,wvx,wvy,wvz)
rc=[0.5*(redges[i]+redges[i+1]) for i in 1:length(redges)-1]
@printf("\n%-9s %-12s %-12s %-8s\n","r[Mpc/h]","v12_ours","v12_wsky","ratio")
for i in eachindex(rc); @printf("%-9.1f %-12.2f %-12.2f %-8.3f\n",rc[i],vo[i],vw[i],vo[i]/(vw[i]==0 ? NaN : vw[i])); end
sel=findall(r->5<=r<=30,rc)
@printf("\nmean ratio (5-30 Mpc/h)=%.3f. v12<0=infall; match -> pairwise-kSZ-relevant velocities reproduce Websky.\n",
        mean(vo[sel]./vw[sel]))
