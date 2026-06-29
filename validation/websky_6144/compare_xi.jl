#!/usr/bin/env julia
# Rigorous 3D two-point correlation ξ(r) via Landy-Szalay, MATCHED GEOMETRY.
# Earlier full-octant-vs-10°patch comparison agreed at r<10 Mpc/h (local pairs) but diverged at large r
# — an apples-to-oranges FOOTPRINT mismatch (octant has large-scale modes the tiny patch can't), NOT a
# catalog difference and NOT a randoms artifact (shuffled-r randoms didn't change it). Fix: carve a
# Websky-SIZED angular cap out of our octant and measure ξ(r) in identical geometry for both.
# Randoms: uniform-in-cap direction (same cap radius for both) × radius drawn from each catalog's own
# radial coords (shuffled-r). ξ_LS = (DD - 2DR + RR)/RR.
using Printf, Random, Statistics, LinearAlgebra
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams
const rho_m=2.775e11*0.31; const hub=0.68
const r1=1600.0; const r2=2000.0; const Mcut=1e13            # thin shell, M>1e13 (Mpc/h)
const BOX=5236.0                                             # our box (Mpc/h)
const redges=collect(10.0 .^ range(log10(1.0), log10(50.0); length=11)); const RMAX=redges[end]
cosmo=CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo;z_max=6.0)
Dr="/home/yguan/projects/aip-aspuru-ab/yguan/websky"
seedrng()=MersenneTwister(12345)   # deterministic (Math.random/argless not allowed)

# ---------- angular cap helpers ----------
ortho_basis(a)=(t=abs(a[1])<0.9 ? [1.0,0.0,0.0] : [0.0,1.0,0.0]; e1=normalize(cross(a,t)); e2=cross(a,e1); (e1,e2))
incap(nx,ny,nz, ax,ay,az, cosw)= (nx*ax+ny*ay+nz*az) >= cosw

# ---------- chaining-mesh pair counts ----------
function pairhist(qx,qy,qz, tx,ty,tz, auto::Bool)
    nb=length(redges)-1; L=RMAX
    xmin=minimum(tx);ymin=minimum(ty);zmin=minimum(tz)
    grid=Dict{NTuple{3,Int},Vector{Int}}()
    @inbounds for j in eachindex(tx)
        c=(floor(Int,(tx[j]-xmin)/L),floor(Int,(ty[j]-ymin)/L),floor(Int,(tz[j]-zmin)/L))
        push!(get!(grid,c,Int[]),j)
    end
    nw=max(Threads.nthreads(),1); rmax2=L^2; le=log10(redges[1]); lr=log10(redges[end])
    function worker(i0,i1)
        h=zeros(Float64,nb)
        @inbounds for i in i0:i1
            cx=floor(Int,(qx[i]-xmin)/L);cy=floor(Int,(qy[i]-ymin)/L);cz=floor(Int,(qz[i]-zmin)/L)
            for dx=-1:1,dy=-1:1,dz=-1:1
                cc=(cx+dx,cy+dy,cz+dz); haskey(grid,cc)||continue
                for j in grid[cc]
                    auto && j<=i && continue
                    d2=(qx[i]-tx[j])^2+(qy[i]-ty[j])^2+(qz[i]-tz[j])^2
                    (d2<rmax2 && d2>0.0)||continue
                    bin=floor(Int,(0.5*log10(d2)-le)/((lr-le)/nb))+1
                    (1<=bin<=nb)&&(h[bin]+=1.0)
                end
            end
        end
        h
    end
    nq=length(qx); bnds=[1+div(nq*(w-1),nw) for w in 1:nw+1]; bnds[end]=nq+1
    ts=[Threads.@spawn worker(bnds[w],bnds[w+1]-1) for w in 1:nw]
    reduce(+, fetch.(ts))
end
function xi_LS(dx,dy,dz, rx,ry,rz)
    Nd=length(dx); Nr=length(rx)
    DD=pairhist(dx,dy,dz, dx,dy,dz, true)
    RR=pairhist(rx,ry,rz, rx,ry,rz, true)
    DR=pairhist(dx,dy,dz, rx,ry,rz, false)
    dn=Nd*(Nd-1)/2; rn=Nr*(Nr-1)/2; xn=Float64(Nd)*Nr
    ddn=DD./dn; rrn=RR./rn; drn=DR./xn
    @. (ddn - 2*drn + rrn)/max(rrn,1e-30)
end

# ---------- load data (cap-filtered) ----------
# ours: stream 26GB, keep shell ∩ M>cut ∩ within cap(axis_o,cosw). Eulerian = q + d1 - d2 (old-conv sign).
function load_ours(cat,obs, ax,ay,az,cosw)
    x=Float64[];y=Float64[];z=Float64[]
    io=open(cat);N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33;chunk=4_000_000;buf=Vector{Float32}(undef,chunk*nf);ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone);read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            b=(k-1)*nf;R=Float64(buf[b+7]);R>0||continue;(4/3*pi*rho_m*R^3)>Mcut||continue
            q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3])
            rq=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2);(r1<=rq<=r2)||continue
            a=1.0/(1.0+chi_to_z(chi2z,rq))
            d11=Float64(buf[b+4])*a;d12=Float64(buf[b+5])*a;d13=Float64(buf[b+6])*a
            d21=Float64(buf[b+8])*a^2;d22=Float64(buf[b+9])*a^2;d23=Float64(buf[b+10])*a^2
            ex=q1+d11-d21;ey=q2+d12-d22;ez=q3+d13-d23                      # Eulerian
            dx=ex-obs;dy=ey-obs;dz=ez-obs;rr=sqrt(dx^2+dy^2+dz^2)
            incap(dx/rr,dy/rr,dz/rr, ax,ay,az, cosw)||continue            # angular cap
            push!(x,ex);push!(y,ey);push!(z,ez)
        end
        ndone+=m
    end
    close(io);(x,y,z)
end
function load_wsky()
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc");read(io,Int32);read(io,Int32);read(io,Int32)
    seekstart(io);Nw=read(io,Int32);read(io,Int32);read(io,Int32);nf=10
    buf=Vector{Float32}(undef,Int(Nw)*nf);read!(io,buf);close(io);x=Float64[];y=Float64[];z=Float64[]
    for i in 1:Nw
        b=(i-1)*nf;R=Float64(buf[b+7])*hub;(4/3*pi*rho_m*R^3)>Mcut||continue
        X=Float64(buf[b+1])*hub;Y=Float64(buf[b+2])*hub;Z=Float64(buf[b+3])*hub
        r=sqrt(X^2+Y^2+Z^2);(r1<=r<=r2)||continue;push!(x,X);push!(y,Y);push!(z,Z)
    end
    (x,y,z)
end

# ---------- cap-uniform randoms (radius from data's own radial coords) ----------
function rand_cap(nr, obs, ax,ay,az,cosw, rsamp, rng)
    e1,e2=ortho_basis([ax,ay,az]); nrs=length(rsamp)
    x=Float64[];y=Float64[];z=Float64[]
    for _ in 1:nr
        ct=cosw+rand(rng)*(1-cosw); st=sqrt(max(1-ct^2,0.0)); ph=2pi*rand(rng)  # uniform in cap
        nx=ct*ax+st*cos(ph)*e1[1]+st*sin(ph)*e2[1]
        ny=ct*ay+st*cos(ph)*e1[2]+st*sin(ph)*e2[2]
        nz=ct*az+st*cos(ph)*e1[3]+st*sin(ph)*e2[3]
        r=rsamp[rand(rng,1:nrs)]
        push!(x,obs+r*nx);push!(y,obs+r*ny);push!(z,obs+r*nz)
    end
    (x,y,z)
end

# ---- Websky first: defines the matched cap (center axis + INSCRIBED circular radius) ----
@info "loading websky..."; wx0,wy0,wz0=load_wsky()
wax=mean(wx0);way=mean(wy0);waz=mean(wz0);wn=sqrt(wax^2+way^2+waz^2);wax/=wn;way/=wn;waz/=wn   # patch center axis
cosmin=minimum((wx0[i]*wax+wy0[i]*way+wz0[i]*waz)/sqrt(wx0[i]^2+wy0[i]^2+wz0[i]^2) for i in eachindex(wx0))  # to corner
cosw=cos(acos(cosmin)/sqrt(2))        # INSCRIBED circle of the square patch (corner θ → θ/√2), fits all-data region
keep=[(wx0[i]*wax+wy0[i]*way+wz0[i]*waz)/sqrt(wx0[i]^2+wy0[i]^2+wz0[i]^2) >= cosw for i in eachindex(wx0)]
wx=wx0[keep];wy=wy0[keep];wz=wz0[keep]                                   # websky data = circular cap (matches randoms)
@info "matched cap" corner_deg=round(acosd(cosmin);digits=2) cap_radius_deg=round(acosd(cosw);digits=2) n_wsky=length(wx)

# ---- ours, carved into a SAME-SIZE cap around the octant center direction (1,1,1) ----
ax=1/sqrt(3);ay=1/sqrt(3);az=1/sqrt(3)
@info "loading ours (Eulerian, capped)..."; ox,oy,oz=load_ours(joinpath(Dr,"catalog_websky_6144_oct000_finecell_AM.pksc"),-2618.0, ax,ay,az,cosw)
rng=seedrng()
if length(ox)>length(wx)*3 && length(ox)>60_000   # keep counts comparable, cap pair cost
    nk=max(length(wx)*3, 60_000); idx=randperm(rng,length(ox))[1:nk]; ox=ox[idx];oy=oy[idx];oz=oz[idx]
end
@info "counts (matched cap)" n_ours=length(ox) n_wsky=length(wx)

orq=[sqrt((ox[i]+2618.0)^2+(oy[i]+2618.0)^2+(oz[i]+2618.0)^2) for i in eachindex(ox)]
wrq=[sqrt(wx[i]^2+wy[i]^2+wz[i]^2) for i in eachindex(wx)]
orx,ory,orz=rand_cap(3*length(ox),-2618.0, ax,ay,az,cosw, orq, rng)
wrx,wry,wrz=rand_cap(3*length(wx), 0.0,    wax,way,waz,cosw, wrq, rng)
@info "computing ξ(r)..."
xo=xi_LS(ox,oy,oz, orx,ory,orz)
xw=xi_LS(wx,wy,wz, wrx,wry,wrz)
rc=[sqrt(redges[i]*redges[i+1]) for i in 1:length(redges)-1]
@printf("\n%-9s %-10s %-10s %-8s\n","r[Mpc/h]","ξ_ours","ξ_wsky","ratio")
for i in eachindex(rc); @printf("%-9.2f %-10.4f %-10.4f %-8.3f\n",rc[i],xo[i],xw[i],xo[i]/(xw[i]==0 ? NaN : xw[i])); end
sel=findall(r->3<=r<=15,rc); bo=sqrt(max(mean(xo[sel]),0)); bw=sqrt(max(mean(xw[sel]),0))
@printf("\nmatched-geometry bias (3-15 Mpc/h): <ξ_ours>=%.4f <ξ_wsky>=%.4f -> b_ours/b_wsky=%.3f\n",mean(xo[sel]),mean(xw[sel]),bo/max(bw,1e-9))
@printf("Same angular cap for both -> apples-to-apples. Match here = clustering reproduced.\n")
