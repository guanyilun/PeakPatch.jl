#!/usr/bin/env julia
# Does the Eulerian shift (finalize_eulerian) fix the clustering deficit (was 0.27x Lagrangian)?
# Counts-in-cells variance σ²_cl = (Var(N)-<N>)/<N>² in a THIN shell (fixes earlier thick-shell
# projection), 3D cubic cells. Compare OUR Lagrangian (q) vs OUR Eulerian (q+disp1+disp2) vs Websky.
#   Lagrangian << Eulerian ~ Websky  ->  the missing Eulerian shift was the clustering deficit.
# Existing catalog is OLD convention: disp1=stored·a, disp2(physical)=-(stored·a²)  (2LPT sign-flip).
using Printf, Statistics
import PeakPatch.Cosmology: build_chi_to_z, chi_to_z, CosmologyParams
const rho_m=2.775e11*0.31; const hub=0.68
const r1=1600.0; const r2=2000.0; const Mcut=1e13; const L=30.0   # thin shell, M>1e13, 30 Mpc/h cells
cosmo=CosmologyParams(0.31,0.049,0.69,0.68,0.965,0.81); chi2z=build_chi_to_z(cosmo;z_max=6.0)
D="/home/yguan/projects/aip-aspuru-ab/yguan/websky"

cellkey(x,y,z)=(floor(Int,x/L),floor(Int,y/L),floor(Int,z/L))
mergeD!(a,b)=(for (k,v) in b; a[k]=get(a,k,0)+v; end; a)
# σ²_cl over interior cells (trim 1-cell margin of the occupied bounding box; include interior zeros)
function clvar(cells)
    isempty(cells) && return (0.0,0.0,0)
    is=[k[1] for k in keys(cells)];js=[k[2] for k in keys(cells)];ks=[k[3] for k in keys(cells)]
    i0,i1=minimum(is)+1,maximum(is)-1;j0,j1=minimum(js)+1,maximum(js)-1;k0,k1=minimum(ks)+1,maximum(ks)-1
    (i1<i0||j1<j0||k1<k0) && return (0.0,0.0,0)
    vals=Int[]
    for i in i0:i1,j in j0:j1,k in k0:k1; c=get(cells,(i,j,k),0); push!(vals,c); end
    m=mean(vals);v=var(vals); (m, m>0 ? (v-m)/m^2 : 0.0, length(vals))
end

# our catalog: q Lagrangian; Eulerian = q + d1 - d2 (old-convention sign). nf=33.
function ours(cat,obs)
    nw=max(Threads.nthreads(),1)
    function rng(buf,k0,k1,nf)
        cl=Dict{NTuple{3,Int},Int}(); ce=Dict{NTuple{3,Int},Int}()
        for k in k0:k1
            b=(k-1)*nf;R=Float64(buf[b+7]);R>0||continue;(4/3*pi*rho_m*R^3)>Mcut||continue
            q1=Float64(buf[b+1]);q2=Float64(buf[b+2]);q3=Float64(buf[b+3])
            rq=sqrt((q1-obs)^2+(q2-obs)^2+(q3-obs)^2); (r1<=rq<=r2)||continue
            z=chi_to_z(chi2z,rq);a=1.0/(1.0+z)
            d11=Float64(buf[b+4])*a;d12=Float64(buf[b+5])*a;d13=Float64(buf[b+6])*a
            d21=Float64(buf[b+8])*a^2;d22=Float64(buf[b+9])*a^2;d23=Float64(buf[b+10])*a^2
            ex=q1+d11-d21;ey=q2+d12-d22;ez=q3+d13-d23   # Eulerian (old-convention 2LPT sign-flip -> -d2)
            kl=cellkey(q1,q2,q3); cl[kl]=get(cl,kl,0)+1
            ke=cellkey(ex,ey,ez); ce[ke]=get(ce,ke,0)+1
        end
        (cl,ce)
    end
    io=open(cat);N=read(io,Int32);read(io,Float32);read(io,Float32)
    nf=33;chunk=4_000_000;buf=Vector{Float32}(undef,chunk*nf);ndone=0
    CL=Dict{NTuple{3,Int},Int}();CE=Dict{NTuple{3,Int},Int}()
    @info "streaming ours ($nw workers)..."
    while ndone<N
        m=min(chunk,Int(N)-ndone);read!(io,view(buf,1:m*nf))
        bnds=[1+div(m*(w-1),nw) for w in 1:nw+1];bnds[end]=m+1
        ts=[Threads.@spawn rng(buf,bnds[w],bnds[w+1]-1,nf) for w in 1:nw]
        for t in ts;cl,ce=fetch(t);mergeD!(CL,cl);mergeD!(CE,ce);end
        ndone+=m
    end
    close(io);(CL,CE)
end
# Websky: pos Mpc ->×h Mpc/h, obs origin; already Eulerian km/s catalog
function wsky()
    io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc");read(io,Int32);read(io,Int32);read(io,Int32)
    seekstart(io);Nw=read(io,Int32);read(io,Int32);read(io,Int32);nf=10
    buf=Vector{Float32}(undef,Int(Nw)*nf);read!(io,buf);close(io);C=Dict{NTuple{3,Int},Int}()
    for i in 1:Nw
        b=(i-1)*nf;R=Float64(buf[b+7])*hub;(4/3*pi*rho_m*R^3)>Mcut||continue
        x=Float64(buf[b+1])*hub;y=Float64(buf[b+2])*hub;z=Float64(buf[b+3])*hub
        r=sqrt(x^2+y^2+z^2);(r1<=r<=r2)||continue;k=cellkey(x,y,z);C[k]=get(C,k,0)+1
    end
    C
end

@info "Websky..."; W=wsky(); mw,s2w,ncw=clvar(W)
CL,CE=ours(joinpath(D,"catalog_websky_6144_oct000_finecell_AM.pksc"),-2618.0)
ml,s2l,ncl=clvar(CL); me,s2e,nce=clvar(CE)
@printf("\nthin shell r=%.0f-%.0f Mpc/h, M>%.0e, cell=%.0f Mpc/h\n",r1,r2,Mcut,L)
@printf("%-18s %-10s %-12s %-7s\n","catalog","<N/cell>","σ²_cl","ncells")
@printf("%-18s %-10.2f %-12.4e %-7d\n","Websky",mw,s2w,ncw)
@printf("%-18s %-10.2f %-12.4e %-7d\n","ours LAGRANGIAN",ml,s2l,ncl)
@printf("%-18s %-10.2f %-12.4e %-7d\n","ours EULERIAN",me,s2e,nce)
@printf("\nσ²_cl ratios vs Websky:  Lagrangian=%.3f   Eulerian=%.3f\n",s2l/max(s2w,1e-30),s2e/max(s2w,1e-30))
@printf("Eulerian/Lagrangian = %.2f (clustering boost from the shift). Eulerian~1 vs Websky -> shift fixes clustering.\n",s2e/max(s2l,1e-30))
