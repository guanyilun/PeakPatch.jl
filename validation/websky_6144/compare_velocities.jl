#!/usr/bin/env julia
# Velocity validation vs Websky — the kSZ-relevant quantity our mass-function match does NOT cover.
# Compares the peculiar-velocity statistics (radial v_r and |v|) as a function of mass:
#   sigma_v(M) = std of v_radial,  rms|v|(M).
# If ours/Websky ~ 1 across mass bins -> our 2LPT velocity normalization matches Websky.
# A constant offset would flag a velocity-normalization mismatch (would bias kSZ).
using Printf
rho_m = 2.775e11*0.31; hub = 0.68
D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
Medges = [10.0^l for l in 12.0:0.5:14.5]
Mc = [sqrt(Medges[i]*Medges[i+1]) for i in 1:length(Medges)-1]
binM(M) = (i=searchsortedfirst(Medges,M)-1; (1<=i<=length(Mc)) ? i : 0)

# accumulators per mass bin: n, sum vr, sum vr^2, sum |v|^2
mutable struct Acc; n::Int; svr::Float64; svr2::Float64; sv2::Float64; end
newaccs() = [Acc(0,0.0,0.0,0.0) for _ in 1:length(Mc)]
function add!(a, vr, v2)
    a.n += 1; a.svr += vr; a.svr2 += vr^2; a.sv2 += v2
end

# Websky patch: pos Mpc (1-3), vel km/s (4-6), R Mpc (7); obs at origin
function wsky(wf)
    io=open(wf); Nw=read(io,Int32); read(io,Int32); read(io,Int32)
    nf=10; buf=Vector{Float32}(undef,Int(Nw)*nf); read!(io,buf); close(io); A=newaccs()
    for i in 1:Nw
        b=(i-1)*nf
        R=Float64(buf[b+7])*hub; M=4/3*pi*rho_m*R^3; ib=binM(M); ib==0 && continue
        x=Float64(buf[b+1]); y=Float64(buf[b+2]); z=Float64(buf[b+3])
        vx=Float64(buf[b+4]); vy=Float64(buf[b+5]); vz=Float64(buf[b+6])
        r=sqrt(x^2+y^2+z^2); vr = r>0 ? (x*vx+y*vy+z*vz)/r : 0.0
        add!(A[ib], vr, vx^2+vy^2+vz^2)
    end
    A
end

# our octant (stream): pos Mpc/h (1-3) rel to obs, vel (4-6), RTHL (7); nf=33
function ours(cat, obs)
    io=open(cat); N=read(io,Int32); read(io,Float32); read(io,Float32)
    nf=33; chunk=1_000_000; buf=Vector{Float32}(undef,chunk*nf); A=newaccs(); ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            b=(k-1)*nf
            R=Float64(buf[b+7]); R>0 || continue
            M=4/3*pi*rho_m*R^3; ib=binM(M); ib==0 && continue
            dx=Float64(buf[b+1])-obs; dy=Float64(buf[b+2])-obs; dz=Float64(buf[b+3])-obs
            vx=Float64(buf[b+4]); vy=Float64(buf[b+5]); vz=Float64(buf[b+6])
            r=sqrt(dx^2+dy^2+dz^2); vr = r>0 ? (dx*vx+dy*vy+dz*vz)/r : 0.0
            add!(A[ib], vr, vx^2+vy^2+vz^2)
        end
        ndone+=m
    end
    close(io); A
end

sig(a) = a.n>1 ? sqrt(max(a.svr2/a.n - (a.svr/a.n)^2, 0.0)) : 0.0
rmsv(a) = a.n>0 ? sqrt(a.sv2/a.n) : 0.0

@info "Websky velocities..."; W = wsky("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc")
@info "our finecell_AM velocities..."; O = ours(joinpath(D,"catalog_websky_6144_oct000_finecell_AM.pksc"), -2618.0)

@printf("\n%-10s %-9s | %-10s %-10s %-7s | %-10s %-10s %-7s\n",
        "M(Msun/h)","Nwsky","σvr_wsky","σvr_ours","ratio","rms|v|_w","rms|v|_o","ratio")
for i in 1:length(Mc)
    sw=sig(W[i]); so=sig(O[i]); rw=rmsv(W[i]); ro=rmsv(O[i])
    @printf("%-10.2e %-9d | %-10.1f %-10.1f %-7.3f | %-10.1f %-10.1f %-7.3f\n",
            Mc[i], W[i].n, sw, so, so/max(sw,1e-9), rw, ro, ro/max(rw,1e-9))
end
@printf("\nUnits km/s. ratio~1 across mass -> velocity normalization matches Websky (kSZ-relevant). Constant offset -> norm mismatch.\n")
