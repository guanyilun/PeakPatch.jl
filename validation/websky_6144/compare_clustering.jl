#!/usr/bin/env julia
# First-pass CLUSTERING amplitude check vs Websky — the other thing our mass-function match
# doesn't cover (the 2-halo / bias term that sets the large-scale power spectra).
# Counts-in-cells in a thin z-shell: equal-solid-angle cells (bin in φ and cosθ → equal area),
# edge-trimmed. Normalized clustering variance  σ²_cl = (Var(N) - <N>) / <N>²  ∝ bias²·σ²_m(cell).
# Compare ours/Websky at fixed PHYSICAL cell size L = r̄·δ. ratio ~1 -> same clustering amplitude.
# CAVEAT: first-pass — approximate edge handling, single patch (cosmic variance), no randoms.
# A rigorous ξ(r)/w(θ) with Landy-Szalay randoms is the proper follow-up.
using Printf, Statistics
rho_m = 2.775e11*0.31; hub = 0.68
D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
r1, r2 = 1400.0, 2300.0          # Mpc/h shell (z~0.55-1.0, both complete)
rbar = (r1+r2)/2
Mcut = 1e13                       # massive halos: strong bias signal, enough numbers
Lcell = 35.0                      # physical cell size Mpc/h
δ = Lcell / rbar                  # angular cell size (rad), equal-area in (φ, cosθ)

# accumulate cell counts into Dict{(iφ,icosθ)=>count} for halos in shell+mass cut
function cic(stream_dirs)  # stream_dirs yields (φ, cosθ) for selected halos
    cells = Dict{Tuple{Int,Int},Int}()
    n=0
    stream_dirs() do φ, ct
        i = floor(Int, φ/δ); j = floor(Int, (ct+1.0)/δ)  # ct in [-1,1]
        cells[(i,j)] = get(cells,(i,j),0) + 1; n+=1
    end
    cells, n
end
# normalized clustering variance over interior cells (trim 1-cell margin, include interior zeros)
function clvar(cells)
    isempty(cells) && return (0.0,0.0,0)
    is=[k[1] for k in keys(cells)]; js=[k[2] for k in keys(cells)]
    i0,i1=minimum(is)+1,maximum(is)-1; j0,j1=minimum(js)+1,maximum(js)-1
    (i1<i0 || j1<j0) && return (0.0,0.0,0)
    vals=Int[]
    for i in i0:i1, j in j0:j1; push!(vals, get(cells,(i,j),0)); end
    m=mean(vals); v=var(vals)
    sig2 = m>0 ? (v - m)/m^2 : 0.0   # shot-noise-subtracted normalized variance
    (m, sig2, length(vals))
end

# Websky patch: pos Mpc (1-3), R Mpc (7); obs origin
function wsky_stream(f)
    return function(cb)
        io=open("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc")
        Nw=read(io,Int32); read(io,Int32); read(io,Int32)
        nf=10; buf=Vector{Float32}(undef,Int(Nw)*nf); read!(io,buf); close(io)
        for i in 1:Nw
            b=(i-1)*nf; R=Float64(buf[b+7])*hub; (4/3*pi*rho_m*R^3)>Mcut || continue
            x=Float64(buf[b+1])*hub; y=Float64(buf[b+2])*hub; z=Float64(buf[b+3])*hub
            r=sqrt(x^2+y^2+z^2); (r1<=r<=r2) || continue
            cb(atan(y,x), z/r)
        end
    end
end
# ours: pos Mpc/h rel to obs; nf=33
function ours_stream(cat, obs)
    return function(cb)
        io=open(cat); N=read(io,Int32); read(io,Float32); read(io,Float32)
        nf=33; chunk=1_000_000; buf=Vector{Float32}(undef,chunk*nf); ndone=0
        while ndone<N
            m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
            @inbounds for k in 1:m
                b=(k-1)*nf; R=Float64(buf[b+7]); R>0 || continue
                (4/3*pi*rho_m*R^3)>Mcut || continue
                x=Float64(buf[b+1])-obs; y=Float64(buf[b+2])-obs; z=Float64(buf[b+3])-obs
                r=sqrt(x^2+y^2+z^2); (r1<=r<=r2) || continue
                cb(atan(y,x), z/r)
            end
            ndone+=m
        end
        close(io)
    end
end

@info "Websky clustering (shell $r1-$r2 Mpc/h, M>$Mcut)..."
cw,nw = cic(wsky_stream(""));  mw,s2w,ncw = clvar(cw)
@info "our finecell_AM clustering..."
co,no = cic(ours_stream(joinpath(D,"catalog_websky_6144_oct000_finecell_AM.pksc"), -2618.0)); mo,s2o,nco = clvar(co)

@printf("\nz-shell r=%.0f-%.0f Mpc/h (rbar=%.0f), M>%.0e, cell L=%.0f Mpc/h (δ=%.4f rad)\n", r1,r2,rbar,Mcut,Lcell,δ)
@printf("%-10s %-9s %-10s %-12s %-8s\n","catalog","Nhalo","<N/cell>","σ²_cl","ncells")
@printf("%-10s %-9d %-10.1f %-12.4e %-8d\n","Websky", nw, mw, s2w, ncw)
@printf("%-10s %-9d %-10.1f %-12.4e %-8d\n","ours",   no, mo, s2o, nco)
@printf("\nσ²_cl(ours)/σ²_cl(Websky) = %.3f   -> relative bias² (clustering amplitude). ~1 = same clustering.\n", s2o/max(s2w,1e-30))
@printf("CAVEAT: first-pass counts-in-cells (approx edges, 1 patch=cosmic variance, no randoms). Rigorous = ξ(r)/w(θ) with LS randoms.\n")
