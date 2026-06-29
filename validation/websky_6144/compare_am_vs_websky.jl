#!/usr/bin/env julia
# DECISIVE TEST of the abundance-matching hypothesis.
# Bin-by-bin N(>M) per deg²: our RAW octant vs our AM'd octant vs Websky public patch.
# Same mass def everywhere: M=(4π/3)ρm R³ (top-hat Lagrangian = Websky M200m).
# Hypothesis: the low-mass deficit (raw ramp 0.41→1.0) is the missing abundance matching.
# If AM'd/Websky ≈ 1 across the range (down to completeness floor), AM was the cause.
using Printf

rho_m = 2.775e11*0.31          # Msun/h per (Mpc/h)^3
hub   = 0.68
masses = [1.23e12, 1.69e12, 3e12, 5e12, 1e13, 3e13, 1e14, 3e14]
D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"

# --- Websky public 10x10 patch (R in Mpc -> Mpc/h: *h) ---
function websky_masses(wf)
    io=open(wf); Nw=read(io,Int32); read(io,Int32); read(io,Int32)
    nf=10; buf=Vector{Float32}(undef,Int(Nw)*nf); read!(io,buf); close(io)
    M=Float64[]
    for i in 1:Nw
        R=Float64(buf[(i-1)*nf+7])*hub
        push!(M, 4/3*pi*rho_m*R^3)
    end
    M
end

# --- our octant catalog (stream, nf=33 ExtHaloRecord; RTHL already Mpc/h) ---
function stream_masses(cat)
    io=open(cat); N=read(io,Int32); read(io,Float32); read(io,Float32)
    nf=33; chunk=2_000_000; buf=Vector{Float32}(undef,chunk*nf)
    M=Float64[]; ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
        @inbounds for h in 1:m
            R=Float64(buf[(h-1)*nf+7]); R>0 && push!(M, 4/3*pi*rho_m*R^3)
        end
        ndone+=m
    end
    close(io); M
end

@info "reading Websky patch..."
Mw = websky_masses("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc")
patch_deg2 = 100.0

@info "streaming RAW octant..."
# ARGS: [1]=raw octant catalog, [2]=AM octant catalog (default = old coarse pkfix pair)
raw_cat = length(ARGS) >= 1 ? ARGS[1] : joinpath(D,"catalog_websky_6144_oct000_pkfix.pksc")
am_cat  = length(ARGS) >= 2 ? ARGS[2] : joinpath(D,"catalog_websky_6144_oct000_pkfix_AM.pksc")
Mraw = stream_masses(raw_cat)
@info "streaming AM octant..."
Mam  = stream_masses(am_cat)
oct_deg2 = 41253.0/8

NgtM(M,M0)=count(>(M0),M)
@printf("\n%-10s %-12s %-12s %-12s %-10s %-10s\n",
        "M(Msun/h)","Websky/dg²","RAW/dg²","AM/dg²","RAW/Wsky","AM/Wsky")
for M0 in masses
    nw=NgtM(Mw,M0)/patch_deg2
    nr=NgtM(Mraw,M0)/oct_deg2
    na=NgtM(Mam,M0)/oct_deg2
    @printf("%-10.2e %-12.3f %-12.3f %-12.3f %-10.3f %-10.3f\n",
            M0,nw,nr,na,nr/max(nw,1e-9),na/max(nw,1e-9))
end
@printf("\nTotals: Websky %d/%.0fdeg²  RAW %d/%.0fdeg²  AM %d/%.0fdeg²\n",
        length(Mw),patch_deg2,length(Mraw),oct_deg2,length(Mam),oct_deg2)
@printf("If AM/Wsky -> 1 across the range, the low-mass deficit was the missing abundance matching.\n")
