#!/usr/bin/env julia
# Localize the deficit: compare N(>M) per deg² for OUR production octant vs the Websky
# public 10x10 patch. SAME mass definition (both M=(4π/3)ρm R³, top-hat Lagrangian).
# If we match at high M and fall short only at low M → low-mass completeness (filter/floor).
# If offset at all M → normalization/config.
using Printf

rho_m = 2.775e11*0.31      # Msun/h per (Mpc/h)^3
hub = 0.68
Medges = [10.0^l for l in 12.0:0.25:15.0]
masses = [1.23e12, 1.69e12, 3e12, 5e12, 1e13, 3e13, 1e14, 3e14]

# --- Websky patch (R in Mpc, convert to Mpc/h: *h) ---
wf = "/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc"
io = open(wf); Nw = read(io, Int32); read(io, Int32); read(io, Int32)
nf = 10
bufw = Vector{Float32}(undef, Int(Nw)*nf); read!(io, bufw); close(io)
Mw = Float64[]
for i in 1:Nw
    R = Float64(bufw[(i-1)*nf+7]) * hub   # Mpc -> Mpc/h
    push!(Mw, 4/3*pi*rho_m*R^3)
end
patch_deg2 = 100.0   # 10x10 deg

# --- our production octant (stream; M=(4/3)π ρm RTHL^3) ---
cat = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_pkfix.pksc"
function stream_masses(cat, rho_m)
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
@info "streaming production..."
Mo = stream_masses(cat, rho_m)
oct_deg2 = 41253.0/8   # octant solid angle

NgtM(M, M0) = count(>(M0), M)
@printf("\n%-10s %-14s %-14s %-8s\n","M(Msun/h)","Websky/deg²","Ours/deg²","Ours/Websky")
for M0 in masses
    nw = NgtM(Mw,M0)/patch_deg2
    no = NgtM(Mo,M0)/oct_deg2
    @printf("%-10.2e %-14.2f %-14.2f %.3f\n", M0, nw, no, no/max(nw,1e-9))
end
@printf("\nWebsky patch: %d halos / %.0f deg²  ;  ours: %d / %.0f deg²\n",
        length(Mw), patch_deg2, length(Mo), oct_deg2)
