#!/usr/bin/env julia
# Discriminate WHY we're short below 3e12. Fine differential mass function (raw, pre-AM)
# ours vs Websky around the smallest-filter mass M(Rf_min=2.068)=3.18e12.
#   - sharp cutoff at ~3.18e12 with little below -> pure FILTER floor (need finer Rf_min)
#   - smooth tail below 3.18e12 but ~0.4x Websky -> shrinkage/completeness, not a hard cutoff
using Printf
rho_m = 2.775e11*0.31; hub = 0.68
D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
# ARGS: [1]=catalog path (default raw pkfix), [2]=cellsize Mpc/h (for Rf_min mark, default 1.25326)
catpath = get(ARGS, 1, joinpath(D,"catalog_websky_6144_oct000_pkfix.pksc"))
cellsz  = parse(Float64, get(ARGS, 2, "1.25326"))
Mfilt = 4/3*pi*rho_m*(1.65*cellsz)^3   # smallest-filter top-hat mass

edges = 10.0 .^ (11.8:0.1:13.4)          # fine log bins 6e11 .. 2.5e13
ctr = [sqrt(edges[i]*edges[i+1]) for i in 1:length(edges)-1]
binM!(h,M) = (i=searchsortedfirst(edges,M)-1; (1<=i<=length(ctr)) && (h[i]+=1.0))

function stream_raw(cat)
    io=open(cat); N=read(io,Int32); read(io,Float32); read(io,Float32)
    nf=33; chunk=2_000_000; buf=Vector{Float32}(undef,chunk*nf); h=zeros(length(ctr)); ndone=0
    while ndone<N
        m=min(chunk,Int(N)-ndone); read!(io,view(buf,1:m*nf))
        @inbounds for k in 1:m
            R=Float64(buf[(k-1)*nf+7]); R>0 && binM!(h, 4/3*pi*rho_m*R^3)
        end
        ndone+=m
    end
    close(io); h
end
function wsky_mf(wf)
    io=open(wf); Nw=read(io,Int32); read(io,Int32); read(io,Int32)
    nf=10; buf=Vector{Float32}(undef,Int(Nw)*nf); read!(io,buf); close(io); h=zeros(length(ctr))
    for i in 1:Nw
        R=Float64(buf[(i-1)*nf+7])*hub; binM!(h, 4/3*pi*rho_m*R^3)
    end
    h
end

@info "raw octant MF..." catpath; ho = stream_raw(catpath)
@info "websky patch MF..."; hw = wsky_mf("/home/yguan/scratch/websky_6144/websky_ref/halos_10x10.pksc")
oct=41253.0/8; patch=100.0
@printf("\nsmallest-filter mass M(Rf_min=2.068) = %.3e Msun/h  (marked * below)\n", Mfilt)
@printf("%-10s %-12s %-12s %-9s\n","M_center","ours/deg²","wsky/deg²","ratio")
for i in 1:length(ctr)
    no=ho[i]/oct; nw=hw[i]/patch
    mark = (edges[i]<=Mfilt<edges[i+1]) ? " *<-Rf_min" : ""
    @printf("%-10.2e %-12.3f %-12.3f %-9.3f%s\n", ctr[i], no, nw, no/max(nw,1e-9), mark)
end
@printf("\nIf ours drops sharply at *<-Rf_min with little below -> pure filter floor.\n")
@printf("If ours has a smooth tail below * but ~0.4x wsky -> shrinkage/completeness, not hard cutoff.\n")
