#!/usr/bin/env julia
# Generate Julia's exact field (256^3, seed 13579, the A/B field) and write it in
# Fortran's Fvec format (n^3 Float32, x-fastest, no header) for ireadfield=1.
# Also count Julia's own local maxima > fcrit on THIS field (nomask + masked,
# summed over the 21 filter scales) so we compare apples-to-apples vs Fortran's
# "peaks found" on the identical field.
using PeakPatch, FFTW, Printf
import PeakPatch.PowerSpectrum: load_pk
import PeakPatch.Filters: smooth_field

n=256; box=320.8333; seed=13579; fcrit=1.686f0
pk=load_pk(joinpath(@__DIR__,"data","pk_websky.dat"))
filters=PeakPatch.read_filterbank(joinpath(@__DIR__,"data","filters_websky.dat"))
sort!(filters; by=f->-f[3])

@info "generating field..."
delta=generate_grf(n,pk,box,seed; fortran_compat=true)   # Float32, column-major (x fastest)
outdir="/home/yguan/scratch/websky_6144/fortran_floortest/fields"
mkpath(outdir)
fn=joinpath(outdir,"Fvec_jfield")
open(fn,"w") do io
    write(io, delta)   # column-major: x fastest, matches Fortran (((d(i,j,k),i),j),k)
end
@printf("wrote %s : %d bytes (expect %d)\n", fn, filesize(fn), n^3*4)
@printf("field mean=%.5f sigma=%.5f\n", sum(delta)/n^3, sqrt(sum(x->x^2,delta)/n^3 - (sum(delta)/n^3)^2))

# Count Julia local maxima > fcrit on this field (global, periodic), nomask+masked
delta_k=rfft(delta)
function is_lmax(ds,i,j,k,n)
    ff=ds[i,j,k]
    @inbounds for kk in -1:1, jj in -1:1, ii in -1:1
        ff < ds[mod1(i+ii,n),mod1(j+jj,n),mod1(k+kk,n)] && return false
    end
    true
end
function count_all(delta_k,filters,n,box,fcrit)
    mask=falses(n,n,n); tm=0; tn=0
    for f in filters
        ds=smooth_field(delta_k,n,box,f[3],1; fortran_compat=true)
        for k in 1:n, j in 1:n, i in 1:n
            @inbounds ds[i,j,k]<fcrit && continue
            is_lmax(ds,i,j,k,n) || continue
            tn+=1
            if !mask[i,j,k]; tm+=1; mask[i,j,k]=true; end
        end
    end
    tm,tn
end
tm,tn=count_all(delta_k,filters,n,box,fcrit)
@printf("\nJULIA local maxima > fcrit on THIS field (global, 256^3):\n")
@printf("  masked (cross-scale) = %d\n  nomask               = %d\n", tm, tn)
@printf("  (Fortran-on-its-own-field found 265858 over ~212^3 cores)\n")
