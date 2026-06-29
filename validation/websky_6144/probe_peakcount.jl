#!/usr/bin/env julia
# Count local maxima > fcrit per filter scale on the GLOBAL 256^3 field, with and
# without cross-scale masking. Isolates peak-finding LOGIC from tiling.
# Fortran total peaks found (8 boxes) = 265858.
using PeakPatch, FFTW, Printf, Statistics
import PeakPatch.PowerSpectrum: load_pk
import PeakPatch.Filters: smooth_field

n = 256; box = 320.8333; seed = 13579; fcrit = 1.686f0
pk = load_pk(joinpath(@__DIR__, "data", "pk_websky.dat"))
filters = PeakPatch.read_filterbank(joinpath(@__DIR__, "data", "filters_websky.dat"))
sort!(filters; by=f->-f[3])   # large Rf first (peak-finding order)

delta = generate_grf(n, pk, box, seed; fortran_compat=true)
delta_k = rfft(delta)

# local maximum test over 3x3x3 (periodic), strict >= all neighbours
function is_localmax(ds, i, j, k, n)
    ff = ds[i,j,k]
    @inbounds for kk in -1:1, jj in -1:1, ii in -1:1
        iii = mod1(i+ii, n); jjj = mod1(j+jj, n); kkk = mod1(k+kk, n)
        ff < ds[iii,jjj,kkk] && return false
    end
    return true
end

function run(delta_k, filters, n, box, fcrit)
    mask = falses(n,n,n)
    tot_masked = 0; tot_nomask = 0
    @printf("%-8s %-12s %-12s\n","Rf","Npk_masked","Npk_nomask")
    for f in filters
        Rf = f[3]
        ds = smooth_field(delta_k, n, box, Rf, 1; fortran_compat=true)
        cnt_m = 0; cnt_n = 0
        for k in 1:n, j in 1:n, i in 1:n
            @inbounds ds[i,j,k] < fcrit && continue
            is_localmax(ds, i, j, k, n) || continue
            cnt_n += 1
            if !mask[i,j,k]
                cnt_m += 1
                mask[i,j,k] = true
            end
        end
        tot_masked += cnt_m; tot_nomask += cnt_n
        Rf < 8.0 && @printf("%-8.3f %-12d %-12d\n", Rf, cnt_m, cnt_n)
    end
    @printf("\nTOTAL  masked=%d  nomask=%d   [Fortran found=265858]\n", tot_masked, tot_nomask)
end
run(delta_k, filters, n, box, fcrit)
