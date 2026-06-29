#!/usr/bin/env julia
# Measure Julia's smoothed-field sigma(Rf) and #cells>fcrit per filter scale on the
# SAME 256^3 box=320.833 config, compare to Fortran's per-scale table.
# Fortran (box-3 of 8): R=2.068 sig=1.787 Npk=8875 ; 2.378 1.668 6719 ; 2.735 1.554 5183
#   3.145 1.447 3733 ; 3.617 1.341 2603 ; 4.159 1.237 1778 ; ... ; 14.632 0.518 1
using PeakPatch, FFTW, Printf, Statistics
import PeakPatch.PowerSpectrum: load_pk
import PeakPatch.Filters: smooth_field
import PeakPatch.HaloFinder

n = 256; box = 320.8333; seed = 13579; fcrit = 1.686
pk = load_pk(joinpath(@__DIR__, "data", "pk_websky.dat"))
filters = PeakPatch.read_filterbank(joinpath(@__DIR__, "data", "filters_websky.dat"))
sort!(filters; by=f->-f[3])

@info "generating GRF (fortran_compat)..."
delta = generate_grf(n, pk, box, seed; fortran_compat=true)
@printf("unsmoothed field: mean=%.4f sigma=%.4f\n", mean(delta), std(delta))
delta_k = rfft(delta)

FORT = Dict(2.068=>(1.787,8875),2.378=>(1.668,6719),2.735=>(1.554,5183),
            3.145=>(1.447,3733),3.617=>(1.341,2603),4.159=>(1.237,1778),
            5.501=>(1.047,772),7.274=>(0.872,273),9.620=>(0.717,75),14.632=>(0.518,1))
@printf("\n%-8s %-10s %-10s %-12s\n","Rf","sig_julia","sig_fort","ratio J/F")
for f in reverse(filters)   # small Rf first
    Rf = f[3]
    ds = smooth_field(delta_k, n, box, Rf, 1; fortran_compat=true)
    sig = sqrt(mean(ds.^2))
    # nearest Fortran R
    kf = argmin([abs(rf-Rf) for rf in keys(FORT)] |> x->x)
    fk = collect(keys(FORT))[argmin([abs(rf-Rf) for rf in collect(keys(FORT))])]
    sf = abs(fk-Rf) < 0.15 ? FORT[fk][1] : NaN
    @printf("%-8.3f %-10.4f %-10s %-12s\n", Rf, sig,
            isnan(sf) ? "-" : @sprintf("%.4f",sf),
            isnan(sf) ? "-" : @sprintf("%.3f",sig/sf))
end
