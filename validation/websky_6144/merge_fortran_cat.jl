#!/usr/bin/env julia
# Apply Julia's binary exclusion (merge_catalog) to the FORTRAN raw catalog (224K).
# Fortran cols: 0-2 Eulerian xyz, 3-5 1LPT disp, 6 RTHL, 7-9 2LPT disp.
# Lagrangian pos = Eulerian - (1LPT + 2LPT). Exclusion uses Lagrangian center + RTHL.
# If Fortran-merged ~ 69K (Julia A/B) -> masking is a valid shortcut.
# If Fortran-merged >> 69K -> Julia under-finds (masking discards real halos).
using PeakPatch, Printf
import PeakPatch.Catalog: HaloRecord
import PeakPatch.Merger: merge_catalog

fn = "/home/yguan/scratch/websky_6144/fortran_floortest/output/catalog_floortest.pksc.13579"
raw = read(fn)
N = reinterpret(Int32, raw[1:4])[1] |> Int
nf = 11
flt = reinterpret(Float32, raw[13:13+N*nf*4-1])
rho_m = 2.775e11 * 0.31
get(h,c) = Float64(flt[(h-1)*nf + c + 1])   # c is 0-based column

halos = Vector{HaloRecord}(undef, N)
for h in 1:N
    # cols 0-2 are xpk = LAGRANGIAN grid position (xbx+alatt*(i-cen)); use directly.
    xl,yl,zl = get(h,0),get(h,1),get(h,2)
    R = get(h,6)
    halos[h] = HaloRecord(Float32(xl),Float32(yl),Float32(zl),
                          0f0,0f0,0f0, Float32(R), 0f0,0f0,0f0, 0f0)
end

NgtM(hs,M0) = count(h -> (4/3*pi*Float64(h.RTHL)^3*rho_m) > M0, hs)
function report(lbl, hs)
    @printf("%-22s N=%8d | >1.69e12=%8d >1e13=%7d >1e14=%6d\n", lbl, length(hs),
            NgtM(hs,1.69e12), NgtM(hs,1e13), NgtM(hs,1e14))
end

report("FORTRAN raw (pre-merge)", halos)
merged = merge_catalog(halos; verbose=true)
report("FORTRAN + Julia merge", merged)
@printf("\nJulia A/B raw (pre-merge) total = 69538 (already masked during find)\n")
