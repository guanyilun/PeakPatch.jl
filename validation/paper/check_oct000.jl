#!/usr/bin/env julia
# Production-painter sanity check on octant 000 alone, against the released Websky maps on
# the same apodized octant footprint — should land near the Tier-B scorecard
# (TIER_B_SUMMARY_2026-07.md: tSZ mean-y 1.10, C_ℓ 1.18–1.28 with the ref pixwin
# uncorrected; kSZ Wc 1.05–1.22 at ℓ≤320, 1.3–1.5 at ℓ=450–1750, 1.09–1.17 at ℓ≥2400).
include(joinpath(@__DIR__, "spectra.jl"))
const D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const WREF = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"
const TCMB_UK = 2.7255e6

function masked_ratio(label, mo, mr, lmax)
    ns = mo.resolution.nside
    w = octant_mask(ns, "000"); w2 = mean(abs2, w.pixels)
    edges = lbins(lmax; lmin=20)
    f(m) = (x = HealpixMap{Float64,RingOrder}(ns); x.pixels .= m.pixels .* w.pixels; x)
    co, le, nm = binned(clof(almof(f(mo), lmax)) ./ w2, edges)
    cr, _, _ = binned(clof(almof(f(mr), lmax)) ./ w2, edges)
    inoct = w.pixels .> 0.999
    @printf("\n== %s (octant 000, apodized): mean ours %.4e ref %.4e ratio %.4f\n", label,
            mean(mo.pixels[inoct]), mean(mr.pixels[inoct]), mean(mo.pixels[inoct]) / mean(mr.pixels[inoct]))
    @printf("%-8s %-12s %-12s %-8s\n", "ell", "C_ours", "C_ref", "ratio")
    for b in eachindex(co)
        @printf("%-8.0f %-12.4e %-12.4e %-8.3f\n", le[b], co[b], cr[b], co[b] / cr[b])
    end
end

# tSZ at the reference's Nside 2048
y = loadmap(joinpath(D, "halomaps_prod", "tsz_y_prod_oct000_AMv2_nside4096.fits"); nside=2048)
yr = loadmap(joinpath(WREF, "tsz_2048.fits"))
masked_ratio("tSZ y", y, yr, 4000)
y = nothing; yr = nothing; GC.gc()

# kSZ total = T_CMB·field + Wc halos, μK
k = loadmap(joinpath(D, "fieldmaps_prod", "ksz_prod_oct000_nside4096.fits"); scale=TCMB_UK)
k.pixels .+= loadmap(joinpath(D, "halomaps_prod", "ksz_halo_Wc_prod_oct000_AMv2_nside4096.fits")).pixels
kr = loadmap(joinpath(WREF, "ksz.fits"))
masked_ratio("kSZ total (field + Wc)", k, kr, 8000)
