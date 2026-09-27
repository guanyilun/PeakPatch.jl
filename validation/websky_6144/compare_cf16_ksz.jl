#!/usr/bin/env julia
# kSZ low-ell Doppler diagnosis, step 2: whole-octant pseudo-C_l of the
# coarse_factor=16 full-matter field maps vs the production coarse_factor=4 maps
# vs ksz.fits. Also kappa as a CONTROL (density projection — should be unchanged).
#
# Hypothesis under test (KSZ_COMPOSITE_2026-07-19.md): the 64^3 coarse grid
# (Nyquist 0.038 h/Mpc) leaves the k~0.03-0.1 band to per-tile fine noise ->
# excess/decoherent large-scale VELOCITY power -> 2-3x kSZ Doppler excess at
# ell<350. If cf16 (Nyquist 0.154) collapses the excess, mechanism confirmed.
using Healpix, Printf, Statistics

const FM = "/home/yguan/scratch/websky_6144/fieldmaps"
const TCMB_UK = 2.7255e6
w2(l, nside) = exp(-l * (l + 1) * (sqrt(π / 3) / nside)^2 / 12)

function degraded(path, scale)
    m = Healpix.readMapFromFITS(path, 1, Float32)
    ns = m.resolution.nside
    md = Healpix.udgrade(m, 1024); m = nothing; GC.gc()
    md.pixels .*= scale
    md, ns
end
function smooth_map(m, fwhm_rad)
    alm = Healpix.map2alm(m; lmax=600)
    σ = fwhm_rad / sqrt(8 * log(2))
    Healpix.almxfl!(alm, [exp(-0.5 * l * (l + 1) * σ^2) for l in 0:600])
    Healpix.alm2map(alm, 1024)
end

# apodized octant mask from the first map's footprint
function octant_weights(md)
    occ = HealpixMap{Float64,RingOrder}(1024)
    occ.pixels .= Float64.(md.pixels .!= 0.0)
    occ_s = smooth_map(occ, deg2rad(2.0))
    apod = HealpixMap{Float64,RingOrder}(1024)
    apod.pixels .= ifelse.(occ_s.pixels .> 0.97, 1.0, 0.0) .* occ.pixels
    apod_s = smooth_map(apod, deg2rad(2.0))
    apod_s.pixels .*= occ.pixels
    apod_s
end
function pseudo_cl(md, wts)
    w2m = mean(wts.pixels .^ 2)
    mu = sum(md.pixels .* wts.pixels) / sum(wts.pixels)
    masked = HealpixMap{Float64,RingOrder}(1024)
    masked.pixels .= (md.pixels .- mu) .* wts.pixels
    Healpix.alm2cl(Healpix.map2alm(masked; lmax=2000)) ./ w2m
end

maps = [
    ("ksz_cf16", joinpath(FM, "ksz_websky_6144_oct000_finecell_cf16_nside2048.fits"), TCMB_UK),
    ("ksz_cf4",  joinpath(FM, "ksz_websky_6144_oct000_finecell_nside4096.fits"), TCMB_UK),
    ("kap_cf16", joinpath(FM, "kappa_websky_6144_oct000_finecell_cf16_nside2048.fits"), 1.0),
    ("kap_cf4",  joinpath(FM, "kappa_websky_6144_oct000_finecell_nside4096.fits"), 1.0),
]
cls = Dict{String,Vector{Float64}}(); nss = Dict{String,Int}()
wts = nothing
for (nm, path, sc) in maps
    md, ns = degraded(path, sc)
    global wts
    wts === nothing && (wts = octant_weights(md))
    cls[nm] = pseudo_cl(md, wts); nss[nm] = ns
    @info "done" nm ns
end
ref, _ = degraded("/home/yguan/scratch/websky_6144/websky_ref/ksz.fits", 1.0)
cls["ksz.fits"] = Healpix.alm2cl(Healpix.map2alm(ref; lmax=2000)); nss["ksz.fits"] = 4096

bands = [(100,140),(140,190),(190,270),(270,380),(380,530),(530,740),(740,1030),(1030,1440),(1440,1990)]
bp(nm, l0, l1) = mean(cls[nm][l0+1:l1+1] ./ [w2(l,1024)*w2(l,nss[nm]) for l in l0:l1])
@printf("%-6s %-11s %-11s %-11s | %-9s %-9s | %-9s\n",
        "ell", "ksz_cf16", "ksz_cf4", "ksz.fits", "cf16/ref", "cf4/ref", "kap16/4")
for (l0, l1) in bands
    k16 = bp("ksz_cf16", l0, l1); k4 = bp("ksz_cf4", l0, l1); kr = bp("ksz.fits", l0, l1)
    kap = bp("kap_cf16", l0, l1) / bp("kap_cf4", l0, l1)
    @printf("%-6.0f %-11.3e %-11.3e %-11.3e | %-9.3f %-9.3f | %-9.3f\n",
            sqrt(l0 * l1), k16, k4, kr, k16 / kr, k4 / kr, kap)
end
println("\ncf16/ref -> 1 at ell<350 while kap16/4 ~ 1 everywhere = coarse-grid velocity")
println("mechanism CONFIRMED (kappa is the density control; realization differs only in")
println("the coarse/fine noise split, so kap16/4 far from 1 would flag a different issue).")
