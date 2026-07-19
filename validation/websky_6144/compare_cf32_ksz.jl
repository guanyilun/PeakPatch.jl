#!/usr/bin/env julia
# cf32 convergence check: whole-octant pseudo-C_l of ksz/kappa for coarse_factor
# 4 / 16 / 32 vs ksz.fits. See compare_cf16_ksz.jl for the mechanism test.
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
    ("ksz_cf32", joinpath(FM, "ksz_websky_6144_oct000_finecell_cf32_nside2048.fits"), TCMB_UK),
    ("ksz_cf16", joinpath(FM, "ksz_websky_6144_oct000_finecell_cf16_nside2048.fits"), TCMB_UK),
    ("ksz_cf4",  joinpath(FM, "ksz_websky_6144_oct000_finecell_nside4096.fits"), TCMB_UK),
    ("kap_cf32", joinpath(FM, "kappa_websky_6144_oct000_finecell_cf32_nside2048.fits"), 1.0),
    ("kap_cf16", joinpath(FM, "kappa_websky_6144_oct000_finecell_cf16_nside2048.fits"), 1.0),
]
cls = Dict{String,Vector{Float64}}(); nss = Dict{String,Int}()
wts = nothing
for (nm, path, sc) in maps
    md, ns = degraded(path, sc)
    global wts
    wts === nothing && (wts = octant_weights(md))
    cls[nm] = pseudo_cl(md, wts); nss[nm] = ns
end
ref, _ = degraded("/home/yguan/scratch/websky_6144/websky_ref/ksz.fits", 1.0)
cls["ksz.fits"] = Healpix.alm2cl(Healpix.map2alm(ref; lmax=2000)); nss["ksz.fits"] = 4096
bands = [(100,140),(140,190),(190,270),(270,380),(380,530),(530,740),(740,1030),(1030,1440),(1440,1990)]
bp(nm, l0, l1) = mean(cls[nm][l0+1:l1+1] ./ [w2(l,1024)*w2(l,nss[nm]) for l in l0:l1])
@printf("%-6s | %-9s %-9s %-9s | %-9s\n", "ell", "cf4/ref", "cf16/ref", "cf32/ref", "kap32/16")
for (l0, l1) in bands
    kr = bp("ksz.fits", l0, l1)
    @printf("%-6.0f | %-9.3f %-9.3f %-9.3f | %-9.3f\n", sqrt(l0*l1),
            bp("ksz_cf4", l0, l1) / kr, bp("ksz_cf16", l0, l1) / kr, bp("ksz_cf32", l0, l1) / kr,
            bp("kap_cf32", l0, l1) / bp("kap_cf16", l0, l1))
end
