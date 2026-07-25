#!/usr/bin/env julia
# ISW validation: whole-octant pseudo-C_l of our :isw octant map (cf32, Nside 1024,
# DeltaT/T) vs Websky's isw.fits (linear-potential LOS integral, halo catalogs unused).
# Realization variance is LARGE at ISW scales for one octant — the test is amplitude
# + shape over the overlap band, not per-band equality.
using Healpix, Printf, Statistics
const FM = "/home/yguan/scratch/websky_6144/fieldmaps"
w2(l, nside) = exp(-l * (l + 1) * (sqrt(π / 3) / nside)^2 / 12)
function degraded(path)
    m = Healpix.readMapFromFITS(path, 1, Float32)
    ns = m.resolution.nside
    md = ns > 512 ? Healpix.udgrade(m, 512) : m
    md, ns
end
function smooth_map(m, fwhm_rad, nside)
    alm = Healpix.map2alm(m; lmax=400)
    σ = fwhm_rad / sqrt(8 * log(2))
    Healpix.almxfl!(alm, [exp(-0.5 * l * (l + 1) * σ^2) for l in 0:400])
    Healpix.alm2map(alm, nside)
end
ours, ns_o = degraded(joinpath(FM, "isw_websky_6144_oct000_finecell_cf32_nside1024.fits"))
ref, ns_r = degraded("/home/yguan/scratch/websky_6144/websky_ref/isw.fits")
@printf("ours: nside %d->512, std=%.3e mean=%.3e (DeltaT/T)\n", ns_o, std(ours.pixels[ours.pixels .!= 0]), mean(ours.pixels[ours.pixels .!= 0]))
@printf("ref:  nside %d->512, std=%.3e mean=%.3e\n", ns_r, std(ref.pixels), mean(ref.pixels))
# octant mask, apodized
occ = HealpixMap{Float64,RingOrder}(512)
occ.pixels .= Float64.(ours.pixels .!= 0.0)
occ_s = smooth_map(occ, deg2rad(3.0), 512)
apod = HealpixMap{Float64,RingOrder}(512)
apod.pixels .= ifelse.(occ_s.pixels .> 0.97, 1.0, 0.0) .* occ.pixels
apod_s = smooth_map(apod, deg2rad(3.0), 512)
apod_s.pixels .*= occ.pixels
w2m = mean(apod_s.pixels .^ 2)
masked = HealpixMap{Float64,RingOrder}(512)
mu = sum(ours.pixels .* apod_s.pixels) / sum(apod_s.pixels)
masked.pixels .= (Float64.(ours.pixels) .- mu) .* apod_s.pixels
cl_o = Healpix.alm2cl(Healpix.map2alm(masked; lmax=700)) ./ w2m
cl_r = Healpix.alm2cl(Healpix.map2alm(ref; lmax=700))
# same-mask reference: if the ratio vs full-sky ref explodes at high l but the
# same-mask ratio ~1, the "excess" is pseudo-C_l leakage (steep spectrum), not map power
maskedr = HealpixMap{Float64,RingOrder}(512)
mur = sum(ref.pixels .* apod_s.pixels) / sum(apod_s.pixels)
maskedr.pixels .= (Float64.(ref.pixels) .- mur) .* apod_s.pixels
cl_rm = Healpix.alm2cl(Healpix.map2alm(maskedr; lmax=700)) ./ w2m
@printf("%-8s %-12s %-12s %-9s %-9s\n", "ell", "ours", "isw.fits", "r_full", "r_samemask")
for (l0, l1) in [(8,15),(15,25),(25,40),(40,70),(70,120),(120,200),(200,330),(330,540)]
    co = mean(cl_o[l0+1:l1+1] ./ [w2(l,512)*w2(l,ns_o) for l in l0:l1])
    cr = mean(cl_r[l0+1:l1+1] ./ [w2(l,512)*w2(l,ns_r) for l in l0:l1])
    crm = mean(cl_rm[l0+1:l1+1] ./ [w2(l,512)*w2(l,ns_r) for l in l0:l1])
    @printf("%-8.0f %-12.3e %-12.3e %-9.3f %-9.3f\n", sqrt(l0*l1), co, cr, co/cr, co/crm)
end
println("\nratio in (DeltaT/T)^2 if ref is DeltaT/T; if ref is in K, ratio scales by (2.7255)^2.")
println("Low-l bands: octant mask + realization -> factor-2 scatter expected; look for O(1) level + shape.")
