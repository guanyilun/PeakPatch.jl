#!/usr/bin/env julia
# Painter test for the tSZ C_ℓ deficit (V2_RESULTS_2026-09-27.md, open item 1): Websky's OWN
# halos (halos_10x10.pksc, converted to our finalized format) painted by our production
# painter (paint_octant.jl, tsz), vs the released post-2022-fix tsz_2048.fits on the SAME
# pixels (apodized disc inside the 10°×10° patch, centred on the patch axis ≈ +x).
# Same sky, same halos → any ratio ≠ 1 is painter/map-making, not catalog.
include(joinpath(@__DIR__, "spectra.jl"))
const W = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"
ours = loadmap("/home/yguan/scratch/websky_6144/wsky_patch_paint/tsz_y_wskypatch_nside2048.fits")
ref = loadmap(joinpath(W, "tsz_2048.fits"))
ns = 2048; res = Healpix.Resolution(ns); npix = 12ns^2
ax = (1.0, -4.1358e-4, -1.4055e-4); axn = sqrt(sum(abs2, ax)); ax = ax ./ axn
for (θr, θap) in ((4.5, 0.5), (3.0, 0.5))
    w = zeros(npix)
    for p in 1:npix
        v = Healpix.pix2vecRing(res, p)
        θ = rad2deg(acos(clamp(v[1] * ax[1] + v[2] * ax[2] + v[3] * ax[3], -1.0, 1.0)))
        θ >= θr && continue
        d = θr - θ
        w[p] = d >= θap ? 1.0 : 0.5 - 0.5cos(π * d / θap)
    end
    inner = w .> 0.999; w2 = mean(abs2, w)
    @printf("\n== disc %.1f° (apod %.1f°): mean y ours %.4e  ref %.4e  ratio %.4f\n", θr, θap,
            mean(ours.pixels[inner]), mean(ref.pixels[inner]), mean(ours.pixels[inner]) / mean(ref.pixels[inner]))
    edges = [150, 250, 400, 650, 1000, 1500, 2200, 3200, 4096]
    f(m) = (x = HealpixMap{Float64,RingOrder}(ns); μ = sum(m.pixels .* w) / sum(w);
            x.pixels .= (m.pixels .- μ) .* w; x)
    co, le, _ = binned(clof(almof(f(ours), 4095)) ./ w2, edges)
    cr, _, _ = binned(clof(almof(f(ref), 4095)) ./ w2, edges)
    cx, _, _ = binned(clof(almof(f(ours), 4095), almof(f(ref), 4095)) ./ w2, edges)
    @printf("%-7s %-11s %-11s %-8s %-8s\n", "ell", "C_ours", "C_ref", "ratio", "r_cross")
    for b in eachindex(le)
        @printf("%-7.0f %-11.3e %-11.3e %-8.3f %-8.3f\n", le[b], co[b], cr[b], co[b] / cr[b], cx[b] / sqrt(co[b] * cr[b]))
    end
end
