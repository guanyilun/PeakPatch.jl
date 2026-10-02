#!/usr/bin/env julia
# Websky's OWN halos (halos_10x10.pksc) painted three ways vs the released tsz_2048.fits on the
# same 4.5° (and 3°) apodized disc, as in check_tsz_painter_wsky.jl:
#   ours      — our port of the Fortran pks2map physics (paint_octant.jl)
#   xg_asis   — upstream XGPaint Battaglia16ThermalSZProfile, Websky mass passed as M200c
#   xg_nfw7   — the same with M200c from NFW c=7 (paint_tsz_xgpaint_wsky.jl)
include(joinpath(@__DIR__, "spectra.jl"))
const W = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"
const PD = "/home/yguan/scratch/websky_6144/wsky_patch_paint"
maps = [("ours", joinpath(PD, "tsz_y_wskypatch_nside2048.fits")),
        ("xg_asis", joinpath(PD, "tsz_y_xgpaint_asis_wskypatch_nside2048.fits")),
        ("xg_nfw7", joinpath(PD, "tsz_y_xgpaint_nfw7_wskypatch_nside2048.fits"))]
# extra maps as label=path arguments (e.g. Fortran-emulation variants of our painter)
for a in ARGS; lab, path = split(a, "="; limit=2); push!(maps, (String(lab), String(path))); end
const OUTF = isempty(ARGS) ? "tsz_three_painters_wskypatch.txt" : "tsz_painter_variants_wskypatch.txt"
ref = loadmap(joinpath(W, "tsz_2048.fits"))
ns = 2048; res = Healpix.Resolution(ns); npix = 12ns^2
ax = (1.0, -4.1358e-4, -1.4055e-4); ax = ax ./ sqrt(sum(abs2, ax))
out = open(joinpath(@__DIR__, "results", OUTF), "w")
say(a...) = begin s = string(a...); println(s); println(out, s) end
for (θr, θap) in ((4.5, 0.5), (3.0, 0.5))
    w = zeros(npix)
    for p in 1:npix
        v = Healpix.pix2vecRing(res, p)
        θ = rad2deg(acos(clamp(v[1] * ax[1] + v[2] * ax[2] + v[3] * ax[3], -1.0, 1.0)))
        θ >= θr && continue
        d = θr - θ; w[p] = d >= θap ? 1.0 : 0.5 - 0.5cos(π * d / θap)
    end
    inner = w .> 0.999; w2 = mean(abs2, w)
    edges = [150, 250, 400, 650, 1000, 1500, 2200, 3200, 4096]
    f(m) = begin x = HealpixMap{Float64,RingOrder}(ns); μ = sum(m.pixels .* w) / sum(w); x.pixels .= (m.pixels .- μ) .* w; x end
    ar = almof(f(ref), 4095); cr, le, _ = binned(clof(ar) ./ w2, edges)
    say(@sprintf("\n== disc %.1f° (apod %.1f°); released tsz_2048 mean y %.4e", θr, θap, mean(ref.pixels[inner])))
    hdr = @sprintf("%-8s", "ell"); for (lab, _) in maps; hdr *= @sprintf(" %-18s", "$(lab): C/ref r"); end; say(hdr)
    rows = [Float64[] for _ in eachindex(le)]
    for (lab, path) in maps
        isfile(path) || (say("missing $path"); continue)
        m = loadmap(path); am = almof(f(m), 4095)
        co, _, _ = binned(clof(am) ./ w2, edges); cx, _, _ = binned(clof(am, ar) ./ w2, edges)
        say(@sprintf("  %-8s mean y / released %.4f", lab, mean(m.pixels[inner]) / mean(ref.pixels[inner])))
        for b in eachindex(le); push!(rows[b], co[b] / cr[b]); push!(rows[b], cx[b] / sqrt(co[b] * cr[b])); end
    end
    for b in eachindex(le)
        s = @sprintf("%-8.0f", le[b]); for k in 1:2:length(rows[b]); s *= @sprintf(" %8.3f %8.3f  ", rows[b][k], rows[b][k+1]); end; say(s)
    end
end
close(out)
