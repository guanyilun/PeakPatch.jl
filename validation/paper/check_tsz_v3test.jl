#!/usr/bin/env julia
# nbuff mass-cap test (V2_RESULTS_2026-09-27.md, update item 3): tSZ on octant 000 alone,
# v2 (nbuff 16, M ≤ 7.5e14) and v3test (nbuff 25) against the released tsz_2048 on the same
# apodized octant footprint, at the reference Nside 2048.
include(joinpath(@__DIR__, "spectra.jl"))
const D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const WREF = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"
ns = 2048; w = octant_mask(ns, "000"); w2 = mean(abs2, w.pixels); inoct = w.pixels .> 0.999
edges = lbins(4000; lmin=20)
function spec(m)
    x = HealpixMap{Float64,RingOrder}(ns); x.pixels .= m.pixels .* w.pixels
    c, le, _ = binned(clof(almof(x, 4000)) ./ w2, edges)
    (c, le, mean(m.pixels[inoct]))
end
cr, le, mr = spec(loadmap(joinpath(WREF, "tsz_2048.fits")))
c2, _, m2 = spec(loadmap(joinpath(D, "halomaps_v2", "tsz_y_v2_oct000_AM_nside4096.fits"); nside=ns))
c3, _, m3 = spec(loadmap("/home/yguan/scratch/websky_6144/v3test_paint/tsz_y_v3test_oct000_AM_nside4096.fits"; nside=ns))
@printf("mean y / websky: v2 %.4f  v3test %.4f\n", m2 / mr, m3 / mr)
@printf("%-8s %-12s %-10s %-10s %-10s\n", "ell", "C_websky", "v2/wsky", "v3/wsky", "v3/v2")
for b in eachindex(le)
    @printf("%-8.0f %-12.4e %-10.3f %-10.3f %-10.3f\n", le[b], cr[b], c2[b] / cr[b], c3[b] / cr[b], c3[b] / c2[b])
end
