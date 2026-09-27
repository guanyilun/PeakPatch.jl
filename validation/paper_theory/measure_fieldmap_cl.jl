#!/usr/bin/env julia
# Absolute pseudo-C_ℓ of the frozen-campaign FIELD maps (κ and kSZ, one octant) for the
# theory-anchor overlay (THEORY_ANCHORS_2026-09-26.md). Estimator as in
# ../websky_6144/compare_cf32_ksz.jl (apodized octant mask, 1/⟨w²⟩ correction,
# Gaussian pixel-window approximation), here at Nside 2048 to reach ℓ≈3000.
#   julia --project=validation validation/paper_theory/measure_fieldmap_cl.jl [oct]
using Healpix, Printf, Statistics
const FM  = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/fieldmaps_prod"
const OUT = joinpath(@__DIR__, "results")
const NS, LMAX = 2048, 3200
const TCMB_UK = 2.7255e6
oct = length(ARGS) >= 1 ? ARGS[1] : "000"
w2(l, nside) = exp(-l * (l + 1) * (sqrt(π / 3) / nside)^2 / 12)
function loadmap(path, scale)
    m = Healpix.readMapFromFITS(path, 1, Float64)
    ns = m.resolution.nside
    md = Healpix.udgrade(m, NS); m = nothing; GC.gc()
    md.pixels .*= scale
    md, ns
end
function smooth_map(m, fwhm_rad; lmax=800)
    alm = Healpix.map2alm(m; lmax=lmax)
    σ = fwhm_rad / sqrt(8 * log(2))
    Healpix.almxfl!(alm, [exp(-0.5 * l * (l + 1) * σ^2) for l in 0:lmax])
    Healpix.alm2map(alm, NS)
end
# Geometric octant mask from pixel directions (NOT map != 0): near-observer cells displaced
# across the octant planes paint signal outside the octant (physical — it lands in the
# neighbouring octant's pixels of an assembled full sky — but it corrupts single-octant
# estimates). octZYX bit = 1 → observer at +2618 on that axis → octant looks toward −axis.
function octant_weights(md, oct)
    sgn = [oct[4-d] == '1' ? -1.0 : 1.0 for d in 1:3]          # x ← oct[3], y ← oct[2], z ← oct[1]
    res = md.resolution
    occ = HealpixMap{Float64,RingOrder}(NS)
    Threads.@threads for p in 1:length(occ.pixels)
        v = Healpix.pix2vecRing(res, p)
        occ.pixels[p] = (sgn[1] * v[1] > 0 && sgn[2] * v[2] > 0 && sgn[3] * v[3] > 0) ? 1.0 : 0.0
    end
    @printf("spill: fraction of map power outside the geometric octant = %.3e\n",
            sum(abs2, md.pixels .* (1 .- occ.pixels)) / sum(abs2, md.pixels))
    occ_s = smooth_map(occ, deg2rad(2.0))
    apod = HealpixMap{Float64,RingOrder}(NS)
    apod.pixels .= ifelse.(occ_s.pixels .> 0.97, 1.0, 0.0) .* occ.pixels
    apod_s = smooth_map(apod, deg2rad(2.0)); apod_s.pixels .*= occ.pixels
    apod_s
end
function pseudo_cl(md, wts)
    w2m = mean(wts.pixels .^ 2)
    mu = sum(md.pixels .* wts.pixels) / sum(wts.pixels)
    masked = HealpixMap{Float64,RingOrder}(NS)
    masked.pixels .= (md.pixels .- mu) .* wts.pixels
    Healpix.alm2cl(Healpix.map2alm(masked; lmax=LMAX)) ./ w2m
end
kap, nsk = loadmap(joinpath(FM, "kappa_prod_oct$(oct)_nside4096.fits"), 1.0)
wts = octant_weights(kap, oct)
@printf("fsky(w>0)=%.4f  <w^2>=%.4f\n", mean(wts.pixels .> 0), mean(wts.pixels .^ 2))
clk = pseudo_cl(kap, wts); kap = nothing; GC.gc()
ksz, nsz = loadmap(joinpath(FM, "ksz_prod_oct$(oct)_nside4096.fits"), TCMB_UK)
clz = pseudo_cl(ksz, wts); ksz = nothing; GC.gc()
edges = unique(round.(Int, exp.(range(log(20), log(3000), length=26))))
mkpath(OUT)
open(joinpath(OUT, "fieldmap_cl_prod_oct$(oct).txt"), "w") do io
    println(io, "# frozen-campaign field maps oct$(oct) (Nside 4096 -> $(NS)), apodized geometric-octant pseudo-C_l / <w^2>, pixwin-corrected (gaussian approx)")
    println(io, "# l_lo l_hi l_eff  Cl_kappa  Dl_ksz[uK^2]")
    for i in 1:length(edges)-1
        l0, l1 = edges[i], edges[i+1] - 1
        ls = l0:l1
        pw = [w2(l, NS) * w2(l, 4096) for l in ls]
        ck = mean(clk[ls .+ 1] ./ pw)
        dz = mean(clz[ls .+ 1] ./ pw .* ls .* (ls .+ 1) ./ (2π))
        @printf(io, "%d %d %.1f %.6e %.6e\n", l0, l1, sqrt(l0 * l1), ck, dz)
    end
end
println("wrote ", joinpath(OUT, "fieldmap_cl_prod_oct$(oct).txt"))
