#!/usr/bin/env julia
# Diagnostics: (1) where are the exactly-zero pixels of the assembled field maps (octant
# seams? near-observer holes?); (2) per-octant tSZ ours/ref C_ℓ ratios.
include(joinpath(@__DIR__, "spectra.jl"))
const F = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/fullsky_prod"
const FM = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/fieldmaps_prod"
const W = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"

function holes(path)
    m = loadmap(path); res = m.resolution
    z = findall(iszero, m.pixels)
    dmin = [minimum(abs.(Healpix.pix2vecRing(res, p))) for p in z]   # sin(dist to nearest plane)
    th = asind.(dmin)
    @printf("%s: %d zero pixels; distance to nearest octant plane: <0.02° %d, 0.02-0.1° %d, 0.1-1° %d, >1° %d\n",
            basename(path), length(z), count(<(0.02), th), count(x -> 0.02 <= x < 0.1, th),
            count(x -> 0.1 <= x < 1, th), count(>=(1), th))
end
holes(joinpath(F, "kappa_field_prod_fullsky_nside4096.fits"))
holes(joinpath(F, "isw_uK_prod_fullsky_nside4096.fits"))

# seam profile: rms (and mean) of a map vs angular distance to the nearest octant plane
function seam_profile(label, path; scale=1.0)
    m = loadmap(path; scale=scale); res = m.resolution; npix = length(m.pixels)
    edges = [0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
    s1 = zeros(7); s2 = zeros(7); nn = zeros(Int, 7)
    for p in 1:7:npix                                  # 1/7 subsample is plenty
        v = Healpix.pix2vecRing(res, p)
        th = asind(minimum(abs.(v))); b = searchsortedlast(edges, th)
        (1 <= b <= 7) || continue
        x = m.pixels[p]; s1[b] += x; s2[b] += x^2; nn[b] += 1
    end
    @printf("%-28s", label)
    for b in 1:7
        @printf(" [%.2f-%.2f°] mean %.3e rms %.3e |", edges[b], edges[b+1], s1[b] / nn[b], sqrt(s2[b] / nn[b]))
    end
    println()
end
seam_profile("ours tSZ y", joinpath(F, "tsz_y_prod_fullsky_nside4096.fits"))
seam_profile("websky tsz_2048", joinpath(W, "tsz_2048.fits"))
seam_profile("ours kSZ total uK", joinpath(F, "ksz_total_uK_prod_fullsky_nside4096.fits"))
seam_profile("websky ksz", joinpath(W, "ksz.fits"))
seam_profile("websky cib545", joinpath(W, "cib_nu0545.fits"))
seam_profile("ours kappa_lt4.5", joinpath(F, "kappa_lt4.5_prod_fullsky_nside4096.fits"))

# per-octant tSZ at Nside 2048
mo = loadmap(joinpath(F, "tsz_y_prod_fullsky_nside4096.fits"); nside=2048)
mr = loadmap(joinpath(W, "tsz_2048.fits"))
edges = [100, 200, 400, 800, 1600, 3200, 4097]
S = octant_spectra(Dict("o" => mo, "r" => mr), [("o", "o"), ("r", "r")], 4096, edges)
@printf("\nper-octant tSZ ours/ref (bands %s)\n", string(edges))
for (i, oct) in enumerate(OCTS)
    @printf("oct%s: %s\n", oct, join((@sprintf("%.3f", S[("o","o")][b, i] / S[("r","r")][b, i]) for b in 1:length(edges)-1), " "))
end
@printf("ours  C(400-800) per octant: %s\n", join((@sprintf("%.2e", S[("o","o")][3, i]) for i in 1:8), " "))
@printf("ref   C(400-800) per octant: %s\n", join((@sprintf("%.2e", S[("r","r")][3, i]) for i in 1:8), " "))
