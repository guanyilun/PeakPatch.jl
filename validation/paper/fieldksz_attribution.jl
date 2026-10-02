#!/usr/bin/env julia
# Which v2→v3 change removed the field-kSZ excess at ℓ = 70–200 (v2 1.13–1.21 × linear theory,
# v3 1.01–1.10)? Two changes coincide: the Ω_m(a) 2LPT fix and the tile layout (n 416→434,
# nbuff 16→25). Map A = oct000 field map with the v2 layout (nbuff 16) and the CURRENT code (fix in);
# map B = the v3 oct000 field map (nbuff 25, fix in). If A/B ≈ 1, the layout is not responsible
# (→ the fix is); if A/B ≈ v2/v3 (~1.1 at ℓ 70–200), the layout is. Octant-000 apodized mask.
include(joinpath(@__DIR__, "spectra.jl"))
const A = get(ENV, "FA_A", "/home/yguan/scratch/websky_6144/fieldattr/ksz_v2_oct000_nside4096.fits")
const B = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/fieldmaps_v3/ksz_v3_oct000_nside4096.fits"
const KA = replace(A, "ksz_" => "kappa_"); const KB = replace(B, "ksz_" => "kappa_")
ns = 4096; w = octant_mask(ns, "000"); w2 = mean(abs2, w.pixels); edges = lbins(1000; lmin=20)
spec(p) = begin m = loadmap(p); x = HealpixMap{Float64,RingOrder}(ns); x.pixels .= m.pixels .* w.pixels
                 c, le, _ = binned(clof(almof(x, 1000)) ./ w2, edges); (c, le) end
ca, le = spec(A); cb, _ = spec(B); ka, _ = spec(KA); kb, _ = spec(KB)
# full-sky v2/v3 field-kSZ ratio for reference (theory tables: col 2 = C_campaign)
rd(f) = [parse.(Float64, split(l)) for l in readlines(joinpath(@__DIR__, "results", f)) if !startswith(l, "#") && !isempty(strip(l))]
t2 = rd("theory_v2_ksz_field.txt"); t3 = rd("theory_v3_ksz_field.txt")
out = open(joinpath(@__DIR__, "results", get(ENV, "FA_OUT", "fieldksz_attribution.txt")), "w")
say(a...) = begin s = string(a...); println(s); println(out, s) end
say("oct000 field maps: A = ", basename(A), "; B = v3 (nbuff 25) + fix")
say(@sprintf("%-7s %-10s %-10s %-12s", "ell", "kSZ A/B", "κ A/B", "[full-sky v2/v3 kSZ]"))
for b in eachindex(le)
    j = argmin(abs.([r[1] for r in t2] .- le[b])); ref = abs(t2[j][1] - le[b]) / le[b] < 0.05 ? t2[j][2] / t3[j][2] : NaN
    say(@sprintf("%-7.0f %-10.4f %-10.4f %-12.4f", le[b], ca[b] / cb[b], ka[b] / kb[b], ref))
end
close(out)
