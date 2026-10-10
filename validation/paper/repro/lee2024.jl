#!/usr/bin/env julia
# Reproduce the unlensed-CIB statistics of Lee, Bond, Motloch, van Engelen & Stein 2024 (MNRAS 529, 2543;
# arXiv:2304.07283), Fig. 5, on both Websky and Huntian (our v5 / public v0.1) with identical code:
#   band-filtered CIB maps (top-hat ℓ bands, Δℓ = 128, as for the Planck comparison) at 217/353/545 GHz,
#   after replacing pixels brighter than the Planck flux cuts (225/315/350 mJy) with the map mean (their §4),
#   -> band variance S2, skewness S3, kurtosis S4 (eqs 6-8), and from them
#      Ĉ_ℓ = S2 / Σ(2ℓ+1)/4π            (eq 11, Jy²/sr)
#      b̂_ℓ ≈ 2√3 π³ S3 Δℓ⁻³ ℓc⁻¹         (eq 18, equilateral bispectrum, Jy³/sr)
#      t̂_ℓ = S4 / ℓc²                    (the paper's trispectrum proxy; arbitrary normalization)
# The lensed statistics need the paper's shell-by-shell lensing pipeline and are not reproduced.
# Pixel cut at the native Nside 4096; bands are synthesized at Nside 2048 (all bands are below ℓ = 2048).
#   env: WS_REF, HUNTIAN, OUT; LMAXB (top band edge, default 2000)
using Healpix, Printf, Statistics

const WREF = get(ENV, "WS_REF", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/websky_6144/websky_ref")
const HT = get(ENV, "HUNTIAN", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/huntian/v0.1")
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "..", "results", "repro")); mkpath(OUT)
const LMAXB = parse(Int, get(ENV, "LMAXB", "2000"))
const DL = 128; const EDGES = collect(50:DL:LMAXB); const LMAX = EDGES[end]
const NS_BAND = 2048
const FREQS = [("0217", 0.225), ("0353", 0.315), ("0545", 0.350)]   # GHz tag, Planck flux cut [Jy]

# load, cut bright pixels (pixel flux > cut -> map mean), convert MJy/sr -> Jy/sr
function prepared(path, cut)
    m = Healpix.readMapFromFITS(path, 1, Float64)
    Ωpix = 4π / length(m.pixels); μ = mean(m.pixels)
    bright = m.pixels .* 1e6 .* Ωpix .> cut
    m.pixels[bright] .= μ
    m.pixels .*= 1e6
    m, count(bright)
end

# ℓ of every a_ℓm in Healpix.jl's m-major storage
lvec(lmax) = [l for m in 0:lmax for l in m:lmax]

function bandstats(m)
    a = Healpix.map2alm(m; lmax=LMAX, mmax=LMAX, niter=0)
    L = lvec(LMAX); @assert length(L) == length(a.alm)
    rows = NTuple{6,Float64}[]
    for b in 1:length(EDGES)-1
        lo, hi = EDGES[b], EDGES[b+1] - 1
        ab = deepcopy(a); ab.alm[(L .< lo) .| (L .> hi)] .= 0
        f = Healpix.alm2map(ab, NS_BAND).pixels
        d = f .- mean(f)
        S2 = mean(d .^ 2); S3 = mean(d .^ 3); S4 = mean(d .^ 4) - 3S2^2
        lc = (lo + hi) / 2
        C = S2 / sum((2l + 1) / 4π for l in lo:hi)
        bq = 2sqrt(3) * π^3 * S3 / DL^3 / lc
        push!(rows, (lc, C, bq, S4 / lc^2, S3 / S2^1.5, S4 / S2^2))
    end
    rows
end

function main()
    io = open(joinpath(OUT, "lee2024_websky_vs_huntian.txt"), "w")
    say(a...) = (local line = string(a...); println(line); println(io, line))
    say("# Lee et al. 2024 (arXiv:2304.07283) Fig. 5 unlensed-CIB band statistics, Websky vs Huntian v0.1")
    say("# top-hat bands Δℓ = $DL from ℓ = $(EDGES[1]) to $LMAX; bright pixels above the Planck cuts set to the map mean")
    say("# C [Jy²/sr] (eq 11), b_equil [Jy³/sr] (eq 18), t = S4/ℓc² [Jy⁴/sr⁴], and the dimensionless skewness/kurtosis")
    for (nu, cut) in FREQS
        mW, nW = prepared(joinpath(WREF, "cib_nu$(nu).fits"), cut)
        mH, nH = prepared(joinpath(HT, "maps", "cib_nu$(nu).fits"), cut)
        say(@sprintf("\n== %s GHz, cut %.0f mJy: pixels replaced Websky %d, Huntian %d", nu, 1e3cut, nW, nH))
        rW = bandstats(mW); rH = bandstats(mH)
        say(@sprintf("%-7s %11s %11s %6s | %11s %11s %6s | %11s %11s %6s | %7s %7s",
                     "ell_c", "C_W", "C_H", "H/W", "b_W", "b_H", "H/W", "t_W", "t_H", "H/W", "skew_W", "skew_H"))
        for (w, h) in zip(rW, rH)
            say(@sprintf("%-7.0f %11.4e %11.4e %6.3f | %11.4e %11.4e %6.3f | %11.4e %11.4e %6.3f | %7.3f %7.3f",
                         w[1], w[2], h[2], h[2] / w[2], w[3], h[3], h[3] / w[3], w[4], h[4], h[4] / w[4], w[5], h[5]))
        end
    end
    close(io)
end
main()
