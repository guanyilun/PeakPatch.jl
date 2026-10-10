#!/usr/bin/env julia
# Calibrate the Huntian CIB normalization to Planck, the way Websky did (Stein+2020; Lee+2024 §2): one overall
# factor s on every galaxy flux, chosen so that the 545 GHz auto power in the Planck ℓ = 411–592 bin, with
# Planck's 350 mJy point-source cut, equals Planck 2013 XXX Table D.2 (1.91e4 Jy²/sr).
# The painter (XGPaint CIB_Planck2013) gives flux ∝ 1/shang_I0, so this is exactly a rescaling of the maps;
# the effective shang_I0 is 92/s. The same s is applied to all six frequencies.
# Then every map (Websky, Huntian v0.1, calibrated) is compared with Planck Table D.2 in Planck's bins.
# Caveat: Planck 2013 absolute calibration at 545/857 GHz is 10% in intensity (Table 1), ~20% in power.
#   env: WS_REF, HUNTIAN (input release), OUTMAPS (dir for calibrated maps)
using Healpix, Printf, Statistics

const WREF = get(ENV, "WS_REF", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/websky_6144/websky_ref")
const HT = get(ENV, "HUNTIAN", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/huntian/v0.1")
const OUTMAPS = get(ENV, "OUTMAPS", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/websky_6144/cib_calibrated_v5"); mkpath(OUTMAPS)
const OUT = joinpath(@__DIR__, "..", "..", "paper", "results", "cib_calibration_v5.txt")
const I0_PAINTED = 92.0
# Planck 2013 XXX: bins (Table D.1), flux cuts [Jy] (Table 1), CIB auto spectra [Jy²/sr] (Table D.2)
const BINS = [(145, 229), (229, 411), (411, 592), (592, 774), (774, 1006), (1006, 1308), (1308, 1701), (1701, 2211), (2211, 3085)]
const CUT = Dict("0143" => 0.350, "0217" => 0.225, "0353" => 0.315, "0545" => 0.350, "0857" => 0.710)
const PLANCK = Dict(
    "0857" => [2.87e5, 1.34e5, 7.20e4, 4.38e4, 3.23e4, 2.40e4, 1.83e4, 1.46e4, 1.16e4],
    "0545" => [6.63e4, 3.34e4, 1.91e4, 1.25e4, 9.17e3, 6.83e3, 5.34e3, 4.24e3, 3.42e3],
    "0353" => [7.88e3, 4.35e3, 2.60e3, 1.74e3, 1.29e3, 9.35e2, 7.45e2, 6.08e2, NaN],
    "0217" => [4.17e2, 2.62e2, 1.75e2, 1.17e2, 8.82e1, 6.42e1, 3.34e1, 4.74e1, NaN],
    "0143" => [3.64e1, 3.23e1, 2.81e1, 2.27e1, 1.84e1, 1.58e1, 1.25e1, NaN, NaN])
const LMAX = 3085

readm(p) = Healpix.readMapFromFITS(p, 1, Float64)

# binned auto spectrum [Jy²/sr] of s·map after the Planck cut (bright pixels -> map mean), pixel window removed
function binned(m, s, cut; lmax=LMAX, bins=BINS)
    x = s .* m.pixels; Ωpix = 4π / length(x); μ = mean(x)
    x[x .* 1e6 .* Ωpix .> cut] .= μ
    hm = HealpixMap{Float64,RingOrder}(x .* 1e6)                 # MJy/sr -> Jy/sr
    cl = Healpix.alm2cl(Healpix.map2alm(hm; lmax=lmax, mmax=lmax, niter=0))
    pw = Healpix.pixwin(m.resolution.nside); pw isa Tuple && (pw = pw[1])
    [mean(cl[l+1] / pw[l+1]^2 for l in lo:hi) for (lo, hi) in bins if hi <= lmax]
end

function main()
    io = open(OUT, "w"); say(a...) = (local line = string(a...); println(line); println(io, line))
    say("# Huntian (v5) CIB calibration to Planck 2013 XXX, 545 GHz, ℓ = 411–592, 350 mJy cut")
    m545 = readm(joinpath(HT, "maps", "cib_nu0545.fits")); target = PLANCK["0545"][3]
    s = 1.0
    for it in 1:4                                                 # the cut makes C(s) slightly non-quadratic
        c = binned(m545, s, CUT["0545"]; lmax=592, bins=[(411, 592)])[1]
        say(@sprintf("iter %d: s = %.5f  C = %.5e  C/Planck = %.5f", it, s, c, c / target))
        abs(c / target - 1) < 1e-4 && break
        s *= sqrt(target / c)
    end
    say(@sprintf("\nscale s = %.5f  ->  effective shang_I0 = %.2f (painted with %.0f)", s, I0_PAINTED / s, I0_PAINTED))
    for nu in ("0100", "0143", "0217", "0353", "0545", "0857")
        m = readm(joinpath(HT, "maps", "cib_nu$(nu).fits")); m.pixels .*= s
        hm = HealpixMap{Float32,RingOrder}(Float32.(m.pixels))
        Healpix.saveToFITS(hm, "!" * joinpath(OUTMAPS, "cib_nu$(nu).fits"), typechar="E")
    end
    say("calibrated maps written to $OUTMAPS")
    say("\n== auto spectra / Planck 2013 Table D.2, per Planck bin (ℓ_lo-ℓ_hi), Planck flux cuts applied")
    for nu in ("0143", "0217", "0353", "0545", "0857")
        mW = readm(joinpath(WREF, "cib_nu$(nu).fits")); mH = readm(joinpath(HT, "maps", "cib_nu$(nu).fits"))
        cW = binned(mW, 1.0, CUT[nu]); cH = binned(mH, 1.0, CUT[nu]); cC = binned(mH, s, CUT[nu])
        say(@sprintf("\n-- %s GHz (cut %.0f mJy)\n%-11s %11s %9s %9s %9s", nu, 1e3CUT[nu], "bin", "Planck", "Websky", "v0.1", "calib"))
        for (k, (lo, hi)) in enumerate(BINS)
            p = PLANCK[nu][k]; isnan(p) && continue
            say(@sprintf("%4d-%-6d %11.4e %9.3f %9.3f %9.3f", lo, hi, p, cW[k] / p, cH[k] / p, cC[k] / p))
        end
    end
    close(io)
end
main()
