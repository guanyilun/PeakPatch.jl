#!/usr/bin/env julia
# Reproduce the Websky-only parts of Rodriguez, Kusiak, Pandey & Hill 2026 (JCAP 06, 004; arXiv:2509.03458)
# on both Websky and Huntian (our v5 / public v0.1), with identical code for the two skies:
#   - Table III: halo counts, 6e13 < M < 5e15 Msun/h, in 0.2<z<0.5, 0.5<z<1, 1<z<2.5
#   - Fig. 2:    v_rms^2(z) of those halos (3D rms peculiar velocity, per Δz = 0.25)
#   - Fig. 3:    tSZ × halo-overdensity cross-spectra in the three z ranges, Δℓ = 250 linear bins
# The paper's kSZ results use their own pasted maps (not Websky products) and are not reproduced here.
# Maps at Nside 2048, because the released Websky y map is Nside 2048; our Nside 4096 y map is degraded to it.
# Masses: M = (4π/3) ρ̄_m R_TH³ in Msun/h for both catalogs (Websky R_TH in Mpc is converted ×h).
#   env: WS_REF (Websky files), HUNTIAN (release dir), OUT (results dir, default results/repro)
using PeakPatch, Healpix, Printf, Statistics
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z

const HUB = 0.68; const RHO_M = 2.775e11 * 0.31
const CHI2Z = build_chi_to_z(CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81); z_max=6.0)
const WREF = get(ENV, "WS_REF", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/websky_6144/websky_ref")
const HT = get(ENV, "HUNTIAN", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/huntian/v0.1")
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "..", "results", "repro")); mkpath(OUT)
const NSIDE = 2048; const RES = Healpix.Resolution(NSIDE); const NPIX = 12NSIDE^2
const MLO, MHI = 6e13, 5e15
const ZBINS = [(0.2, 0.5), (0.5, 1.0), (1.0, 2.5)]
const VZ = collect(0.0:0.25:4.0)                      # v_rms^2(z) bins (Fig. 2 range)
const LMAX = 4000; const DL = 250
const MAXH = parse(Int, get(ENV, "MAXH", "0"))           # smoke test: first MAXH halos per file
const TAG = MAXH > 0 ? "_smoke" : ""
mass(R) = 4 / 3 * π * RHO_M * R^3
obs_of(o) = (o[3] == '1' ? 2618.0 : -2618.0, o[2] == '1' ? 2618.0 : -2618.0, o[1] == '1' ? 2618.0 : -2618.0)

mutable struct Acc
    maps::Vector{Vector{Float64}}                    # halo counts per z range
    v2::Vector{Float64}; nv::Vector{Int}              # Σ|v|², N per VZ bin
end
Acc() = Acc([zeros(NPIX) for _ in ZBINS], zeros(length(VZ) - 1), zeros(Int, length(VZ) - 1))

# one halo in observer-centred Mpc/h, km/s, Msun/h
@inline function add!(a::Acc, dx, dy, dz, vx, vy, vz, M)
    (MLO < M < MHI) || return
    r = sqrt(dx^2 + dy^2 + dz^2); r > 0 || return
    z = chi_to_z(CHI2Z, r)
    iv = searchsortedlast(VZ, z)
    if 1 <= iv < length(VZ); a.v2[iv] += vx^2 + vy^2 + vz^2; a.nv[iv] += 1; end
    for (k, (lo, hi)) in enumerate(ZBINS)
        lo < z < hi || continue
        p = Healpix.ang2pixRing(RES, acos(clamp(dz / r, -1, 1)), mod2pi(atan(dy, dx)))
        a.maps[k][p] += 1
    end
end

function stream(path, ncol, f)
    nh = open(io -> Int(read(io, Int32)), path); MAXH > 0 && (nh = min(nh, MAXH))
    open(path) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        ch = 4_000_000; b = Vector{Float32}(undef, ch * ncol); nd = 0
        while nd < nh
            m = min(ch, nh - nd); read!(io, view(b, 1:m*ncol))
            @inbounds for k in 1:m; f(b, (k - 1) * ncol); end
            nd += m
        end
    end
    nh
end

function catalogs()
    W = Acc(); H = Acc()
    # Websky: 10 columns, Mpc and km/s, observer at the origin
    nw = stream(joinpath(WREF, "halos.pksc"), 10, (b, o) -> add!(W, HUB * b[o+1], HUB * b[o+2], HUB * b[o+3],
                b[o+4], b[o+5], b[o+6], mass(HUB * Float64(b[o+7]))))
    @info "Websky streamed" nw
    for oc in ("000", "001", "010", "011", "100", "101", "110", "111")
        ob = obs_of(oc)
        stream(joinpath(HT, "catalogs", "halos_oct$(oc).pksc"), 33, (b, o) -> add!(H, b[o+1] - ob[1], b[o+2] - ob[2],
               b[o+3] - ob[3], b[o+4], b[o+5], b[o+6], mass(Float64(b[o+7]))))
        @info "Huntian streamed" oc
    end
    W, H
end

delta(m) = (n̄ = sum(m) / NPIX; HealpixMap{Float64,RingOrder}(m ./ n̄ .- 1))
readmap(p) = (m = Healpix.readMapFromFITS(p, 1, Float64); m.resolution.nside == NSIDE ? m : Healpix.udgrade(m, NSIDE))
alm(m) = Healpix.map2alm(m; lmax=LMAX, mmax=LMAX, niter=0)

function main()
    W, H = catalogs()
    io = open(joinpath(OUT, "rodriguez2026_websky_vs_huntian$(TAG).txt"), "w")
    say(a...) = (local line = string(a...); println(line); println(io, line))
    say("# Rodriguez et al. 2026 (arXiv:2509.03458) Websky measurements, repeated on Websky and Huntian v0.1")
    say("# halos 6e13 < M < 5e15 Msun/h (M = 4/3 π ρm R_TH³); maps Nside $NSIDE")
    say("\n== Table III: halo counts (paper, Websky: 332434 / 707105 / 384444)")
    for (k, (lo, hi)) in enumerate(ZBINS)
        nw, nh = Int(sum(W.maps[k])), Int(sum(H.maps[k]))
        say(@sprintf("%.1f < z < %.1f   Websky %9d   Huntian %9d   Huntian/Websky %.3f", lo, hi, nw, nh, nh / nw))
    end
    say("\n== Fig. 2: v_rms^2 = <|v|^2> of the selected halos [(m/s)^2]")
    say(@sprintf("%-11s %14s %14s %8s", "z", "Websky", "Huntian", "H/W"))
    for i in 1:length(VZ)-1
        (W.nv[i] > 0 && H.nv[i] > 0) || continue
        w = 1e6 * W.v2[i] / W.nv[i]; h = 1e6 * H.v2[i] / H.nv[i]
        say(@sprintf("%.2f-%.2f  %14.4e %14.4e %8.3f", VZ[i], VZ[i+1], w, h, h / w))
    end
    yW = readmap(joinpath(WREF, "tsz_2048.fits")); yH = readmap(joinpath(HT, "maps", "tsz_y.fits"))
    aW = alm(yW); aH = alm(yH)
    pw = Healpix.pixwin(NSIDE); pw isa Tuple && (pw = pw[1]); pw = pw[1:LMAX+1]   # one pixel window, as the paper does for Websky
    edges = collect(0:DL:LMAX)
    say("\n== Fig. 3: tSZ × halo, D_ℓ = ℓ(ℓ+1)C_ℓ/2π ×1e6, one pixel window removed, Δℓ = $DL")
    for (k, (lo, hi)) in enumerate(ZBINS)
        cW = Healpix.alm2cl(aW, alm(delta(W.maps[k]))) ./ pw .^ 2
        cH = Healpix.alm2cl(aH, alm(delta(H.maps[k]))) ./ pw .^ 2
        say(@sprintf("\n-- %.1f < z < %.1f\n%-8s %12s %12s %8s", lo, hi, "ell", "Websky", "Huntian", "H/W"))
        for b in 1:length(edges)-1
            ls = max(edges[b], 2):edges[b+1]-1
            le = mean(ls); dw = mean(l * (l + 1) * cW[l+1] for l in ls) / 2π; dh = mean(l * (l + 1) * cH[l+1] for l in ls) / 2π
            say(@sprintf("%-8.0f %12.4e %12.4e %8.3f", le, 1e6 * dw, 1e6 * dh, dh / dw))
        end
    end
    close(io)
end
main()
