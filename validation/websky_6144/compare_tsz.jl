#!/usr/bin/env julia
# tSZ (Compton-y) halo painting vs tsz_2048.fits — PURE catalog product (no field
# component), replicated VERBATIM from Fortran pks2map (the code that made the map):
#
#   mass:    mh = sqrt(deltacrit(z)/200) * M_RTH   (SIS M200c proxy, Bryan-Norman)
#   profile: Battaglia 2012 AGN Δ=200 pressure (bbps_ptilde, model 1):
#            P̃ = P0 (x/xc)^-0.3 (1+x/xc)^-β, α=1
#            P0 = 18.1  m14^0.154   (1+z)^-0.758
#            xc = 0.497 m14^-0.00865 (1+z)^0.731
#            β  = 4.35  m14^0.0393  (1+z)^0.415ъ    x = r/rvir(200c comoving)
#   trunc:   spherical x <= 4
#   y(θ):    y0 · mh · E²(z) · Σ̃P(b/rvir)          (maptable.f90)
#            y0 = fb·200·h²·ypermass/2, ypermass = 3σT·H100²/(8π·me·c²)·Msun/1.932
#
# This independently validates the Battaglia machinery (mass conversion, rvir, table
# projection, amplitude constants) used by the kSZ halo painter. No mass cut (the
# released tSZ map is painted from the full catalog). Caveat: released tsz maps were
# REPLACED 2022 (7% normalization + high-mass center fixes) — the local Fortran clone
# postdates those fixes, so ratios should be ~1 if our catalog+painter are right.
#
# Usage: julia --project=validation -t N compare_tsz.jl [ref] [--selftest]
using Healpix, FFTW, Printf, Statistics, LinearAlgebra
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z

const WREF = "/home/yguan/scratch/websky_6144/websky_ref"
const Dr = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const CAT = joinpath(Dr, "catalog_websky_6144_oct000_finecell_AM.pksc")
const CAPDEG = 8.0
const NPIXFLAT = 2048
const NCAPS = 6
const hub = 0.68
const Om = 0.31; const OL = 0.69; const OmB = 0.049
const fb = OmB / Om
const rho_mh = 2.775e11 * Om
const OBS = -2618.0

# Fortran constants (cgs)
const SIGMAT = 6.65246e-25
const ME_G = 9.109383e-28
const C_CMS = 2.99792458e10
const MSUN_G = 1.9891e33
const H100_S = 1e2 / 3.08567758e19
const YPERMASS = 3 * SIGMAT * H100_S^2 / (8π * ME_G * C_CMS^2) * MSUN_G / 1.932
const Y0 = fb * 200 * hub^2 * YPERMASS / 2      # per Msun (maptable.f90)

Ez2(z) = Om * (1 + z)^3 + OL
function deltacrit(z)
    omz = Om * (1 + z)^3 / Ez2(z)
    x = omz - 1
    18π^2 + 82x - 39x^2
end
rhocrit_com(z) = 2.775e11 * hub^2 * Ez2(z) / (1 + z)^3
rvir_com_mpc(mh, z) = cbrt(3 * mh / (4π * 200 * rhocrit_com(z)))

# Battaglia 2012 AGN pressure (bbps_ptilde, model 1)
@inline function ptilde(x, m14, zp1)
    P0 = 18.100 * m14^0.154 * zp1^-0.758
    xc = 0.497 * m14^-0.00865 * zp1^0.731
    β = 4.35 * m14^0.0393 * zp1^0.415
    P0 * (x / xc)^-0.3 * (1 + (x / xc))^-β
end

# Σ̃P table (projected pressure, rvir units, spherical x<=4)
const XB_N = 160; const XB_LOG = range(log(1e-3), log(4.0), length=XB_N)
const MH_N = 80;  const MH_LOG = range(log(1e11), log(2e16), length=MH_N)
const ZT_N = 48;  const ZT_LOG = range(log(1.0), log(5.8), length=ZT_N)
const SIG_TAB = Array{Float64,3}(undef, XB_N, MH_N, ZT_N)
function build_sigma_table!()
    Threads.@threads for k in 1:ZT_N
        zp1 = exp(ZT_LOG[k])
        for j in 1:MH_N
            m14 = exp(MH_LOG[j]) / 1e14
            for i in 1:XB_N
                b = exp(XB_LOG[i])
                s0 = sqrt(max(16.0 - b^2, 0.0))
                n = max(Int(ceil(s0 / 1e-2)), 1); ds = s0 / n
                acc = 0.0
                for is in 1:n
                    s = (is - 0.5) * ds
                    acc += ptilde(sqrt(s^2 + b^2), m14, zp1)
                end
                SIG_TAB[i, j, k] = 2 * acc * ds
            end
        end
    end
end
@inline function _lerp3(A, i, j, k, fx, fy, fz)
    c00 = A[i, j, k] * (1 - fx) + A[i+1, j, k] * fx
    c10 = A[i, j+1, k] * (1 - fx) + A[i+1, j+1, k] * fx
    c01 = A[i, j, k+1] * (1 - fx) + A[i+1, j, k+1] * fx
    c11 = A[i, j+1, k+1] * (1 - fx) + A[i+1, j+1, k+1] * fx
    (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz
end
function sigma_interp(xb, mh, z)
    xb >= 4.0 && return 0.0
    tx = clamp((log(max(xb, 1.1e-3)) - XB_LOG[1]) / step(XB_LOG) + 1, 1.0, XB_N - 1e-6)
    ty = clamp((log(mh) - MH_LOG[1]) / step(MH_LOG) + 1, 1.0, MH_N - 1e-6)
    tz = clamp((log(1 + z) - ZT_LOG[1]) / step(ZT_LOG) + 1, 1.0, ZT_N - 1e-6)
    i = floor(Int, tx); j = floor(Int, ty); k = floor(Int, tz)
    _lerp3(SIG_TAB, i, j, k, tx - i, ty - j, tz - k)
end
y_of(mh, z, xb) = Y0 * mh * Ez2(z) * sigma_interp(xb, mh, z)

@info "building Battaglia pressure Σ̃ table ($(XB_N)x$(MH_N)x$(ZT_N))..."
build_sigma_table!()

# self-tests: central y of a massive cluster ~1e-4; y(1e14) ~ few e-6
let
    y15 = y_of(1e15, 0.2, 1.2e-3)
    y14 = y_of(1e14, 0.5, 1.2e-3)
    @printf("y(b≈0; 1e15 Msun, z=0.2) = %.3e (expect 2e-5..5e-4)\n", y15)
    @printf("y(b≈0; 1e14 Msun, z=0.5) = %.3e (expect 5e-7..5e-5)\n", y14)
    (2e-5 < y15 < 5e-4 && 5e-7 < y14 < 5e-5) || error("y amplitude anchor failed")
end
length(ARGS) >= 1 && ARGS[end] == "--selftest" && (println("self-tests passed"); exit(0))

# ---------- caps machinery ----------
cosmo = CosmologyParams(Om, OmB, OL, hub, 0.965, 0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0)
ortho_basis(a) = (t = abs(a[1]) < 0.9 ? [1.0,0,0] : [0.0,1,0]; e1 = normalize(cross(a,t)); e2 = cross(a,e1); (e1,e2))
function cap_axes()
    diag = [1.0, 1.0, 1.0] ./ sqrt(3.0)
    e1, e2 = ortho_basis(diag)
    axes = Vector{Vector{Float64}}()
    push!(axes, diag)
    for k in 0:NCAPS-2
        ϕ = 2π * k / (NCAPS - 1)
        push!(axes, normalize(cos(deg2rad(16)) .* diag .+ sin(deg2rad(16)) .* (cos(ϕ) .* e1 .+ sin(ϕ) .* e2)))
    end
    axes
end
function flat_patch(m::HealpixMap, ax, s, npixf)
    e1, e2 = ortho_basis(ax)
    res = m.resolution
    out = zeros(npixf, npixf)
    Threads.@threads for j in 1:npixf
        for i in 1:npixf
            gx = (-s + (i - 0.5) * 2s / npixf); gy = (-s + (j - 0.5) * 2s / npixf)
            out[i, j] = m.pixels[Healpix.vec2pixRing(res, ax[1] + gx*e1[1] + gy*e2[1],
                                                     ax[2] + gx*e1[2] + gy*e2[2],
                                                     ax[3] + gx*e1[3] + gy*e2[3])]
        end
    end
    out
end
function cl_flat(A, L, ledges)
    m = A .- mean(A)
    F = fft(m); n = size(A, 1)
    f = fftfreq(n, n / L) .* 2π
    nb = length(ledges) - 1; S = zeros(nb); Nm = zeros(nb)
    for j in 1:n, i in 1:n
        l = sqrt(f[i]^2 + f[j]^2); (ledges[1] <= l < ledges[end]) || continue
        b = searchsortedlast(ledges, l); S[b] += abs2(F[i, j]); Nm[b] += 1
    end
    [Nm[b] > 0 ? S[b] / Nm[b] * L^2 / n^4 : 0.0 for b in 1:nb]
end

# stream catalog: per-cap (gx, gy, r, mh) — Eulerian positions, no velocity needed
function select_cap_halos(path, chi2z, axes, s)
    ncap = length(axes)
    E1 = [ortho_basis(ax)[1] for ax in axes]
    E2 = [ortho_basis(ax)[2] for ax in axes]
    out = [(gx=Float32[], gy=Float32[], r=Float32[], mh=Float32[]) for _ in 1:ncap]
    nh_tot = open(path) do io Int(read(io, Int32)) end
    open(path) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        nf = 33; chunk = 1_000_000; buf = Vector{Float32}(undef, chunk * nf); ndone = 0
        while ndone < nh_tot
            m = min(chunk, nh_tot - ndone); read!(io, view(buf, 1:m*nf))
            @inbounds for k in 1:m
                b = (k - 1) * nf
                R = Float64(buf[b+7])
                q1 = Float64(buf[b+1]); q2 = Float64(buf[b+2]); q3 = Float64(buf[b+3])
                rq = sqrt((q1 - OBS)^2 + (q2 - OBS)^2 + (q3 - OBS)^2)
                a = 1.0 / (1.0 + chi_to_z(chi2z, rq))
                vx = q1 + Float64(buf[b+4]) * a - Float64(buf[b+8]) * a^2 - OBS
                vy = q2 + Float64(buf[b+5]) * a - Float64(buf[b+9]) * a^2 - OBS
                vz = q3 + Float64(buf[b+6]) * a - Float64(buf[b+10]) * a^2 - OBS
                r = sqrt(vx^2 + vy^2 + vz^2)
                (30.0 <= r <= 5170.0) || continue
                z = chi_to_z(chi2z, r)
                M_RTH = 4 / 3 * π * rho_mh * R^3 / hub          # Msun
                mh = sqrt(deltacrit(z) / 200) * M_RTH
                rh = cbrt(3 * mh / (4π * 200 * 2.775e11 * hub^2 * Om)) * hub  # Mpc/h
                gmax = 4 * rh / r * 1.3
                for c in 1:ncap
                    ax = axes[c]
                    na = (vx*ax[1] + vy*ax[2] + vz*ax[3]) / r
                    na > 0.9 || continue
                    e1 = E1[c]; e2 = E2[c]
                    gxh = (vx*e1[1] + vy*e1[2] + vz*e1[3]) / (r * na)
                    gyh = (vx*e2[1] + vy*e2[2] + vz*e2[3]) / (r * na)
                    (abs(gxh) - gmax <= s && abs(gyh) - gmax <= s) || continue
                    push!(out[c].gx, Float32(gxh)); push!(out[c].gy, Float32(gyh))
                    push!(out[c].r, Float32(r));    push!(out[c].mh, Float32(mh))
                end
            end
            ndone += m
        end
    end
    return out
end

function paint_halos_cap(h, chi2z, ax, s, npixf)
    e1, e2 = ortho_basis(ax)
    hy = zeros(npixf, npixf)
    dpix = 2s / npixf
    npaint = 0
    @inbounds for n in eachindex(h.mh)
        gxh = Float64(h.gx[n]); gyh = Float64(h.gy[n])
        r = Float64(h.r[n]); mh = Float64(h.mh[n])
        z = chi_to_z(chi2z, r)
        amp = Y0 * mh * Ez2(z)
        θv = rvir_com_mpc(mh, z) * hub / r
        θmax = 4 * θv
        vx = ax[1] + gxh*e1[1] + gyh*e2[1]
        vy = ax[2] + gxh*e1[2] + gyh*e2[2]
        vz = ax[3] + gxh*e1[3] + gyh*e2[3]
        vn = sqrt(vx^2 + vy^2 + vz^2)
        gmax = θmax * 1.3 * vn
        ilo = max(1, floor(Int, (gxh - gmax + s) / dpix) + 1)
        ihi = min(npixf, ceil(Int, (gxh + gmax + s) / dpix))
        jlo = max(1, floor(Int, (gyh - gmax + s) / dpix) + 1)
        jhi = min(npixf, ceil(Int, (gyh + gmax + s) / dpix))
        (ilo <= ihi && jlo <= jhi) || continue
        npaint += 1
        for j in jlo:jhi
            gy = -s + (j - 0.5) * dpix
            for i in ilo:ihi
                gx = -s + (i - 0.5) * dpix
                px = ax[1] + gx*e1[1] + gy*e2[1]
                py = ax[2] + gx*e1[2] + gy*e2[2]
                pz = ax[3] + gx*e1[3] + gy*e2[3]
                cosang = (px*vx + py*vy + pz*vz) / (sqrt(px^2 + py^2 + pz^2) * vn)
                θ = acos(clamp(cosang, -1.0, 1.0))
                θ > θmax && continue
                hy[i, j] += amp * sigma_interp(θ / θv, mh, z)
            end
        end
    end
    (hy, npaint)
end

# ---------- run ----------
ref_path = length(ARGS) >= 1 && !startswith(ARGS[1], "--") ? ARGS[1] : joinpath(WREF, "tsz_2048.fits")
s = tan(deg2rad(CAPDEG)) / sqrt(2)
L = 2s
ledges = [10.0^l for l in range(log10(100.0), log10(4000.0); length=12)]
lc = [sqrt(ledges[i] * ledges[i+1]) for i in 1:length(ledges)-1]
axes = cap_axes()

@info "extracting reference cap patches..." ref_path
refmap = Healpix.readMapFromFITS(ref_path, 1, Float32)
PK = [flat_patch(refmap, ax, s, NPIXFLAT) for ax in axes]
refmap = nothing; GC.gc()

@info "streaming catalog (all halos, no cut)..." CAT
caphalos = select_cap_halos(CAT, chi2z, axes, s)
@info "selected" nper=[length(h.mh) for h in caphalos]

nb = length(lc)
cY = zeros(nb, NCAPS); cK = zeros(nb, NCAPS)
ymeans = zeros(NCAPS); yrefmeans = zeros(NCAPS)
for (ic, ax) in enumerate(axes)
    hy, npaint = paint_halos_cap(caphalos[ic], chi2z, ax, s, NPIXFLAT)
    cY[:, ic] = cl_flat(hy, L, ledges)
    cK[:, ic] = cl_flat(PK[ic], L, ledges)
    ymeans[ic] = mean(hy); yrefmeans[ic] = mean(PK[ic])
    @info "cap $ic done" npaint mean_y=round(mean(hy); sigdigits=3) mean_ref=round(mean(PK[ic]); sigdigits=3)
end

@printf("\nmean y: ours %.3e  ref %.3e  ratio %.3f\n",
        mean(ymeans), mean(yrefmeans), mean(ymeans) / mean(yrefmeans))
@printf("%-7s %-11s %-8s %-8s\n", "ell", "Cl_ref", "ours/ref", "+-")
for i in eachindex(lc)
    rY = [cY[i, c] / cK[i, c] for c in 1:NCAPS]
    @printf("%-7.0f %-11.3e %-8.3f %-8.3f\n", lc[i], mean(cK[i, :]), mean(rY), std(rY))
end
sel = findall(l -> 150 <= l <= 2500, lc)
@printf("\nband mean 150<=ell<=2500: ours/ref = %.3f\n",
        mean([mean(cY[i, :]) / mean(cK[i, :]) for i in sel]))
@printf("Caveats: different octant/realization (cap scatter = error); tSZ C_l is\n")
@printf("shot/1-halo dominated -> more realization-sensitive than kSZ; mean-y is the\n")
@printf("robust absolute-normalization check. Ref pixel window (Nside 2048) uncorrected.\n")
