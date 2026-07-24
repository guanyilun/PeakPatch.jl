#!/usr/bin/env julia
# Composite kSZ v2: Websky's ACTUAL halo model, replicated VERBATIM from the Fortran
# pks2map source (~/work/peakpatch/src/pks2map/{pks2map,bbps_profile,maptable,
# integrate_profiles,cosmology}.f90) — the code that made ksz.fits:
#
#   mass:    mh = sqrt(deltacrit(z)/200) * M_RTH   (SIS conversion, Bryan-Norman
#            deltacrit; M_RTH = 4π/3 ρ̄₀ RTH³ in Msun, NO h)
#   cut:     mh > 1e13 Msun (paper §4.4.2: smaller halos live in the 2LPT field)
#   profile: Battaglia AGN gas density (bbps_rhotilde, model 1 = AGN Δ=200):
#            ρ̃ = P0 (x/0.5)^-0.2 (1+(x/0.5)^αρ)^-β × (1+ΩΛ/(Ωm(1+z)³))  [ρ̄_m units]
#            P0 = 4e3 m14^0.29 zp^-0.66, αρ = 0.88 m14^-0.03 zp^0.19,
#            β = 3.83 m14^0.04 zp^-0.025;  x = r/rvir, rvir = R200c (comoving crit)
#   trunc:   spherical cut at x = 4 (integrate_profiles.f90: rmax=4, s0=sqrt(16-b²))
#   kSZ:     ΔT/T = -(v_r/c) · tau0persigma · (Ωm h²)^(2/3) (1+z)² (mh/200)^(1/3)
#                   · fb · Σ̃(b/rvir)                        (pks2map.f90:270)
#            tau0persigma = σT/(μe·mp) · (3Msun/4π)^(1/3) · rho100^(2/3), μe=1.136
#
# Halo map variants: W  = plain painting (what the Fortran does — NO compensation in
# the kSZ path, unlike κ's (gas+DM-1)); Wc = minus a uniform sphere of the same total
# painted τ, radius 4·rvir (zero net electrons per halo) — brackets the paper's §3.2.1
# compensation language. The halo term is small at low ℓ, so W vs Wc mostly matters
# for the cross term; ksz.fits decides.
#
# Field component: FULL-matter kSZ map (compensated construction), coarse_factor=32
# (velocity-coherence fix, KSZ_COMPOSITE_2026-07-19.md).
#
# Usage: julia --project=. -t N compare_composite_ksz_battaglia.jl <ksz_field_all.fits> [ref] [--selftest]
using Healpix, FFTW, Printf, Statistics, LinearAlgebra
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z, chi,
                            Dlinear_tables, Dlinear_ab

const WREF = "/home/yguan/scratch/websky_6144/websky_ref"
const Dr = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const CAT = joinpath(Dr, "catalog_websky_6144_oct000_finecell_AM.pksc")
const CAPDEG = 8.0
const NPIXFLAT = 2048
const NCAPS = 6
const hub = 0.68
const Om = 0.31; const OL = 0.69; const OmB = 0.049
const fb = OmB / Om
const rho_mh = 2.775e11 * Om            # Msun/h per (Mpc/h)^3 (our catalog units)
const OBS = -2618.0
const C_KMS = 299792.458
const TCMB_UK = 2.7255e6
const MMIN_MSUN = 1e13                  # Websky halo-map cut (M200c proxy, Msun)

# ---- Fortran cosmology.f90 constants (cgs) ----
const SIGMAT = 6.65246e-25
const MUE = 1.136
const MP_G = 1.67e-24
const MSUN_G = 1.9891e33
const RHO100 = 1.88e-29                 # g/cm^3, rho_crit for H=100
const TAU0PERSIGMA = SIGMAT / (MUE * MP_G) * (3 * MSUN_G / (4π))^(1 / 3) * RHO100^(2 / 3)

Ez2(z) = Om * (1 + z)^3 + OL
function deltacrit(z)                   # Bryan-Norman, crit-density units
    omz = Om * (1 + z)^3 / Ez2(z)
    x = omz - 1
    18π^2 + 82x - 39x^2
end
rhocrit_com(z) = 2.775e11 * hub^2 * Ez2(z) / (1 + z)^3    # Msun/Mpc^3 comoving

# Battaglia AGN gas density (bbps_rhotilde, model 1), in mean-matter units
@inline function rhotilde(x, m14, zp1, rhofac)
    P0 = 4.0e3 * m14^0.29 * zp1^-0.66
    αρ = 0.88 * m14^-0.03 * zp1^0.19
    β = 3.83 * m14^0.04 * zp1^-0.025
    P0 * (x / 0.5)^-0.2 * (1 + (x / 0.5)^αρ)^-β * rhofac
end

# ---- Σ̃(x_b; mh, z) table: projected profile in rvir units, spherical cut x<=4 ----
const XB_N = 160; const XB_LOG = range(log(1e-3), log(4.0), length=XB_N)
const MH_N = 60;  const MH_LOG = range(log(1e13), log(2e16), length=MH_N)
const ZT_N = 48;  const ZT_LOG = range(log(1.0), log(5.8), length=ZT_N)   # 1+z
const SIG_TAB = Array{Float64,3}(undef, XB_N, MH_N, ZT_N)
const ITOT_TAB = Array{Float64,2}(undef, MH_N, ZT_N)      # ∫Σ̃ 2πx dx (for comp/deposits)
function build_sigma_table!()
    Threads.@threads for k in 1:ZT_N
        zp1 = exp(ZT_LOG[k]); z = zp1 - 1
        rhofac = 1 + OL / Om / zp1^3
        for j in 1:MH_N
            m14 = exp(MH_LOG[j]) / 1e14
            for i in 1:XB_N
                b = exp(XB_LOG[i])
                s0 = sqrt(max(16.0 - b^2, 0.0))
                n = max(Int(ceil(s0 / 1e-2)), 1); ds = s0 / n
                acc = 0.0
                for is in 1:n
                    s = (is - 0.5) * ds
                    acc += rhotilde(sqrt(s^2 + b^2), m14, zp1, rhofac)
                end
                SIG_TAB[i, j, k] = 2 * acc * ds
            end
            it = 0.0
            for i in 1:XB_N-1
                b0 = exp(XB_LOG[i]); b1 = exp(XB_LOG[i+1])
                it += 0.5 * (SIG_TAB[i, j, k] * b0 + SIG_TAB[i+1, j, k] * b1) * (b1 - b0) * 2π
            end
            # add the inner disk (0..b_min, Σ̃≈Σ̃(b_min))
            ITOT_TAB[j, k] = it + SIG_TAB[1, j, k] * π * exp(XB_LOG[1])^2
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
function itot_interp(mh, z)
    ty = clamp((log(mh) - MH_LOG[1]) / step(MH_LOG) + 1, 1.0, MH_N - 1e-6)
    tz = clamp((log(1 + z) - ZT_LOG[1]) / step(ZT_LOG) + 1, 1.0, ZT_N - 1e-6)
    j = floor(Int, ty); k = floor(Int, tz)
    fy = ty - j; fz = tz - k
    (ITOT_TAB[j, k] * (1 - fy) + ITOT_TAB[j+1, k] * fy) * (1 - fz) +
    (ITOT_TAB[j, k+1] * (1 - fy) + ITOT_TAB[j+1, k+1] * fy) * fz
end

# per-halo kSZ amplitude (pks2map.f90:270): ΔT/T per unit Σ̃, × v_r/c separately
tau_amp(mh, z) = TAU0PERSIGMA * (Om * hub^2)^(2 / 3) * (1 + z)^2 * (mh / 200)^(1 / 3) * fb
rvir_com_mpc(mh, z) = cbrt(3 * mh / (4π * 200 * rhocrit_com(z)))   # Mpc comoving

@info "building Battaglia Σ̃ table ($(XB_N)x$(MH_N)x$(ZT_N))..."
build_sigma_table!()

# ---- self-tests ----
let
    # 1) central optical depth of a massive cluster: τ(b=0) ~ few e-3 (Battaglia 2016)
    mh = 3e14; z = 0.55
    τ0 = tau_amp(mh, z) * sigma_interp(1.2e-3, mh, z)
    @printf("τ(b≈0; 3e14 Msun, z=0.55) = %.3e  (expect ~5e-4..1.2e-2)\n", τ0)
    5e-4 < τ0 < 1.2e-2 || error("central-τ anchor failed")
    # 2) implied gas mass = ∫τ dΩ χ² μe mp / σT vs fb·mh — NOTE this is gas within
    #    the 4·rvir truncation sphere, and M(<4·R200c) ≈ 2.2×M200c (NFW c~5), so
    #    cluster-mass ratios land ~2-3; groups are feedback-depleted (<~1.3)
    for (mh2, lo, hi) in ((1e13, 0.3, 1.5), (1e15, 1.0, 3.5))
        rv = rvir_com_mpc(mh2, z)                                 # comoving Mpc
        tot = tau_amp(mh2, z) * itot_interp(mh2, z) * rv^2 / (1 + z)^2  # ∫τdΩ·χ², proper area
        mgas = tot * MUE * MP_G / SIGMAT * (3.0857e24)^2 / MSUN_G  # Msun
        ratio = mgas / (fb * mh2)
        @printf("gas-mass check M=%.0e: Mgas/(fb·M) = %.3f (expect %.1f-%.1f)\n", mh2, ratio, lo, hi)
        lo < ratio < hi || error("gas-mass sanity failed")
    end
    # 3) kSZ of a 1e14 halo at 300 km/s: ~1-20 μK central
    ΔT = TCMB_UK * (300 / C_KMS) * tau_amp(1e14, 0.5) * sigma_interp(1.2e-3, 1e14, 0.5)
    @printf("kSZ(1e14, z=0.5, b≈0, 300 km/s) = %.2f μK (expect 0.2-20)\n", ΔT)
    0.2 < ΔT < 20 || error("kSZ amplitude anchor failed")
end
length(ARGS) >= 1 && ARGS[end] == "--selftest" && (println("self-tests passed"); exit(0))

# ---------- cosmology / caps / spectra (same machinery as compare_composite_ksz.jl) ----------
cosmo = CosmologyParams(Om, OmB, OL, hub, 0.965, 0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0)
gt = Dlinear_tables(cosmo)

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

# ---- stream catalog: per-cap lists of (gx, gy, r, mh[Msun], vr) with the 1e13 cut ----
function select_cap_halos(path, chi2z, gt, axes, s)
    ncap = length(axes)
    E1 = [ortho_basis(ax)[1] for ax in axes]
    E2 = [ortho_basis(ax)[2] for ax in axes]
    out = [(gx=Float32[], gy=Float32[], r=Float32[], mh=Float32[], vr=Float32[]) for _ in 1:ncap]
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
                d11 = Float64(buf[b+4]) * a;  d12 = Float64(buf[b+5]) * a;  d13 = Float64(buf[b+6]) * a
                d21 = -Float64(buf[b+8]) * a^2; d22 = -Float64(buf[b+9]) * a^2; d23 = -Float64(buf[b+10]) * a^2
                vx = q1 + d11 + d21 - OBS
                vy = q2 + d12 + d22 - OBS
                vz = q3 + d13 + d23 - OBS
                r = sqrt(vx^2 + vy^2 + vz^2)
                (30.0 <= r <= 5170.0) || continue
                z = chi_to_z(chi2z, r)
                # Fortran mass pipeline: M_RTH[Msun] (our Msun/h ÷ h) → SIS M200c proxy
                M_RTH = 4 / 3 * π * rho_mh * R^3 / hub
                mh = sqrt(deltacrit(z) / 200) * M_RTH
                mh > MMIN_MSUN || continue
                aE = 1.0 / (1.0 + z)
                fE = Dlinear_ab(aE, gt)[2]
                vfac = aE * 100.0 * sqrt(Om * aE^-3 + OL) * fE
                vr = vfac * ((d11 + 2d21) * vx + (d12 + 2d22) * vy + (d13 + 2d23) * vz) / r
                for c in 1:ncap
                    ax = axes[c]
                    na = (vx*ax[1] + vy*ax[2] + vz*ax[3]) / r
                    na > 0.9 || continue
                    # paint window: 4·asin(rh/χ), rh = mean-density-200 radius (Fortran)
                    rh = cbrt(3 * mh / (4π * 200 * 2.775e11 * hub^2 * Om)) * hub  # Mpc/h
                    gmax = 4 * rh / r * 1.3
                    e1 = E1[c]; e2 = E2[c]
                    gxh = (vx*e1[1] + vy*e1[2] + vz*e1[3]) / (r * na)
                    gyh = (vx*e2[1] + vy*e2[2] + vz*e2[3]) / (r * na)
                    (abs(gxh) - gmax <= s && abs(gyh) - gmax <= s) || continue
                    push!(out[c].gx, Float32(gxh)); push!(out[c].gy, Float32(gyh))
                    push!(out[c].r, Float32(r));    push!(out[c].mh, Float32(mh))
                    push!(out[c].vr, Float32(vr))
                end
            end
            ndone += m
        end
    end
    return out
end

# paint one cap: W (plain, Fortran-literal) and Wc (uniform-sphere compensated, R=4rvir)
function paint_halos_cap(h, chi2z, ax, s, npixf)
    e1, e2 = ortho_basis(ax)
    hw = zeros(npixf, npixf); hc = zeros(npixf, npixf)
    dpix = 2s / npixf
    npaint = 0
    @inbounds for n in eachindex(h.mh)
        gxh = Float64(h.gx[n]); gyh = Float64(h.gy[n])
        r = Float64(h.r[n]); mh = Float64(h.mh[n])
        z = chi_to_z(chi2z, r)
        kfac = -TCMB_UK * Float64(h.vr[n]) / C_KMS
        amp = tau_amp(mh, z)
        rvir_h = rvir_com_mpc(mh, z) * hub          # Mpc/h comoving (catalog units)
        θv = rvir_h / r
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
        itot = itot_interp(mh, z)
        ucomp = itot * 3 / (256π)                   # uniform sphere: Σc = ucomp·2√(16-xb²)
        sumw = 0.0; sumc = 0.0
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
                xb = θ / θv
                Σ = sigma_interp(xb, mh, z)
                w = kfac * amp * Σ
                wc = kfac * amp * (Σ - ucomp * 2 * sqrt(max(16.0 - xb^2, 0.0)))
                hw[i, j] += w; hc[i, j] += wc
                sumw += w; sumc += wc
            end
        end
        # per-halo exactness deposits (plain: analytic total; comp: zero)
        icen = floor(Int, (gxh + s) / dpix) + 1
        jcen = floor(Int, (gyh + s) / dpix) + 1
        if 1 <= icen <= npixf && 1 <= jcen <= npixf &&
           gxh - gmax > -s && gxh + gmax < s && gyh - gmax > -s && gyh + gmax < s
            ki = kfac * amp * itot * θv^2 / dpix^2
            hw[icen, jcen] += ki - sumw
            hc[icen, jcen] += -sumc
        end
    end
    (hw, hc, npaint)
end

# ---------- run ----------
field_all_path = ARGS[1]
ref_path = length(ARGS) >= 2 && !startswith(ARGS[2], "--") ? ARGS[2] : joinpath(WREF, "ksz.fits")

s = tan(deg2rad(CAPDEG)) / sqrt(2)
L = 2s
ledges = [10.0^l for l in range(log10(100.0), log10(4000.0); length=12)]
lc = [sqrt(ledges[i] * ledges[i+1]) for i in 1:length(ledges)-1]
axes = cap_axes()

function patches_of(path, axes, s; scale=1.0)
    m = Healpix.readMapFromFITS(path, 1, Float32)
    ns = m.resolution.nside
    p = [flat_patch(m, ax, s, NPIXFLAT) .* scale for ax in axes]
    m = nothing; GC.gc()
    p, ns
end
@info "extracting cap patches..."
PFA, ns_fa = patches_of(field_all_path, axes, s; scale=TCMB_UK)
PK,  ns_k  = patches_of(ref_path, axes, s)
@info "maps" field=field_all_path ref=ref_path ns_field=ns_fa ns_ref=ns_k

@info "streaming catalog (mh>1e13 Msun cut)..." CAT
caphalos = select_cap_halos(CAT, chi2z, gt, axes, s)
@info "selected" nper=[length(h.mh) for h in caphalos]

nb = length(lc)
cW = zeros(nb, NCAPS); cWc = zeros(nb, NCAPS); cK = zeros(nb, NCAPS)
cFA = zeros(nb, NCAPS); cH = zeros(nb, NCAPS)
for (ic, ax) in enumerate(axes)
    hw, hc, npaint = paint_halos_cap(caphalos[ic], chi2z, ax, s, NPIXFLAT)
    pW = PFA[ic] .+ hw
    pWc = PFA[ic] .+ hc
    cW[:, ic] = cl_flat(pW, L, ledges); cWc[:, ic] = cl_flat(pWc, L, ledges)
    cK[:, ic] = cl_flat(PK[ic], L, ledges)
    cFA[:, ic] = cl_flat(PFA[ic], L, ledges); cH[:, ic] = cl_flat(hw, L, ledges)
    @info "cap $ic done" npaint std_W=round(std(pW); digits=3) std_ref=round(std(PK[ic]); digits=3)
end

@printf("\n%-7s %-11s %-8s %-8s %-8s %-8s %-9s %-9s\n",
        "ell", "Cl_ref[uK2]", "W/ref", "+-", "Wc/ref", "+-", "fld/ref", "halo/ref")
for i in eachindex(lc)
    rW = [cW[i, c] / cK[i, c] for c in 1:NCAPS]
    rWc = [cWc[i, c] / cK[i, c] for c in 1:NCAPS]
    @printf("%-7.0f %-11.3e %-8.3f %-8.3f %-8.3f %-8.3f %-9.3f %-9.3f\n",
            lc[i], mean(cK[i, :]), mean(rW), std(rW), mean(rWc), std(rWc),
            mean(cFA[i, :]) / mean(cK[i, :]), mean(cH[i, :]) / mean(cK[i, :]))
end
sel = findall(l -> 150 <= l <= 2500, lc)
mW = mean([mean(cW[i, :]) / mean(cK[i, :]) for i in sel])
mWc = mean([mean(cWc[i, :]) / mean(cK[i, :]) for i in sel])
@printf("\nband mean 150<=ell<=2500:  W/ref = %.3f   Wc/ref = %.3f\n", mW, mWc)
@printf("W = Fortran-literal (uncompensated) Battaglia halos; Wc = zero-net compensated.\n")
