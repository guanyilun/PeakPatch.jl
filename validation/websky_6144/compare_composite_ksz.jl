#!/usr/bin/env julia
# Composite kSZ vs ksz.fits — same construction as compare_composite_kappa.jl, which
# settled the scheme question: Websky composites are B (full-matter field + halo
# profile compensated by a uniform Δ=3 sphere of the full painted mass, zero net).
# ksz.fits is z<4.5 websky only (no high-z tail issue; ksz_patchy.fits is a separate
# uncorrelated reionization model) — so ksz.fits is apples-to-apples from the start.
#
#   A (exclusion reading, for contrast): halo-excluded field kSZ + plain NFW halo kSZ.
#   B (Websky's construction):           full-matter field kSZ + Δ=3-compensated halo kSZ.
#
# Halo kSZ: ΔT = −T_CMB·(v_r/c)·τ_halo(θ), τ_halo = σ_T·n_e,0·x_e(z)·(1+z)²·Σ̃(θ)
# with Σ̃ = Σ/ρ̄_m [Mpc/h comoving] from the same truncated NFW (c=7, xmax=2) as κ —
# constants copied VERBATIM from src/FieldMap.jl (:tau kernel, Stein+2020 eq 3.22)
# so halo and field components share the normalization exactly. This is the first
# ABSOLUTE-units test of the pipeline (κ ratios were dimensionless; ksz.fits is μK).
#
# Halo v_r from the old-convention catalog (same reconstruction as the κ script /
# FieldMap comment): stored ψ₁ = Sbar·(D/a), stored ψ₂ carries +3/7 and (D/a)² →
# true d1 = stored·a, d2 = −stored·a²; x_E = q + d1 + d2;
# v = a_E·100·E(a_E)·f(a_E)·(d1 + 2·d2) km/s (finalize_eulerian, f₂≈2f).
#
# Usage: julia --project=validation -t N compare_composite_ksz.jl <ksz_field_excl.fits> <ksz_field_all.fits> [ref]
#        julia --project=validation compare_composite_ksz.jl --selftest
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
const OmB = 0.049
const rho_mh = 2.775e11 * 0.31          # Msun/h per (Mpc/h)^3
const OBS = -2618.0
const CNFW = 7.0; const XMAX = 2.0; const DCOMP = 3.0
const C_KMS = 299792.458
const TCMB_UK = 2.7255e6

# Thomson-depth constant, VERBATIM from FieldMap.jl (σ_T·n_e,0 per Mpc/h of ρ̄-column)
const F_E = 0.9; const Y_HE = 0.245
const SIGT_NE0 = 6.65246e-29 * 11.2299 * F_E * OmB * hub^2 * 3.0857e22 / hub
x_e_of_z(z) = z < 3 ? (1 - Y_HE / 2) : (1 - 3 * Y_HE / 4)   # He doubly/once ionized

f_nfw(c) = log(1 + c) - c / (1 + c)
const MTOT_FAC = 1 + CNFW^2 * (XMAX - 1) / ((1 + CNFW)^2 * f_nfw(CNFW))   # 1.636

# ---------- truncated-NFW projected profile (same as compare_composite_kappa.jl) ----------
_rho_dimless(u) = u < CNFW ? 1 / (u * (1 + u)^2) :
                  (u < XMAX * CNFW ? CNFW / ((1 + CNFW)^2 * u^2) : 0.0)
function _g_of_x(x)
    umax = XMAX * CNFW
    x >= umax && return 0.0
    lmax = sqrt(umax^2 - x^2); n = 512; h = lmax / n; acc = 0.0
    for i in 0:n-1
        l0 = i * h; lm = l0 + h / 2; l1 = l0 + h
        acc += h / 6 * (_rho_dimless(sqrt(l0^2 + x^2)) + 4 * _rho_dimless(sqrt(lm^2 + x^2)) +
                        _rho_dimless(sqrt(l1^2 + x^2)))
    end
    2acc
end
const GX_N = 2048
const GX_LOGX = range(log(1e-8), log(XMAX * CNFW), length=GX_N)
const GX_TAB = [_g_of_x(exp(lx)) for lx in GX_LOGX]
function g_of_x(x)
    lx = log(max(x, 1e-8))
    lx >= GX_LOGX[end] && return 0.0
    t = (lx - GX_LOGX[1]) / step(GX_LOGX) + 1
    i = clamp(floor(Int, t), 1, GX_N - 1)
    GX_TAB[i] * (1 - (t - i)) + GX_TAB[i+1] * (t - i)
end
const RHOS_OVER_RHOM = 200 * CNFW^3 / (3 * f_nfw(CNFW))

# τ of one halo at angle θ (comoving Mpc/h throughout).
# comp: 0 = plain truncated NFW; 1 = minus uniform Δ=3 sphere of the TOTAL painted
# mass 1.636·M200m (zero net electrons — Websky's construction).
@inline function tau_halo(θ, M, z, χ; comp::Int=0)
    r200 = cbrt(3 * M / (800π * rho_mh))
    rs = r200 / CNFW
    b = θ * χ
    Σ = b / rs < XMAX * CNFW ? RHOS_OVER_RHOM * rs * g_of_x(b / rs) : 0.0
    if comp > 0
        Rc = cbrt(3 * MTOT_FAC * M / (4π * DCOMP * rho_mh))
        b < Rc && (Σ -= DCOMP * 2 * sqrt(Rc^2 - b^2))
    end
    SIGT_NE0 * x_e_of_z(z) * (1 + z)^2 * Σ
end
# analytic ∫τ dΩ of one plain halo
tau_halo_integral(M, z, χ) =
    SIGT_NE0 * x_e_of_z(z) * (1 + z)^2 * MTOT_FAC * M / (rho_mh * χ^2)
function paint_radius(M; comp::Bool)
    r200 = cbrt(3 * M / (800π * rho_mh))
    fac = comp ? max(XMAX, cbrt(MTOT_FAC * 200 / DCOMP)) : XMAX
    return fac * r200
end

# self-tests: ∫τdΩ vs analytic (plain / compensated), plus an independent amplitude
# anchor — a 1e15 Msun/h cluster at z=0.5 has τ(1′) ~ 6e-3 (real-cluster scale, breaks
# circularity), and its kSZ at v_r=300 km/s must be ~10 μK.
let M = 3e14, z = 0.7
    cosmo0 = CosmologyParams(0.31, 0.049, 0.69, hub, 0.965, 0.81)
    χ = chi(z, cosmo0)
    for comp in (0, 1)
        θmax = paint_radius(M; comp=comp > 0) / χ
        n = 40000; h = θmax / n; acc = 0.0
        for i in 1:n
            θ = (i - 0.5) * h
            acc += tau_halo(θ, M, z, χ; comp=comp) * 2π * θ * h
        end
        want = (comp == 0 ? 1.0 : 0.0) * tau_halo_integral(M, z, χ)
        @printf("NFW τ self-test comp=%d: ∫τdΩ = %+.4e (expect %+.4e)\n", comp, acc, want)
        abs(acc - want) < 0.01 * tau_halo_integral(M, z, χ) || error("τ self-test failed")
    end
    χ5 = chi(0.5, cosmo0)
    τ1 = tau_halo(deg2rad(1 / 60), 1e15, 0.5, χ5; comp=0)
    ΔT = TCMB_UK * (300.0 / C_KMS) * τ1
    @printf("amplitude anchor: τ(1e15, z=0.5, 1′) = %.3e (expect 1e-3..3e-2); kSZ@300km/s = %.1f μK\n", τ1, ΔT)
    1e-3 < τ1 < 3e-2 || error("halo τ amplitude anchor failed — normalization wrong")
end
length(ARGS) >= 1 && ARGS[1] == "--selftest" && (println("self-tests passed"); exit(0))

# ---------- cosmology ----------
cosmo = CosmologyParams(0.31, 0.049, 0.69, hub, 0.965, 0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0)
gt = Dlinear_tables(cosmo)

# ---------- caps and flat-sky machinery (same as compare_composite_kappa.jl) ----------
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

# ---------- stream the 26 GB catalog ONCE: compact per-cap halo lists ----------
# Same as the κ script, plus per-halo radial velocity v_r [km/s] reconstructed from
# the stored 2LPT displacements (see header). Kept streaming/serial for the same
# login-node-memory and Julia-1.12 task-ownership reasons.
function select_cap_halos(path, chi2z, gt, axes, s)
    ncap = length(axes)
    E1 = [ortho_basis(ax)[1] for ax in axes]
    E2 = [ortho_basis(ax)[2] for ax in axes]
    out = [(gx=Float32[], gy=Float32[], r=Float32[], M=Float32[], vr=Float32[]) for _ in 1:ncap]
    nh_tot = open(path) do io Int(read(io, Int32)) end
    Om, OL = cosmo.Om, cosmo.OL
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
                # true displacements: d1 = stored·a, d2 = −stored·a² (stored ψ₂ carries +3/7)
                d11 = Float64(buf[b+4]) * a;  d12 = Float64(buf[b+5]) * a;  d13 = Float64(buf[b+6]) * a
                d21 = -Float64(buf[b+8]) * a^2; d22 = -Float64(buf[b+9]) * a^2; d23 = -Float64(buf[b+10]) * a^2
                vx = q1 + d11 + d21 - OBS
                vy = q2 + d12 + d22 - OBS
                vz = q3 + d13 + d23 - OBS
                r = sqrt(vx^2 + vy^2 + vz^2)
                (30.0 <= r <= 5170.0) || continue
                aE = 1.0 / (1.0 + chi_to_z(chi2z, r))
                fE = Dlinear_ab(aE, gt)[2]
                vfac = aE * 100.0 * sqrt(Om * aE^-3 + OL) * fE
                vr = vfac * ((d11 + 2d21) * vx + (d12 + 2d22) * vy + (d13 + 2d23) * vz) / r
                M = 4 / 3 * π * rho_mh * R^3
                gmax = paint_radius(M; comp=true) / r * 1.3
                for c in 1:ncap
                    ax = axes[c]
                    na = (vx*ax[1] + vy*ax[2] + vz*ax[3]) / r
                    na > 0.9 || continue
                    e1 = E1[c]; e2 = E2[c]
                    gxh = (vx*e1[1] + vy*e1[2] + vz*e1[3]) / (r * na)
                    gyh = (vx*e2[1] + vy*e2[2] + vz*e2[3]) / (r * na)
                    (abs(gxh) - gmax <= s && abs(gyh) - gmax <= s) || continue
                    push!(out[c].gx, Float32(gxh)); push!(out[c].gy, Float32(gyh))
                    push!(out[c].r, Float32(r));    push!(out[c].M, Float32(M))
                    push!(out[c].vr, Float32(vr))
                end
            end
            ndone += m
        end
    end
    return out
end

# paint one cap's halos: ΔT[μK] = −T_CMB·(v_r/c)·τ(θ); (plain, compensated) in one pass.
# Per-halo mass-exactness residual deposits as in the κ script (plain: analytic total;
# compensated: zero) — pks2map has no such correction, so at high ℓ we may sit ABOVE
# ksz.fits for the same reason we sit above kap_lt4.5.
function paint_halos_cap(h, chi2z, ax, s, npixf)
    e1, e2 = ortho_basis(ax)
    h0 = zeros(npixf, npixf); hc = zeros(npixf, npixf)
    dpix = 2s / npixf
    npaint = 0
    @inbounds for n in eachindex(h.M)
        gxh = Float64(h.gx[n]); gyh = Float64(h.gy[n])
        r = Float64(h.r[n]); M = Float64(h.M[n])
        kfac = -TCMB_UK * Float64(h.vr[n]) / C_KMS
        vx = ax[1] + gxh*e1[1] + gyh*e2[1]
        vy = ax[2] + gxh*e1[2] + gyh*e2[2]
        vz = ax[3] + gxh*e1[3] + gyh*e2[3]
        vn = sqrt(vx^2 + vy^2 + vz^2)
        θmaxc = paint_radius(M; comp=true) / r
        gmax = θmaxc * 1.3 * vn
        z = chi_to_z(chi2z, r)
        ilo = max(1, floor(Int, (gxh - gmax + s) / dpix) + 1)
        ihi = min(npixf, ceil(Int, (gxh + gmax + s) / dpix))
        jlo = max(1, floor(Int, (gyh - gmax + s) / dpix) + 1)
        jhi = min(npixf, ceil(Int, (gyh + gmax + s) / dpix))
        (ilo <= ihi && jlo <= jhi) || continue
        npaint += 1
        sum0 = 0.0; sumc = 0.0
        for j in jlo:jhi
            gy = -s + (j - 0.5) * dpix
            for i in ilo:ihi
                gx = -s + (i - 0.5) * dpix
                px = ax[1] + gx*e1[1] + gy*e2[1]
                py = ax[2] + gx*e1[2] + gy*e2[2]
                pz = ax[3] + gx*e1[3] + gy*e2[3]
                cosang = (px*vx + py*vy + pz*vz) /
                         (sqrt(px^2 + py^2 + pz^2) * vn)
                θ = acos(clamp(cosang, -1.0, 1.0))
                θ > θmaxc && continue
                w0 = kfac * tau_halo(θ, M, z, r; comp=0)
                wc = kfac * tau_halo(θ, M, z, r; comp=1)
                h0[i, j] += w0; hc[i, j] += wc
                sum0 += w0; sumc += wc
            end
        end
        icen = floor(Int, (gxh + s) / dpix) + 1
        jcen = floor(Int, (gyh + s) / dpix) + 1
        if 1 <= icen <= npixf && 1 <= jcen <= npixf &&
           gxh - gmax > -s && gxh + gmax < s && gyh - gmax > -s && gyh + gmax < s
            ki = kfac * tau_halo_integral(M, z, r) / dpix^2
            h0[icen, jcen] += ki - sum0
            hc[icen, jcen] += -sumc
        end
    end
    (h0, hc, npaint)
end

# ---------- run ----------
field_excl_path = ARGS[1]
field_all_path = ARGS[2]
ref_path = length(ARGS) >= 3 ? ARGS[3] : joinpath(WREF, "ksz.fits")

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
@info "extracting cap patches (maps read one at a time)..."
# our field maps store ΔT/T_CMB (dimensionless) — convert to μK; ksz.fits is already μK
PFE, ns_fe = patches_of(field_excl_path, axes, s; scale=TCMB_UK)
PFA, ns_fa = patches_of(field_all_path, axes, s; scale=TCMB_UK)
PK,  ns_k  = patches_of(ref_path, axes, s)
@info "reference map" ref_path
@info "map nsides" field_excl=ns_fe field_all=ns_fa ref=ns_k

w2_hp(l, nside) = exp(-l * (l + 1) * (sqrt(π / 3) / nside)^2 / 12)
w2_flat(l, dpix) = exp(-l * (l + 1) * dpix^2 / 12)

@info "streaming catalog (26 GB) for per-cap halo lists (+v_r)..." CAT
caphalos = select_cap_halos(CAT, chi2z, gt, axes, s)
@info "selected" nper=[length(h.M) for h in caphalos] vr_rms=[round(sqrt(mean(abs2, h.vr)); digits=1) for h in caphalos]

nb = length(lc)
cA = zeros(nb, NCAPS); cB = zeros(nb, NCAPS); cK = zeros(nb, NCAPS)
cFE = zeros(nb, NCAPS); cFA = zeros(nb, NCAPS); cH0 = zeros(nb, NCAPS); cHC = zeros(nb, NCAPS)
for (ic, ax) in enumerate(axes)
    pfe = PFE[ic]; pfa = PFA[ic]; pk_ = PK[ic]
    h0, hc, npaint = paint_halos_cap(caphalos[ic], chi2z, ax, s, NPIXFLAT)
    pA = pfe .+ h0                     # A: excluded field + plain NFW kSZ
    pB = pfa .+ hc                     # B: full field + Δ=3 comp (Websky construction)
    cA[:, ic] = cl_flat(pA, L, ledges);  cB[:, ic] = cl_flat(pB, L, ledges)
    cK[:, ic] = cl_flat(pk_, L, ledges)
    cFE[:, ic] = cl_flat(pfe, L, ledges); cFA[:, ic] = cl_flat(pfa, L, ledges)
    cH0[:, ic] = cl_flat(h0, L, ledges);  cHC[:, ic] = cl_flat(hc, L, ledges)
    @info "cap $ic done" npaint std_A=round(std(pA); digits=3) std_B=round(std(pB); digits=3) std_ref=round(std(pk_); digits=3)
end

@printf("\n%-7s %-11s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
        "ell", "Cl_ref[uK2]", "A/ref", "+-", "B/ref", "+-", "B/ref*", "fldA/ref", "halo0/ref")
dpixf = 2s / NPIXFLAT
for i in eachindex(lc)
    rA = [cA[i, c] / cK[i, c] for c in 1:NCAPS]
    rB = [cB[i, c] / cK[i, c] for c in 1:NCAPS]
    # window-corrected B/ref from the measured per-band decomposition
    wfa = w2_hp(lc[i], ns_fa); wh = w2_flat(lc[i], dpixf); wk = w2_hp(lc[i], ns_k)
    fldA = mean(cFA[i, :]); haloC = mean(cHC[i, :]); crossB = mean(cB[i, :]) - fldA - haloC
    Bstar = (fldA / wfa + haloC / wh + crossB / sqrt(wfa * wh)) / (mean(cK[i, :]) / wk)
    @printf("%-7.0f %-11.3e %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",
            lc[i], mean(cK[i, :]), mean(rA), std(rA), mean(rB), std(rB), Bstar,
            fldA / mean(cK[i, :]), mean(cH0[i, :]) / mean(cK[i, :]))
end
sel = findall(l -> 150 <= l <= 2500, lc)
mA = mean([mean(cA[i, :]) / mean(cK[i, :]) for i in sel])
mB = mean([mean(cB[i, :]) / mean(cK[i, :]) for i in sel])
@printf("\nband mean 150<=ell<=2500:  A/ref = %.3f   B/ref = %.3f\n", mA, mB)
@printf("Different realizations: cap scatter is the error bar. ksz.fits is z<4.5 websky\n")
@printf("(late-time only) — apples-to-apples with our z<4.6 octant, no tail correction.\n")
