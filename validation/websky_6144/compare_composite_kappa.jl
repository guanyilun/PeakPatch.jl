#!/usr/bin/env julia
# Phase-C composite test (docs/field_lightcone_plan.md): assemble full κ maps from our
# field + halo components and compare BOTH double-counting schemes against kap.fits:
#
#   A (paper §3.1.3 literal): field map with halo Lagrangian spheres EXCLUDED
#                             + plain truncated-NFW halo κ.
#   B (paper §3.2.4 literal): FULL-matter field map
#                             + NFW halo κ compensated by a uniform Δ=3 sphere of the
#                               same total mass (halo map adds zero net mass).
#
# Both conserve mass by construction; they differ in how the collapsed-region power is
# split between components. kap.fits decides empirically which reproduces Websky.
#
# Usage: julia --project=validation -t 8 compare_composite_kappa.jl <kappa_field_excl.fits> <kappa_field_all.fits>
#
# The halo κ is painted here standalone (same truncated-NFW c=7/xmax=2 as the XGPaint
# fork, validated in-script against the analytic mass integral) on the SAME gnomonic
# pixel grid as the field patches, with exact per-pixel angles — so components add
# pixel-wise with no projection mismatch. Everything in Mpc/h and Msun/h.
using Healpix, FFTW, Printf, Statistics, LinearAlgebra, DelimitedFiles
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z, chi

const WREF = "/home/yguan/scratch/websky_6144/websky_ref"
const Dr = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const CAT = joinpath(Dr, "catalog_websky_6144_oct000_finecell_AM.pksc")
const CAPDEG = 8.0
const NPIXFLAT = 2048
const NCAPS = 6
const hub = 0.68
const rho_mh = 2.775e11 * 0.31          # Msun/h per (Mpc/h)^3
const chistar = 14200.0 * hub           # Websky hardwired CMB source plane, Mpc/h
const OBS = -2618.0
const CNFW = 7.0; const XMAX = 2.0; const DCOMP = 3.0
const H0C = 1.0 / 2997.92458            # H0/c in h/Mpc

f_nfw(c) = log(1 + c) - c / (1 + c)
const MTOT_FAC = 1 + CNFW^2 * (XMAX - 1) / ((1 + CNFW)^2 * f_nfw(CNFW))   # 1.636

# ---------- truncated-NFW projected profile: g(x) table (Simpson), x = R⊥/r_s ----------
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

# κ of one halo at angle θ: W_κ(χ)·(Σ/ρ̄)(θχ), comoving Mpc/h throughout.
# comp=true subtracts the uniform Δ=3 sphere of the same TOTAL painted mass.
# Kernel W_κ = (3/2)Ω_m(H0/c)²(1+z)·χ·(1−χ/χ*) — the χ MULTIPLIES (Born kernel);
# an earlier version divided, suppressing halos by χ² ≈ 10⁷ (and its self-test
# used the same wrong W, so it "passed" — hence the independent aggregate check
# against ∫W f_coll dχ below).
# comp: 0 = plain NFW; 1 = subtract sphere of the TOTAL painted mass 1.636·M200m
# (zero net mass); 2 = paper-literal "same mass as the halo" = M200m (net +0.636·M).
@inline function kappa_halo(θ, M, z, χ; comp::Int=0)
    r200 = cbrt(3 * M / (800π * rho_mh))
    rs = r200 / CNFW
    b = θ * χ
    Σ = b / rs < XMAX * CNFW ? RHOS_OVER_RHOM * rs * g_of_x(b / rs) : 0.0
    if comp > 0
        Msph = comp == 1 ? MTOT_FAC * M : M
        Rc = cbrt(3 * Msph / (4π * DCOMP * rho_mh))
        b < Rc && (Σ -= DCOMP * 2 * sqrt(Rc^2 - b^2))
    end
    1.5 * 0.31 * H0C^2 * (1 + z) * χ * (1 - χ / chistar) * Σ
end

# analytic solid-angle integral of one plain halo: ∫κ dΩ = W_κ·M_tot/(ρ̄χ²)
kappa_halo_integral(M, z, χ) =
    1.5 * 0.31 * H0C^2 * (1 + z) * χ * (1 - χ / chistar) * MTOT_FAC * M / (rho_mh * χ^2)
function paint_radius(M; comp::Bool)
    r200 = cbrt(3 * M / (800π * rho_mh))
    fac = comp ? max(XMAX, cbrt(MTOT_FAC * 200 / DCOMP)) : XMAX   # comp sphere: 4.78·r200
    return fac * r200
end

# self-test: ∫κ 2πθ dθ == kappa_halo_integral (plain) and ≈0 (compensated), PLUS an
# independent amplitude anchor: κ of an M=1e15 cluster at z=0.5, θ=1′ must be O(0.1-1)
# (breaks the circularity of testing the kernel against itself).
let M = 3e14, z = 0.7
    cosmo0 = CosmologyParams(0.31, 0.049, 0.69, hub, 0.965, 0.81)
    χ = chi(z, cosmo0)
    for comp in (0, 1, 2)
        θmax = paint_radius(M; comp=comp > 0) / χ
        n = 40000; h = θmax / n; acc = 0.0
        for i in 1:n
            θ = (i - 0.5) * h
            acc += kappa_halo(θ, M, z, χ; comp=comp) * 2π * θ * h
        end
        want = (comp == 0 ? 1.0 : comp == 1 ? 0.0 : 1 - 1 / MTOT_FAC) * kappa_halo_integral(M, z, χ)
        @printf("NFW self-test comp=%d: ∫κdΩ = %+.4e (expect %+.4e)\n", comp, acc, want)
        abs(acc - want) < 0.01 * kappa_halo_integral(M, z, χ) || error("NFW self-test failed")
    end
    χ5 = chi(0.5, cosmo0)
    κ1 = kappa_halo(deg2rad(1 / 60), 1e15, 0.5, χ5; comp=0)
    @printf("amplitude anchor: κ(1e15 Msun/h, z=0.5, θ=1′) = %.3f (expect 0.05-2)\n", κ1)
    0.05 < κ1 < 2 || error("halo κ amplitude anchor failed — kernel wrong")
end

# ---------- cosmology ----------
cosmo = CosmologyParams(0.31, 0.049, 0.69, hub, 0.965, 0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0)

# ---------- caps and flat-sky machinery (same as compare_fieldmap_kappa.jl) ----------
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
# The login node silently kills processes between 4 and 6 GB RSS, so we never hold
# the full catalog: only in-cap halos are kept as (gnomonic gx, gy, distance r, M)
# — 16 B each, ~5M per cap. Selection is serial: the pass is disk-bound, and pushing
# to main-task vectors from threads is a Julia 1.12 ConcurrencyViolationError anyway.
function select_cap_halos(path, chi2z, axes, s)
    ncap = length(axes)
    E1 = [ortho_basis(ax)[1] for ax in axes]
    E2 = [ortho_basis(ax)[2] for ax in axes]
    out = [(gx=Float32[], gy=Float32[], r=Float32[], M=Float32[]) for _ in 1:ncap]
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
                # old-convention catalog: Eulerian = q + d1·a − d2·a² (stored ψ₂ carries +3/7)
                vx = q1 + Float64(buf[b+4]) * a - Float64(buf[b+8]) * a^2 - OBS
                vy = q2 + Float64(buf[b+5]) * a - Float64(buf[b+9]) * a^2 - OBS
                vz = q3 + Float64(buf[b+6]) * a - Float64(buf[b+10]) * a^2 - OBS
                r = sqrt(vx^2 + vy^2 + vz^2)
                (30.0 <= r <= 5170.0) || continue
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
                end
            end
            ndone += m
        end
    end
    return out
end

# paint one cap's (pre-selected) halos on the SAME gnomonic grid, exact per-pixel
# angles. Halo direction reconstructed exactly from its gnomonic coords.
# Returns (plain, compensated) patches in one pass.
function paint_halos_cap(h, chi2z, ax, s, npixf)
    e1, e2 = ortho_basis(ax)
    h0 = zeros(npixf, npixf); hc = zeros(npixf, npixf); hp = zeros(npixf, npixf)
    dpix = 2s / npixf
    npaint = 0
    @inbounds for n in eachindex(h.M)
        gxh = Float64(h.gx[n]); gyh = Float64(h.gy[n])
        r = Float64(h.r[n]); M = Float64(h.M[n])
        vx = ax[1] + gxh*e1[1] + gyh*e2[1]
        vy = ax[2] + gxh*e1[2] + gyh*e2[2]
        vz = ax[3] + gxh*e1[3] + gyh*e2[3]
        vn = sqrt(vx^2 + vy^2 + vz^2)
        θmaxc = paint_radius(M; comp=true) / r
        gmax = θmaxc * 1.3 * vn           # conservative gnomonic-plane radius
        z = chi_to_z(chi2z, r)
        ilo = max(1, floor(Int, (gxh - gmax + s) / dpix) + 1)
        ihi = min(npixf, ceil(Int, (gxh + gmax + s) / dpix))
        jlo = max(1, floor(Int, (gyh - gmax + s) / dpix) + 1)
        jhi = min(npixf, ceil(Int, (gyh + gmax + s) / dpix))
        (ilo <= ihi && jlo <= jhi) || continue
        npaint += 1
        sum0 = 0.0; sumc = 0.0; sump = 0.0
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
                w0 = kappa_halo(θ, M, z, r; comp=0)
                wc = kappa_halo(θ, M, z, r; comp=1)
                wp = kappa_halo(θ, M, z, r; comp=2)
                h0[i, j] += w0; hc[i, j] += wc; hp[i, j] += wp
                sum0 += w0; sumc += wc; sump += wp
            end
        end
        # per-halo exactness: pixel-center sampling misses the NFW cusp for halos
        # near/below the pixel scale — deposit the residual vs the analytic totals
        # (plain: W·M_tot/ρ̄χ²; full comp: 0; partial comp: (1−1/1.636)·plain) into
        # the nearest pixel. Only when the full paint window fits the grid.
        icen = floor(Int, (gxh + s) / dpix) + 1
        jcen = floor(Int, (gyh + s) / dpix) + 1
        if 1 <= icen <= npixf && 1 <= jcen <= npixf &&
           gxh - gmax > -s && gxh + gmax < s && gyh - gmax > -s && gyh + gmax < s
            ki = kappa_halo_integral(M, z, r) / dpix^2
            h0[icen, jcen] += ki - sum0
            hc[icen, jcen] += -sumc
            hp[icen, jcen] += ki * (1 - 1 / MTOT_FAC) - sump
        end
    end
    (h0, hc, hp, npaint)
end

# ---------- run ----------
field_excl_path = ARGS[1]
field_all_path = ARGS[2]
ref_path = length(ARGS) >= 3 ? ARGS[3] : joinpath(WREF, "kap.fits")   # e.g. kap_lt4.5.fits

s = tan(deg2rad(CAPDEG)) / sqrt(2)
L = 2s
ledges = [10.0^l for l in range(log10(100.0), log10(4000.0); length=12)]
lc = [sqrt(ledges[i] * ledges[i+1]) for i in 1:length(ledges)-1]
axes = cap_axes()

# maps one at a time (kap.fits alone is 0.8 GB): extract all cap patches, then free
function patches_of(path, axes, s)
    m = Healpix.readMapFromFITS(path, 1, Float32)
    ns = m.resolution.nside
    p = [flat_patch(m, ax, s, NPIXFLAT) for ax in axes]
    m = nothing; GC.gc()
    p, ns
end
@info "extracting cap patches (maps read one at a time)..."
PFE, ns_fe = patches_of(field_excl_path, axes, s)
PFA, ns_fa = patches_of(field_all_path, axes, s)
PK,  ns_k  = patches_of(ref_path, axes, s)
@info "reference map" ref_path
@info "map nsides" field_excl=ns_fe field_all=ns_fa kap=ns_k

# Gaussian pixel-window approximations (good to <1% at these ell):
# HEALPix pixel size sqrt(pi/3)/nside; flat halo grid dpix; sigma = size/sqrt(12)
w2_hp(l, nside) = exp(-l * (l + 1) * (sqrt(π / 3) / nside)^2 / 12)
w2_flat(l, dpix) = exp(-l * (l + 1) * dpix^2 / 12)

@info "streaming catalog (26 GB) for per-cap halo lists..." CAT
caphalos = select_cap_halos(CAT, chi2z, axes, s)
@info "selected" nper=[length(h.M) for h in caphalos]

nb = length(lc)
cA = zeros(nb, NCAPS); cB = zeros(nb, NCAPS); cB2 = zeros(nb, NCAPS); cK = zeros(nb, NCAPS)
cFE = zeros(nb, NCAPS); cFA = zeros(nb, NCAPS); cH0 = zeros(nb, NCAPS); cHC = zeros(nb, NCAPS)
for (ic, ax) in enumerate(axes)
    pfe = PFE[ic]; pfa = PFA[ic]; pk_ = PK[ic]
    h0, hc, hp, npaint = paint_halos_cap(caphalos[ic], chi2z, ax, s, NPIXFLAT)
    pA = pfe .+ h0                     # scheme A: excluded field + plain NFW
    pB = pfa .+ hc                     # scheme B: full field + Δ=3 comp (mass 1.636·M)
    pB2 = pfa .+ hp                    # scheme B2: full field + Δ=3 comp (mass M, paper literal)
    cA[:, ic] = cl_flat(pA, L, ledges);  cB[:, ic] = cl_flat(pB, L, ledges)
    cB2[:, ic] = cl_flat(pB2, L, ledges)
    cK[:, ic] = cl_flat(pk_, L, ledges)
    cFE[:, ic] = cl_flat(pfe, L, ledges); cFA[:, ic] = cl_flat(pfa, L, ledges)
    cH0[:, ic] = cl_flat(h0, L, ledges);  cHC[:, ic] = cl_flat(hc, L, ledges)
    @info "cap $ic done" npaint mean_A=round(mean(pA); digits=4) mean_B=round(mean(pB); digits=4) mean_B2=round(mean(pB2); digits=4) mean_kap=round(mean(pk_); digits=4)
end

@printf("\n%-7s %-11s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
        "ell", "Cl_ref", "A/ref", "+-", "A/ref*", "B/ref", "+-", "B2/ref", "+-", "fldA/ref", "halo/ref")
dpixf = 2s / NPIXFLAT
for i in eachindex(lc)
    rA = [cA[i, c] / cK[i, c] for c in 1:NCAPS]
    rB = [cB[i, c] / cK[i, c] for c in 1:NCAPS]
    rB2 = [cB2[i, c] / cK[i, c] for c in 1:NCAPS]
    # window-corrected A/ref from the measured per-band decomposition:
    # A = fldE (field-map window) + halo (flat-grid window) + cross (geometric mean)
    wfe = w2_hp(lc[i], ns_fe); wh = w2_flat(lc[i], dpixf); wk = w2_hp(lc[i], ns_k)
    fldE = mean(cFE[i, :]); halo = mean(cH0[i, :]); crossA = mean(cA[i, :]) - fldE - halo
    Astar = (fldE / wfe + halo / wh + crossA / sqrt(wfe * wh)) / (mean(cK[i, :]) / wk)
    @printf("%-7.0f %-11.3e %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",
            lc[i], mean(cK[i, :]), mean(rA), std(rA), Astar, mean(rB), std(rB),
            mean(rB2), std(rB2),
            mean(cFA[i, :]) / mean(cK[i, :]), halo / mean(cK[i, :]))
end
sel = findall(l -> 150 <= l <= 2500, lc)
mA = mean([mean(cA[i, :]) / mean(cK[i, :]) for i in sel])
mB = mean([mean(cB[i, :]) / mean(cK[i, :]) for i in sel])
mB2 = mean([mean(cB2[i, :]) / mean(cK[i, :]) for i in sel])
@printf("\nband mean 150<=ell<=2500:  A/ref = %.3f   B/ref = %.3f   B2/ref = %.3f\n", mA, mB, mB2)
@printf("Caveats: different realizations (cap scatter is the error bar). Reference map\n")
@printf("choice matters: kap.fits INCLUDES the z>4.5 Gaussian tail (15-45%% of C_l,\n")
@printf("measured from kap_gt4.5.fits); use kap_lt4.5.fits for apples-to-apples.\n")
