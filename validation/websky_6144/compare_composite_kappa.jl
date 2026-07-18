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
# Usage: julia --project=. -t 8 compare_composite_kappa.jl <kappa_field_excl.fits> <kappa_field_all.fits>
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
@inline function kappa_halo(θ, M, z, χ; comp::Bool)
    r200 = cbrt(3 * M / (800π * rho_mh))
    rs = r200 / CNFW
    b = θ * χ
    Σ = b / rs < XMAX * CNFW ? RHOS_OVER_RHOM * rs * g_of_x(b / rs) : 0.0
    if comp
        Rc = cbrt(3 * MTOT_FAC * M / (4π * DCOMP * rho_mh))
        b < Rc && (Σ -= DCOMP * 2 * sqrt(Rc^2 - b^2))
    end
    1.5 * 0.31 * H0C^2 * (1 + z) * (1 - χ / chistar) / χ * Σ
end
function paint_radius(M; comp::Bool)
    r200 = cbrt(3 * M / (800π * rho_mh))
    fac = comp ? max(XMAX, cbrt(MTOT_FAC * 200 / DCOMP)) : XMAX   # comp sphere: 4.78·r200
    return fac * r200
end

# self-test: ∫κ 2πθ dθ · χ²ρ̄/W == M_tot (plain) and ≈0 (compensated)
let M = 3e14, z = 0.7
    cosmo0 = CosmologyParams(0.31, 0.049, 0.69, hub, 0.965, 0.81)
    χ = chi(z, cosmo0)
    W = 1.5 * 0.31 * H0C^2 * (1 + z) * (1 - χ / chistar) / χ
    for comp in (false, true)
        θmax = paint_radius(M; comp=comp) / χ
        n = 40000; h = θmax / n; acc = 0.0
        for i in 1:n
            θ = (i - 0.5) * h
            acc += kappa_halo(θ, M, z, χ; comp=comp) * 2π * θ * h
        end
        got = acc * χ^2 * rho_mh / W / M
        want = comp ? 0.0 : MTOT_FAC
        @printf("NFW self-test comp=%-5s: mass integral / M200m = %+.4f (expect %+.4f)\n",
                comp, got, want)
        abs(got - want) < 0.01 * MTOT_FAC || error("NFW painting self-test failed")
    end
end

# ---------- load catalog once: Eulerian positions + masses ----------
# NOTE: all hot loops live in functions taking the arrays as ARGUMENTS — non-const
# globals in Julia dynamically dispatch every operation (100-1000× slowdown).
function load_catalog(path, chi2z)
    nh_tot = open(path) do io Int(read(io, Int32)) end
    EX = Vector{Float32}(undef, nh_tot); EY = similar(EX); EZ = similar(EX)
    MM = similar(EX)
    open(path) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        nf = 33; chunk = 1_000_000; buf = Vector{Float32}(undef, chunk * nf); ndone = 0
        while ndone < nh_tot
            m = min(chunk, nh_tot - ndone); read!(io, view(buf, 1:m*nf))
            @inbounds Threads.@threads for k in 1:m
                b = (k - 1) * nf
                R = Float64(buf[b+7])
                q1 = Float64(buf[b+1]); q2 = Float64(buf[b+2]); q3 = Float64(buf[b+3])
                rq = sqrt((q1 - OBS)^2 + (q2 - OBS)^2 + (q3 - OBS)^2)
                zq = chi_to_z(chi2z, rq)
                a = 1.0 / (1.0 + zq)
                # old-convention catalog: Eulerian = q + d1·a − d2·a² (stored ψ₂ carries +3/7)
                j = ndone + k
                EX[j] = Float32(q1 + Float64(buf[b+4]) * a - Float64(buf[b+8]) * a^2 - OBS)
                EY[j] = Float32(q2 + Float64(buf[b+5]) * a - Float64(buf[b+9]) * a^2 - OBS)
                EZ[j] = Float32(q3 + Float64(buf[b+6]) * a - Float64(buf[b+10]) * a^2 - OBS)
                MM[j] = Float32(4 / 3 * π * rho_mh * R^3)
            end
            ndone += m
        end
    end
    return EX, EY, EZ, MM
end
@info "loading catalog (26 GB single pass)..." CAT
cosmo = CosmologyParams(0.31, 0.049, 0.69, hub, 0.965, 0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0)
EX, EY, EZ, MM = load_catalog(CAT, chi2z)
@info "catalog loaded" nh=length(MM)

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

# paint the cap's halos on the SAME gnomonic grid, exact per-pixel angles.
# Returns (plain, compensated) patches in one pass. Arrays passed as args (see NOTE).
function paint_halos_cap(EX, EY, EZ, MM, chi2z, ax, s, npixf)
    e1, e2 = ortho_basis(ax)
    h0 = zeros(npixf, npixf); hc = zeros(npixf, npixf)
    dpix = 2s / npixf
    npaint = 0
    @inbounds for n in eachindex(MM)
        vx = Float64(EX[n]); vy = Float64(EY[n]); vz = Float64(EZ[n])
        r = sqrt(vx^2 + vy^2 + vz^2)
        (30.0 <= r <= 5170.0) || continue
        na = (vx*ax[1] + vy*ax[2] + vz*ax[3]) / r
        na > 0.9 || continue
        gxh = (vx*e1[1] + vy*e1[2] + vz*e1[3]) / (r * na)
        gyh = (vx*e2[1] + vy*e2[2] + vz*e2[3]) / (r * na)
        M = Float64(MM[n])
        θmaxc = paint_radius(M; comp=true) / r
        gmax = θmaxc * 1.3 / na           # conservative gnomonic-plane radius
        (abs(gxh) - gmax <= s && abs(gyh) - gmax <= s) || continue
        z = chi_to_z(chi2z, r)
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
                cosang = (px*vx + py*vy + pz*vz) /
                         (sqrt(px^2 + py^2 + pz^2) * r)
                θ = acos(clamp(cosang, -1.0, 1.0))
                θ > θmaxc && continue
                h0[i, j] += kappa_halo(θ, M, z, r; comp=false)
                hc[i, j] += kappa_halo(θ, M, z, r; comp=true)
            end
        end
    end
    (h0, hc, npaint)
end

# ---------- run ----------
field_excl_path = ARGS[1]
field_all_path = ARGS[2]
@info "reading maps..." field_excl_path field_all_path
mfe = Healpix.readMapFromFITS(field_excl_path, 1, Float32)
mfa = Healpix.readMapFromFITS(field_all_path, 1, Float32)
mk = Healpix.readMapFromFITS(joinpath(WREF, "kap.fits"), 1, Float32)

s = tan(deg2rad(CAPDEG)) / sqrt(2)
L = 2s
ledges = [10.0^l for l in range(log10(100.0), log10(4000.0); length=12)]
lc = [sqrt(ledges[i] * ledges[i+1]) for i in 1:length(ledges)-1]
axes = cap_axes()

nb = length(lc)
cA = zeros(nb, NCAPS); cB = zeros(nb, NCAPS); cK = zeros(nb, NCAPS)
cFE = zeros(nb, NCAPS); cFA = zeros(nb, NCAPS); cH0 = zeros(nb, NCAPS); cHC = zeros(nb, NCAPS)
for (ic, ax) in enumerate(axes)
    pfe = flat_patch(mfe, ax, s, NPIXFLAT)
    pfa = flat_patch(mfa, ax, s, NPIXFLAT)
    pk_ = flat_patch(mk, ax, s, NPIXFLAT)
    h0, hc, npaint = paint_halos_cap(EX, EY, EZ, MM, chi2z, ax, s, NPIXFLAT)
    pA = pfe .+ h0                     # scheme A: excluded field + plain NFW
    pB = pfa .+ hc                     # scheme B: full field + Δ=3-compensated NFW
    cA[:, ic] = cl_flat(pA, L, ledges);  cB[:, ic] = cl_flat(pB, L, ledges)
    cK[:, ic] = cl_flat(pk_, L, ledges)
    cFE[:, ic] = cl_flat(pfe, L, ledges); cFA[:, ic] = cl_flat(pfa, L, ledges)
    cH0[:, ic] = cl_flat(h0, L, ledges);  cHC[:, ic] = cl_flat(hc, L, ledges)
    @info "cap $ic done" npaint mean_A=round(mean(pA); digits=4) mean_B=round(mean(pB); digits=4) mean_kap=round(mean(pk_); digits=4)
end

@printf("\n%-7s %-11s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
        "ell", "Cl_kap", "A/kap", "+-", "B/kap", "+-", "fldE/kap", "fldA/kap", "halo/kap")
for i in eachindex(lc)
    rA = [cA[i, c] / cK[i, c] for c in 1:NCAPS]
    rB = [cB[i, c] / cK[i, c] for c in 1:NCAPS]
    @printf("%-7.0f %-11.3e %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",
            lc[i], mean(cK[i, :]), mean(rA), std(rA), mean(rB), std(rB),
            mean(cFE[i, :]) / mean(cK[i, :]), mean(cFA[i, :]) / mean(cK[i, :]),
            mean(cH0[i, :]) / mean(cK[i, :]))
end
sel = findall(l -> 150 <= l <= 2500, lc)
mA = mean([mean(cA[i, :]) / mean(cK[i, :]) for i in sel])
mB = mean([mean(cB[i, :]) / mean(cK[i, :]) for i in sel])
@printf("\nband mean 150<=ell<=2500:  A/kap = %.3f   B/kap = %.3f\n", mA, mB)
@printf("Caveats: different realizations (cap scatter is the error bar); kap.fits includes\n")
@printf("the z>4.5 Gaussian tail (small at ell>~100) and its own Nside-4096 pixel window.\n")
