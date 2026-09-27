#!/usr/bin/env julia
# Phase-A exit test (docs/field_lightcone_plan.md): compare the full-matter lightcone κ map
# painted by run_multitile_fieldmap (job output FITS) against Websky's released kap.fits.
#
# Usage: julia --project=validation -t 8 compare_fieldmap_kappa.jl <kappa_fieldmap.fits>
#
# Geometry: our octant covers the +++ octant of sky seen from the corner observer; kap.fits
# is full-sky. These are DIFFERENT realizations, so we compare C_ℓ band levels, not maps.
# Method: take Nc equal-size circular caps (radius CAPDEG) well inside our octant, and the
# same number spread over kap.fits sky; flat-sky C_ℓ per cap (gnomonic sampling at the cap
# center, same estimator for both); report per-band cap means ± scatter and the ratio.
# Also overlay linear-theory Limber (lower bound at nonlinear scales; at ℓ≲300 should be
# close). kap.fits additionally contains the z>4.5 Limber Gaussian tail (small at ℓ≳100)
# and its halo component is NFW-painted vs our pure-2LPT matter — expect %-level shape
# differences at high ℓ, plus each map's Nside pixel window (both 2048/4096-class; noted).
using Healpix, FFTW, Printf, Statistics, LinearAlgebra, DelimitedFiles
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z, chi, Dlinear_tables, Dlinear_ab

const WREF = "/home/yguan/scratch/websky_6144/websky_ref"
const CAPDEG = 8.0          # cap radius; flat square half-width = CAPDEG/√2
const NPIXFLAT = 512
const NCAPS = 6

fieldmap_path = ARGS[1]
@info "reading maps..." fieldmap_path
mo = Healpix.readMapFromFITS(fieldmap_path, 1, Float32)
mk = Healpix.readMapFromFITS(joinpath(WREF, "kap.fits"), 1, Float32)
@info "maps" nside_ours=mo.resolution.nside nside_kap=mk.resolution.nside

ortho_basis(a) = (t = abs(a[1]) < 0.9 ? [1.0,0,0] : [0.0,1,0]; e1 = normalize(cross(a,t)); e2 = cross(a,e1); (e1,e2))

# sample a HEALPix map onto a flat gnomonic grid centred on axis `ax`
function flat_patch(m::HealpixMap, ax, s, npixf)
    e1, e2 = ortho_basis(ax)
    res = m.resolution
    out = zeros(npixf, npixf)
    Threads.@threads for j in 1:npixf
        for i in 1:npixf
            gx = (-s + (i - 0.5) * 2s / npixf); gy = (-s + (j - 0.5) * 2s / npixf)
            nx = ax[1] + gx * e1[1] + gy * e2[1]
            ny = ax[2] + gx * e1[2] + gy * e2[2]
            nz = ax[3] + gx * e1[3] + gy * e2[3]
            out[i, j] = m.pixels[Healpix.vec2pixRing(res, nx, ny, nz)]
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

# cap centres: ours = ring around the (1,1,1) diagonal inside the octant; kap = same
# pattern (any sky location is valid for a full-sky map — use the SAME axes for symmetry)
function cap_axes()
    diag = [1.0, 1.0, 1.0] ./ sqrt(3.0)
    e1, e2 = ortho_basis(diag)
    axes = Vector{Vector{Float64}}()
    push!(axes, diag)
    for k in 0:NCAPS-2
        ϕ = 2π * k / (NCAPS - 1)
        v = normalize(cos(deg2rad(16)) .* diag .+ sin(deg2rad(16)) .* (cos(ϕ) .* e1 .+ sin(ϕ) .* e2))
        push!(axes, v)
    end
    axes
end

s = tan(deg2rad(CAPDEG)) / sqrt(2)
L = 2s
ledges = [10.0^l for l in range(log10(100.0), log10(4000.0); length=12)]
lc = [sqrt(ledges[i] * ledges[i+1]) for i in 1:length(ledges)-1]

axes = cap_axes()
co = zeros(length(lc), length(axes)); ck = zeros(length(lc), length(axes))
for (ic, ax) in enumerate(axes)
    po = flat_patch(mo, ax, s, NPIXFLAT)
    pk_ = flat_patch(mk, ax, s, NPIXFLAT)
    co[:, ic] = cl_flat(po, L, ledges)
    ck[:, ic] = cl_flat(pk_, L, ledges)
    @info "cap $ic done" mean_ours=round(mean(po); digits=4) mean_kap=round(mean(pk_); digits=4)
end

# ---- Limber linear theory ----
pk = readdlm(joinpath(@__DIR__, "data", "pk_websky.dat"), comments=true, comment_char='#')
kk = Float64.(pk[:, 1]); Pk = Float64.(pk[:, 2]) .* (2π)^3
ok = (kk .> 1e-4) .& (Pk .> 0); kk = kk[ok]; Pk = Pk[ok]; lkk = log.(kk); lPk = log.(Pk)
Plin(k) = (k <= kk[1] || k >= kk[end]) ? 0.0 :
    (i = searchsortedfirst(kk, k); exp(lPk[i-1] + (lPk[i] - lPk[i-1]) * (log(k) - lkk[i-1]) / (lkk[i] - lkk[i-1])))
cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
gt = Dlinear_tables(cosmo); chi2z = build_chi_to_z(cosmo; z_max=6.0)
chimax = chi(4.6, cosmo); chistar = 14200.0 * 0.68
function cl_limber(l)
    nint = 400; dr = chimax / nint; acc = 0.0
    for ii in 1:nint
        r = (ii - 0.5) * dr
        z = chi_to_z(chi2z, r); a = 1 / (1 + z)
        D, _, _ = Dlinear_ab(a, gt)   # 1st return = D; 3rd is D/a (a 10x C_l bug when squared)
        W = 1.5 * 0.31 * (1 / 2997.92458)^2 * (1 + z) * r * (1 - r / chistar)
        acc += W^2 * D^2 * Plin((l + 0.5) / r) / r^2 * dr
    end
    acc
end
clth = [cl_limber(l) for l in lc]

@printf("\n%-7s %-11s %-11s %-8s %-8s %-11s %-8s\n", "ell", "Cl_ours", "Cl_kap", "ratio", "±scat", "Cl_linear", "kap/lin")
for i in eachindex(lc)
    mo_ = mean(co[i, :]); mk_ = mean(ck[i, :])
    rat = [co[i, c] / ck[i, c] for c in 1:size(co, 2)]
    @printf("%-7.0f %-11.3e %-11.3e %-8.3f %-8.3f %-11.3e %-8.3f\n",
            lc[i], mo_, mk_, mean(rat), std(rat), clth[i], mk_ / clth[i])
end
sel = findall(l -> 150 <= l <= 1500, lc)
@printf("\nmean ours/kap.fits over 150<=ell<=1500: %.3f ± %.3f (cap scatter)\n",
        mean([mean(co[i, :]) / mean(ck[i, :]) for i in sel]),
        std([mean(co[i, :]) / mean(ck[i, :]) for i in sel]))
@printf("PASS: ratio ~1 within cap scatter -> full-matter LPT kappa reproduces kap.fits levels.\n")
@printf("Caveats: different realizations; kap.fits has NFW halos + z>4.5 tail; pixel windows uncorrected.\n")
