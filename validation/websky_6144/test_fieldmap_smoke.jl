#!/usr/bin/env julia
# Smoke test for run_multitile_fieldmap (Phase A): small CPU lightcone box.
# Checks (1) mass conservation: Σ :mass map == ρ̄·a³latt·N_cells(rmin ≤ r ≤ χmax),
#        (2) mean κ: map mean == (Ω_octant/4π)·∫W_κ dχ analytically (grid-quantized),
#        (3) pixel bookkeeping sanity (all mass lands in the observer's octant of sky).
# Usage: julia --project=. -t 4 validation/websky_6144/test_fieldmap_smoke.jl [--gpu]
using PeakPatch, Healpix, Printf, Statistics, Random
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z, chi

const USE_GPU = "--gpu" in ARGS
if USE_GPU
    using CUDA
end

datadir = joinpath(@__DIR__, "data")
alatt = 8.0; nmesh = 64; nbuff = 8; ntile = 2
nsub = nmesh - 2nbuff; N = nsub * ntile + 2nbuff   # 112
boxfull = N * alatt
config = Dict{String,Any}(
    "cosmology" => Dict{String,Any}("Om" => 0.31, "OB" => 0.049, "OL" => 0.69, "h" => 0.68),
    "grid" => Dict{String,Any}("n" => nmesh, "boxsize" => nmesh * alatt, "nbuff" => nbuff,
                    "cenx" => -384.0, "ceny" => -384.0, "cenz" => -384.0),
    "run" => Dict{String,Any}("ievol" => 1, "z_max" => 0.1, "z_out" => 0.0, "ilpt" => 2,
                   "ioutshear" => 0),
    "files" => Dict{String,Any}("pk" => joinpath(datadir, "pk_websky.dat"),
                     "filterbank" => joinpath(datadir, "filters_websky.dat"),
                     "homeltab" => joinpath(datadir, "HomelTab_websky.dat"),
                     "output" => "/tmp/fieldmap_smoke.pksc"),
)
cfg = PipelineConfig(config)

nside = 64; res = Resolution(nside); npix = nside2npix(nside)
v2p = (x, y, z) -> Healpix.vec2pixRing(res, x, y, z)

# ---- ang2pix_ring (own implementation, used by the GPU painter) vs Healpix.jl ----
let a2p = PeakPatch.MultiResolution.ang2pix_ring, rng = Random.MersenneTwister(7), nbad = 0
    for _ in 1:200_000
        x = randn(rng); y = randn(rng); z = randn(rng)
        a2p(nside, x, y, z) == Healpix.vec2pixRing(res, x, y, z) || (nbad += 1)
    end
    @printf("ang2pix_ring vs Healpix.jl: %d/200000 mismatches (exact-boundary ties only)\n", nbad)
end

KERNELS = [:kappa, :mass, :tau, :ksz]
@info "running fieldmap (subdiv_max=1: exact bookkeeping)..." N ntile use_gpu=USE_GPU
maps = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix, vec2pix=v2p,
                              kernels=KERNELS, subdiv_max=1,
                              use_gpu=USE_GPU, devices=USE_GPU ? [0] : nothing,
                              verbose=false)
@info "running fieldmap (subdiv_max=3: exercises splitting)..."
maps3 = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix, vec2pix=v2p,
                               kernels=KERNELS, subdiv_max=3,
                               use_gpu=USE_GPU, devices=USE_GPU ? [0] : nothing,
                               verbose=false)

# ---- Multi-worker dispatch: 4 CPU workers must reproduce the single-worker maps ----
# (Painted totals are psi-independent, so sums must match to summation-order roundoff.
#  Guards against shared/aliased worker accumulators: job 4280770 painted 2^(n-1)=8x.)
@info "running fieldmap (cpu_workers=4: multi-worker dispatch)..."
maps4 = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix, vec2pix=v2p,
                               kernels=KERNELS, subdiv_max=3,
                               use_gpu=false, cpu_workers=4, verbose=false)
@printf("cpu_workers=4 vs 1: mass ratio=%.8f  kappa ratio=%.8f (both MUST be 1.00000000)\n",
        sum(maps4[:mass]) / sum(maps3[:mass]), sum(maps4[:kappa]) / sum(maps3[:kappa]))

# ---- Phase B cross-check: on-device painting (gpu_paint) vs CPU pixelization ----
# Non-fatal: reports and continues so a Phase-B regression never blocks the octant job.
if USE_GPU
    @info "running fieldmap (gpu_paint=true: device RING pixelization)..."
    try
        mapsg = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix,
                                       kernels=KERNELS, subdiv_max=3,
                                       use_gpu=true, devices=[0],
                                       gpu_paint=true, nside=nside, verbose=false)
        mr = sum(mapsg[:mass]) / sum(maps3[:mass])
        kr = sum(mapsg[:kappa]) / sum(maps3[:kappa])
        dk = maximum(abs.(mapsg[:kappa] .- maps3[:kappa])) / maximum(maps3[:kappa])
        ndiff = count(abs.(mapsg[:mass] .- maps3[:mass]) .> 1e-6 .* maximum(maps3[:mass]))
        @printf("gpu_paint vs cpu: mass ratio=%.8f  kappa ratio=%.8f  max|dkappa|/max=%.2e  npix-diff=%d\n",
                mr, kr, dk, ndiff)
        dt = maximum(abs.(mapsg[:tau] .- maps3[:tau])) / maximum(maps3[:tau])
        ksz_rms = sqrt(sum(abs2, maps3[:ksz]) / count(!iszero, maps3[:ksz]))
        dz = maximum(abs.(mapsg[:ksz] .- maps3[:ksz])) / ksz_rms
        @printf("gpu_paint vs cpu: max|dtau|/max=%.2e  max|dksz|/rms=%.2e\n", dt, dz)
    catch err
        @error "gpu_paint cross-check FAILED (non-fatal)" exception=(err, catch_backtrace())
    end
end

# ---- expected cell count / mass ----
cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.808)
chimax = chi(0.1, cosmo)
rho_m = 2.775e11 * 0.31
obs = (-384.0, -384.0, -384.0)
rmin_eff = 2 * alatt
ncell = 0
Nc = nsub * ntile
for k in 1:Nc, j in 1:Nc, i in 1:Nc
    x = (i - (Nc + 1) / 2) * alatt; y = (j - (Nc + 1) / 2) * alatt; zz = (k - (Nc + 1) / 2) * alatt
    r = sqrt((x - obs[1])^2 + (y - obs[2])^2 + (zz - obs[3])^2)
    (rmin_eff <= r <= chimax) && (global ncell += 1)
end
Mexp = rho_m * alatt^3 * ncell
Mgot = sum(maps[:mass])
@printf("mass conservation: painted=%.6e  expected=%.6e  ratio=%.6f\n", Mgot, Mexp, Mgot/Mexp)

# ---- mean kappa vs analytic ----
chi2z = build_chi_to_z(cosmo; z_max=1.0)
chistar = chi(1089.0, cosmo)
nint = 2000; dr = (chimax - rmin_eff) / nint; acc = 0.0
for ii in 0:nint-1
    r = rmin_eff + (ii + 0.5) * dr
    z = chi_to_z(chi2z, r)
    global acc += 1.5 * 0.31 * (1 / 2997.92458)^2 * (1 + z) * (1 - r / chistar) * r * dr
end
kexp = acc / 8                       # octant of sky
kgot = sum(maps[:kappa]) / npix
@printf("mean kappa:        painted=%.6e  expected=%.6e  ratio=%.6f\n", kgot, kexp, kgot/kexp)
@printf("subdiv=3 vs 1:     mass ratio=%.6f (boundary sub-cell truncation, expect ~0.99)\n", sum(maps3[:mass])/Mgot)
@printf("subdiv=3 vs 1:     kappa-mean ratio=%.6f\n", (sum(maps3[:kappa])/npix)/kgot)

# ---- mean tau vs analytic (same construction as mean kappa) ----
sigT_ne0 = 6.65246e-29 * 11.2299 * 0.9 * 0.049 * 0.68^2 * 3.0857e22 / 0.68
tacc = 0.0
for ii in 0:nint-1
    r = rmin_eff + (ii + 0.5) * dr
    z = chi_to_z(chi2z, r)
    x_e = z < 3 ? (1 - 0.245 / 2) : (1 - 3 * 0.245 / 4)
    global tacc += sigT_ne0 * x_e * (1 + z)^2 * dr
end
texp = tacc / 8
tgot = sum(maps[:tau]) / npix
@printf("mean tau:          painted=%.6e  expected=%.6e  ratio=%.6f\n", tgot, texp, tgot/texp)

# ---- kSZ: signed map, mean must cancel against rms; rms at the tau*v/c scale ----
cov = findall(!iszero, maps[:tau])
kszm = sum(maps[:ksz][cov]) / length(cov)
kszr = sqrt(sum(abs2, maps[:ksz][cov]) / length(cov))
vr_eff = abs(kszm) / (sum(maps[:tau][cov]) / length(cov)) * 299792.458
@printf("ksz (covered pix): mean=%.3e  rms=%.3e  |mean|/rms=%.3f (should be <<1)\n",
        kszm, kszr, abs(kszm) / kszr)
@printf("ksz implied bulk v_r: %.1f km/s (few-hundred km/s coherent flow OK at this tiny volume)\n", vr_eff)

# ---- octant containment: pixels with mass should have direction in +++ octant (from obs at corner) ----
pix_on = findall(>(0), maps[:mass])
bad = 0
for p in pix_on
    v = pix2vecRing(res, p)
    (v[1] > -0.15 && v[2] > -0.15 && v[3] > -0.15) || (global bad += 1)
end
@printf("octant containment: %d/%d nonzero pixels outside +++ octant (should be ~0)\n", bad, length(pix_on))
@printf("nonzero pixels: %d/%d (octant fraction=%.3f, expect ~0.125 of sky covered by chi_max cone)\n",
        length(pix_on), npix, length(pix_on)/npix)

# ---- exclude_halos: mass deficit must equal the exact count of excluded cells ----
@info "running fieldmap (exclude_halos: 2 synthetic Lagrangian spheres)..."
halos = (x=[-250.0, -150.0], y=[-250.0, -320.0], z=[-250.0, -300.0], R=[30.0, 20.0])
maps_ex = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix, vec2pix=v2p,
                                 kernels=[:mass], subdiv_max=1, exclude_halos=halos,
                                 use_gpu=USE_GPU, devices=USE_GPU ? [0] : nothing,
                                 verbose=false)
n_excl = 0
for n in 1:2
    hx, hy, hz, R = halos.x[n], halos.y[n], halos.z[n], halos.R[n]
    for k in 1:Nc, j in 1:Nc, i in 1:Nc
        x = (i - (Nc + 1) / 2) * alatt; y = (j - (Nc + 1) / 2) * alatt; zz = (k - (Nc + 1) / 2) * alatt
        (x - hx)^2 + (y - hy)^2 + (zz - hz)^2 <= R^2 || continue
        r = sqrt((x - obs[1])^2 + (y - obs[2])^2 + (zz - obs[3])^2)
        (rmin_eff <= r <= chimax) && (global n_excl += 1)
    end
end
deficit = Mgot - sum(maps_ex[:mass])
@printf("exclude_halos: deficit=%.6e  expected=%.6e (%d cells)  ratio=%.6f\n",
        deficit, rho_m * alatt^3 * n_excl, n_excl, deficit / (rho_m * alatt^3 * n_excl))

# gpu_paint × exclusion (the combination the field-4096 production job uses):
# mass is psi-independent, so the sum must match the CPU-pixelized excluded run exactly
if USE_GPU
    maps_exg = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix,
                                      kernels=[:mass], subdiv_max=1, exclude_halos=halos,
                                      use_gpu=true, devices=[0],
                                      gpu_paint=true, nside=nside, verbose=false)
    @printf("gpu_paint exclusion: mass ratio=%.8f (MUST be 1.00000000)\n",
            sum(maps_exg[:mass]) / sum(maps_ex[:mass]))
end

# ---- cross-tile mask rasterization: sphere straddling the x=0 tile boundary ----
let bmask = PeakPatch.MultiResolution._build_exclusion_mask,
    tc = PeakPatch.MultiResolution.tile_center
    hb = (x=[0.0], y=[-200.0], z=[-200.0], R=[30.0])
    dcore = nsub * alatt
    x0 = -(ntile / 2) * dcore
    bins = Dict{NTuple{3,Int},Vector{Int}}()
    bt = (clamp(floor(Int, (hb.x[1] - x0) / dcore) + 1, 1, ntile),
          clamp(floor(Int, (hb.y[1] - x0) / dcore) + 1, 1, ntile),
          clamp(floor(Int, (hb.z[1] - x0) / dcore) + 1, 1, ntile))
    bins[bt] = [1]
    nmask = 0
    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        xbx, ybx, zbx = tc(it, jt, kt, ntile, dcore)
        m = bmask(hb, bins, it, jt, kt, ntile, dcore, nmesh, nbuff, alatt, xbx, ybx, zbx)
        nmask += count(m)
    end
    nbrute = 0
    for k in 1:Nc, j in 1:Nc, i in 1:Nc
        x = (i - (Nc + 1) / 2) * alatt; y = (j - (Nc + 1) / 2) * alatt; zz = (k - (Nc + 1) / 2) * alatt
        (x - hb.x[1])^2 + (y - hb.y[1])^2 + (zz - hb.z[1])^2 <= hb.R[1]^2 && (nbrute += 1)
    end
    @printf("cross-tile mask: masked=%d  brute-force=%d  (must be equal)\n", nmask, nbrute)
end
