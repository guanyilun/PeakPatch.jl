#!/usr/bin/env julia
# Smoke test for run_multitile_fieldmap (Phase A): small CPU lightcone box.
# Checks (1) mass conservation: Σ :mass map == ρ̄·a³latt·N_cells(rmin ≤ r ≤ χmax),
#        (2) mean κ: map mean == (Ω_octant/4π)·∫W_κ dχ analytically (grid-quantized),
#        (3) pixel bookkeeping sanity (all mass lands in the observer's octant of sky).
# Usage: julia --project=. -t 4 validation/websky_6144/test_fieldmap_smoke.jl [--gpu]
using PeakPatch, Healpix, Printf, Statistics
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

@info "running fieldmap (subdiv_max=1: exact bookkeeping)..." N ntile use_gpu=USE_GPU
maps = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix, vec2pix=v2p,
                              kernels=[:kappa, :mass], subdiv_max=1,
                              use_gpu=USE_GPU, devices=USE_GPU ? [0] : nothing,
                              verbose=false)
@info "running fieldmap (subdiv_max=3: exercises splitting)..."
maps3 = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix, vec2pix=v2p,
                               kernels=[:kappa, :mass], subdiv_max=3,
                               use_gpu=USE_GPU, devices=USE_GPU ? [0] : nothing,
                               verbose=false)

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
