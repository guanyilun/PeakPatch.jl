#!/usr/bin/env julia
# Numerical probe for the multi-resolution 1.6x over-production.
# Build the field BOTH ways from the SAME Threefry noise (same seed) and compare
# sigma(R) per scale, to localize where the coarse+residual split has excess power.
#
#   global = generate_grf(N) (full FFT, Fortran-equivalent, correct)
#   split  = interpolate(coarse_field) + isolated_conv(residual)  for a central tile
#
# Usage: julia --project=validation -t 8 validation/websky_6144/probe_field_sigma.jl

using PeakPatch, FFTW, Statistics, Printf
const MR = PeakPatch.MultiResolution

# --- geometry: same as the A/B (N=256, cellsize 1.25326) ---
N = 256; ntile = 4; nbuff = 24; cf = 4
nsub  = (N - 2*nbuff) ÷ ntile          # 52
nmesh = nsub + 2*nbuff                  # 100
M     = ntile * cf                      # 16
alatt = 1.25326
box_full  = N * alatt                   # 320.8
box_local = nmesh * alatt               # 125.3
seed = 13579
pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "data", "pk_websky.dat"))

@info "geometry" N nmesh nsub M box_full

# --- A. global field (correct reference) ---
delta_global = PeakPatch.RandomField.generate_grf(N, pk, box_full, seed)

# --- B. split field for central tile (it,jt,kt) = (2,2,2) ---
it = jt = kt = 2
coarse_noise = MR._downsample_noise(N, M, seed)
ck = rfft(coarse_noise)
MR._periodic_convolve!(ck, pk, M, box_full)
delta_coarse = irfft(ck, M)
delta_long = MR._interpolate_to_tile(delta_coarse, it, jt, kt, nsub, nmesh, N, M)
residual = MR._generate_extended_residual(it, jt, kt, nsub, nmesh, N, seed, coarse_noise, M, 0)
delta_self = MR._isolated_convolve(residual, pk, box_local, nmesh)
delta_tile = delta_self .+ delta_long

# core slices: tile core in tile-local coords, and the matching global region
clo, chi = nbuff+1, nmesh-nbuff                       # 25..76 (52 cells)
gi0 = (it-1)*nsub                                     # global offset of tile core
gslice = (gi0+1):(gi0+nsub)                           # 53..104

core_tile(a)   = @view a[clo:chi, clo:chi, clo:chi]
core_global(a) = @view a[gslice, gslice, gslice]

@printf("\n%-8s | %-12s | %-12s | %-12s | %-12s\n", "R[Mpc/h]", "sigma_global", "sigma_split", "split/global", "(also: coarse,resid)")
# unsmoothed (R=0)
sg = std(core_global(delta_global)); ss = std(core_tile(delta_tile))
sc = std(core_tile(delta_long)); sr = std(core_tile(delta_self))
@printf("%-8s | %-12.4f | %-12.4f | %-12.3f | coarse=%.3f resid=%.3f\n", "raw", sg, ss, ss/sg, sc, sr)

gk_full = rfft(delta_global)
tk_full = rfft(delta_tile)
lk_full = rfft(delta_long); rk_full = rfft(delta_self)
for R in (2.0, 4.0, 8.0, 16.0, 32.0)
    g = PeakPatch.Filters.smooth_field(gk_full, N, box_full, R, 1)
    t = PeakPatch.Filters.smooth_field(tk_full, nmesh, box_local, R, 1)
    l = PeakPatch.Filters.smooth_field(lk_full, nmesh, box_local, R, 1)
    r = PeakPatch.Filters.smooth_field(rk_full, nmesh, box_local, R, 1)
    sg = std(core_global(g)); ss = std(core_tile(t))
    @printf("%-8.1f | %-12.4f | %-12.4f | %-12.3f | coarse=%.3f resid=%.3f\n",
            R, sg, ss, ss/sg, std(core_tile(l)), std(core_tile(r)))
end
@info "done — split/global > 1 at small R localizes the excess to coarse vs resid"
