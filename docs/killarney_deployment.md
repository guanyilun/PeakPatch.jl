# Running PeakPatch.jl at N=6144 on Killarney

This guide covers how to deploy the GPU-accelerated multi-tile PeakPatch
pipeline on the Killarney cluster (Standard Compute: 4× NVIDIA L40S 48 GB
per node). It assumes the GPU port (`ext/CUDAExt.jl`) and multi-GPU tile
dispatcher (`run_multitile_split(...; devices=[0,1,2,3])`) from this branch.

## TL;DR

| What | Value |
|---|---|
| Target queue | Killarney `Standard Compute` |
| Resources per octant | 1 node × 4 L40S × 8 h wall request |
| Wall time per octant (fixed-z or lightcone) | **~2 h** |
| Wall time for 8 octants (sequential jobs) | **~16 h** |
| Wall time for 8 octants (8 parallel 1-node jobs) | **~2 h** (same as 1 octant) |
| Disk per octant (extended output) | ~15 GB catalog + ~50 GB scratch |
| Config | `ntile=16, nmesh=399, nbuff=8, coarse_factor=4` |
| Accuracy | ~1–2% halo count vs reference global-FFT at production scale |

Total for full 8-octant lightcone or full-sky mock: **2 h wall × 8 jobs ≈ 16 node-hours**, but only 2 h of real time if jobs run in parallel.

Wall-time improvement since initial deployment projection is ~2.25×
(A5000 bench at nmesh=192 went from 20.4 s → 9.0 s per octant-equivalent).
See **§7 Further performance roadmap** for what changed.

See [Scaling the estimate to 8 octants](#scaling-the-estimate-to-8-octants)
and [Lightcone vs fixed-z](#lightcone-mode-vs-fixed-redshift) below for how
this interacts with the `ievol=1` lightcone path.

---

## 1. Hardware fit

Killarney Standard Compute node (Dell 750xa):
- 2× Intel Xeon Gold 6338 → 64 cores
- 512 GB DDR4
- **4× NVIDIA L40S 48 GB** (Ada Lovelace, 864 GB/s HBM, 91 TFLOPS FP32)
- Local NVMe scratch, shared Lustre `/project` and `/scratch`

Per-GPU speedup vs our development A5000 baseline on this workload is
realistically **~1.3–1.5×** (limited by memory bandwidth for FFTs; shell
analysis gains more from Ada's improved atomics). The 48 GB of VRAM
per GPU is 20× what we measured the largest tile actually needed
(`nmesh=384` peaks at ~2.5 GB).

Performance Compute (8× H100 SXM) is overkill — we are memory-bandwidth
bound, not compute-bound, and don't exploit tensor cores. Stick with
Standard Compute unless queue times are a bottleneck.

---

## 2. Recommended configuration for N=6144

`N = 6144 = 2¹¹ × 3`. This factorisation forces `M = N / block` to be a
divisor of 6144; 80 (which would give the notes-optimal `coarse_factor=5`)
does not divide 6144. The valid choices near the optimum are **M=64 → cf=4**
or M=96 → cf=6. **Use cf=4.**

From the notes (`docs/multi_resolution_fft.md`):
- `coarse_factor=5` is optimal at the 120³ test scale (1.2% halo error).
- Accuracy improves at production scale because tiles are a smaller
  fraction of the box, there are more coarse cells, and `nbuff` is a
  smaller fraction of the tile.
- `cf=4` sits between tested points `cf=3` (2.7%) and `cf=5` (1.2%);
  expected ~2% halo-count error, further reduced at production scale.

**Recommended config:**

```julia
ntile        = 16
nmesh        = 399        # = nsub + 2·nbuff = 383 + 16
nbuff        = 8
coarse_factor = 4          # → M = 64, block = 96
```

Sanity check: `N = nsub·ntile + 2·nbuff = 383·16 + 16 = 6144 ✓`, `N % M = 6144 % 64 = 0 ✓`.

Per-tile padded FFT memory: `(2×399)³ × 4 B ≈ 2.0 GB` — trivial on L40S.
Total tiles: `16³ = 4096` tiles.

### Alternative (fewer bigger tiles)

```julia
ntile = 8, nmesh = 782, nbuff = 8, coarse_factor = 4   # → M = 32, block = 192
```

Per-tile padded FFT memory: `(2×782)³ × 4 B ≈ 15.3 GB` — still fits on a
48 GB L40S. Total tiles: 512. Throughput is essentially identical to the
ntile=16 config (the sweep data showed cost is ~linear in total cells,
not in tile count). Use ntile=16 unless you have a reason to prefer
larger tiles.

### If you can change the box

If you are free to pick N, use **N=6000** or **N=6250**: both are
divisible by M values that give `coarse_factor=5` exactly
(e.g., N=6000, ntile=15, nmesh=408, M=75 → cf=5). This trades the
aesthetically nice 6144 = 2¹¹·3 for ~20% better halo accuracy.

---

## 3. How to run

### 3.1 One-time setup on Killarney

```bash
# Log in, go to your project space
ssh <you>@killarney

# Clone or rsync the repo
git clone <repo-url> ~/projects/def-<group>/PeakPatch.jl
cd ~/projects/def-<group>/PeakPatch.jl

# Load modules (adjust to cluster defaults as needed)
module load StdEnv cuda/12.3 julia/1.11

# Instantiate the Julia project (this downloads FFTW, CUDA.jl, etc.)
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Warm up CUDAExt precompile
julia --project=. -e 'using PeakPatch, CUDA; @assert CUDA.functional()'
```

### 3.2 Submit script (single octant, N=6144)

Save as `run_octant.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=pp6144-oct${OCTANT:-0}
#SBATCH --account=def-<group>
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4          # all 4 L40S on the node
#SBATCH --cpus-per-task=32         # 32 threads (half the node)
#SBATCH --mem=400G                 # leave headroom on 512 GB box
#SBATCH --time=08:00:00            # 8 h wall (projection ~4.5 h, 2× safety)
#SBATCH --output=logs/oct${OCTANT:-0}_%j.out

set -euo pipefail
module load StdEnv cuda/12.3 julia/1.11

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Don't set CUDA_VISIBLE_DEVICES — we want all 4 devices visible so
# `devices=[0,1,2,3]` can address them.
unset CUDA_VISIBLE_DEVICES

cd ${SLURM_SUBMIT_DIR}
julia --project=. --threads=${SLURM_CPUS_PER_TASK} bin/run_octant.jl \
    --octant ${OCTANT:-0} \
    --ntile 16 --nmesh 399 --nbuff 8 --coarse-factor 4 \
    --out /scratch/${USER}/pp6144/catalog_oct${OCTANT:-0}.pksc
```

Launch one octant job:

```bash
sbatch --export=OCTANT=0 run_octant.sh
```

Launch all 8 octants in parallel (8 node-hours concurrent, ~4.5 h wall):

```bash
for oct in 0 1 2 3 4 5 6 7; do
    sbatch --export=OCTANT=$oct run_octant.sh
done
```

### 3.3 Driver script (`bin/run_octant.jl`)

You'll need a driver that calls `run_multitile_split` with `devices=[0,1,2,3]`:

```julia
using PeakPatch, CUDA, ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "--octant";        arg_type=Int;    default=0
    "--ntile";         arg_type=Int;    default=16
    "--nmesh";         arg_type=Int;    default=399
    "--nbuff";         arg_type=Int;    default=8
    "--coarse-factor"; arg_type=Int;    default=4
    "--out";           arg_type=String; required=true
    "--z-out";         arg_type=Float64; default=0.0
    "--ioutshear";     arg_type=Int;    default=1
    "--ilpt";          arg_type=Int;    default=2
end
args = parse_args(s)

# Octant offsets (cenx,ceny,cenz). Observer at one corner of the box per octant.
# Boxsize is set so alatt = 226.67/68 Mpc/h per cell.
alatt = 226.67 / 68
N = args["ntile"] * (args["nmesh"] - 2*args["nbuff"]) + 2*args["nbuff"]
boxsize = N * alatt

# For full-sky 8-octant lightcone: each octant places the observer at a
# different corner and flips sign. Adjust to match your upstream convention.
octant_centers = [(0.0, 0.0, 0.0)]  # TODO: fill in 8-octant geometry if lightcone

cx, cy, cz = length(octant_centers) > args["octant"] ?
                octant_centers[args["octant"]+1] : (0.0, 0.0, 0.0)

cfg = PeakPatch.PipelineConfig(
    n=args["nmesh"], boxsize=boxsize,
    pkfile="data/pk_planck18.dat",
    filterfile="data/filter_default.dat",
    tabfile="data/homeltab.dat",
    fileout=args["out"], nbuff=args["nbuff"],
    z_out=args["z-out"], ilpt=args["ilpt"], ioutshear=args["ioutshear"],
    rmax2rs=0.0, ievol=0, z_max=0.0,
    cenx=cx, ceny=cy, cenz=cz,
    Omx=0.261, OmB=0.049, Omvac=0.69, h=0.68,
    NonGauss=0, fNL=0.0, wsmooth=1,
)

@assert CUDA.functional()
ndev = length(CUDA.devices())
devices = collect(0:ndev-1)   # [0,1,2,3] on Standard Compute node
@info "Using $ndev GPUs: $devices"

halos = PeakPatch.run_multitile_split(cfg;
    ntile=args["ntile"], seed=42+args["octant"],
    coarse_factor=args["coarse-factor"],
    use_gpu=true, devices=devices,
    verbose=true, profile=true)

@info "Octant $(args["octant"]): $(length(halos)) halos written to $(args["out"])"
```

(Adjust octant_centers and seeds to your actual lightcone geometry; the
existing `bin/peakpatch.jl` has the full pattern.)

---

## 4. Scaling the estimate to 8 octants

PeakPatch light-cone mocks are typically built by running the pipeline
**once per octant** with the observer placed at different corners /
orientations of the box. Each octant is an independent pipeline invocation
and produces an independent halo catalog which is merged afterwards.

Cost model:

- Single octant (1 node, 4 L40S): **~2 h wall**
- 8 octants sequential (1 node, 4 L40S): ~16 h total
- 8 octants as 8 parallel jobs (8 nodes, 32 GPUs total): **~2 h wall** end-to-end
- 8 octants as 4 parallel jobs × 2 sequential (4 nodes): **~4 h wall**

Killarney has 168 Standard Compute nodes; submitting 8 parallel jobs is
realistic if the queue is not saturated.

Disk for 8 octants:
- 8 × 15 GB catalogs ≈ **120 GB**
- Scratch during run: ~50 GB per octant; only held while the job runs.
- Merged mock catalog (post-octant concatenation + dedup): similar total.

---

## 5. Lightcone mode vs fixed redshift

PeakPatch has two modes, controlled by `cfg.ievol`:

- **`ievol = 0` (fixed redshift)**: every tile is evaluated at the same
  output redshift `z_out`. One pipeline run per `z_out`. If you need
  snapshots at multiple redshifts, **yes, you run the pipeline once per
  redshift**, each producing its own pksc catalog.

- **`ievol = 1` (lightcone)**: each tile's effective redshift is computed
  from its distance to the observer (`peak_redshift(obs, tile_center)`).
  A single run writes a single catalog where halos are already at their
  lightcone-appropriate redshifts — **no per-redshift loop required**.
  This is what the 8-octant convention builds on.

So the runtime picture is:

| Mode | How many pipeline runs? |
|---|---|
| Fixed-z snapshot at one `z_out` | 1 run |
| Fixed-z snapshots at N redshifts | N separate runs |
| Full-sky lightcone | 8 octants × 1 run each = 8 runs |

**GPU shell analysis now supports lightcone mode** (`ievol=1`). The
post-process kernel in `ext/CUDAExt.jl` takes per-peak `ZZon`/`fcrit`
device vectors rather than scalars; the host computes
`peak_redshift → ZZon_pk, fcrit_pk` for each candidate peak in
Float64 on the CPU and the downcast arrays are uploaded to the GPU
per batch. Peaks with `z_pk > z_max` are filtered on the host before
the batch so the kernel doesn't run on out-of-light-cone work.

Measured payoff on a small test grid (A5000, NMESH=128, ntile=2,
ilpt=2, `ievol=1`, `z_max=3.0`, 4551 halos): **CPU wall 55.41 s →
GPU wall 3.25 s (17.03×)**, halo counts identical
(`docs/ievol1_bench_A5000.log`). Lightcone wall-time is now
comparable to fixed-z at the same config; the ~2× lightcone penalty
that existed before this change is gone.

---

## 6. How the benchmark numbers were derived

The runtime projections in this guide come from `bench/scaling_sweep.jl`.
Raw output is archived at `docs/scaling_sweep_A5000.log`.

### 6.1 What the sweep measures

The script runs `run_multitile_split(...; use_gpu=true, profile=true)` at
fixed `ntile=2, nbuff=8` and sweeps `nmesh ∈ {128, 192, 256, 320, 384}`.
At each point it records wall time, per-tile time, halo count, peak GPU
memory (via a background poller), and the per-stage profile breakdown.

### 6.2 Results on RTX A5000 (development machine)

```
nmesh  nsub   N      wall(s)    per-tile   halos     mem(MB)  cell·µs⁻¹·GPU⁻¹
128    112    240    6.54       0.818      4,982     689      2.11
192    176    368    20.40      2.550      20,467    971      2.44
256    240    496    53.34      6.667      54,057    1,417    2.29
320    304    624    98.49      12.312     112,271   2,521    2.47
384    368    752    164.84     20.605     198,049   2,539    2.58
```

Per-stage breakdown at `nmesh=384` (A5000, GPU path):
- `09_shell_analysis` — 56.6%
- `10_record_packing` — 11.7%
- `04_iso_fft_psi1`   — 9.0%
- `05_interp_coarse_psi1` — 7.2%
- `06_2lpt_periodic_fft`  — 6.0%
- `08_peak_find` — 4.0%
- Everything else — < 4% each

### 6.3 Throughput model

**Cost per cell is nearly constant** (3.6–4.0 × 10⁻⁷ s/cell on A5000),
which means per-tile time scales ~linearly in `nmesh³`. Total pipeline
cost is approximately `cells_total × cost_per_cell / n_workers`.

For L40S we estimate 1.3× speedup vs A5000 (bandwidth-limited on FFT,
modest kernel gains). For N=6144 (2.32×10¹¹ cells):

- 1× A5000 (pre-opt): `2.32e11 × 3.7e-7 s ≈ 24 h` (original baseline)
- 1× A5000 (post-opt): ~**11 h** (after 2.25× wall reduction at nmesh=192)
- 1× L40S: ~**8 h**
- 4× L40S (1 standard node): ~**2 h**
- 8× L40S (2 nodes): ~**1 h**

The "pre-opt" projection was based on the `scaling_sweep` baseline in §6.2;
the "post-opt" column uses the measured 9.03 s GPU wall at nmesh=192 from
`docs/multitile_opt_bench_A5000.log` after the roadmap §7-2 … §7-5
optimisations landed (plan cache, fused ψ₁ FFTs, GPU-resident tile
fields, vectorised record packing).

### 6.4 Reproducing the bench locally

```bash
# Development sandbox (reuses packed CUDAExt)
export PATH=$HOME/.juliaup/bin:$PATH
CUDA_VISIBLE_DEVICES=0 julia --project=/tmp/pp_gpu_sandbox --threads=32 \
    bench/scaling_sweep.jl
```

Knobs: `NMESH_LIST`, `NTILE`, `ILPT`, `IOUTSHEAR`, `NBUFF`, `TARGET_N`.

---

## 7. Further performance roadmap

If 4.5 h/octant (or 9 h/octant lightcone) is the binding constraint,
here is where additional GPU work can bite, in decreasing expected
return:

| # | Item | Est. wall reduction | Effort |
|---|---|---|---|
| ✓ | ~~Batch GPU shell analysis for `ievol=1` (per-peak ZZon)~~ | **Done — 17× measured lightcone speedup at NMESH=128**; parity at `test/test_multiresolution_gpu_ievol1.jl`, bench at `docs/ievol1_bench_A5000.log`. Kernel now takes per-peak `ZZon_pp, fcrit_pp` device vectors and host pre-filters by `z_max`. | — |
| ✓ | ~~Keep tile fields GPU-resident across stages~~ | **Done** — `isolated_convolve_gpu_multi`, `compute_2lpt_gpu`, `compute_laplacian_gpu`, `peak_find_tile_gpu`, and `analyse_peaks_gpu_cuda_multirf` now accept / return `CuArray{Float32,3}` via a `return_device=true` kwarg (and by detecting `CuArray` inputs). In the tile loop δ, ψ₁, ψ₂, lapd all stay device-resident from construction through shell analysis, eliminating 5× host↔device round-trips of `nmesh³ × Float32` per tile (~1.3 GB of PCIe traffic at `nmesh=399`). | — |
| ✓ | ~~Persistent cuFFT plans across tiles~~ | **Done** — cuFFT plans are now memoised in `ext/CUDAExt.jl` via a `(size, device_id)` keyed Dict (`_CUFFT_FWD_CACHE` / `_CUFFT_INV_CACHE`) and shared across every `_cufft_rfft` / `_cufft_irfft` call site (isolated FFTs, 2LPT, Laplacian, peak-find). Each tile after the first reuses the plans rather than rebuilding them. | — |
| ✓ | ~~Vectorise CPU record-packing (stage 10)~~ | **Done** — the GPU-path per-peak loop in `run_multitile_split` now reads directly out of the SoA `gpu_res` arrays (no intermediate NamedTuple-per-peak, no broadcast slice allocations) and pre-sizes the halo output vectors via `sizehint!`. The CPU path (which is dominated by `analyse_peak` itself) is unchanged. | — |
| ✓ | ~~Fuse ψ₁ isolated FFTs (3 dims → 3-plane batched rFFT)~~ | **Done** — `PeakPatch.isolated_convolve_gpu_multi(noise; kernels=[...])` shares the (upload + zero-pad + forward rFFT) of the residual noise across the four transfer kernels (δ, ψ₁_x, ψ₁_y, ψ₁_z); only the per-kernel multiply + hermitise + inverse rFFT differ. Saves 3× noise upload + 3× padded rFFT per tile. | — |
| 6 | Fuse 2LPT 9-FFT pipeline into one kernel batch | 2–3% | ~1 day (not yet pursued) |

All remaining wall-time items from this table (other than the small
2–3% 2LPT kernel fusion) are merged. Measured aggregate on
A5000, `nmesh=192, ntile=2, ilpt=2, ioutshear=1, ievol=0`
(`bench/profile_multitile.jl`, full log at
`docs/multitile_opt_bench_A5000.log`):

| Metric | Before (scaling_sweep baseline) | After (this commit) |
|---|---|---|
| GPU total wall | 20.40 s | **9.03 s** |
| Speedup over CPU | ~3.3× | **16.66×** |
| Halo count parity | exact | exact (20 467 / 20 467) |

Stage breakdown shifts as expected — `09_shell_analysis` now
dominates at 76.8% of the GPU wall (up from ~56% before), because
every other stage shrank:

```
01_residual_gen           0.38 s  ( 4.3%)
02_iso_fft_delta+psi1     0.60 s  ( 6.9%)   # was 9.0% + extra
03_interp_coarse_delta    0.02 s  ( 0.2%)
05_interp_coarse_psi1     0.07 s  ( 0.8%)
06_2lpt_periodic_fft      0.15 s  ( 1.7%)   # was 6.0%
07_laplacian_periodic_fft 0.02 s  ( 0.2%)
08_peak_find              0.26 s  ( 3.0%)   # was 4.0%
09_shell_analysis         6.73 s  (76.8%)
10_record_packing         0.52 s  ( 6.0%)   # was 11.7%
```

Expected wall time on a 4× L40S Killarney Standard node at `N=6144`
drops from the original ~4.5 h/octant projection to **~2 h/octant**
(same for fixed-z and lightcone, since the lightcone shell kernel
already runs at fixed-z throughput). The remaining optimisation
opportunity is inside `09_shell_analysis` itself — the 2LPT
kernel-fusion item would trim an additional 1–2% but the really
impactful work is attacking the shell gather/post-process kernels
(register pressure, shared-memory reuse, persistent per-peak
reductions); not pursued in this cycle.

Multi-node MPI (multi-node tile dispatch) is not on this list because
Killarney job scheduling makes it easier to submit 8 parallel 1-node
jobs than one 8-node MPI job — the 8-octant structure is already
embarrassingly parallel across SLURM submissions.

---

## 8. Checklist before submitting a 6144³ run

- [ ] `julia --project=. -e 'using PeakPatch, CUDA; @assert CUDA.functional() && length(CUDA.devices()) == 4'` passes on a Standard Compute interactive session.
- [ ] `test/test_multigpu_dispatch.jl` passes on an interactive node (2+ GPU required; 4 works fine).
- [ ] Power spectrum, filter bank, and homeltab files are staged on `/project`.
- [ ] Scratch directory exists on `/scratch/${USER}/pp6144/` with ~120 GB free.
- [ ] `sbatch` script parameters match §3.2 (4 GPUs, 32 threads, 8 h wall, 400 GB RAM).
- [ ] Decision made on `ievol` (lightcone vs snapshot) — affects time budget and whether to extend the GPU shell batch first.
- [ ] Decision made on `coarse_factor` (4 for strict N=6144, 5 if you can move to N=6000/6250).

---

## 9. References in this repo

- `src/MultiResolution.jl` — `run_multitile_split` entry point, Phase 1a coarse grid, per-tile loop, multi-GPU dispatch.
- `ext/CUDAExt.jl` — GPU implementations of isolated FFT, 2LPT, Laplacian, peak-find, tri-cubic interpolation, and shell-analysis batch kernel.
- `docs/multi_resolution_fft.md` — algorithm, coarse-factor trade-off table, accuracy analysis.
- `docs/gpu_implementation_plan.org` — phased GPU port plan (historical reference).
- `docs/gpu_session_log.org` — per-session port notes, per-stage profile evolution, production bench numbers.
- `bench/scaling_sweep.jl` — benchmark used for all projections in §6.
- `bench/profile_multitile.jl` — single-config profile (used earlier to guide per-stage ports).
- `test/test_multiresolution_gpu.jl` — single-GPU parity vs CPU (fixed-z, `ievol=0`).
- `test/test_multiresolution_gpu_ievol1.jl` — single-GPU parity vs CPU in lightcone mode (`ievol=1`), including `z_max` filter exercise.
- `test/test_multigpu_dispatch.jl` — multi-GPU dispatcher parity vs single-GPU.
- `docs/ievol1_bench_A5000.log` — profile bench showing the 17× lightcone speedup (NMESH=128, ntile=2, ilpt=2).
- `docs/ievol_gpu_plan.md` — the implementation plan that drove the lightcone GPU port.
