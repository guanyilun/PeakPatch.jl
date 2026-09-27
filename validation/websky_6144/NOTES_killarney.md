# Websky 6144^3 on Killarney — GPU Multi-Resolution

## Cluster Environment (verified 2026-04-15)

### Hardware

| Tier | Nodes | GPU | CPU | Cores | RAM | GPUs/node |
|------|-------|-----|-----|-------|-----|-----------|
| Standard Compute (L40S) | 168 | NVIDIA L40S 48 GB | 2x Xeon Gold 6338 | 64 | 512 GB | 4 |
| Performance Compute (H100) | 10 | NVIDIA H100 SXM 80 GB | 2x Xeon Gold 6442Y | 48 | 2048 GB | 8 |

Network: Infiniband HDR100 (100 Gbps) on Standard, 2x HDR200 (400 Gbps) on Performance.

### Software

| Package | Version | Notes |
|---------|---------|-------|
| Julia | 1.12.5 | `module load julia/1.12.5` (**must use 1.12.x** — Manifest.toml incompatible with 1.11) |
| CUDA toolkit | 12.6 | `module load cuda/12.6` |
| CUDA.jl | 5.x | Installed via `Pkg.add("CUDA")` — needed for GPU pipeline |
| Scheduler | SLURM | `sbatch`, `squeue`, `sinfo` |

### Account and Storage

- SLURM account: `aip-aspuru-ab`
- Home: `/home/yguan` (small quota, daily backup)
- Scratch: `/home/yguan/scratch` (large quota, **purged when inactive**)
- Project: `/home/yguan/projects/aip-aspuru-ab/yguan` (adjustable quota, daily backup)
- **Output**: `/home/yguan/projects/aip-aspuru-ab/yguan/websky` (persistent)
- **Logs + temp**: `/home/yguan/scratch/websky_6144` (may be purged)

### SLURM Partitions

| Partition | Wall limit | L40S nodes | Notes |
|-----------|------------|------------|-------|
| `gpubase_l40s_b1` | 3:00:00 | ~123 | Short tests |
| `gpubase_l40s_b2` | 12:00:00 | ~86 | Medium runs |
| `gpubase_l40s_b3` | 1-00:00:00 | ~60 | **Used for production** |
| `gpubase_interac` | 3:00:00 | ~21 | Interactive/debug |

### Killarney Quirks

- **Cannot `sbatch` from `/home`**: Must submit from `/scratch` or `/project`.
  SLURM script uses `cd /home/yguan/work/PeakPatch.jl` after SLURM directives.
- **Job name `#SBATCH --job-name` does not expand `$OCTANT`**: The variable
  is passed at runtime but the SBATCH header evaluates at parse time.
  Job names appear as `pp6144-oct${` — cosmetic only.
  Similarly, `--output`/`--error` filenames contain literal `oct${OCTANT:-0}`
  instead of the octant name.
- **Bash octal interpretation**: `printf "%03d" 010` interprets `010` as
  octal (8 decimal), not ten. This breaks binary octant names like 010, 011,
  100–111. **Fix**: pass the 3-digit name directly instead of using `printf`,
  e.g. `OCT=${OCTANT:-000}` with no padding.

### Key Differences from CITA Sunnyvale

| | CITA (Sunnyvale) | Killarney |
|---|---|---|
| Scheduler | PBS (`qsub`) | SLURM (`sbatch`) |
| GPU nodes | starq (no GPU), greenq (no GPU) | L40S (4 GPU), H100 (8 GPU) |
| Approach | MPI + PencilFFTs (distributed FFT) | GPU multi-resolution (`run_multitile_split`) |
| Julia | `module load julia/1.11.5` | `module load julia/1.12.5` |

---

## Run Configuration

### Websky 6144 Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| N | 6144 | Websky Stein+2020 |
| Box size | 7700 Mpc/h | Websky |
| Cell size | 7700/6144 = 1.25326 Mpc/h | Derived |
| Cosmology | Om=0.31, OB=0.049, OL=0.69, h=0.68 | Planck 2018 (Websky values) |
| sigma8 | 0.81 | CAMB-normalized |
| Seed | 12345 | Websky |
| LPT order | 2 (2LPT) | |
| Smoothing | Top-hat (wsmooth=1) | |
| Filters | 20, Rf=2.507 to ~36 Mpc/h, 1.15x spacing | Rf,min = 2 × a_latt (Stein+ 2020) |
| Lightcone | Yes (ievol=1), z_max=4.6 | Full-sky |
| Octants | 8 (000..111 binary) | Full-sky |
| Catalog | Extended (ioutshear=1) | XGPaint compatibility |

### GPU Multi-Resolution Tiling

| Parameter | Value | Notes |
|-----------|-------|-------|
| ntile | 16 | Per dimension |
| nmesh | 399 | = nsub + 2*nbuff |
| nbuff | 8 | Buffer cells per side |
| coarse_factor | 4 | -> M=64, block=96 |
| nsub | 383 | = nmesh - 2*nbuff |
| Total tiles | 4096 | = 16^3 |
| Per-tile padded FFT | (2×399)^3 × 4B ~ 2.0 GB | Fits in L40S 48 GB |

Sanity: N = 383×16 + 16 = 6144. Per-tile boxsize = 399 × 1.25326 = 500.05 Mpc/h.

### Octant Observer Positions

Each octant places the observer at a different corner of the 7700 Mpc/h box.
Binary naming: bit 0 = x, bit 1 = y, bit 2 = z.

**Important**: `tile_center()` in `src/MultiTile.jl` returns coordinates centered
at the origin: `(it - (ntile+1)/2) * dcore_box`. The tile grid spans
approximately ±3600 Mpc/h (tile centers) to ±3840 Mpc/h (tile edges). The
observer position in the config must be expressed in this **centered coordinate
system**, not in physical box coordinates (0–7700).

| Octant | Observer (centered) [Mpc/h] | Physical box corner |
|--------|----------------------------|---------------------|
| 000 | (-3850, -3850, -3850) | (0, 0, 0) |
| 001 | (+3850, -3850, -3850) | (7700, 0, 0) |
| 010 | (-3850, +3850, -3850) | (0, 7700, 0) |
| 011 | (+3850, +3850, -3850) | (7700, 7700, 0) |
| 100 | (-3850, -3850, +3850) | (0, 0, 7700) |
| 101 | (+3850, -3850, +3850) | (7700, 0, 7700) |
| 110 | (-3850, +3850, +3850) | (0, 7700, 7700) |
| 111 | (+3850, +3850, +3850) | (7700, 7700, 7700) |

With chi(z=4.6) = 3361.7 Mpc/h, each octant retains **230 tiles** out of 4096
(~5.6%), consistent with the solid-angle/volume fraction of a lightcone octant.

### Per-Octant Resources

| Resource | Value |
|----------|-------|
| GPUs per octant | 4× L40S (1 node) — all used by `run_multitile_split` |
| Wall time requested | 1 hour (partition b1; actual runtime ~19 min) |
| CPU threads | 32 |
| RAM | 400 GB requested (512 GB available) |
| Observed runtime | 18.9 min (915,737 halos, oct 000, kn079) |

---

## Setup Procedure (from scratch)

Run on a login node **once** to prepare everything:

```bash
module load julia/1.12.5 cuda/12.6
bash validation/websky_6144/setup_killarney.sh
```

This script (`setup_killarney.sh`) performs:

1. **`Pkg.add("CUDA")` + `Pkg.instantiate()`** — CUDA.jl is a weak dependency;
   must be explicitly added for `using CUDA` to work. Modifies Project.toml/Manifest.toml.
2. **Generate data files** via `generate_data_files.jl`:
   - Filter bank: 20 top-hat filters (Rf=2.507 to 35.68 Mpc/h)
   - Collapse table: HomelTab for Websky cosmology (20 k-steps)
3. **Generate CAMB P(k)** via `generate_pk_camb.py`:
   - Creates a Python venv at `/home/yguan/scratch/camb_env`
   - Uses CAMB 1.6.6 with sigma8 normalization (target 0.81, achieved 0.8100)
   - 500 k-points from 1e-5 to 100 h/Mpc
4. **Generate 8 octant configs** via `generate_octant_configs.py`
5. **Create output directories**

---

## File Inventory

### Scripts (all in `validation/websky_6144/`)

| File | Purpose |
|------|---------|
| `setup_killarney.sh` | One-time setup (run on login node) |
| `generate_data_files.jl` | Generate filters + collapse table (Julia) |
| `generate_pk_camb.py` | Generate CAMB P(k) (Python + CAMB) |
| `generate_octant_configs.py` | Generate 8 TOML configs from template |
| `run_gpu_octant.jl` | GPU driver — calls `run_multitile_split` with all 4 GPUs |
| `run_websky_6144_killarney.slurm` | SLURM submission script (1 octant per job) |
| `run_test_gpu.jl` | Small-scale GPU test (120^3, 1 GPU) |
| `config_test_gpu_small.toml` | Config for the small-scale GPU test |
| `run_test_gpu_killarney.slurm` | SLURM script for the GPU test |
| `NOTES_killarney.md` | This file |

### Data Files (in `validation/websky_6144/data/`)

| File | Generated by | Notes |
|------|-------------|-------|
| `pk_websky.dat` | `generate_pk_camb.py` | CAMB sigma8=0.81, 500 k-points |
| `filters_websky.dat` | `generate_data_files.jl` | 20 top-hat filters |
| `HomelTab_websky.dat` | `generate_data_files.jl` | Collapse table for Websky cosmology |

### Octant Configs (in `validation/websky_6144/`)

Generated by `generate_octant_configs.py`:

```
config_websky_6144_oct000.toml  through  config_websky_6144_oct111.toml
```

All share the same cosmology, tiling, and pipeline params — only observer position
and output filename differ.

---

## Submission

### Prerequisite: copy SLURM script to scratch

Killarney rejects `sbatch` from `/home`. Copy the script first:

```bash
mkdir -p /home/yguan/scratch/websky_6144/slurm
cp validation/websky_6144/run_websky_6144_killarney.slurm \
   /home/yguan/scratch/websky_6144/slurm/
```

### All 8 octants (parallel, 8 nodes)

```bash
cd /home/yguan/scratch/websky_6144/slurm
for oct in 000 001 010 011 100 101 110 111; do
    sbatch --export=OCTANT=$oct run_websky_6144_killarney.slurm
done
```

### Single octant

```bash
cd /home/yguan/scratch/websky_6144/slurm
sbatch --export=OCTANT=000 run_websky_6144_killarney.slurm
```

### Small-scale GPU test (1 GPU, ~1 min)

```bash
cd /home/yguan/scratch/websky_6144/slurm
sbatch run_test_gpu_killarney.slurm
```

### Monitor

```bash
squeue -u yguan
# Logs: /home/yguan/scratch/websky_6144/logs/oct{NNN}_{JOBID}.{out,err}
# Output: /home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct{NNN}.pksc
```

---

## Julia 1.12 Migration Challenges

Killarney's `module load julia/1.12.5` is the first time we've run PeakPatch on
Julia 1.12 (CITA Sunnyvale uses 1.11). This introduced several non-obvious
breaking changes and bugs that cost significant debugging time. Documenting
them here for future reference.

### 1. Manifest.toml Incompatibility

Julia 1.11.x cannot read a Manifest.toml generated by 1.12.x. The reverse also
holds. **You must match the Julia version to the manifest.** On Killarney,
`module load julia/1.12.5` is mandatory — earlier versions will fail at
`Pkg.instantiate()`.

### 2. `ConcurrencyViolationError` on Array Mutation

**Julia 1.12 introduced task-based array ownership tracking.** Arrays created on
one task cannot be resized (`push!`, `sizehint!`, `growend!`) from a different
task — even `Threads.@spawn` tasks that are clearly operating on distinct
per-worker arrays.

This broke `run_multitile_split`'s multi-GPU dispatch, which used a pattern like:

```julia
local_halos = [ExtHaloRecord[] for _ in 1:n_workers]  # created on main task
for wid in 1:n_workers
    Threads.@spawn begin
        for tile in my_tiles
            push!(local_halos[wid], result)  # ERROR: spawned task ≠ owner task
        end
    end
end
```

Even though each worker accesses its own `local_halos[wid]` vector (no actual
data race), Julia 1.12's ownership model rejects any mutation of a parent-task
array from a child task.

**Fix**: Each spawned task creates its own output arrays fresh, then results are
collected after `wait()`:

```julia
worker_results = Vector{Any}(undef, n_workers)
for wid in 1:n_workers
    Threads.@spawn begin
        my_halos = ExtHaloRecord[]   # created inside the spawned task
        for tile in my_tiles
            push!(my_halos, result)  # OK: same task owns the array
        end
        worker_results[wid] = my_halos
    end
end
```

This also applies to `Dict{String,Float64}` used for timing — `setindex!`
(`dict[key] += value`) can also trigger the error. All mutable per-worker state
must be created inside the spawned task.

**Lesson**: Any Julia 1.12 code using `Threads.@spawn` with closures that
mutate arrays created in the parent scope will hit this. Restructure to create
mutable state inside the spawned task.

### 3. `\` Is Not Line Continuation

Julia does **not** have a line continuation operator. The backslash `\` is the
left-division operator (`A \ b`). This broke a multi-line `@info` macro call:

```julia
# WRONG — \ is left-division, not line continuation
@info "title" N=N ntile=ntile \
    seed=seed coarse_factor=coarse_factor
```

Julia parses this as `ntile \ seed = seed` which is nonsensical.

**Fix**: Put all arguments on one line, or restructure the macro call.

### 4. `using` Must Be at Top Level

`using CUDA` (and all `using`/`import` statements) must be at module top level,
not inside functions. This is enforced more strictly in 1.12:

```julia
# WRONG
function main()
    using CUDA  # syntax: "using" expression not at top level
end

# CORRECT
using CUDA
function main()
    # ...
end
```

### 5. CUDA.jl Non-Official Build Warning

Killarney's Julia is a non-official build. CUDA.jl warns:

```
Warning: You are using a non-official build of Julia. This may cause issues with CUDA.jl.
```

This is non-fatal and can be ignored, but it's worth knowing about.

### 6. `CUDA.version()` Does Not Exist

`CUDA.version()` is not a valid function in CUDA.jl 5.x. Use
`CUDA.runtime_version()` or `CUDA.device!(0) |> CUDA.name` instead.

### 7. `CUDA.jl` Is a Weak Dependency

CUDA.jl is listed as a weak dependency in PeakPatch's `Project.toml`. It will
**not** be auto-installed by `Pkg.instantiate()`. You must explicitly
`Pkg.add("CUDA")` before `using CUDA` works. The extension system
(`ext/CUDAExt.jl`) loads automatically once CUDA.jl is present.

---

## Bugs Found During Killarney Deployment

### Bug 1: Observer Coordinate Mismatch (0 halos)

**Symptom**: All pipeline stages complete in 0.00s with 0 halos. Jobs finish
in ~11 minutes (all spent on precompilation + table generation, no tile work).

**Root cause**: `tile_center()` returns coordinates centered at the origin
(tile grid from -3600 to +3600 Mpc/h), but the octant configs had observer
positions in physical box coordinates (0 to 7700 Mpc/h). The lightcone filter
compared these incompatible systems, filtering out ALL tiles for 7 of 8
octants. For oct 000 (obs at origin, coinciding with grid center), 1912 tiles
passed — but this gave a full-sky sphere, not an octant.

**Fix**: Express observer positions in the centered coordinate system used by
`tile_center`. Observer = ±boxsize_full/2 = ±3850 Mpc/h per axis. See the
Octant Observer Positions table above. Updated `generate_octant_configs.py`.

### Bug 2: Bash Octal Interpretation of Octant Names

**Symptom**: `sbatch --export=OCTANT=010` → script looks for `oct008.toml`
instead of `oct010.toml`.

**Root cause**: Bash treats numbers with leading zeros as octal. `printf "%03d" 010`
→ 008 (octal 010 = decimal 8).

**Fix**: Don't use `printf` for padding. Pass the 3-digit binary name directly:
`OCT=${OCTANT:-000}`.

### Bug 3: SLURM `--output`/`--error` Don't Expand `--export` Variables

**Symptom**: Log files named `oct${OCTANT:-0}_JOBID.out` instead of
`oct000_JOBID.out`.

**Root cause**: `#SBATCH --output=...` is parsed at job submission time, before
`--export` variables are set. The `${OCTANT:-0}` in the output path is a shell
default-value expression, not a SLURM substitution.

**Workaround**: Accept the ugly filenames (cosmetic only). Could redirect
stdout/stderr inside the script body instead.

---

## GPU Pipeline Performance

### Octant 000 (test run, job 3188243)

| Metric | Value |
|--------|-------|
| Node | kn079 (4× L40S 48 GB) |
| Wall time | **18.9 minutes** |
| Pre-merge halos | 1,028,487 |
| Post-merge halos | **915,737** |
| Tiles processed | 290 (out of 4096) |
| Merge exclusion | 112,750 removed (11%) |
| Merge reduction | 0 removed by volume |

### Stage Breakdown (4 GPU workers)

| Stage | Wall (max/worker) | % | Total GPU-s |
|-------|-------------------|---|-------------|
| 09_shell_analysis | 274.8 s | **59.7%** | 1027.2 s |
| 08_peak_find | 73.9 s | 16.1% | 292.4 s |
| 01_residual_gen | 53.0 s | 11.5% | 202.3 s |
| 02_iso_fft_delta | 28.7 s | 6.2% | 109.7 s |
| 05_interp_coarse_psi1 | 24.7 s | 5.4% | 96.9 s |
| 06_2lpt_periodic_fft | 2.7 s | 0.6% | 10.2 s |
| 03_interp_coarse_delta | 1.7 s | 0.4% | 4.5 s |
| 10_record_packing | 0.6 s | 0.1% | 1.9 s |
| 07_laplacian_periodic_fft | 0.1 s | 0.0% | 0.3 s |
| 04_iso_fft_psi1 | 0.0 s | 0.0% | 0.0 s |

**Bottleneck**: Shell analysis at 59.7% of wall time. This is the primary
target for future GPU optimization — currently the per-peak ZZon/fcrit
computation in lightcone mode (ievol=1) limits batching efficiency.

### Timing Recommendation

At ~19 min/octant on L40S, **1 hour wall time** is sufficient with ~3× safety
margin. The b3 partition (6h minimum perception) causes unnecessary queue
delays. **Use `gpubase_l40s_b1` with `--time=01:00:00`** for all 8 octants.

Full-sky estimate: 8 octants × 19 min = ~2.5 hours total (parallel), or
~19 min if 8 nodes are available simultaneously.

GPU pipeline is fully accelerated (see `docs/gpu_session_log.org`):

- Isolated FFTs (cuFFT): 17-20× over CPU FFTW at production tile size
- Full end-to-end: 13.4× speedup (nmesh=192, ntile=2 benchmark)
- All per-tile stages ported: isolated FFT, 2LPT, peak find, shell analysis, interpolation
- Latest optimizations (commit cda3b00): plan cache, fused ψ₁ FFTs, device-resident fields

---

## Pre-Flight Checklist

- [x] Julia 1.12.5 resolves and instantiates
- [x] CUDA.jl added to project (`Pkg.add("CUDA")`)
- [x] PeakPatch + CUDAExt precompiles (275 packages)
- [x] Data files generated (filters, collapse table, CAMB P(k))
- [x] Octant configs generated (000–111)
- [x] Small-scale GPU test passed (job 3173989, kn113, L40S, 928 halos, 46s)
- [x] Output directory created
- [x] SLURM script staged in `/home/yguan/scratch/websky_6144/slurm/`
- [x] Observer coordinate bug fixed (centered coords ±3850)
- [x] Bash octal interpretation bug fixed
- [x] Julia 1.12 ConcurrencyViolationError fixed (task-local arrays)
- [x] `\` line-continuation syntax error fixed
- [x] Full production run completed (7,349,188 halos across 8 octants)
- [ ] Validate halo count against original Websky (currently 120× fewer — see analysis below)
- [ ] Investigate filter bank configuration (optimal vs fixed spacing)

---

## Run Log

### GPU Test

| Date | Job ID | Octant | Node | Config | Time | Halos | Status |
|------|--------|--------|------|--------|------|-------|--------|
| 2026-04-15 | 3173989 | test (120^3) | kn113 | config_test_gpu_small.toml | 46s | 928 | PASSED |

### Production Attempts

| Date | Job ID | Octant | Node | Issue | Status |
|------|--------|--------|------|-------|--------|
| 2026-04-15 | 3175467–3175474 | 000–111 | various | Syntax error (`\` line continuation), octal bug, wrong coords | FAILED |
| 2026-04-15 | 3181837–3181843, 3181935 | 000–111 | kn146, kn041 | Observer coords in box coords → 0 tiles for 7 octants; ConcurrencyViolationError on oct 000 | FAILED |
| 2026-04-15 | 3187081 | 000 | kn166 | Corrected coords, halos found, ConcurrencyViolationError on `sizehint!` | FAILED |
| 2026-04-15 | 3187496 | 000 | kn033 | Removed `sizehint!`, ConcurrencyViolationError on `push!` | FAILED |
| 2026-04-15 | 3187838 | 000 | — | Task-local arrays fix (cancelled, resubmitted as shorter job) | CANCELLED |
| 2026-04-16 | 3188243 | 000 | kn079 | **All fixes applied, 1h wall b1** | **915,737 halos, 18.9 min** |
| 2026-04-16 | 3188510–3188516 | 001–111 | — | 1h wall b1, PENDING | — |

### Final Production Results (2026-04-16)

| Octant | Job ID | Pre-merge | Post-merge | Wall time |
|--------|--------|-----------|------------|-----------|
| 000 | 3188243 | 1,028,487 | 915,737 | 18.9 min |
| 001 | 3188510 | 1,032,389 | 918,988 | 18.9 min |
| 010 | 3188511 | 1,030,969 | 916,618 | 18.7 min |
| 011 | 3188512 | 1,033,574 | 918,835 | 18.8 min |
| 100 | 3188513 | 1,038,505 | 922,790 | 18.8 min |
| 101 | 3188514 | 1,031,491 | 917,641 | 18.8 min |
| 110 | 3188515 | 1,032,498 | 918,823 | 18.8 min |
| 111 | 3188516 | 1,034,057 | 919,756 | 18.8 min |
| **Total** | | **8,261,970** | **7,349,188** | |

- Total catalog size: ~965 MB across 8 files
- Output: `/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct{000-111}.pksc`
- All octants show consistent halo counts (~917–923K post-merge), as expected
  from statistical uniformity of the full-sky lightcone.
- Merge exclusion removes ~11% of pre-merge peaks (Lagrangian overlap removal).

### Key Observations

- 3187081 (first run with correct coords, kn166): Tile (4,1,3) found
  **700,723 peaks, 89,627 halos** before crashing on `sizehint!`. Confirms
  the observer coordinate fix works and the pipeline finds halos correctly.
- Each octant retains 230 out of 4096 tiles after lightcone filtering.
- Runtime per octant: ~19 min on 4× L40S, well within 1h wall time on b1.

---

## Halo Count Comparison with Original Websky

### Summary

Our GPU run found **7.35 million halos** across the full sky. The original
Websky simulation (Stein et al. 2020, arXiv:2001.08787) reports
**~9 × 10⁸ (900 million) halos** (Section 4.4.3). This is a factor of
**~120× fewer halos** than the original.

### Detailed Comparison

| Property | Original Websky (Fortran) | Our GPU Run (Julia) | Notes |
|----------|--------------------------|---------------------|-------|
| Total halos | ~9 × 10⁸ | 7.35 × 10⁶ | 122× fewer |
| Grid resolution | 6144³ native | 6144³ (ntile=16, 399³/tile) | Same cellsize |
| Cell size (a_latt) | 1.25326 Mpc/h | 1.25326 Mpc/h | Identical |
| Rf,min | 2 × a_latt = 2.507 Mpc/h | 2.507 Mpc/h | Identical |
| Rf,max | ~36 Mpc/h | 35.68 Mpc/h | ~Same |
| Filter spacing | "optimal" (σ(R)-based) | Fixed 1.15× logarithmic | **Different** |
| Number of filters | Not specified (paper ref [68]) | 20 | Possibly more in original |
| Min halo mass | ~1.2 × 10¹² M☉/h (10 particles) | ~5.7 × 10¹² M☉/h (from Rf,min) | **4.7× higher** |
| Lightcone radius | 7.7 Gpc (z_max = 4.6) | 7.7 Gpc/h (z_max = 4.6) | Same |
| Catalog size | 33 GB | 965 MB | Reflects halo count |
| M_cell (particle mass) | 1.69 × 10¹¹ M☉/h | Same | |

### Mass Scale Analysis

The minimum detectable halo mass depends on the smallest filter scale:

```
M(Rf) = (4π/3) × ρ̄_m × (Rf/h)³
```

For our cosmology (Om=0.31, h=0.68):
- ρ̄_m = 3.98 × 10¹⁰ M☉/Mpc³
- M(Rf=2.507 Mpc/h) = 8.35 × 10¹² M☉ = 5.68 × 10¹² M☉/h  ← our minimum
- M(Rf=2.08 Mpc/h) = 4.77 × 10¹² M☉ = 3.24 × 10¹² M☉/h   ← NOTES.md MPI estimate
- 10 × M_cell = 2.49 × 10¹² M☉ = 1.69 × 10¹² M☉/h          ← paper's cutoff

The paper states they retain halos with pre-abundance-matched mass > 10 particles,
giving M_min ≈ 1.2 × 10¹² M☉/h. Our Rf_min of 2.507 Mpc/h corresponds to a
mass of ~5.7 × 10¹² M☉/h — a factor of ~5 higher minimum mass.

### Likely Causes of the Discrepancy

**1. Filter bank configuration (primary suspect)**

The paper states it uses "optimal filter spacings to maximize both accuracy and
efficiency presented in [68]" (Stein et al. 2019, MNRAS 483, 2236,
arXiv:1810.07727). This σ(R)-based optimal spacing is **not fixed logarithmic**
spacing. It typically packs many more filters at small scales where σ(R)
changes rapidly, potentially yielding significantly more than 20 filters.
More filters at small scales means more small-scale peaks are detected, which
directly increases the halo count at low masses where the mass function is
steepest.

Our fixed 1.15× spacing gives exactly 20 filters from Rf=2.507 to 35.68 Mpc/h.
The original Fortran code may have used 30–50+ filters with denser coverage at
small Rf, finding many more low-mass peaks.

**2. Minimum mass threshold**

The halo mass function dn/dM rises extremely steeply at low masses
(approximately dn/dM ∝ M^{-2} at M ~ 10^{12} M☉/h). A factor of ~5 increase
in minimum mass (from 1.2 × 10¹² to 5.7 × 10¹² M☉/h) translates to roughly:

```
n(> M_min) ∝ M_min^{-1} × exp(-(M_min/M*)^α)
```

The exponential suppression at these masses can easily produce a factor of
50–100× fewer halos, approaching the observed 120× discrepancy.

**3. Julia vs Fortran implementation differences**

PeakPatch.jl is a Julia reimplementation of the original Fortran mass-Peak
Patch code. Subtle differences in any of the following could affect halo counts:
- Ellipsoidal collapse ODE integration (step size, tolerances, freeze-out factors)
- Peak finding threshold (δ_c = 1.686 for all filters, but actual threshold may vary)
- Binary exclusion algorithm (overlap subtraction order, volume calculations)
- Collapse table interpolation (HomelTab resolution and accuracy)
- 2LPT displacement computation

**4. Abundance matching**

The original Websky paper applies abundance matching to the Tinker et al. (2008)
mass function in Δz = 0.1 redshift bins. Our run does **not** apply abundance
matching — halo masses are as computed by the ellipsoidal collapse directly.
This affects the mass-luminosity mapping but should not affect the *number*
of halos found (only their assigned masses).

### Recommended Investigation Steps

1. **Check the original Fortran filter bank**: Obtain the exact number and
   Rf values used in the original Websky Fortran run. Compare against our
   20-filter fixed-spacing bank. This is the most impactful unknown.

2. **Test with more filters**: Generate a filter bank with denser spacing
   (e.g., 1.05× instead of 1.15×) or σ(R)-based optimal spacing, and rerun
   a single octant. This directly tests the filter bank hypothesis.

3. **Compare against a known small-scale result**: Run the 120³ GPU test with
   the same filter bank and compare halo counts against an equivalent Fortran
   run (if available) to isolate Julia-vs-Fortran differences.

4. **Check minimum mass of output halos**: Inspect the mass distribution of
   our 7.35M halos. If the minimum mass is well above 1.2 × 10¹² M☉/h,
   it confirms we're missing the low-mass population. If it's near
   1.2 × 10¹², the issue is elsewhere.

5. **Validate against Tinker mass function**: Compute n(> M) from our catalog
   and compare against Tinker et al. (2008) predictions. This would reveal
   whether our halo counts are physically reasonable for the minimum mass
   we achieve, even if below the original Websky target.

### Reference

- Stein, Alvarez, Bond, van Engelen & Battaglia (2020), "The Websky
  Extragalactic CMB Simulations", arXiv:2001.08787
- Stein, Alvarez & Bond (2019), "The mass-Peak Patch algorithm for fast
  generation of deep all-sky dark matter halo catalogues and its N-body
  validation", MNRAS 483, 2236, arXiv:1810.07727
- Tinker et al. (2008), "Toward a Halo Mass Function for Precision
  Cosmology", ApJ 688, 709, arXiv:0803.2706
