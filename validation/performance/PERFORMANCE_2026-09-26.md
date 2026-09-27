# Performance: frozen campaign + controlled benchmarks (2026-09-26)

Material for paper §9 / Table T3 (plan item A7 in `docs/paper_comparison_plan_2026-09.md`).
Hardware: Killarney L40S nodes (4× L40S 48 GB, 64 cores, 512 GB); one H100 80 GB point.
Baseline: Websky (Stein+2020 §4.1): full-sky 8-octant catalog in **3.84 h on 1128 Skylake
cores = 4336 core-h**, **peak memory 7.67 TB**, **5.9 TB of stored initial conditions**,
33 GB catalog.

## 1. Frozen campaign (8 octants, N=6144, box 5236 Mpc/h, cf=32, z_max=4.5)

Reproduce: `bash validation/performance/harvest_production.sh`. The catalog logs have no
per-stage timestamps, so the post-pipeline split comes from file mtimes (raw catalog written
→ merge+finalize+write done; AM catalog written → AM done).

Catalog jobs (1 node, 4× L40S, 32 cores):

| octant | job wall (h) | pipeline (min) | merge+finalize+write (h) | AM (min) | GPU-s (stages) | host MaxRSS (GiB) |
|---|---|---|---|---|---|---|
| 000 | 7.13 | 121.5 | 4.85 | 15.4 | 26 248 | 135 |
| 001 | 6.90 | 121.8 | 4.61 | 15.2 | 26 313 | 135 |
| 010 | 7.08 | 124.5 | 4.75 | 15.5 | 26 565 | 134 |
| 011 | 6.95 | 123.4 | 4.64 | 15.3 | 26 527 | 135 |
| 100 | 10.01 (TIMEOUT) | 126.0 | 7.66 (slow node kn103) | rerun, 16 min | 27 122 | 135 |
| 101 | 7.18 | 124.4 | 4.84 | 15.8 | 26 801 | 134 |
| 110 | 7.06 | 123.5 | 4.74 | 15.4 | 26 731 | 135 |
| 111 | 7.08 | 122.8 | 4.78 | 14.9 | 26 498 | 134 |

Pipeline stage profile (oct000, `profile=true`; per-worker max, 4 workers, 2413 active tiles of 4096):

| stage | max/worker (s) | share |
|---|---|---|
| 01 residual noise generation | 306 | 4.6% |
| 02 isolated FFT δ | 67 | 1.0% |
| 03 coarse→fine δ interpolation | 250 | 3.8% |
| 05 coarse→fine ψ1 interpolation | 118 | 1.8% |
| 06–07 2LPT + Laplacian FFTs | 3 | 0.05% |
| 08 peak finding | 582 | 8.8% |
| **09 shell analysis (GPU batch)** | **5192** | **78.4%** |
| 10 record packing | 108 | 1.6% |
| stages total | 6626 (sum 26 248 GPU-s) | |

Pipeline `elapsed_min` (121.5 min = 7290 s) minus the stage wall (6626 s) ≈ 11 min is
startup/JIT plus the global coarse phase (Phase 1a: block-averaging all 6144³ Threefry
draws + M=512 coarse FFTs, on the CPU).

Fieldmap jobs (1× L40S, 16 cores; 5 kernels κ/mass/τ/kSZ/ISW painted together, Nside 4096,
2413 tiles): **76.9–82.2 min painting**, 1.30–1.38 h job wall, host MaxRSS 21–31 GiB. There is
no per-kernel breakdown, because all kernels are painted in one pass per tile (1.9–2.0 s/tile).

### Where the ~7 h goes (and why "19 min at cf4" is not a comparison)
- **Halo finding is ~2.0 h**: 78% of it is GPU shell analysis (~540k peaks/tile × 20
  filters). The comment at `src/MultiResolution.jl:532-534` saying shell analysis runs on the
  CPU is stale: the production path calls `_gpu_analyse_peaks_batch`.
- **Merging is ~4.6 h, which is 2/3 of the job, with all 4 GPUs idle**: `merge_catalog` →
  `lagrangian_exclusion!` is a serial greedy loop over the 333.7M pre-merge halos
  (→195.5M). Its spatial hash is a fixed `NC = 256` cells per axis
  (`src/Merger/Exclusion.jl:22`). At a 5236 Mpc/h box that makes ~20.9 Mpc/h cells holding
  ~20 halos each, so every halo walks ≥27 linked-list cells scattered over ~8 GB of
  coordinates, which is cache-miss bound. See §3.
- **cf32 is not the cause.** The production-cellsize run at cf=4 (job 3968163, 2026-06-17)
  already took 6h39m against 7h13m at cf=32 (job 4387351). That puts the cf32 overhead at
  ~0.5 h, consistent with stages 03/05 (coarse interpolation, ~6% of the pipeline).
- The "19 min/octant" figure (NOTES_killarney.md, 2026-04) came from the pre-fix
  configuration: 1/h-coarser cells (box 7700 Mpc/h), the `fsc_of_z` threshold bug and
  ~0.9M halos/octant. That is a 200× smaller catalog; do not quote it.
- Also note: `run_gpu_octant.jl`'s last line "Done: … in 121.5 minutes" re-prints the
  *pipeline* time, not the job total.

## 2. Controlled benchmarks (job IDs 5692900, 5692991–5692994, 5693005)

Scripts: `gen_bench_configs.jl` (configs in `configs/`), `bench_octant.jl` (timed driver,
warm-up run first so JIT is excluded), `bench.slurm`, `mem_sweep.slurm`,
`submit_bench.sh`. Logs are in `/home/yguan/scratch/websky_6144/logs/bench_*`, and nvidia-smi
traces in `/home/yguan/scratch/websky_6144/bench/gpumem_*.csv`.

**Strong scaling config `scale`**: production tile (n=414, nbuff=16) and production physics
(cell 0.852 Mpc/h, 2LPT, ioutshear=1, finecell filters, lightcone) with ntile=4 → N=1560
(1330 Mpc/h), 64 tiles, cf=30. z_max=1.016 puts the horizon beyond the far corner, so every
tile is active and the work is fixed. Its per-tile cost is **11.1 GPU-s against 10.9 in
production**, so the benchmark is representative. 8 CPU cores per GPU.

| hardware | GPUs | pipeline (s) | speedup | efficiency | merge (s) | finalize (s) | write (s) | peak GPU mem (GiB) | mean GPU util |
|---|---|---|---|---|---|---|---|---|---|
| L40S | 1 | 709.7 | 1.00 | — | 51.2 | 1.3 | 15.0 | 23.5 | 34% |
| L40S | 2 | 374.6 | 1.89 | 95% | 56.6 | 1.6 | 15.4 | 21.2 | 29% |
| L40S | 4 | 225.0 | 3.15 | 79% | 50.3 | 1.5 | 14.8 | 21.3 | 19–20% |
| H100 80GB | 1 | 419.1 | **1.69× one L40S** | — | 50.9 | 1.2 | 11.3 | 23.7 | 21% |

(25.7M pre-merge → 10.0M halos; H100 job had 6 CPU cores.) The 79% at 4 GPUs is tail
imbalance with only 16 tiles per worker. Production, with ~600 tiles per worker, is 99%
balanced: stage wall 6626 s against sum/4 = 6562 s.

**GPU memory vs tile size** (1× L40S, ntile=2, production physics; device-level
nvidia-smi peak, which includes the CUDA.jl pool cache and the ~0.5 GB context):

| tile n | 256 | 320 | 384 | 414 (prod) | 448 | 512 | 576 |
|---|---|---|---|---|---|---|---|
| peak GPU mem (GiB) | 4.5 | 7.9 | 13.1 | 21.2–23.7 | 20.6 | 30.3 | 44.4 (of 45.0) |
| buffer overhead 1−(nsub/n)³ | 33% | 27% | 23% | 21% | 20% | 18% | 16% |

Memory scales ≈ n³. Production n=414 uses about half an L40S. n=512 fits with headroom,
and n=576 fits only with no headroom; pool-limited, so treat it as the ceiling. The ~21–24
GiB at n=414 is higher than the n=448 sweep point because the `scale` tiles span a wider z
range, so more halos are retained per tile. Moving production to n=512 would trim the
buffer overhead from 21% to 18%, a small gain.

## 3. The merge bottleneck: measured on real data (job 5693005)

`exclusion_scaling.jl` runs Lagrangian exclusion on the 9.81M oct000 halos inside the
1309 Mpc/h observer-corner sub-box, the densest region. It times the `src` loop unchanged,
with the hash resolution varied, against a counting-sorted CSR cell list. Every variant
returns **identical survivor sets** (6 038 367; halos are re-excluded because the file
holds Eulerian positions, which is irrelevant for timing).

| cell (Mpc/h) | src linked-list (s) | CSR cell list (s) |
|---|---|---|
| **20.85 (= production NC=256)** | **346.2** | 54.1 |
| 10.0 | 55.9 | 16.4 |
| 5.0 | 13.1 | 10.1 |
| 2.5 | 8.1 | 8.5 |

At production resolution this is 35 µs per halo, against the ~50 µs/halo implied by
production (16.6 ks / 333.7M; larger arrays mean more cache misses), so the mechanism
accounts for the production time. With 5 Mpc/h cells the unchanged algorithm is **26× faster**,
which extrapolates the production merge from **~4.6 h to ~10–15 min**. A full-box 5 Mpc/h
hash is 1071³ Int32 = 4.9 GB, which is fine. That brings the octant job to ~2.6 h.

## 4. Comparison with Websky (full sky = 8 octants)

| | Websky (Stein+2020) | PeakPatch.jl as run | PeakPatch.jl, merge fix (projected) |
|---|---|---|---|
| hardware | 1128 Skylake cores (~28 Niagara nodes) | 8 × (4× L40S node) | same |
| catalog time-to-solution | 3.84 h (8 octants serial) | 7.0 h (8 nodes in parallel) / 56 h (one node) | ~2.6 h / ~21 h |
| catalog compute | 4336 core-h ≈ 108 node-h | 56 node-h = 225 L40S-GPU-h (of which ~59 GPU-h GPU-active) | ~21 node-h ≈ 83 GPU-h |
| field maps (5 LOS kernels, Nside 4096) | from stored ICs; no timing given | 8 × 1.35 h on 1 L40S = 10.8 GPU-h | same |
| peak memory | 7.67 TB (whole 6144³ box) | 135 GiB host + 4 × ~22 GiB GPU ≈ 0.24 TB per octant (~30× less) | same |
| stored ICs | 5.9 TB | **0** (fields regenerated from the seed by the counter RNG) | 0 |
| catalog on disk | 33 GB (10 floats, ~9e8 halos) | 25.8 GB/octant raw (33 floats, 195M halos) + same for AM | 11-float basic format would be 8.6 GB/octant |

Framing for the paper: the node-hour comparison is **not like-for-like** (CPU nodes of 2019
vs 2024-era GPU nodes, and a different halo count; we keep all 195M halos per octant against
Websky's pre-AM 10-particle cut). The robust claims are:
- (i) the memory footprint is ~30× smaller, because no global fine-grid field is ever held;
- (ii) zero IC storage;
- (iii) a single 4-GPU node produces a Websky-equivalent octant in ~7 h as run, or ~2.6 h
  with the merge fix, and the full sky in ~2–3 days on one node.

## Caveats
- The production stage split uses mtimes. JIT/startup (~minutes) is folded into "pipeline".
- nvidia-smi memory is device-level and includes pool caching, so it is an upper bound on
  what is needed. The n=576 point is pool-limited.
- H100 used 6 CPU threads against 8 for the L40S (node core/GPU ratio). CPU-side stages
  (01, 08, 10, and the whole merge) are therefore not like-for-like.
- **Reproducibility across worker count and hardware.** Post-merge counts are 9 976 347 /
  9 976 343 / 9 976 353 for 1/2/4 L40S and 9 976 357 on H100. The H100 *pre-merge* count is
  also +1 (25 729 744), from floating-point kernel differences across architectures. The
  ~1e-6 worker-count dependence is consistent with tie-breaking in the merge:
  `sortperm(r; rev=true)` is stable and the input order follows tile completion order.
  Catalogs are therefore bitwise-reproducible only for a fixed hardware and worker count;
  otherwise they agree statistically. This matters for how plan item A5 (tiling
  independence) is worded.
- Mean GPU utilization is only 19–34%, with CPU stages between GPU kernels. That is
  headroom, not a measured limit.

## Suggested code changes (not made; src/ untouched)
1. `src/Merger/Exclusion.jl`: pick the hash resolution from the data (cell ≈ 5 Mpc/h or
   ~1 halo/cell) instead of fixed `NC = 256`, optionally with CSR storage. This is 4.6 h →
   ~10–15 min at production, with identical results. Longer term, parallelize over spatial
   slabs.
2. Make the merge order deterministic: break `RTHL` ties by (x, y, z), or sort the input by
   global tile index before merging. The catalog then no longer depends on the number of GPUs.
3. Fix the stale comment at `src/MultiResolution.jl:532-534` (shell analysis is on the GPU).
4. `run_gpu_octant.jl`: log timestamps per stage (merge/finalize/write) and the true total.
5. Throughput: overlap the CPU stages (01 residual generation, 08 peak finding, 10 packing)
   with GPU work, e.g. 2 workers per GPU.
