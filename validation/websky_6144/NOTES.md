# Websky 6144³ Full Reproduction

## Goal
Reproduce the Websky extragalactic CMB simulation (Stein et al. arXiv:2001.08787)
using PeakPatch.jl on CITA's Sunnyvale cluster.

## Configuration

### Grid
- N = 6144, cellsize = 7700/6144 = 1.25326 Mpc/h (Websky exact)
- nmesh = 768, nbuff = 48, nsub = 672, ntile = 9
- Total tiles: 729, active per octant (z_max=4.6): ~341
- Box: 7700 Mpc/h, core: 7579.7 Mpc/h

### Lightcone
- 8 octants (observer at each corner of the core region, ±3789.84 Mpc/h)
- z_max = 4.6, chi_max ~ 6700 Mpc/h
- Expected ~100-120M halos per octant, ~900M total (full sky)

### Cosmology (Planck 2018, Websky values)
- Om = 0.31, OB = 0.049, OL = 0.69, h = 0.68

### Physics
- seed = 12345, 2LPT, top-hat smoothing
- 21 filters (Rf = 2.08 to 34 Mpc/h, spacing 1.15×)
- Extended catalog (ioutshear=1) for XGPaint compatibility

## Implementation: MPI Lowmem Batched Mode

Standard PeakPatch MPI extracts all displacement tiles during Phase 1 and holds
them in memory for Phases 2-4. At 6144³ this requires ~5 TB of tile data across
nodes — too much even for 12 × 1 TB nodes.

The **lowmem batched mode** (`lowmem=true`) restructures the computation:

### Phase 1: Field Generation
Only computes `delta_k` and `src2_k` (distributed via PencilFFTs). No tile
extraction. Peak memory during incremental src2 computation: 4 distributed
fields simultaneously (delta_k + src2 + 2×phi).

### Phase 2: Peak Finding
Same as standard — iterates over 21 filter scales, does distributed IFFT of
smoothed delta_k, extracts smoothed tiles temporarily for peak detection.
Only needs delta_k + 1 smoothed field at a time.

### Phase 3-4: Batched Shell Analysis
Instead of using pre-extracted tiles, processes tiles in batches:
1. Determine batch size from available memory (auto-computed from `available_gb`)
2. For each batch of ~26 tiles:
   - 8 distributed IFFTs (delta, psi×3, psi2×3, lapd) → extract batch tiles
   - Run shell analysis on batch tiles
   - Free all tile data
3. Total: 3 batches × 8 IFFTs = 24 distributed IFFTs (vs 8 standard)
4. Overhead: ~16 extra IFFTs ≈ 10-15 minutes

### Float32 Distributed Fields
Changed PencilFFT plan from Float64 to Float32 (ComplexF32 k-space).
Halves all distributed array sizes. This was essential for fitting on 6 nodes:
- Float64: 2LPT peak = 4 × 309 GB = 1237 GB/node (exceeds 1 TB)
- Float32: 2LPT peak = 4 × 155 GB = 619 GB/node (fits comfortably)

K-space helper functions compute using Float64 k-vectors (small 1D arrays)
but store results as ComplexF32. Tile extraction is Float32 (unchanged from
standard). Roundtrip FFT error ~3e-7 — acceptable for Float32.

### `_extract_my_tiles` Batch Support
Added optional `global_tile_list` kwarg. When provided, only those tiles
participate in MPI communication (sends/receives). All ranks must pass the
same list — deterministic batch decomposition ensures this.

## Memory Budget (6 starq nodes, Float32)

| Phase | Peak memory/node | Distributed arrays |
|-------|-----------------|-------------------|
| 1e (src2, incremental) | 619 GB | delta_k + src2 + 2×phi |
| 2 (peak finding) | ~465 GB | delta_k + src2_k + smoothed temp |
| 3-4 (batched shells) | ~619 GB | delta_k + src2_k + IFFT temp + batch tiles |

Batch size: 26 tiles/batch, 3 batches per octant.

## Timing Estimate (per octant)

| Phase | Estimated time |
|-------|---------------|
| Precompilation | ~10 min (one-time) |
| Phase 1 (noise + delta_k + src2_k) | ~2 hours |
| Phase 2 (21 filters × distributed smooth) | ~1.5 hours |
| Phase 3-4 (shells, 3 batches × 8 IFFTs) | ~2 hours |
| Phase 5 (gather + merge + write) | ~15 min |
| **Total** | **~6 hours** |

Full sky (8 octants): ~48 hours sequential.

## Differences from Websky Production

| Parameter | Websky (Stein+2020) | This run |
|-----------|---------------------|----------|
| Code | Fortran peakpatch + Python post | PeakPatch.jl |
| Grid | 6144³ | 6144³ (same) |
| Box | 7.7 Gpc/h | 7.7 Gpc/h (same) |
| Octants | 8 (full sky) | 8 (same) |
| z_max | 4.6 | 4.6 (same) |
| RNG | Fortran LCG | Threefry |
| Field precision | Float32 | Float32 distributed, Float32 tiles |
| Abundance matching | Sheth-Tormen | Tinker (2008) |
| Map projection | pks2map (Fortran) | XGPaint.jl |

## Cluster Deployment

### Resource Request
- Queue: starq (1 TB RAM, 128 cores per node, 24h max walltime)
- Nodes: 6 × ppn=96 = 576 cpus (fits within 1536 queue limit)
- Walltime: 12 hours per octant
- Output: `/fs/lustre/scratch/yguan/work/peakpatch/websky_6144/output/`

### Deployment Issues Encountered
1. **"Unable to find compatible target in cached code image"**: Cross-node CPU
   target mismatch in Julia pkgimages. Fix: precompile on the job's head node
   before launching mpiexec (added single-rank `using PeakPatch` before mpiexec).

2. **"Disk quota exceeded"**: `--pkgimages=no` causes 6 ranks to simultaneously
   write compilation artifacts to ~/.julia, exceeding NFS quota. Fix: don't use
   `--pkgimages=no`, instead precompile once and clear stale cache before runs.

3. **Queue resource limits**: starq max 1536 ncpus total. With other users'
   jobs occupying 5 nodes (360 cpus), our 10-node request (960 cpus) wouldn't
   start despite fitting the cpu limit — the scheduler needs exclusive nodes.
   Fix: reduced to 6 nodes.

4. **Node availability**: 5 of 12 starq nodes occupied by long-running jobs
   from another user. Only 6 free nodes available. Float32 fields made 6 nodes
   feasible (Float64 needed 8+ nodes).

## Performance Fixes Applied (2026-04-14)

1. **Float32 distributed fields** (`ext/MPIExt.jl`): PencilFFT plan uses Float32
   instead of Float64. Halves all distributed array sizes.

2. **Extract-then-convert** (`src/MultiTile.jl:241,261,579,597`): Changed
   `extract_tile(Tf.(full_grid), ...)` to `Tf.(extract_tile(full_grid, ...))`.
   Avoids allocating a full N³ temporary array per tile per filter.

3. **Laplacian Float32 extraction** (`ext/MPIExt.jl`): Extract tiles as Float32
   directly with in-place scaling, instead of Float64 extraction + conversion.

4. **MPI lowmem batched mode** (`ext/MPIExt.jl`): New `lowmem=true` kwarg for
   `run_multitile_mpi`. Defers tile extraction to batched Phase 3-4.

### Remaining Performance Opportunities
- `smooth_field` (Filters.jl:89): copies entire delta_k each filter call.
  Could multiply into a pre-allocated buffer instead.
- `displacements_2lpt` (LPT.jl): holds 6 phi_ij arrays simultaneously in
  non-lowmem serial path. Could use `compute_src2_k` incremental approach.
- Missing `GC.gc()` calls in serial `run_multitile` after filter loop.

## Files
- `config_websky_6144_oct{000..111}.toml` — 8 octant configs
- `filters_websky.dat` — 21 Websky-style filters (same as 1760³ run)
- `run_websky_6144.pbs` — PBS script (pass `OCTANT=oct000` etc.)
- `generate_octant_configs.py` — generates all 8 configs from oct0 template

## Submission

```bash
# Single octant:
qsub -v OCTANT=oct000 validation/websky_6144/run_websky_6144.pbs

# All 8 octants (run one at a time, check output between):
for oct in 000 001 010 011 100 101 110 111; do
    qsub -v OCTANT=oct${oct} validation/websky_6144/run_websky_6144.pbs
done
```

## Run Log

| Job ID | Date | Octant | Nodes | Status | Notes |
|--------|------|--------|-------|--------|-------|
| 596545 | 2026-04-14 | oct000 | 10×96 | Q (never started) | Blocked: queue ncpus limit |
| 596547 | 2026-04-14 | oct000 | 8×96 | Q (never started) | Only 6 nodes free |
| 596548 | 2026-04-14 | oct000 | 6×96 | Crashed (5 min) | "Unable to find compatible target" — stale pkgimages |
| 596549 | 2026-04-14 | oct000 | 6×96 | Crashed (20 min) | "Disk quota exceeded" — --pkgimages=no fills home |
| 596551 | 2026-04-14 | oct000 | 6×96 | Crashed (20 min) | MPI.Initialized() MethodError — JULIA_DEPOT_PATH broke MPI prefs |
| 596552 | 2026-04-14 | oct000 | 6×96 | Crashed (OOM) | Out-of-place RFFT needs ~1108 GB, exceeds 1007 GB |
| 596561 | 2026-04-14 | oct000 | 6×96 | Crashed (OOM) | Same OOM during Phase 1b forward FFT |
| 596571 | 2026-04-14 | oct000 | 6×96 | Killed | In-place RFFT! fix, but FFTW 1-thread bug, ran 9h on 1 core |
| 596572 | 2026-04-14 | oct000 | 6×96 | Running | FFTW threading fix (commit 4976a6b) |
| 596573 | 2026-04-14 | oct000 (multires) | 6×32 greenq | Running | Multi-resolution comparison, no MPI |

## Multi-Resolution Comparison Run (2026-04-14)

Submitted a parallel run using `run_multitile_split` (MUSIC-style decomposition)
on greenq nodes to compare against the MPI distributed FFT approach on starq.

### Why
The multi-resolution approach eliminates the need for distributed FFT entirely.
Each tile generates its fields independently using:
1. A cheap coarse-grid periodic FFT (M=24, negligible memory)
2. Per-tile isolated FFT of the Hoffman-Ribak residual noise (2×768³ zero-padded)

This reduces peak memory from ~620 GB/node (MPI) to ~40 GB/tile, enabling
the use of much smaller nodes (128 GB greenq instead of 1 TB starq).

### Configuration
- Same physics: 6144³, ntile=9, 2LPT, lightcone oct000, seed=12345
- Coarse grid: auto-selected M=24, block=256 (optimal for halo counts)
- Queue: greenq (6 nodes × 128 GB × 32 cores)
- Script: `run_websky_6144_multires.pbs`
- Output: `/fs/lustre/scratch/yguan/work/peakpatch/websky_6144_multires/`

### Memory Budget (per tile, greenq node)
| Component | Memory |
|-----------|--------|
| Padded FFT array (1536³ Float32) | 14.5 GB |
| Complex k-space (1536³ ComplexF32) | 14.5 GB |
| Residual noise (768³) | 1.8 GB |
| Working arrays (delta, psi, etc.) | ~9 GB |
| **Total per tile** | **~40 GB** |

Fits comfortably in 128 GB with room for 2 tiles simultaneously.

### Accuracy (validated on smaller grids)
| Grid | Halo count Δ vs global FFT |
|------|---------------------------|
| N=120 (ntile=2, cf=5) | 1.2% |
| N=260 (ntile=2, cf=2) | 0.52% |
| N=6144 (ntile=9, M=24) | Expected <0.5% (larger tile:box ratio) |

### Performance Comparison
| | MPI (starq) | Multi-resolution (greenq) |
|---|-------------|--------------------------|
| Nodes | 6 × 1 TB × 96 cores | 6 × 128 GB × 32 cores |
| Peak RAM/node | ~620 GB | ~40 GB |
| MPI required | Yes (PencilFFTs) | No |
| Parallelism | Distributed FFT | Independent tiles |
| Expected time | ~6 hours | TBD |

### Run Log

| Job ID | Date | Approach | Queue | Nodes | Status | Notes |
|--------|------|----------|-------|-------|--------|-------|
| 596571 | 2026-04-14 | MPI | starq | 6×96 | Killed | FFTW single-threaded (bug: missing set_num_threads) |
| 596572 | 2026-04-14 | MPI | starq | 6×96 | Running | Fixed FFTW threading, resubmitted |
| 596573 | 2026-04-14 | Multi-res | greenq | 6×32 | Running | First multi-resolution 6144³ run |

### FFTW Threading Bug (found 2026-04-14)
`FFTW.set_num_threads()` was called in `Pipeline.jl` (single-tile) but NOT
in `MPIExt.jl` (MPI driver). This caused all 6144³ MPI runs to use 1 FFTW
thread per rank despite launching with `-t 96`. The FFT phase that should
take minutes was taking hours. Fixed in commit `4976a6b`.

## Troubleshooting Playbook

When a job fails, check:

```bash
# 1. What happened?
grep -E "ERROR|signal|Bus error|Killed|quota" websky_6144.out | head -10
```

### "Unable to find compatible target in cached code image"
Stale pkgimages compiled on a different CPU architecture.
```bash
rm -rf ~/.julia/compiled/v1.11/PeakPatch* ~/.julia/compiled/v1.11/MPIExt*
qsub -v OCTANT=oct000 validation/websky_6144/run_websky_6144.pbs
```
The PBS script already precompiles on the head node before launching mpiexec.

### "Disk quota exceeded"
Home NFS quota hit during compilation. Do NOT use `--pkgimages=no` (causes
6 ranks to compile simultaneously). Instead:
- Clear old compiled cache: `rm -rf ~/.julia/compiled/v1.11/PeakPatch*`
- The PBS script precompiles once on head node before mpiexec, avoiding parallel writes

### "MethodError: no method matching Initialized()"
MPI preferences not found. Usually caused by `JULIA_DEPOT_PATH` not including
`~/.julia`. The PBS script should NOT set `JULIA_DEPOT_PATH`.

### Bus error / signal 7 / OOM
Out of memory. Check:
```bash
# Which phase crashed?
grep -E "Phase|Filter|Batch" websky_6144.out | tail -20
```
- **During Phase 1e (2LPT)**: Peak is 4 distributed fields. With Float32 at
  6 nodes: 619 GB/node. If this crashes, need more nodes or there's a leak.
- **During Phase 2 (filters)**: Need delta_k + src2_k + smoothed temp.
  ~465 GB/node. Should not OOM.
- **During Phase 3-4 (batched shells)**: Batch too large. Reduce `available_gb`
  in config or code (default 900). Check batch size in log output.

### Job stuck in Q state
```bash
qstat -f <jobid> | grep comment
qstat -a | grep starq                        # see who's using nodes
pbsnodes -a | awk '/^tpb/{n=$1} /starq/{f=1} /state =/{s=$3} f{print n,s;f=0}'
```
- "exceed overall limit on resource ncpus": too many cpus requested.
  Reduce nodes or ppn.
- No comment, just "Not Running": not enough exclusive nodes free.
  Check which starq nodes are occupied. Need 6 free.

### Resubmit after fixing
```bash
# Always clear PeakPatch compiled cache after code changes:
rm -rf ~/.julia/compiled/v1.11/PeakPatch* ~/.julia/compiled/v1.11/MPIExt*

# Submit:
cd /home/yguan/work/cita/PeakPatch.jl
qsub -v OCTANT=oct000 validation/websky_6144/run_websky_6144.pbs

# Monitor:
qstat -u yguan
tail -f /fs/lustre/scratch/yguan/work/peakpatch/websky_6144/output/oct000.log
```

### Changing node count
If more/fewer nodes become available, edit `run_websky_6144.pbs`:
- `nodes=N:ppn=96` and `mpiexec -np N`
- Minimum: 6 nodes (Float32 fields, 619 GB peak during 2LPT)
- Maximum: 12 nodes (most comfortable, fewer batches)
- Queue limit: N × 96 + other_users_cpus ≤ 1536
- Batch size auto-adjusts from `available_gb` parameter (default 900 GB)
