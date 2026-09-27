# Production-v2 campaign (paper dataset) — submitted 2026-09-26

This reruns the frozen 2026-07 campaign (`../production/`) with the two defects found by the
full-sky comparisons fixed:

1. **Seam gap** (`../../paper/FULLSKY_COMPARISON_2026-09-26.md`). The frozen layout
   N = nsub·ntile + 2nbuff left 32 cells per axis that were never in a tile core, i.e. an
   unsimulated 13.6 Mpc/h slab next to each octant plane. v2 uses
   `periodic_cores = true`, n = 416, nsub = 384 and N = nsub·ntile = 6144, so the
   buffers wrap.
2. **Multires splice P(k) bump** (`../../tiling/SPLICE_COMPENSATION_2026-09-26.md`), +2–11%
   in δ and ψ at k = 0.05–0.25 h/Mpc. v2 uses `coarse_compensation = true`
   (coarse kernels × D/T).

Code fixes are bundled too: deterministic merge order (total order on (−R, x, y, z)), a
merge hash sized from the data (`auto_hash_nc`, about 4.6 h → about 10–15 min), and
tile-local φ_ij zeroing only k=0 on both GPU and CPU (2LPT trace identity).

**Unchanged:** seed 12345, N=6144, cell 5236/6144 = 0.852213 Mpc/h, box 5236 Mpc/h
(= 7700 Mpc), ntile 16, nbuff 16, cf 32 (M = 512, block 12), z_max 4.5, 2LPT,
ioutshear 1, finecell filter bank. The observer is at ±2618 per axis (octZYX bits,
1 → +), now exactly the core corner.

## Jobs (`bash submit_all.sh`; per octant)

| stage | script | partition | notes |
|---|---|---|---|
| catalog + AM | `run_catalog.slurm` | b2, 4×L40S, 8 h | → `catalog_websky_6144_v2_octZYX.pksc`, `_AM.pksc` |
| field maps | `run_fieldmap.slurm` | b1, 1×L40S | → `/project/.../websky/fieldmaps_v2/{kappa,mass,tau,ksz,isw}_v2_octZYX_nside4096.fits` |
| halo maps | `../production/run_paint.slurm` (afterok catalog) | b1 | → `halomaps_v2/*_v2_octZYX_AM_nside4096.fits` |
| CIB | `../production/run_cib_prod.slurm` (afterok catalog) | b1 | → `cibmaps_v2/cib_nuXXXX_wcut_v2_octZYX_AM_nside4096.fits` |

Submitted 2026-09-26:
- catalogs 5696115/119/123/127/131/135/139/143, each followed by paint (+1) and CIB (+2)
- field maps 5696740–747

A first batch of field maps (5696118–146) was cancelled before writing any output,
when the splice compensation was added.

**Afterwards:** `CAMPAIGN=v2` for `../../paper/assemble_fullsky.jl`, `run_fullsky_compare.sh`
and `compare_cross.jl`. Check the seams on every product first: signal against distance to
the octant planes, as in `../../paper/diag_holes_tsz.jl`.

**Pre-flight** (`examples/config_gpu_octant.toml`, which is now periodic + compensated;
jobs 5696038/5696088/5696647):

| layout | halos per 5 Mpc/h of plane distance |
|---|---|
| periodic | 15.5k, 19.3k, 19.4k, … |
| legacy | 0.3k, 2.3k, 8.2k, … |

Compensation changes the small test octant by +1.2% in N, +1.8% in N(>1e13) and −0.9% in σ_vx.
