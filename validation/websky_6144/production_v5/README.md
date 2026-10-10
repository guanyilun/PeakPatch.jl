# Production-v5 campaign ("v5fs"), Flatiron Rusty

v5 is v4 (`../production_v4/`: nbuff 25, volume reduction, per-tile threshold, full-sky AM with
`tail_N` = 10) plus **`gaussian_split = true`**. That is the Gaussian long/short multires handoff
that removes the split ψ aliasing (`../../paper/MATCHED_FORTRAN_2026-10.md` §9–10). The
production-size test on oct000 (v5test) moved the Websky ξ ratio from 0.942 to 0.980, put the bias
at 0.977–0.997, and cost about 9% more pipeline time.

Differences from v4 besides the config flag:
- **The AM table is rebuilt from the v5 raw catalogs.** v5test reused the v4 table, which shifted
  dN/dz by 8–16% at z ≈ 3–4.5.
- **Field maps are regenerated.** v4 reused the v3 maps, but gaussian_split also changes the
  field (`src/FieldMap.jl`, c880c24).
- **Cluster.** This campaign runs on Flatiron Rusty, not Killarney (the account closed 2026-10-08).
  There is no `--account`, and nothing is copied to scratch before `sbatch`.
- **Paths.** They come from `env.sh` (`WS_ROOT`, `WS_CATS`, `WS_REF`). The analysis scripts in
  `validation/paper/` read the same variables and fall back to the Killarney paths, so older
  campaigns still resolve as before.

## Layout (`env.sh` defaults)

```
/mnt/ceph/users/yguan/projects/uoft/peakpatch/websky_6144/
  catalogs_v5/        raw + _AMfs catalogs (15 GB each), am_table_fullsky_tail10.txt, .prov.toml   ~240 GB
  fieldmaps_v5fs/     5 field maps per octant (v5_oct* files + v5fs_oct* symlinks)                 ~31 GB
  halomaps_v5fs/  cibmaps_v5fs/  fullsky_v5fs/                                                     ~81 GB
  websky_ref/         released Websky files (copied from the Killarney backup, md5-checked)        ~41 GB
  fullsky_v4fs/       only ksz_field_uK_v4fs_*: the "before" column of the theory step (from the backup)
  logs/
```

## Jobs (`bash submit_all.sh` from a clean, committed tree)

| step | n | partition | per job | measured, first v5 run on Rusty (2026-10-10) |
|---|---|---|---|---|
| catalog | 8 | gpu, a100-80gb&rocky9 | 2 GPU, 32 CPU, 128G, 3.5 h | 1.75–1.82 h; MaxRSS 90.1–93.5 GiB |
| field maps | 8 | gpu, a100-80gb&rocky9 | 1 GPU, 16 CPU, 64G, 2.5 h | 1.26–1.31 h; MaxRSS 35–40 GiB |
| AM table | 1 | genx | 8 CPU, 96G, 1.5 h | 19 min; MaxRSS 64.0 GiB (hit the old 64G cap) |
| AM apply | 8 | genx | 8 CPU, 96G, 30 min | 8 min; MaxRSS 57.8 GiB |
| halo paint | 8 | genx | 32 CPU, 48G, 45 min | 6.5–7.5 min; MaxRSS 28.5 GiB |
| CIB | 8 | genx | 32 CPU, 48G, 1.5 h | 34–37 min; MaxRSS 27.2 GiB |
| analysis + AM tails | 1 | genx | 32 CPU, 64G, 1.5 h | 23.5 min to the cross spectra; MaxRSS 35.7 GiB |
| Tier-A | 1 | genx | 16 CPU, 128G, 1.5 h | 16.5 min; MaxRSS 96.1 GiB (hit the old 96G cap) |

The gpu QoS allows 16 GPUs per user. Catalogs use 2 GPUs each, so all 8 fit under the cap at once
(field maps then queue behind them). With 4 GPUs per job they waited many hours for a fully free
node, while 2 free GPUs on a node are common. The `rocky9` constraint is required: `env.sh` loads
modules from the Rocky 9 tree. Memory requests are about 1.3-1.7× the measured peaks; per-job MaxRSS is from
`sacct -j <id> -o JobID,MaxRSS` (the `.batch` step).

## Before the first submission
1. Instantiate the environments under julia/1.12.7 (`validation/` and `$XGPAINT`).
2. Commit. `run_gpu_octant.jl` refuses a dirty tree, and `submit_all.sh` checks this up front.
3. Do not edit `src/` until the catalog jobs have started, because each job loads the live tree.

## Rerunning analysis steps
`bash validation/websky_6144/production_v5/rerun_analysis.sh [step ...]` resubmits
`run_analysis.slurm` for some steps (assemble seams autos cross theory tails; default all). It works
from any directory. The first run (job 7206290) failed at `theory` because `fullsky_v4fs/` was missing.

## First run (2026-10-10)
Jobs 7206249–7206291, plus 7211354 (theory + AM tails rerun) and 7211375 (full-sky Tier-A). Results
are in `../../paper/V5_RESULTS_2026-10-10.md`. Against the full Websky catalog, ξ(3–15 Mpc/h)
W/ours is 0.978 ± 0.003 (v4: 0.948 ± 0.002).
