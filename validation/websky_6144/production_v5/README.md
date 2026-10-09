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
  logs/
```

## Jobs (`bash submit_all.sh` from a clean, committed tree)

| step | n | partition | per job | basis (Killarney, v4/v5test) |
|---|---|---|---|---|
| catalog | 8 | gpu, a100-80gb | 4 GPU, 32 CPU, 160G, 4 h | 1.62–1.95 h (+9% for v5); MaxRSS 89.7 GiB / 96 GB |
| field maps | 8 | gpu, a100-80gb | 1 GPU, 16 CPU, 96G, 3 h | ~1.35 h (v3); memory not measured |
| AM table | 1 | genx | 8 CPU, 64G, 1.5 h | 25 min; reads one 15 GB catalog at a time |
| AM apply | 8 | genx | 8 CPU, 64G, 30 min | 7–8 min; raw + AM copy ≈ 30 GB |
| halo paint | 8 | genx | 32 CPU, 48G, 45 min | 7 min; MaxRSS 21–31 GiB (frozen campaign) |
| CIB | 8 | genx | 32 CPU, 128G, 1.5 h | 40 min; memory not measured |
| analysis + AM tails | 1 | genx | 32 CPU, 128G, 1.5 h | ~26 + 9 min; memory not measured |
| Tier-A | 1 | genx | 16 CPU, 96G, 1.5 h | 16 min (v4), 30 min (v5test) |

The gpu QoS allows 16 GPUs per user, so catalogs run 4 at a time. Once the first jobs finish,
check `seff <jobid>` for one job of each type and tighten the requests, especially the rows whose
memory was never measured.

## Before the first submission
1. Instantiate the environments under julia/1.12.7 (`validation/` and `$XGPAINT`).
2. Commit. `run_gpu_octant.jl` refuses a dirty tree, and `submit_all.sh` checks this up front.
3. Do not edit `src/` until the catalog jobs have started, because each job loads the live tree.
