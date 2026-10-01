# Production-v4 campaign ("v4fs"), submitted 2026-09-29

v4 is v3 (`../production_v3/`: nbuff=25, shell early exit, Ω_m(a) fix) plus the two
Fortran-faithful catalog steps found missing
(`../../paper/CLUSTERING_EXCESS_2026-09-28.md`):

1. **`volume_reduction = true`**: the `merge_pkvd` volume-reduction pass after exclusion.
   Without it, our catalogs were ~7% more biased than Websky. A same-raw-catalog A/B test
   gives ξ(1–3) 0.796, against Websky/ours 0.79.
2. **`peak_threshold_per_tile = true`**: the per-tile lightcone candidate threshold from
   `hpkvd`. Its effect is negligible above 1e12 (`../../paper/results/threshold_ab.txt`).

AM uses one full-sky table from all 8 raw octants with `tail_N` = 10, as for v3fs. Field
maps are halo-independent and are symlinked from `fieldmaps_v3`.

**Storage.** All outputs are on scratch:
- catalogs: `scratch/websky_6144/catalogs_v4/`, raw and `_AMfs`;
- maps: `scratch/websky_6144/v4fs/`, reached through the /project symlinks
  `{halomaps,cibmaps,fullsky}_v4fs` and the symlink directory `fieldmaps_v4fs`.

## Jobs (`bash submit_all.sh`; log `scratch/.../logs/v4_submit.txt`)

| step | jobs |
|---|---|
| catalogs oct000…111 | 5752129–5752136 |
| full-sky AM table | 5752137 (after all catalogs) |
| AM / paint / CIB per octant | 5752138 + 3k / 5752139 + 3k / 5752140 + 3k (k = 0…7) |
| full-sky analysis (CAMPAIGN=v4fs) | 5752162 |
| tail check | 5752163 |
| Tier-A (TIERA_CAMP=v4) | 5752164 |

Run from commit 756fa4b. **Do not edit `src/` while catalog jobs are queued**, because each
job loads the live working tree.
