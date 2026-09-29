# Why v3fs is more clustered than Websky: the missing volume-reduction pass (2026-09-28)

**Question.** The Tier-A rerun (`TIERA_V3FS_2026-09-28.md`) found our halos clustered more
than Websky's at fixed number density:
- ξ(3–15 Mpc/h), Websky/ours = 0.877 ± 0.029 (4σ);
- ξ at 1–3 Mpc/h is 1.2–1.33× theirs;
- b(M) is 4–7% higher below 3e13.

## Diagnostics (`tierA_diag.jl`, job 5736375; `results/tierA_diag_v3fs.txt`, log `results/tierA_diag_job5736375.log`)

All diagnostics use the same caps and shell as Tier-A: Websky's patch against the same cap
at each of our 8 octant centres.

1. **AM is not the cause (D4).** Selecting halos by raw-mass rank instead of by AM mass
   picks 99.86–99.91% the same halos, and ξ(3–15) changes by 0.990–1.012.
2. **Ours matches Tinker10; Websky is below it (D1).**

   | M bin | 5–8e12 | 0.8–1.3e13 | 1.3–2e13 | 2–3.2e13 | 3.2–7.9e13 | 0.8–2.5e14 |
   |---|---|---|---|---|---|---|
   | b_ours / Tinker10 | 1.010 | 1.003 | 0.997 | 0.974 | 0.929 | 0.845 |
   | b_Websky / Tinker10 | 0.932 | 0.991 | 0.928 | 0.991 | 0.882 | 0.913 |

   b = √(ξ̄_hh(6–18)/ξ̄_mm,lin), with ξ_mm from our CAMB P(k) (σ8 = 0.810).
3. **Ours has more close neighbours (D2).** Pairs per halo with u = d/(R_i+R_j) (Eulerian d,
   AM R_TH), Websky/ours: 0.827 / 0.838 / 0.894 / 0.928 / 0.933 / 0.954 / 0.975 in
   u = 0–0.25 / 0.25–0.5 / 0.5–0.75 / 0.75–1 / 1–1.5 / 1.5–2 / 2–3. That is 5.9σ and 6.3σ in
   the first two bins.
4. **Removing close pairs removes the excess (D3).** A common post-hoc Eulerian exclusion
   (drop a halo whose centre lies within f·R of a larger kept halo) moves ξ(3–15)
   Websky/ours from 0.877 (f=0) to 0.893 (f=1) and 0.986 (f=2).

## Code audit: two differences from the Fortran

Both were checked against the source by hand.

1. **Volume reduction.**
   - Fortran `merge_pkvd` runs a second pass after exclusion (`merge_pkvd.f90:139-140`,
     `exclusion.f90:120-204`). Every surviving pair with d < r_i + r_j adds each halo's own cap
     beyond the mid-plane to that halo's deficit ('shared'). After all pairs,
     r → (r³ − 3ΔV/4π)^(1/3).
   - Stein+2020 §2.3 describes the same procedure. Only a per-pair re-centring variant is
     commented out (`exclusion.f90:286-295`).
   - **Our `Merger.jl` skipped this pass**, with a comment that wrongly said Fortran had it
     commented out. Every catalog up to v3fs is exclusion-only.
2. **Lightcone peak threshold.**
   - Fortran `get_pks` uses fv = fsc_of_z(z), with z from the **tile centre's** distance for
     ievol=1 (`hpkvd.f90:510-514, 612`; `peakvoidsubs.f90:82`).
   - We use fsc_of_z(z_out = 0) everywhere. The 2026-04 note (`INVESTIGATION_2026-04-16.md`)
     read only the commented-out per-cell block and missed that z is per tile.
   - Not yet tested.

Minor Fortran quirks:
- equal-radius pairs are reduced twice;
- the hash search window is truncated for the smallest halos;
- r³ < ΔV gives a NaN radius, which is never flagged.

## A/B test of volume reduction (`merge_ab.jl`, job 5736480; `results/merge_ab.txt`, log `results/merge_ab_job5736480.log`)

**Setup.** One periodic snapshot at z = 0.7, the Tier-A shell redshift, with v3 physics:
N = 1536, L = 1309 Mpc/h, nbuff 25, cf 32 (`configs/merge_ab_z07.toml`).
- 21.27M raw halos → 9,223,160 after exclusion.
- Reduction then removes 7 halos, and mean R_TH goes from 1.817 to 1.780 Mpc/h.
- Both variants are rank-matched to the Tinker08 abundance of the box.

**Reduction on/off compared with Websky/ours:**

| | B/A (reduction on / off) | Websky/ours (Tier-A) |
|---|---|---|
| ξ(1–3 Mpc/h) | **0.796** | **0.79** |
| ξ(3–15 Mpc/h) | **0.916** | 0.877 ± 0.029 |
| close pairs, u = 0–0.25 … 2–3 | 0.776 0.862 0.935 0.953 0.966 0.977 0.987 | 0.827 0.838 0.894 0.928 0.933 0.954 0.975 |
| b(M), 6 bins 5e12–2.5e14 | 0.975 0.961 0.967 0.958 0.953 0.942 | 0.92–1.08 (noisy) |

**Bias against Tinker10 in the snapshot, by mass bin:**

| | 5–8e12 | 0.8–1.3e13 | 1.3–2e13 | 2–3.2e13 | 3.2–7.9e13 | 0.8–2.5e14 |
|---|---|---|---|---|---|---|
| A: exclusion only | 1.027 | 1.027 | 1.035 | 1.020 | 0.985 | 0.877 |
| B: exclusion + reduction | 1.001 | 0.987 | 1.001 | 0.977 | 0.939 | 0.827 |

**Conclusion.** The missing volume reduction reproduces the small-scale excess exactly:
ξ(1–3) 0.796 against 0.79. It also reproduces the close-pair deficit profile, and about
70% of the ξ(3–15) offset. What remains is 0.877/0.916 = 0.96 ± 0.03, i.e. 1.4σ, which is
consistent with sample variance. The per-tile threshold is the other candidate and has not
been tested.

Two consequences:
- **Websky's clustering is what the Fortran algorithm produces.** It is not a problem with
  their catalog.
- **Our exclusion-only catalogs are more clustered than the reference algorithm makes
  them.** With reduction, our low-mass bias sits at 0.98–1.00 × Tinker10.

## Code added (default off, so existing results are unchanged)

- `merge_catalog(...; volume_reduction=true)`: a Fortran 'shared' reduction. Each pair is
  visited once from its larger member within 2 r_i; ties broken by the merge order. The old
  unused `volume_reduction!` searched r_i + r_max, which is too slow at production scale.
  Tests in `test_merger.jl` cover a two-halo analytic case and a 600-halo brute-force case.
- `[run] peak_threshold_per_tile = true`: the Fortran per-tile lightcone threshold.
  Untested beyond loading.
- CPU suite 1796/1796.

## Decisions needed

1. **Rerun the catalogs with `volume_reduction=true`?** The raw pre-merge halos were not
   saved, so this is a full catalog rerun: 8 × ~2.7 h on 4×L40S, plus AM, paint and CIB.
   Field maps are not affected.
2. **Test `peak_threshold_per_tile` first?** It needs a lightcone A/B test, because the
   snapshot threshold is the same in both codes. If adopted, it goes into the same rerun.
