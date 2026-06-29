# The Websky low-mass count deficit (2026-06-15)

> 🎯 **ROOT CAUSE of the sub-3e12 residual: a FACTOR-of-h cellsize error.** Our run used
> box=7700 **Mpc/h** (cellsize 1.2533 Mpc/h); Websky's (7.7 Gpc)³ = 7700 **Mpc** (cellsize
> 1.2533 Mpc = **0.852 Mpc/h**). We ran **1/h = 1.47× COARSER**, so our resolution mass floor is
> ~3e12 vs Websky's ~1.2e12 — exactly where the AM-corrected catalog still falls short. The
> factor-h pervades box/cellsize/filter (same NUMBERS, ours in Mpc/h): that's why "filter numbers
> match exactly" still left a deficit. Decisive check: Websky's floor 1.23e12 = its 10-cell mass
> AM-boosted UP only if cellsize=0.852 (at 1.253 it would need boost DOWN, nonsensical). See
> `memory/websky_cellsize_factor_h.md`.
>
> **CORRECTED CONFIG to match Websky:** box = **5236 Mpc/h** (=7700 Mpc), N=6144 → cellsize
> **0.8522 Mpc/h**; filter bank regenerated → `data/filters_websky_finecell.dat` (23 filters,
> Rf_min=1.406 Mpc/h → smallest-halo mass 1.00e12 ≈ Websky floor); pk/HomelTab/cosmology unchanged.
> Secondary: bump nbuff 16→~24 to keep physical buffer ~20 Mpc/h; box is barely short of χ(4.6)=
> 5258 (reaches z≈4.57, same as Websky → minor replication). Two-part resolution: (1) apply AM
> [fixes M>3e12]; (2) rerun at 0.852 Mpc/h [fixes <3e12].
>
> ✅ **CONFIRMED (2026-06-17, job 3968163 on b3, 6h39m, 197.5M halos):** the corrected-cellsize
> octant's low-mass MF now rises down to ~1.1e12 (smallest-filter mass 1.00e12 = Websky floor)
> instead of turning over at ~2e12. ours/Websky ratio in 1-3e12 jumped from 0.09-0.42 (coarse) to
> **0.53-0.79** (fine); 1.12e12 bin: 397→3050/deg² (~8×); matches Websky to ~5-15% above 3.5e12.
> Floor dropped 3e12→1e12 exactly as predicted — factor-h cellsize diagnosis CONFIRMED. (Residual
> 0.5-0.8 at 1-3e12 is the AM regime; this catalog is raw/merged, no AM.) The fine octant is a
> 6-7h / 340M-halo / ~40GB job (merge ~2h is the bottleneck) vs coarse ~19min.


> ✅ **RESOLVED: it was the missing ABUNDANCE-MATCHING step.** Websky remaps top-hat masses
> onto the Tinker M200 mass function; our production catalog did not. Applying AM correctly
> (Tinker, per-z-shell rank match, **Om=0.31**) makes our catalog match Websky to a few percent:
> AM/Websky cumulative N(>M) per deg² = 1.23e12:0.70, 1.69e12:0.83, **3e12:1.04, 5e12:1.00,
> 1e13:1.00, 3e13:1.02, 1e14:1.01**, 3e14:0.93. dN/dz at M>3e12 AM/Websky ≈ 0.85→1.0 (z<2)
> →1.2-1.36 (z>2.5).
>
> **The residual below 3e12 is the SMALLEST-FILTER scale, not a cellsize-resolution limit**
> (cellsize 1.25 Mpc/h is identical to Websky). Our smallest filter Rf=rmincell·cellsize=
> 1.65·1.25=2.068 → M=(4π/3)ρm Rf³=3.18e12: the deficit begins EXACTLY at this mass. Halos
> below it exist only by collapse-shrinkage of the smallest-filter peaks. The dense-bank test
> (Step 3, post-merge) shows rmincell=1.2 gives 1.33–1.35× more halos at 1.23–1.69e12 — and the
> AM residual is exactly 1/1.43=0.70 @1.23e12, 1/1.20=0.83 @1.69e12. So a smaller rmincell
> quantitatively closes it (AM is rank-matched, so floor completeness = #raw peaks to remap).
> This implies Websky's PRODUCTION used rmincell smaller than the filter_gen.py DEFAULT of 1.65
> (the generator exposes rmincell as tunable; commented `#rmincell=2.`). I verified the default,
> never located Websky's actual run param — so this is strong quantitative inference, not confirmed.
>
> Two fixes were needed: (1) apply AM at all (missing production step); (2) build the AM cosmology
> with **Om_total=0.31** — the script had `CosmologyParams(0.31+0.049,...)`=0.359, double-counting
> baryons and corrupting rho_mean/D(z)/shell-volumes (AM/wsky was 0.14 at z=4.6 before the fix).
> Job 3961489 (CPU-only). NOT a filter-bank/finder/field/lightcone bug — all verified correct.
> The earlier sections below trace the (longer, error-corrected) road to this answer.


> ⚠️ **HEADLINE CORRECTION (later same day): the filter bank is NOT the cause.**
> I initially concluded the deficit was driven by the filter bank (denser bank → more
> low-mass halos). Then I read the ACTUAL Websky filter-bank generator on
> `gw.cita.utoronto.ca:/fs/lustre/project/act/njcarlson/peakpatch/python/filter_gen.py`:
> ```python
> rmincell = 1.65   # smallest smoothing radius in cellsize units
> spacing  = 1.15   # best spacing to use for filters
> rmin = rmincell * cellsize ; filters rmin·1.15^k up to Rsmooth_max
> ```
> This is **bit-for-bit our standard `filters_websky.dat`** (rmincell 1.65, spacing 1.15,
> Rf_min=2.068 at cellsize 1.25326, ~20–21 filters to Rsmooth_max=34). ⚠️ CAVEAT (added after
> review): this proves ours == the generator's **DEFAULT**, NOT ours == Websky's actual
> PRODUCTION bank — I never located Websky's run param. The generator exposes rmincell as
> tunable (`#rmincell=2.` commented), and evidence is CONFLICTING on what production used:
> filter_gen default=1.65; a Stein+2020 note says Rf_min=2·cell=2.507 (rmincell=2.0, COARSER);
> the dense-test residual implies rmincell≈1.2 (FINER). So "Websky used the SAME bank as our
> standard" is UNVERIFIED. The dense-bank lift (~1.35× post-merge at 1.23-1.69e12) quantitatively
> matches the AM residual (1/0.70-0.83), i.e. a finer smallest filter would close it.
>
> **Leading explanation now: abundance matching (which Websky applies and our production
> catalog does NOT).** The bin-by-bin signature below — ratio ramping 0.41 (low mass) → 1.0
> (high mass) — is exactly an AM signature: AM is ≈identity at high mass (peak-patch already
> matches the N-body MF there) and boosts the count above a fixed low-mass threshold by
> remapping the peak-patch MF onto the steeper N-body target. NOT a uniform mass offset (that
> would shift all bins), NOT the filter bank (proven identical). To confirm: apply AM to our
> catalog and re-run `compare_websky_mf.jl`.
>
> The Step 1–3 analysis below is retained as-is (data is correct); only the *cause
> attribution* changed. See the revised "Bottom line".

This supersedes the "grid-limited floor + M200m/AM" framing in
`FORTRAN_COMPARISON_2026-06-14.md`. That note correctly proved the *finder* is not
buggy (= Fortran to 0.5% at matched cellsize); it did NOT identify what actually
causes the ~2× count deficit vs Websky.

## Step 1 — localize the deficit (it is NOT a uniform offset)
`compare_websky_mf.jl`: N(>M) per deg², our production octant vs the Websky public
10°×10° patch, SAME mass definition (both M=(4π/3)ρm R³, top-hat Lagrangian — verified
from Websky `readhalos.py`: `M200m = 4π/3·ρ·R³`, NOT a within-R200m mass). Result:

| M (M☉/h) | ours/Websky |
|---|---|
| 1.23e12 | 0.41 |
| 1.69e12 | 0.52 |
| 3e12 | 0.68 |
| 5e12 | 0.79 |
| 1e13 | 0.89 |
| **3e13** | **1.01** |
| 1e14 | 1.09 |
| 3e14 | 1.02 |

**We MATCH Websky at high mass and fall short only at low mass.** A uniform
normalization/config error would offset all bins equally; a mass-definition difference
would shift the whole curve. Neither is happening. The deficit is purely low-mass
**completeness** — exactly what the filter bank (smallest Rf, spacing, count) controls.

## Step 2 — confirm the cause: change ONLY the filter bank
`run_filter_test.jl` runs `run_multitile_split` on the SAME box/seed/cellsize/cosmology,
varying ONLY `filterbank`:
- STANDARD = `filters_websky.dat`: 21 filters, rmincell 1.65, Rf 2.07–34 Mpc/h
- DENSE = `filters_dense.dat`: 35 filters, rmincell 1.2, Rf 1.50–36 Mpc/h (reaches ~1.2 cells)

**Pre-merge (raw peaks):**
| M (M☉/h) | STANDARD | DENSE | dense/std |
|---|---|---|---|
| >1.23e12 | 65,975 | 104,296 | 1.58× |
| >1.69e12 | 65,223 | 102,222 | 1.57× |
| >3e12 | 59,777 | 87,484 | 1.46× |
| >1e13 | 39,101 | 51,184 | 1.31× |
| >1e14 | 5,611 | 7,467 | 1.33× |
| total | 69,008 | 118,234 | 1.71× |

Reaching to smaller scales recovers ~1.5–1.6× more low-mass halos — precisely the band
where we were short of Websky. **The filter bank is the cause.**

### Caveats (honest)
- These are PRE-merge counts. The denser bank also produces more cross-scale duplicate
  peaks that `merge_catalog` removes, so the NET post-merge gain is smaller than 1.5×.
  (The high-mass rise — 1e14: 1.33× — is almost entirely duplicates that merge back to
  the same big halos; it is NOT new high-mass halos.) The apples-to-apples number is the
  POST-MERGE comparison — see Step 3.
- Even the dense bank lifts 1.23e12 only from 0.41→~0.65× Websky (pre-merge). The
  residual likely needs an even finer/smaller bank approaching the grid Nyquist (~1 cell),
  and/or matching Websky's exact production bank, which we have not yet reproduced.

## Step 3 — POST-MERGE apples-to-apples (2026-06-15, DONE)
`run_filter_test.jl` now also reports `merge_catalog(h)` for each bank, so the comparison
uses the production exclusion/merge (no cross-scale duplicate inflation):

| M (M☉/h) | STANDARD merged | DENSE merged | dense/std |
|---|---|---|---|
| >1.23e12 | 23,403 | 31,698 | **1.35×** |
| >1.69e12 | 22,841 | 30,381 | 1.33× |
| >3e12 | 18,661 | 21,685 | 1.16× |
| >1e13 | 8,456 | 8,824 | 1.04× |
| >1e14 | 676 | 699 | 1.03× |
| total | 25,898 | 42,555 | 1.64× |

Two things the merge confirms:
1. **High-mass bins converge** (1e13: 1.04×, 1e14: 1.03×). The pre-merge high-mass "rise"
   (1.33× at 1e14) WAS cross-scale duplicates — they merge back to the same giants, as
   expected. So the dense bank adds NO spurious high-mass halos. The merge works correctly.
2. **Low-mass genuinely rises** but by **~1.33–1.35×**, not the 1.58× the raw counts
   suggested — the duplicate inflation is removed. (Total 1.64× is inflated by the dense
   bank also finding many sub-1.23e12 halos *below* Websky's completeness floor: dense has
   ~10.9k halos under 1.23e12 vs standard's ~2.5k — not Websky-relevant.)

## Bottom line (REVISED after reading filter_gen.py — see headline)
The deficit is a low-mass **shape** difference, and:
- **The filter bank is NOT the cause.** Our standard `filters_websky.dat` (rmincell 1.65,
  spacing 1.15) is exactly what Websky's own `filter_gen.py` produces. Proven by reading
  the reference generator, not inferred. The dense-bank test over-resolves relative to
  Websky and is therefore not a "match-Websky" experiment.
- **Leading explanation = abundance matching**, which Websky applies (README: "Due to the
  abundance matching performed on the halo catalogue there is a slight redshift dependent
  minimum halo mass … M_{200,M}") and our production catalog does NOT. The bin-by-bin ramp
  (0.41 low → 1.0 high) is the AM signature: identity where peak-patch already matches the
  N-body MF (high mass), boost above a fixed low threshold where it doesn't (low mass).
- **Caveat (don't overclaim):** AM is count-PRESERVING in total, so it can fix the low-mass
  *shape* but a residual *total*-count question (ours ~53M vs Websky ~110M) may remain and
  could instead be a z-range / completeness-floor / volume-normalization difference in the
  comparison. Needs the AM test to separate the two.

NOT a field/finder/lightcone/GPU/multi-res bug and NOT the filter bank (all verified). 
**NEXT (decisive): apply abundance matching to our catalog and re-run `compare_websky_mf.jl`.**
The AM target / procedure is in the reference repo
(`gw…/njcarlson/peakpatch/python/peakpatchtools/`).
