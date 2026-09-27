# Tiling invariance and coarse/fine splice: field and catalog tests (2026-09-26)

Plan item A5 (`docs/paper_comparison_plan_2026-09.md`), data for the §3 multires figures (F2 and
tiling-invariance). Commit at run time: 91e5a9c (src unchanged by this work).

**Setup (all tests).** Global grid N=480, Websky cellsize a=0.85221 Mpc/h (box 409.06 Mpc/h),
seed 12345, nbuff=16, so the tiled core region is NC=448 cells. Coarse grid fixed at M=N/block,
with block=12 as in production cf32 (coarse cell 10.23 Mpc/h, coarse Nyquist k_Nc=0.307 h/Mpc).
ntile ∈ {1,2,4,8} gives nsub = 448/ntile and nmesh = nsub+32 ∈ {480, 256, 144, 88}. The
production run has nmesh=414 and nsub=382, between the ntile=1 and ntile=2 cases here.
ntile=1 is one isolated tile with no internal faces, so it separates *split* error from
*tiling* error.

## What is bit-identical (and what is not)

| Quantity | Across tilings | Run-to-run (same tiling, same GPU) |
|---|---|---|
| White noise of every cell (Threefry, global index) | **bitwise identical**; also identical to the global-FFT noise (40/40 cells checked across 4 tilings, including buffer cells shared by neighbouring tiles) | bitwise |
| Coarse noise and coarse fields (depend on N, M, seed only) | **bitwise identical** by construction | bitwise |
| Residual noise ξ − spread(ξ_coarse) | **bitwise identical** (0/40 mismatches) | bitwise |
| Fine fields δ, ψ in a tile core | **not identical**: the residual convolution is isolated per tile, so it truncates differently for different tile extents. δ differs by 0.28–0.62% of σ_δ between tilings; ψx by 7–13% of σ_ψ | not compared directly; the raw peak set is bit-identical |
| Raw (pre-merge) peak list | 98.5–99.5% of peak *positions* shared; **0% of records bit-identical** (RTHL, ψ differ slightly; median \|Δln R_TH\| 0.1–0.7%) | **identical as a set; list order differs** (GPU) |
| Merged catalog | ~95–98% of halos matched; N(>M) equal to ≤0.5% below 1e14 | **as returned: 109/225k records differ**; with the canonical tie-break: bit-identical (0 differ) |

So the claim "same seed → same realization at any tiling" holds exactly for the **noise, coarse
field and residual**. The synthesized fine field is tiling-independent only up to the truncation
of the per-tile isolated residual convolution. That residual is small for δ and dominated by
the coarse-Nyquist band.

## (b) Field-level: split vs exact global FFT (`field_split_test.jl`)

Reference: exact periodic N³ FFT of the same noise with the same √P kernel convention.
Jobs 5693060 (block 12), 5693061 (block 6), 5693062 (block 24).

| block | ntile | δ rms err/σ_δ | δ corr | ψx rms err/σ_ψ | ψx corr |
|---|---|---|---|---|---|
| 12 | 1 | 5.10% | 0.99875 | 15.2% | 0.9897 |
| 12 | 2 | 5.11% | 0.99875 | 16.1% | 0.9881 |
| 12 | 4 | 5.12% | 0.99874 | 17.6% | 0.9856 |
| 12 | 8 | 5.13% | 0.99874 | 18.0% | 0.9849 |
| 6 | 1 / 4 | 6.84 / 6.85% | 0.9978 | 10.7 / 12.3% | 0.995 / 0.993 |
| 24 | 1 / 4 | 3.73 / 3.78% | 0.9993 | 20.3 / 23.9% | 0.981 / 0.973 |

Tiling-vs-tiling δ differences (block 12): 1v2 0.28%, 1v4 0.47%, 1v8 0.62%, 2v4 0.38%,
4v8 0.49% of σ_δ. The δ error vs distance to the nearest tile face is **flat**: 4.8–5.6% from 1 to
16 cells for every ntile. There is no seam at tile faces.

**Where the error lives in k** (error-power decomposition of the binned spectra; P ratio and
cross-correlation r(k) of stitched vs global on the same cosine-tapered NC³ region):
- δ: 55% of the error power is in 0.5–2 k_Nc and ~43% above 2 k_Nc (small r deficits at many
  high-k modes). Below 0.5 k_Nc there is ≤2%. The P ratio shows the splice signature
  **+11% at 0.5 k_Nc, −7% at 1.2 k_Nc**, with r ≥ 0.993 everywhere. It returns to 1.000 above
  ~4 k_Nc and matches to <1% below 0.1 k_Nc. **Identical for all ntile**: this is the coarse/fine
  split (tricubic interpolation of the coarse field plus the block-mean residual), not tiling.
- ψx: 85% of the error power sits at 0.5–2 k_Nc (**r(k≈1.2 k_Nc) = 0.82**) for ntile=1.
  Tiling adds a low-k part: 35% of the error power below 0.5 k_Nc at ntile=4, because the 1/k
  kernel is long-range and the isolated tile convolution truncates it. The ψ error grows with block
  (10.7% at block 6, 15% at 12, 20% at 24), which is consistent with the coarse-grid
  velocity-coherence criterion (cf≈32) already used in production.
- σ(R) for the top-hat filter-bank radii (the quantity that sets peak counts), stitched/global:
  0.999–1.000 at R ≤ 2 Mpc/h, 1.011 at 4.9, up to **1.035 at R ≈ 11 Mpc/h** (block 12). This is
  identical to <0.1% for all ntile. The split therefore slightly boosts σ at R ~ one coarse cell.
  Block 24 moves the bump to R ≈ 26 Mpc/h (+3.6%) and block 6 to R ≈ 5 Mpc/h (+3.0%).
- Estimator sanity: global/input P = 0.96–1.04 for 0.1 < k < 1 (the lowest bins are
  window- and cosmic-variance-limited in a 382 Mpc/h region).

**Paper consequence.** The multires field reproduces the exact FFT field with r ≥ 0.993 at all
k for δ. It carries a ±10% P(k) wiggle around the coarse Nyquist and a ≤3.5% σ(R) excess at
R ≈ one coarse cell. These are properties of the split, independent of tiling, and should be
stated in §3 along with the block-size dependence. The abundance-matching step absorbs the
mass-function consequence. The ψ decoherence at k_Nc is the known reason for cf≈32.

F2 data: `results/field_N480_b12/field_slice_ntile4_k50.f32` holds four 448×448 Float32 planes at
global z-index 50: coarse (interpolated), self (isolated residual), total, exact global. It
crosses tile faces at x,y = 112, 224, 336. See `field_slice_README.txt`. The same files exist for
blocks 6 and 24.

## (a) Catalog-level (`catalog_tiling_test.jl`, GPU, job 5695201)

z=0 snapshot, 2LPT, finecell filter bank, `run_multitile_split(...; coarse_grid=40, use_gpu=true)`
on one L40S; 717k raw peaks and 225k merged halos per tiling.

### (i) Pre-merge peak lists (vs ntile=1; peak positions are cell centres, so compared exactly)

| ntile | raw peaks | positions found | …within 0–2 cells of a face | …>16 cells from faces | median \|Δln R_TH\| |
|---|---|---|---|---|---|
| 2 | 716,814 (−0.03%) | 99.55% | 98.65% | 99.7–99.8% | 0.11% |
| 4 | 716,800 (−0.03%) | 99.07% | 98.16% | 99.5–99.6% | 0.40% |
| 8 | 716,482 (−0.07%) | 98.54% | 97.87% | 99.2% | 0.73% |

The mismatch is mildly concentrated at tile faces (about 2× more missing peaks within 2 cells
than >32 cells away) and grows with the number of faces. It is present everywhere, though,
because the fine field differs at the 0.3–0.6% level throughout each tile.

### (ii) Merge-order effects (same raw list, different input order)

`_merge_impl` (src/Merger/Merger.jl) sorts with `sortperm(r; rev=true)`, which is stable, so
**halos with equal R_TH are visited in input order**. Ties are common: **22% of raw peaks
(160k/717k) share their R_TH exactly with another peak**. They sit on the discrete shell radii
(1, √2, √3, 2, √5, √6 cells = 0.852…2.087 Mpc/h, i.e. M ≈ 2e11–3.3e12 M⊙/h) from the
uninterpolated branch `RTHL = rad[m0-1]` in RadialShell.jl. Measured on a fixed raw list:
- 3 random input permutations vs a canonical order (−R_TH, then x, y, z): **6,300–6,900 records
  differ** in the merged catalog (~3% of 225k; halo count changes by up to ±80).
- The as-returned order also differs from canonical by 6,300–6,600 records.
- On GPU the raw list **order changes run to run**, although the set is bit-identical. Repeat test
  (job 5695214, ntile=4 run twice in one process): the as-returned merged catalogs **differ by
  109 records** (224,728 vs 224,741 halos). With the canonical order, **0 records differ**, and
  the canonical ntile=4 merge has the same halo count as in job 5695201 (224,730). This is consistent with the performance fork's few-tens-in-10M
  difference between 1/2/4 GPUs, where only the tile order changes and within-tile order is
  preserved.

**Proposed src fix (not applied):** in `_merge_impl` use a total order,
`order = sortperm(eachindex(r); by = i -> (-r[i], x[i], y[i], z[i]))`. The merged catalog then
depends only on the raw peak *set*: it becomes bit-reproducible run to run, across GPU counts and
across tile completion order. All catalog comparisons below already use this canonical order
(the script pre-sorts the input). The production catalogs were made with the as-returned order;
the effect is ~3% of records swapped among tied, overlapping small halos at fixed total count
(±0.04%). After AM this is a relabelling within tied masses. It is not a bias, but it is a
reproducibility defect.

### (iii) Merged catalogs across tilings (canonical merge, vs ntile=1)

Match: nearest unclaimed halo within max(1 cell, 0.2 R_TH), capped at 3 cells.

| ntile | N ratio | matched 1–3e12 | 3e12–1e13 | 1e13–3e13 | 3e13–1e14 | >1e14 | per-halo ln M rms | ψx diff rms/σ |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.99935 | 0.983 | 0.972 | 0.963 | 0.953 | 0.962 | 2.6–6.6% | 0.077 |
| 4 | 0.99930 | 0.964 | 0.943 | 0.917 | 0.907 | 0.918 | 4.2–7.5% | 0.118 |
| 8 | 0.99939 | 0.948 | 0.915 | 0.882 | 0.863 | 0.874 | 5.4–9.0% | 0.131 |

The median ln(M ratio) is 0 to ±0.002 in every bin, so there is **no systematic mass shift**.
Cumulative N(>M) ratios vs ntile=1 are within ±0.2% for M < 1e13, within ±0.5% up to
1e14, and within Poisson noise above (1–4% with ≤1000 halos). Matched positions are exactly
equal for 99.0–99.7%. The unmatched fraction rises toward tile faces (ntile=8: 8.5% within 2
cells, 4.7% at 16–32 cells). High-mass matching is lower because large, flat peaks move and
merge differently under small field changes. The statistics are unchanged.

Summary for the paper: tiling changes **individual** halos at the few-% level (which halo
survives exclusion, R_TH to ~1–9% rms). It leaves **statistics** unchanged at ≤0.2–0.5%
(counts) with no mean mass shift. The noise and the long-wavelength field are exactly
tiling-independent.

## Files

- `field_split_test.jl`, `catalog_tiling_test.jl`, `run_field_split.slurm`, `run_catalog_tiling.slurm`
- `results/field_N480_b{6,12,24}/`: `field_summary.csv`, `field_pk_ratio.csv` (δ and ψx: P_split,
  P_global, P_input, ratio, transfer, r, nmodes), `field_sigmaR_ratio.csv`,
  `field_edge_profile_ntile*.csv`, and the F2 slices
- `results/repeat_ntile4_N480/`: repeat-determinism run (ntile=4)
- `results/catalog_N480/`:
  - `raw_ntile*.f32` and `merged_ntile*.f32`: 10 Float32 per halo (x y z vx vy vz RTHL vx2 vy2 vz2), canonical order
  - `raw_compare.csv`, `catalog_match_summary.csv`, `catalog_unmatched_vs_face_ntile*.csv`, `catalog_cumcounts.csv`
- Logs: `/home/yguan/scratch/websky_6144/logs/tiling_{field,cat}_<jobid>.out`

Cost: each field job ran 47–93 min, CPU-bound in the serial `_isolated_convolve` k-loop. The
catalog job took 4 min on one L40S.
