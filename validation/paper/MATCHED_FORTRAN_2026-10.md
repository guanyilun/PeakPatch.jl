# Same-field Julia vs Fortran comparison (2026-10-02)

**Question.** The full-sky Tier-A comparison (TIERA_FULLSKY_2026-10-02.md) finds the full Websky catalogue
~5% less clustered than v4fs: ξ(3–15 Mpc/h) Websky/ours = 0.948 ± 0.002, b ≈ 0.95–0.96 below 8e13.
How much of that is a difference between the two halo-finding codes?

**Method.** One linear field, both codes. Julia writes its exact global field; Fortran `hpkvd` reads it
(`ireadfield = 1`); both chains then run on identical input, so cosmic variance cancels in every ratio.

- Box: z = 0.7 snapshot, production cell 0.8522 Mpc/h, N = 1056, ntile 4, nbuff 26, non-periodic cores
  (the Fortran tiling), production filter bank (`filters_websky_finecell.dat`).
  Config: `configs/matched_fortran.toml`.
- Fortran: the local clone `~/work/peakpatch`, built unmodified in scratch at nmesh 303
  (`fortran_matched/src`; link fixes only: `-lstdc++`, `-DOMPI_SKIP_MPICXX`). The clone is not guaranteed to be
  the binary Websky ran in production (see memory note on local patches; `hpkvd` itself was not hand-edited).
- Analysis: 852 Mpc/h interior cube, rank-matched to Tinker08. Cross-bias with the shared linear field at
  k < 0.1 h/Mpc, with 8-subvolume jackknife errors. ξ by pair counts against region randoms.

## 1. Result: Fortran is 2% less clustered on the same field

`results/matched_fortran.txt`. F = Fortran hpkvd + merge_pkvd. FJ = Fortran raw → Julia merge.
JE = Julia exact field.

| | F/FJ | F/JE |
|---|---|---|
| ξ(3–15), M>5e12 | 0.9999 | **0.980** |
| b_E, M>5e12 | 0.9999 ± 0.0001 | 0.992 ± 0.001 |
| b_E, 1.3–2e13 | 1.000 | 0.981 ± 0.003 |
| b_E, 5–8e12 | 0.999 | 1.000 ± 0.002 |

**The merge is exonerated.** The Julia merge on Fortran raw peaks reproduces `merge_pkvd` to 1e-4.

**Fields, displacements and peak selection agree.** Density rms is identical (3.43969).

- Matched peaks agree in Eulerian position to a median 0.002 Mpc/h.
- 99.3% of Fortran peaks have a Julia peak in the same cell (99.2% in the same filter).

**R_TH carries the whole difference.** Give the Julia raw peaks Fortran's R_TH, then merge and abundance-match
them identically (`results/matched_rth_hybrid.txt`). The result reproduces F:

- ξ ratio 0.9995;
- b_E 0.999–1.001 in every mass bin.

## 2. Two bugs in the Fortran `get_homel` (`src/hpkvd/peakvoidsubs.f90`)

Both are coding errors visible in the source. Both make R_TH larger for some peaks in the lowest-mass part of
a mass bin, and smaller in overdense environments. The net effect is less clustered halos at fixed abundance.

### (a) Dead outward search (l. 447–448 vs 497–498)

The gradient-at-R_f loop `do jp=2,npart` reuses `jp`. The outward search `j0=jp+1; do jp=j0,npart` therefore
never executes.

- A peak with Fbar ≥ fcrit at the first tested shell m0 (first r² > ir2min ≈ (1.75 R_f/a)²) cannot grow beyond
  rad(m0).
- In 97% of the cases where Julia is larger, the Fortran R_TH sits exactly on a lattice shell √n·a.
- Affected peaks (Julia's R beyond the Fortran limit): 23% at 5e12–1.3e13, 14% at 1.3–5e13, 10% above 5e13.
- Julia/Fortran R_TH in that class: median 1.30 / 1.18 / 1.11.

### (b) Strain kernel truncated by a clobbered `mupp` (l. 471–476, then 600–605)

The R_f-gradient block sets `mupp=mrf` and then loops `do m1=mrf+1,mupp`, an empty range. That overwrites
the `mupp` built in the first loop. The inward branch then computes `muppnew` over `m0+1..mupp`, also empty, so
`mupp = m0`.

- The shell-averaged strain at each tested radius therefore uses only shells inside it. The intended window is
  two-sided, out to r + 2 cells.
- This underestimates the strain trace and the ellipticity. The collapse redshift comes out higher, so R_TH is
  larger.
- Julia uses separate variables (`mupp_rf`), i.e. it does what the Fortran evidently intended.

**Per-peak harness evidence.** 132 peaks in tile (2,2,2) were traced
(`results/matched_peak_traces/`, `matched/fortran_dbg_patch.py`, `matched/radialshell_dbg_hook.patch`).

- The per-shell sums Σ ηL·xK agree to float32 precision.
- The kernel range differs. For peak (−95.02, −103.54, −108.66), R_f 2.459, Fortran uses shells 3–10 (m0 = 10)
  while Julia uses 3–24 (r + 2).
- As a result, the strain trace is 1.20 (F) vs 2.17 (J), against Fbar = 2.35.
- e_v is 0.163 vs 0.247.
- R_TH is 2.427 vs 2.368 Mpc/h.
- Fortran's kernel is truncated for 129 of the 132 peaks. Of the 44 with |ΔR| > 2%, 38 have Fortran larger,
  mostly R_f ≤ 1.86.
- Catalogue-wide, among peaks unaffected by the cap, Julia's R_TH is more than 2% smaller than Fortran's for
  51% of R_f = 1.41 peaks, 39% at 1.62, and 28% at 1.86. That falls below 1% by R_f ≈ 3.7. The median offset
  is 0.07 cells.

## 3. Decomposition (same field, `results/matched_rth_hybrid.txt`)

Hybrids: Fortran R_TH applied to a subset of the Julia raw peaks.

| Fortran R_TH applied to | F/hybrid ξ(3–15) | share of the 2.0% |
|---|---|---|
| none (plain JE) | 0.980 | — |
| cap class only, bug (a) | 0.985 | ~0.5 point |
| everything else, mainly bug (b) | 0.993 | ~1.4 points |
| all matched peaks | 0.9995 | all |

- Bug (a) carries the bias difference above 1.3e13 (F/JEcap 0.993–1.000 there).
- Bug (b) mostly moves ξ, through which small-R_f peaks cross the abundance-matching threshold and survive
  exclusion and volume reduction.

## 4. What the code differences alone explain

The two bugs account for ~2.0 of the 5.2-point full-sky ξ gap (F/JE), in the right direction. §6–7 locate the
rest.

**Framing (neutral):** these are differences between our code and the local Fortran clone. Whether Websky's
production run had them is not established. No upstream contact without the user's go-ahead.

## 5. Side finding: GPU shell-table overflow at nbuff ≥ 26 (fixed)

`run_multitile_split` returned ~10× too few halos for every nbuff = 26 configuration (periodic or not, cf 22 or
24). The cause: the GPU post-process kernel keeps each peak's radial profile in shared memory sized for
`_MAX_SHELLS_GPU = 512` distinct shells (`ext/CUDAExt.jl`).

- nhunt = nbuff − 1 = 25 gives 523 shells. The overflow is silent, and ~80% of peaks then fail to collapse.
- nbuff 25 (482 shells) fits, so production v1–v4 (nbuff ≤ 25) is unaffected.
- Fix: `_check_shell_capacity` now errors at both kernel launch sites. Tested: nbuff 26 errors; nbuff 25 output
  is identical (6,907,598 raw → 2,997,974 merged).

## 6. Code × filter bank × field construction on the same field (`results/matched_bank.txt`)

Labels:
- F2: Fortran with the bank the Websky paper documents (R_f,min = 2 a_latt → 36 Mpc;
  `filters_paper_2cell.dat`, 1.704 → 24.26 Mpc/h, 20 × 1.15).
- JE2: Julia exact with the same bank.
- JSP: our production path (GPU multires split, periodic cores, nbuff 25, cf 22 → block 12 as in production
  6144/(16·32)) on the same δ. Threefry noise is indexed globally, so the field is identical.

| ratio | ξ(3–15) | b_E M>5e12 |
|---|---|---|
| F/JE (code) | 0.980 | 0.992 ± 0.001 |
| F2/F (bank, Fortran) | 0.995 | 0.997 ± 0.002 |
| JE2/JE (bank, Julia) | 0.984 | 0.992 ± 0.001 |
| JSP/JE (split vs exact field) | **1.037** | **1.017 ± 0.002** |
| **F2/JSP (documented Websky vs our production)** | **0.940** | **0.973 ± 0.002** |
| full sky, Websky/ours v4fs | 0.948 ± 0.002 | ≈ 0.95–0.96 below 8e13 |

**The gap is reproduced in one box on one field.** F2/JSP b_E is 0.956–0.986 below 8e13. Rough budget in ξ
points: split ~3.7, code ~2.0, bank ~0.5. These are not exactly additive, and the box is a snapshot, not the
lightcone.

## 7. The multires split changes halo clustering (`results/matched_split_scan.txt`)

Production split path vs the exact global field, same δ, periodic, nbuff 25. The exact field is the ground
truth, so any departure from 1 is a split error.

| cf | block | ξ(3–15) JS/JE | b_E JS/JE M>5e12 |
|---|---|---|---|
| 4 | 66 | 0.986 | 0.930 ± 0.001 |
| 8 | 33 | 0.968 | 1.004 ± 0.001 |
| 11 | 24 | 0.974 | 1.010 ± 0.001 |
| 22 | 12 (production block) | 1.037 | 1.017 ± 0.002 |
| 33 | 8 | 1.058 | 1.018 ± 0.002 |

**The effect is not a monotone convergence.**
- At block 66 the coarse Nyquist (k ≈ 0.056 h/Mpc) falls inside the k < 0.1 bias band, giving large-scale
  decorrelation (b_E 0.93).
- At small blocks, clustering is enhanced.

**Caveat:** this box has smaller tiles than production (nsub 264 vs 384). The isolated-convolution error
depends on tile size too.

**Next:**
- (a) field-level split error vs k and tile position (`compare_fields_split`);
- (b) split vs exact at the production tile size (N 1536, nsub 384);
- (c) fix, then decide on rerunning production.

## 8. Where the split differs: peaks, radii, displacements (`results/split_vs_exact_n1056_cf22.txt`)

Same field, production split path (block 12) vs exact, raw peaks inside the analysis region.

**The two fields differ at the cell level, far more than Fortran vs Julia did.**
- Only 85% of exact peaks have a split peak in the same cell (76% in the same filter). Fortran vs Julia
  exact: 99.3%.
- R_TH of matched peaks differs by more than 2% for 74% of them (median split/exact 0.993, p05 0.81,
  p95 1.16).
- Eulerian positions differ by a median of 0.84 Mpc/h (p90 1.39), against ~6 Mpc/h displacements.
  Fortran vs Julia exact: 0.002. That is a ~14% displacement error, consistent with Audit C's ~15% ψ rms
  error near the coarse Nyquist.

**Most of the clustering change comes through R_TH.** Exact peaks given the split R_TH (hybrid E_Rsplit)
already reproduce ~70% of the effect:

| | S/E | E_Rsplit/E |
|---|---|---|
| ξ(3–15) | 1.038 | 1.026 |
| b_E M>5e12 | 1.017 ± 0.002 | 1.013 ± 0.002 |

The rest comes from peak selection and positions.

**The error sits in ψ1, in a band just above the coarse Nyquist** (`results/split_field_error_cf{22,8}.txt`;
2³ central tiles, core cells vs the exact global fields).

| block | k_Nyq,coarse | δ rms err | ψ1 rms err | ψ2 rms err | worst ψ1 band: P_err/P_ref, r(k) |
|---|---|---|---|---|---|
| 12 (production) | 0.307 h/Mpc (λ 20.5) | 4.5% | **14.1%** | 11.5% | **0.30 at k 0.38**, r 0.85; 0.11 at 0.54 |
| 33 | 0.112 h/Mpc (λ 56) | 2.9% | **23.0%** | 12.3% | **0.23 at k 0.13**, r 0.88; 0.14 at 0.19 |

- **δ is accurate.** P_err/P_ref ≤ 1.2% at every k, with r ≥ 0.994. So peak finding sees nearly the right
  field.
- **ψ1 is badly wrong in a band just above the coarse Nyquist**, and that band moves with the block size. At
  block 12 it sits at λ ≈ 11–20 Mpc/h, which are the scales of the shell-averaged strain that sets ellipticity
  and R_TH.
- This is the uncompensated aliased part of the coarse/fine splice. The D/T compensation fixes only the
  diagonal part below the coarse Nyquist. ψ ∝ δ/k weights the coarse-dominated band, so the aliasing
  residual that is ~1% in δ becomes ~30% in ψ1.
- **ψ2 errors are largest at low k** (22% at k 0.033). That is expected from the tile-periodic 2LPT, and
  matters less.
- **The missing outer residual is not the cause.** An extended residual shell (nshell 24) changes nothing,
  and the error is flat with distance from the tile edge.

**Interpretation.** The production field's displacements are wrong by ~14% rms, concentrated at ~10–20 Mpc/h.
That changes the strain, hence R_TH (§8), and gives the ~0.8 Mpc/h Eulerian position errors. It is the
dominant cause of the split's +3.7% ξ.

The block-33 ξ/b behaviour (§7) is the same error moved into the k < 0.1 bias band, where it decorrelates
(b_E 0.93 at block 66).

**Fix direction (not yet done).** Make the splice consistent for ψ:
- either subtract the same interpolated coarse noise that is added back, instead of the nearest-neighbour block
  mean;
- or split in k-space with an anti-aliasing filter, as in MUSIC.

Then re-measure the ψ1 band error and the halo ratios.

**Production tile size confirms it** (`results/split_vs_exact_n1536_cf32.txt`). Setup: N 1536, nsub 384,
block 12. The split is `merge_ab_z07.toml`; the exact field is `configs/exact_n1536.toml`, on the same δ.

| | ξ(3–15) S/E | b_E S/E | same cell | R_TH > 2% off | Eulerian Δ (median) |
|---|---|---|---|---|---|
| this box | 1.035 | 1.017 ± 0.002 | 85.4% | 74% | 0.82 Mpc/h |

The result is identical to the nsub 264 box. The error does not depend on tile size, and the v4 catalogues
carry it.

## 9. Fix: Gaussian long/short handoff (`[run] gaussian_split`, opt-in)

Design and literature: `docs/split_psi_fix_literature_2026-10.md`.

**The construction.** G = exp(−k² r_s²), with r_s = 9 fine cells (0.75 block).

- δ = interp(G·δ_coarse) + isolated[√P·residual + (1 − G)·√P·spread(block mean)]
- ψ = interp(G·ψ_coarse) + isolated Poisson(δ_short)

CPU and GPU (`ext/CUDAExt.jl: isolated_poisson_psi_gpu`). The default is off, and the off path is unchanged.

**Field level** (`results/split_fix_proto_cf22_deltafix.txt`): δ rms error 4.5% → 0.3%; ψ1 13.8% → 3.8%. The
remaining ψ error is a low-k (k < k_N/2) bulk term from the isolated solves, at the same level as the original
splice.

**Halo level** (same field, GPU production path vs exact; `results/split_vs_exact_n1056_cf22_gsplit.txt`):

| split vs exact | original | ψ from δ only | **gaussian_split (δ + ψ)** |
|---|---|---|---|
| ξ(3–15), M>5e12 | 1.037 | 1.024 | **1.001** |
| b_E M>5e12 | 1.017 ± 0.002 | 1.007 ± 0.001 | **1.000 ± 0.0005** |
| b_E per bin (5e12–2.5e14) | 0.99–1.02 | 1.00–1.01 | **0.999–1.002** |
| peaks in the same cell | 85.4% | 87.7% | **99.3%** |
| R_TH off by >2% | 74% | 53% | **2.9%** |
| Eulerian Δ median | 0.84 Mpc/h | 0.31 Mpc/h | 0.23 Mpc/h (low-k bulk ψ) |

**Verdict.** With gaussian_split the production path reproduces the exact-field halo catalogue on the same δ:
clustering to 0.1% and bias to 0.05%. Peak selection agrees as well as Fortran vs Julia exact did (99.3%).

**Open:**
- the low-k ψ bulk error (0.23 Mpc/h), which matters for velocity products such as kSZ;
- 2LPT (tile-periodic, unchanged);
- GPU timing;
- the user's decision on a production rerun (v5).

## Files

- `matched/matched_julia.jl` (field / exact / split), `matched/gen_fortran_inputs.py`.
- `matched/matched_compare.jl`, `matched/matched_rth_diag.jl`, `matched/matched_rth_shell.jl`,
  `matched/matched_rth_hybrid.jl`.
- `matched/matched_peak_dbg.jl` + `matched/radialshell_dbg_hook.patch` + `matched/fortran_dbg_patch.py`.
- `matched/matched_bank_compare.jl`, `matched/matched_split_scan.jl`, `matched/split_vs_exact.jl`,
  `matched/split_field_error.jl`.
- Configs: `configs/matched_fortran_bank2.toml`, `configs/split_iso_*.toml`, `configs/split_scan_cf*.toml`
  (`matched_split_periodic.toml` = the nbuff 26 overflow case).
- Run data: `/home/yguan/scratch/websky_6144/fortran_matched/{run,run_dbg}`.
