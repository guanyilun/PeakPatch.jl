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

## 4. What this does and does not explain

- **Explained:** the code differences account for ~2.0 of the 5.2-point full-sky ξ gap, in the right direction.
- **Unexplained:** ~3.2 points come from something outside this box. Candidates:
  - the lightcone and per-tile threshold;
  - periodic production tiling;
  - inputs (filter bank: −1.7 points in the Julia-only bank A/B, not additive with the above; P(k));
  - Websky's production binary differing from this clone.
- **Framing (neutral):** these are differences between our code and the local Fortran clone. Whether Websky's
  production run had them is not established. No upstream contact without the user's go-ahead.

## 5. Side finding: the non-periodic split path

`run_multitile_split` with `periodic_cores = false` (same field geometry, cf 22, block 12) returned 252k merged
halos vs 2.58M exact. Its tile offsets assume the periodic layout (`i0 = (it−1)·nsub + 1 − nbuff`). Production
uses `periodic_cores = true`, so production is unaffected. Under investigation.

## Files

- `matched/matched_julia.jl` (field / exact / split), `matched/gen_fortran_inputs.py`.
- `matched/matched_compare.jl`, `matched/matched_rth_diag.jl`, `matched/matched_rth_shell.jl`,
  `matched/matched_rth_hybrid.jl`.
- `matched/matched_peak_dbg.jl` + `matched/radialshell_dbg_hook.patch` + `matched/fortran_dbg_patch.py`.
- Run data: `/home/yguan/scratch/websky_6144/fortran_matched/{run,run_dbg}`.
