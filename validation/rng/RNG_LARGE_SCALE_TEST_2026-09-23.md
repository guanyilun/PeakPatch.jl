# RNG large-scale null test — 2026-09-23

## Question
J. R. Bond: the Fortran peak-patch generator (`random.f90`, 48-bit LCG) was built so
that drawing a huge number of deviates does not create spurious correlations on large
(angular) scales — `rans()` splits the 2^46 cycle into disjoint streams via jump-ahead.
Is that an issue with the RNG we use?

## What production uses
- **Threefry-2x64-20** (Random123), counter-based: each cell's deviate is a pure function
  of (seed, global linear index) (`RandomField._threefry_gaussian`, Box–Muller on one
  Threefry call per pair of cells along x). No state, no per-tile streams: tiles, buffers
  and the coarse grid *recompute* the same numbers for the same global cells.
- Coarse grid = exact normalized 12³ block-mean of the same fine noise
  (`MultiResolution._downsample_noise(6144, 512, 12345)`); fine residual subtracts it
  exactly → coarse/fine modes consistent.
- `LCG.jl` (Fortran port) is used only with `fortran_compat`.

Structural argument: the classic failure modes (overlapping/correlated streams; period
exhaustion; LCG lattice/lag structure mapping onto rows/planes) do not apply — there are no
streams, the counter space is 2^64 pairs (we use ~1.2e11), and Threefry has no linear
structure (passes BigCrush from ~13 rounds; we use 20).

## Test (`rng_large_scale_test.jl`, `run_rng_test.slurm`)
Everything is compared against independent **Xoshiro** (Julia default) white noise pushed
through identical statistics and geometry, so box replication, pixelization and
mode-counting cancel and only the generator is tested.

Statistics per field: mean, variance, skewness, kurtosis; shell P(k) in log bins; low-k
and axis-aligned power (row/plane-fill artifacts would be axis-aligned); lag-ξ along each
axis; **max |ξ| over every 3D separation**; angular C_ℓ (ℓ bins 2–5, 6–10, 11–20, 21–40,
41–80) of HEALPix maps in radial shells r/L = 0.30–0.50, 0.50–0.75, 0.75–1.00, 0.30–1.00
from an observer at the periodic-box origin (= the 8-octant full-sky production geometry,
including replication).

- **A. Ensemble (generator soundness)** — 40 Threefry seeds vs 40 Xoshiro seeds;
  two-sample z and KS per statistic.
- **B. Production realization** — the exact production coarse noise (seed 12345, 2.3e11
  deviates) ranked within 40 Xoshiro 512³ fields. Says whether our seed is a typical draw
  (cosmic variance), not whether the generator is sound.

## Results

### A — PASSED (job 5635142)
| ensemble | # stats | max \|z\| (chance ≈) | Bonferroni KS p |
|---|---|---|---|
| coarse, N=1536 → M=128 (block 12³), Nside 32, incl. C_ℓ | 61 | 2.78 (≈2.9) | 1.00 |
| fine 256³ subcubes at random offsets in the N=6144 grid (counters to 2.3e11) | 45 | 2.59 (≈2.8) | 1.00 |

7/106 statistics at |z|>2 (≈5 expected), scattered across unrelated statistics/axes, no
pattern. Low-ℓ C_ℓ ratio Threefry/Xoshiro at ℓ 2–10: 0.91–1.12 across shells, no
systematic excess. Full tables: `results/REPORT_A_job5635142.md`, `results/ensemble_*.csv`.

### B — pending (job 5637544)
First attempt (5635142) generated the production coarse noise (529 s, cached at
`scratch/websky_6144/rng_test/coarse_threefry_N6144_M512_seed12345.f32`) then segfaulted in
threaded FFTW (`spawn_apply`) on the 512³ transform under Julia 1.12; rerun with
`RNG_PART=B RNG_FFTW_THREADS=1`.

## Conclusion
No detectable difference between the production Threefry noise and an independent
generator on any tested scale, including large angular scales in the replicated
full-sky geometry. The concern that motivated the Fortran stream partitioning does not
apply to a counter-based generator.

## Lessons / caveats
- A first single-realization-vs-4-nulls smoke run showed a 1.47× low-ℓ "excess" in one
  shell — a test-design artifact: low-ℓ C_ℓ in a replicated box depends on a few
  fundamental 3D modes, so single-realization distributions are heavy-tailed. Generator
  soundness needs ensemble-vs-ensemble.
- **Separate, not an RNG issue, untested:** the full sky is 8 views of ONE periodic
  5236 Mpc/h box, so structures repeat on the sky beyond χ ≈ L/2 = 2618 Mpc/h (z ≳ 1.2).
  Same as Websky (observer at box centre, replicated). Potential large-angle effect worth
  its own check if large-angle fidelity matters.
- Healpix `map2alm` default `niter=3` costs ~2 s/call regardless of size; `niter=0`
  (0.3 s) used — identical treatment for both ensembles, ℓ ≤ 80 ≪ pixel limit.
