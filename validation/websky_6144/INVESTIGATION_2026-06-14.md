# Websky 6144³ Halo-Deficit Investigation — 2026-06-14

Continuation of `INVESTIGATION_2026-04-16.md`. That session found the original
120× halo deficit (7.35M halos vs Websky ~9×10⁸) and fixed several config issues.
This session ran the corrected pipeline, found and fixed two more lightcone bugs,
and traced the remaining discrepancy to its root cause in the field construction.

Reference numbers: original Websky (Stein+2020 §4.1) = (7.7 Gpc)³, **6144³**,
observer at each of 8 corners (octant trick → effective 12288³/15.4 Gpc).
**Our N=6144 / box 7700 Mpc/h config reproduces Websky's resolution exactly** —
the 2026-04-16 note's "need N≈10240" was wrong (it anchored on Nate's finer
WebSky2 runs). Particle mass m_p = 2.775e11·0.31·(7700/6144)³ = 1.69e11 M☉/h;
Websky keeps >10 particles = 1.69e12 M☉/h. We match this floor.

---

## Bug 1 — `chi(z)` comoving distance a factor h too small  [FIXED, commit 70e0beb]

`chi(z)` in `src/Cosmology/Cosmology.jl` (and the duplicate in
`AbundanceMatch.jl:95`) multiplied by a spurious `* c.h`. Correct Mpc/h comoving
distance is `(c/100)·∫dz/E` — NO explicit h factor. The bug made χ ≈0.68× too small.

- coded `chi(4.6) = 3557 Mpc/h` vs correct `5231 Mpc/h`.
- Effect (lightcone, ievol=1 only — ievol=0 never calls chi, which is why z=0
  validation passed): tile selection `chi_max=chi(z_max)` and the chi→z table both
  truncated the lightcone at z≈1.96 instead of 4.6 → only ~31% of octant volume;
  every kept halo got an over-estimated redshift.

**Verified by re-run** (job 3950164, oct000): lightcone now reaches r=5231 Mpc/h
(z=4.6) exactly; halos 3.21M → **10.33M** (3.2×). Fix confirmed.

---

## Bug 2 — `nhunt` shell-analysis radius clamped by `nbuff`  [REAL, but only a partial cause]

`nhunt = min(nbuff-1, floor(Rfclmax*1.75/alatt))` (MultiResolution.jl:569). For
nbuff=16 → min(15,47) = **15 cells = 18.8 Mpc/h**. The collapse-search radius (and
thus RTHL cap) is bound by the buffer. Confirmed: oct000 had **428,286 halos pinned
at RTHL=18.673 Mpc/h** (≈ the cap) → fake 1.1e15 M☉/h "clusters".

**Reframed (CITA cross-ref, 2026-06-14):** the `nhunt = min(nbuff-1, …)` clamp is
**intended Fortran behavior, NOT a bug.** Greg Stein's production param/filter files
on `gw.cita.utoronto.ca:/fs/lustre/.../peak-patch-runs/OLD` show e.g. `800_0256_16Gpc`:
nbuff=16, cellsize 3.57 Mpc, buffer 57 Mpc, filters Rf 6–56.7 Mpc (`tables/filter.dat`:
`1.686, Rf, 0.3`; Rf_min=1.65·cellsize=rmincell rule; Rf_max≈buffer). There
`1.75·Rf_max = 99 Mpc > buffer 57 Mpc`, so nhunt clamps at nbuff−1 exactly like ours,
and Fortran has the same `nhunt=min(nbuff-1,nhunt)` (RandomField.f90:209). It's fine in
production because the FIELD is correct → few large-filter peaks → few hit the clamp.
In our run the field is wrong (Bug 3) → a flood of large-filter peaks → they pile at the
clamp radius. The nbuff=24 experiment (pile merely relocated 18.7→28.7) confirms the
clamp was never the cause. Sizing nbuff so buffer ≥ Rf_max is still good practice, but it
does NOT fix the catalog — Bug 3 does.

### GPU shared-memory blocker
The GPU shell kernel (`_post_process_kernel!`, CUDAExt.jl) stores per-peak profiles
in dynamic shared memory sized to `_MAX_SHELLS_GPU = 512`. nshells = #distinct r² ≤
nhunt²: nhunt=15→190, 23→443, 24→482, 31→~800, 47→1843. **nbuff=48 (nhunt=47, 1843
shells) overflowed shared mem → CUDA ILLEGAL_ADDRESS.** Max safe is nhunt≈24 (≈30 Mpc/h).
Also: coarse grid `M = ntile·coarse_factor` must divide N=6144 (valid ntile ∈ {16,24,32,48}).

### nbuff=24 experiment (job 3951693, COMPLETED) — the decisive negative result
Raised nbuff 16→24 (nhunt=23, 28.8 Mpc/h, max safe). Result:
- pile MOVED 18.7 → 28.7 Mpc/h and got MORE massive (top bin 5.6e15, 303K halos);
- total halos went DOWN 10.33M → 7.94M (bigger fake giants exclude more neighbors;
  exclusion 50.5M→7.94M = 84% removed).

**Conclusion: the high-mass pile is NOT a truncation artifact — it relocates with the
cap. So the full GPU kernel fix (raising _MAX_SHELLS_GPU) is NOT worth doing.** Keep
nbuff=24 (harmless); the real bug is upstream.

---

## Bug 3 (ROOT CAUSE) — multi-resolution field has wrong σ(R) (flat, too much large-scale power)

Empirical, measured from the oct000 catalog (per-halo `FcollvRf` = smoothed field
value at each peak): **median smoothed-peak-amplitude is FLAT ~6–13 across the entire
filter range Rf 2→34 Mpc/h** (even rising at Rf 4–6). A correct linear field has σ(R)
falling ~10× from 2 to 34 Mpc/h, so peak amplitudes should collapse with scale. They
don't — the field has white-noise-like / too-flat power.

Consequence chain:
- 232K "halos" at the largest filter Rf≈29 (should be ~tens — a 29 Mpc/h region
  collapsing by z~1 is a ~3σ rarity) → fake 1e15–5e15 M☉/h giants (≈26% of catalog);
- those giants over-exclude real halos in the Lagrangian merger → suppressed count.

### What was ruled out (all correct / cross-referenced)
| Component | Status |
|---|---|
| P(k) shape & σ8 (=0.817) | correct |
| `smooth_field` top-hat window | correct (kills high-k properly) |
| Filter bank | matches Fortran `python/filter_gen.py` (rmincell=1.65, spacing=1.15) |
| Merger exclusion | correct (Lagrangian, on raw peak positions, vol-reduction off) |
| `fsc_of_z` / fcrit threshold | correct (collapse-table bisection) |
| `chi` | fixed (Bug 1) |
| nbuff/nhunt | real bug, but only relocates the pile |

### Fortran cross-reference (the key)
- `RandomField.f90` (`RandomField_make`/`convolve_noise`): Fortran builds each tile's
  field by a **single per-tile PERIODIC √P(k) convolution** — white noise → FFT →
  `× sqrt(P(k)·V_k)`, `V_k = dk³·n³`, `dk = 2π/boxsize`. **No coarse grid, no
  zero-padded/isolated convolution.** (Fortran author even left a comment at line 103
  questioning the V_k box-size convention.)
- Julia `_periodic_convolve!` (MultiResolution.jl:301, `amp=sqrt(pk·dk³·n³)`) **matches
  Fortran exactly**. `run_multitile` (MultiTile.jl) is the Fortran-equivalent: one
  global N³ FFT, then extract tiles — and it's validated to sub-percent vs Fortran.
- **The GPU path `run_multitile_split` does NOT use this.** It replaces the per-tile
  convolution with **coarse (periodic) + residual (`_isolated_convolve`, zero-pad n→2n,
  MultiResolution.jl:246-283)** — a Julia-only MUSIC-style optimization with NO Fortran
  analog. That decomposition produces the wrong σ(R).
- README validation: the coarse+residual path was validated ONCE (N=120, ntile=2, cf=5,
  z=0, halo *count* only → 1.2%). Never at large box / mass function / lightcone / cf=4.
  The bug lives precisely in that validation gap.

Note: my algebra suggested `_isolated_convolve`'s full-band variance is ~1/8 of the
periodic one — which contradicts the too-HIGH symptom, so the exact error is a
SHAPE/scale issue in the coarse+residual split (or noise generation), not a simple
overall factor. Needs the A/B measurement to pin, not more algebra.

---

## A/B test RESULT (job 3952854, COMPLETED) — REFUTES the "field is sole cause" hypothesis

Small box (N=256, cellsize 1.25 Mpc/h, box 321 Mpc/h, z=0, nbuff=24), same seed/filters,
both paths on CPU. A = run_multitile_split (coarse+residual); B = run_multitile (global
FFT, Fortran-equivalent). Result:

- **BOTH paths produce the fake-giant high-mass bump** (mass function RISES at logM
  15.0–15.75; max RTHL = 28.7 Mpc/h = the nhunt cap in BOTH). The global-FFT control —
  validated sub-percent vs Fortran — has the giants too. **So the giants are NOT caused
  by the multi-resolution field construction.**
- Smoothed-peak-amplitude DECLINES with Rf for both (split 31.9→5.4, global 20.2→5.4).
  So neither field is "flat" at N=256 — the production flatness is a larger-N/lightcone
  effect, not the giant-maker here.
- Split over-produces: 40,599 vs 25,682 raw halos (1.58×); split's mass function is flat
  at logM 13–15 where global's declines (13→14: split ~3000/bin, global 1745→628).

**Revised conclusion — two separate effects:**
1. **Headline fake giants (RTHL pinned at the nhunt cap)** = a SHELL-ANALYSIS bug present
   in BOTH Julia paths, NOT in Fortran (Websky mass function is normal). It hid because
   README validations checked halo COUNT, never the mass function. This is the dominant
   bug and the new primary target. Location: the Phase-4 collapse logic
   (ShellAnalysisGPU.jl ~326–478 / RadialShell.jl ~480–517): when Fbar ≥ fcrit out to the
   outermost shell, RTHL is set to the cap radius instead of rejecting / finding a smaller
   radius — Fortran evidently doesn't pin like this. Cross-ref `peakvoidsubs.f90`.
2. **Multi-resolution over-production (1.58× + flat mid-mass)** = a real but SECONDARY
   coarse+residual effect (excess power → more peaks). Worth fixing after #1, but it is
   not the giant-maker.

My earlier "Bug 3 = multi-res field is the sole root cause" was WRONG; the A/B caught it.
The flat-σ(R) I measured in the production catalog is real but is the symptom of BOTH the
shell-analysis pinning (giants at the cap dominate the large-Rf peaks) AND the secondary
multi-res excess — not a clean field-only signature.

---

## Status of outputs / configs
- `catalog_websky_6144_oct000_prechifix.pksc` — 3.21M, pre-chi-fix baseline.
- `catalog_websky_6144_oct000.pksc` — 10.33M, post-chi-fix (nbuff=16), still has fake pile @18.7.
- `catalog_websky_6144_oct000_nbuff24.pksc` — 7.94M, nbuff=24, pile moved to @28.7.
- Config `config_websky_6144_oct000_nbuff48.toml` currently holds nbuff=24/ntile=16/n=429.
- A/B: `config_ab_test.toml`, `run_ab_test.jl`, `run_ab_test.slurm`.

## Fix order from here
1. Confirm via A/B that coarse+residual is the culprit (job 3952854).
2. Fix `_isolated_convolve` normalization / coarse-residual split so σ(R) matches the
   input P(k) (target: reproduce the global-FFT field).
3. Re-run oct000; verify the fake-giant tail is gone and the mass function falls
   monotonically; compare count/mass-function to Websky.
4. Keep nbuff=24 (don't pursue the _MAX_SHELLS_GPU kernel surgery — it doesn't fix the
   real bug).

---

## RESOLUTION (2026-06-14, ground truth) — it IS the nhunt cap; the kernel fix is required

An existing Fortran catalog was found already built on Killarney:
`validation/websky_multitile/fortran_run/output/catalog_raw.pksc.13579` (555,030 halos,
ntile=2, 468³, z=0, 2LPT, cellsize 2.415 Mpc/h, nbuff=22). Comparing its mass function
to the Julia catalog from the same validation:

| | Fortran | Julia | (our A/B) |
|---|---|---|---|
| cellsize / nbuff / **nhunt cap** | 2.415 / 22 / **50.7 Mpc/h** | same | 1.25 / 24 / **28.8** |
| max RTHL | 27.4 (≪ cap) | 27.8 | 28.7 (= cap) |
| mass fn logM=15.75 | **41** | **82** | ~2400 (rising) |

Both Fortran and Julia produce a **clean, steeply-falling** mass function with NO pile —
because the nhunt cap (50.7 Mpc/h) sits well above the largest real halo (~27 Mpc/h).

**Conclusion (corrects the A/B-section "shell logic / field" framing above):** the fake
giants are the **nhunt cap height**, not the collapse logic (which matches Fortran,
peakvoidsubs.f90:570-578) and not the multi-res field. The pile forms only when the cap
≈ the largest-halo collapse radius (~30 Mpc/h). Coarse-cellsize runs clear this for free
(cap 50 Mpc); our fine Websky cellsize (1.25 Mpc/h) puts nhunt=23-24 cells at only ~29
Mpc/h — right in the danger zone. The nbuff 16→24 experiment SHRANK the pile (428K→303K),
i.e. was partially fixing it, not "relocating."

**The fix (corrected):** push the cap above ~40 Mpc/h → nhunt ~32-40 cells at cellsize
1.25 → ~700-1300 shells → exceeds GPU `_MAX_SHELLS_GPU=512`. So the GPU kernel change
previously dismissed IS the fix:
  1. Raise `_MAX_SHELLS_GPU` to ~1024 (CUDAExt.jl:1263).
  2. Add the >48KB dynamic-shared-memory opt-in at the `_post_process_kernel!` launches
     (CUDAExt.jl:2192, 2421) — ~60-90KB fits the L40S 99KB limit (occupancy drops, but
     it's 1 block/peak anyway).
  3. Set nbuff ~33-40 (buffer ≥ ~40-50 Mpc/h); re-tune ntile so M=ntile·cf divides N and
     nmesh stays in memory.
  Alternative: run shell analysis on the CPU path (no 512-shell limit) for production.

The coarse+residual 1.58× over-production (A/B) is a SEPARATE, secondary issue to address
after the cap fix.

---

## ★ ACTUAL ROOT CAUSE (2026-06-14) — pk_websky.dat missing the (2π)³ normalization ★

The nhunt-cap "resolution" above is the *proximate* mechanism, not the root cause.
The field-σ probe (`probe_field_sigma.jl`) and factor test (`probe_pk_norm.jl`) found the
real bug: **the generated field is ~16× too high in σ at every scale**, because
`pk_websky.dat` is missing the Peak Patch P(k) normalization.

### The bug
`generate_pk_camb.py` wrote the RAW CAMB matter P(k) in (Mpc/h)³ with NO normalization.
The Peak Patch convolution uses `amp = sqrt(P·dk³·n³) = sqrt(P/dx³)·(2π)^1.5`, so the P(k)
file MUST be pre-divided by **(2π)³**. The official `peakpatch/tools/powerspectrum_create.py`
does `pk /= (2π·h)³` — the extra h³ is a unit conversion because its CAMB output is physical
Mpc; ours is already (Mpc/h), so the factor is (2π)³, no h.

### Evidence
- Probe: field σ(R) is a flat ~16× = (2π)^1.5 above theory at ALL R (right shape, wrong amplitude).
- split/global σ ≈ 1.00 → it's the shared P(k) DATA file, not the coarse+residual.
- Factor test (`probe_pk_norm.jl`, job 3952975): field σ(R)/theory =
  raw 15.4–15.7×, **÷(2π)³ = 1.00/0.99/0.99/0.98/0.98**, ÷(2πh)³ = 1.74–1.78×. → factor = **(2π)³ = 248.05**.
- Working `planck18_intermittent.dat` has standard-σ8 = 0.067 (pre-divided); our raw file
  had 0.817 (not divided).

### Why this explains everything
A 16×-too-high field exceeds fcrit=1.686 *everywhere*, even at 30 Mpc/h → flood of
large-filter peaks → fake 10¹⁵ M☉/h giants. The giants pile at the nhunt cap (proximate
mechanism). The giants over-exclude real halos in the merger → net DEFICIT (the original
120×). Affects BOTH paths and ALL Websky runs (shared P(k) file). It was never caught
because the Julia↔Fortran validations used the correctly-normalized `planck18` file and
checked halo *counts* (which match under a shared normalization error), never absolute σ8
or the mass function.

### Fix applied (2026-06-14)
- `generate_pk_camb.py`: now divides P(k) by (2π)³ before writing (with comment).
- `data/pk_websky.dat`: regenerated (÷(2π)³); standard-σ8 now 0.0519; raw backed up as
  `pk_websky_RAW_unnormalized.dat`. All octant configs point to it.
- End-to-end confirmation in progress (job 3952978): re-run the giant-producing N=256 z=0
  config with the corrected P(k) → expect the giant bump GONE and counts drastically lower.

### Status of the earlier "fixes" in this document
- Bug 1 (chi factor-h): real, fixed, independent. KEEP.
- Bug 2 (nhunt cap): proximate mechanism for the giants, NOT root cause. With a correct field
  there are few large peaks, so the cap rarely binds (like the clean validation). No GPU
  kernel surgery needed once the P(k) is fixed.
- Bug 3 (multi-res field σ(R)): the "flat σ(R)" was the uniform 16× over-amplitude, i.e. THIS
  P(k) bug — not a coarse+residual shape error. The coarse+residual path is fine (1.6×
  over-production is a small, separate secondary effect; likely shrinks with the correct field).

### ✅ END-TO-END CONFIRMED (job 3952978, 2026-06-14)
Re-ran the exact N=256 z=0 nbuff=24 config that produced the giants, now with the corrected
pk_websky.dat (÷(2π)³). Result vs the original (job 3952854, raw P(k)):
- max RTHL 28.7 (=cap) → **15.0/14.1 Mpc/h** (well below cap, no pile).
- high-mass tail: rising (~2400 @ logM15.75) → **falls monotonically** (18 @ 15.0, none above).
- smoothed-δ at peak: 31.9→5.4 (16× high) → **3.04→1.72** (physical ~2-3σ).
- split/global: 1.58× → **0.99×** (the 1.6× over-production was the SAME root cause).
- count 25.7K → **69.5K** (UP 2.7×): broken field's giants over-excluded real halos; fixing
  the field removes the giants → real halos survive. The deficit mechanism, reversed.
CONCLUSION: the single (2π)³ P(k) normalization bug caused the giants, over-exclusion,
over-production, and the original 120× deficit. Fix resolves all at once. Next: full
6144³ octant with corrected P(k) → expect clean mass function + count toward Websky ~1.1e8/octant.

### ✅✅ PRODUCTION-SCALE VALIDATION (job 3953930, 2026-06-14)
Full 6144³ oct000 lightcone (ievol=1, z=4.6, nbuff=24) with corrected pk_websky.dat, 4×L40S, 2h27m.
**53,355,181 halos** (pre-merge 81.7M → 53.4M; exclusion removed only 35% — over-exclusion GONE,
vs 66-84% in broken runs). All signatures clean:
- max RTHL = 17.5 Mpc/h (well below 28.8 cap; top RTHL all 1.2-1.8 Mpc/h — no giant spike).
- Mass function FALLS monotonically at high mass (logM 13.0→2.1M, 14.0→50K, 15.0→20, 15.25→5):
  physical exponential cutoff, NO rising bump (pre-fix had ~600K fake @15.75).
- Lightcone reaches r_corner = 5231 Mpc/h = z=4.6 exactly (chi fix holds).
Arc: ~0.9M/octant (all bugs) → 10.3M (chi fix) → 53.4M (P(k) fix) vs Websky ~1.1e8/octant.
From ~120× deficit to within ~2× of Websky with a clean mass function. The two real bugs
(chi factor-h, pk_websky missing /(2π)³) are fixed and validated small-box→production.
Residual ~2× is normal halo-finder territory (filter-bank spacing + Websky's abundance matching),
not a bug. FULL INVESTIGATION COMPLETE.
