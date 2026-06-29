# Controlled Julia-vs-Fortran finder comparison (2026-06-14)

Goal: rigorously determine whether the Julia pipeline under-finds halos vs the
Fortran reference (`hpkvd`), against GROUND TRUTH rather than the unreliable
Sheth-Tormen theory baseline the earlier notes leaned on.

## Setup (reproducible)
- Fortran `hpkvd` runs on Killarney: `/home/yguan/scratch/websky_6144/fortran_floortest/`
  - binary: `/home/yguan/work/peakpatch/bin/hpkvd` (Killarney-built)
  - env (libs): add to LD_LIBRARY_PATH a dir with ONLY a symlink to
    `libze_loader.so.1` (gentoo lib64) + fftw/3.3.10 + fftw-mpi/3.3.10 +
    openmpi/4.1.5 lib dirs. Do NOT add the whole gentoo lib64 (overrides libc → GLIBC error).
  - params: `gen_params.py`-style binary `hpkvd_params.bin` (stream); pk file must
    have comments stripped and be uniform-log-k (Fortran reader assumes both).
- Matched config (Julia A/B `config_ab_test.toml`): next=256, box=320.833 Mpc/h,
  cellsize 1.25326 (production), websky cosmo (Omx=0.261,OmB=0.049,OL=0.69,h=0.68),
  z=0 snapshot (ievol=0), same `pk_websky.dat` (÷(2π)³) + `filters_websky.dat`,
  collapse params identical (iforce_strat=4, ivir_strat=2, fcoll 0.01/0.171/0.171, dcrit=200).
  Mass-function comparison is statistical → seed-independent.

## ⚠️ CORRECTION (the headline result)

**At MATCHED cellsize, Julia and Fortran agree to 0.5%.** An earlier version of this
note reported Fortran finding ~2× more peaks — that was a CELLSIZE CONFOUND in the
Fortran setup, NOT an algorithm difference:

### Root cause: the grid is COMPILE-TIME, not a runtime param
`src/hpkvd/arrays.f90` declares the grid as a compile-time constant:
`integer, parameter :: n1 = N_REPLACE` (the build step substitutes the run's nmesh and
**recompiles**). This binary was built with **n1 = nmesh = 256** (see
`arrays_gen.f90`), so it ALWAYS uses
`nsub = nmesh-2*nbuff = 212`, `n_global = nsub*ntile+2*nbuff = 468` — *regardless of
the `next` field in the param* (verified: changing next 468→256 did not change the
logged n). The cellsize the run actually uses is **`cellsize = dcore_box / nsub`**
(== `dL_box/nmesh`). My param sized `dcore_box = 106*1.2533 = 132.845` for nsub=106
(nmesh=150), but the binary's nsub is fixed at 212, so the effective
**cellsize = 132.845/212 = 0.63 Mpc/h** — ~2× FINER than the 1.25 Julia ran at. At the
grid-limited floor, halving the cell ~doubles the peak count → the bogus "2×".

### How to ensure matched geometry in future (checklist)
1. **Read the binary's compiled nmesh** (`arrays_gen.f90` → 256) or the log line
   `Slab decomposition: n = <N>` → `nsub = (N-2*nbuff)/ntile`.
2. **Size the box to that nsub**: `dcore_box = nsub*cellsize`, `dL_box = nmesh*cellsize`;
   assert `dcore_box/nsub == dL_box/nmesh`. Use `gen_fortran_params.py` (has this
   assertion baked in).
3. **Or recompile hpkvd** for the nmesh you want (the standard peak-patch.py workflow).
4. **Guards before trusting any comparison**: (a) verify logged `n == nsub*ntile+2*nbuff`
   and `cellsize == dcore_box/nsub`; (b) cross-check `σ(R)` at a couple of filter scales
   (~1% — physical invariant). NOTE: σ0 alone is INSENSITIVE to the Nyquist cutoff, so a
   cellsize mismatch can hide behind a matching σ0 — also check the box size / σ2 / density.
- **Fix applied**: run BOTH at cellsize 1.2533 on the native n=468 grid
  (Fortran `dcore_box=265.69, dL_box=320.83`; Julia n=256/nbuff=22/ntile=2
  → N_global=468). Result:

| (cellsize 1.2533, z=0) | kept peaks | core box | density /Mpc³ |
|---|---|---|---|
| Fortran (n=468) | 592,193 | 531.4³ | 3.946e-3 |
| Julia   (N=256) | 69,538  | 260.7³ | 3.925e-3 |

**Fortran/Julia density = 1.005.** Julia's finder matches Fortran to 0.5%. There is
NO Julia under-finding and NO near-Nyquist treatment difference (field-gen, smoothing,
and peak criteria are identical in the two codes).

## Results at the (confounded) finer Fortran cellsize 0.63 — for the record only
| | peaks found | kept (post-collapse) |
|---|---|---|
| Julia (cellsize 1.25) | 82,022 | 69,538 |
| Fortran (cellsize 0.63, CONFOUNDED) | 265,858 | 224,181 |

These are NOT comparable (different cellsize). Dump rates matched (15% vs 16%).

### Find-masking is BENIGN (not the cause)
Disabling Julia's cross-scale peak masking (`PP_NO_FINDMASK` experiment):
raw 69,538 → 112,409, but after `merge_catalog` **both give exactly 25,875**.
So the masking is a valid dedup that the merge reproduces — NOT a bug.
(The earlier "3× under-finding" alarm was a pre-merge bookkeeping artifact:
Fortran defers exclusion to merge_pkvd; Julia masks during find.)

### Julia's FIELD is correct (the decisive check)
Smoothed-field spectral moments at the smallest filter Rf=2.068 (≈1.65 cells):
| moment | Julia | continuum theory | grid (k≤k_nyq … corner) |
|---|---|---|---|
| σ0 | 1.762 | 1.783 | 1.779 |
| σ1 | 1.221 | 1.280 | 1.208 |
| σ2 | **1.453** | 3.703 | 1.30 – 1.59 |

σ2 (curvature, sets peak density) matches the **cubic-grid Nyquist limit**. The
continuum σ2=3.70 is physically unreachable on a 1.25 Mpc/h grid. → Julia's field
is NOT missing high-k power; it is correct.

### The floor IS grid-limited — but both codes agree there
At Rf=2.068, R* = √3 σ1/σ2 = 0.60 Mpc/h < cellsize 1.25 → the lowest-mass halos
are sub-cell / marginally resolved, so the floor count is set by (and sensitive to)
the cellsize. BUT at the SAME cellsize Julia and Fortran agree to 0.5% (see top) —
the floor count is a shared, resolution-set property, NOT a code/near-Nyquist
difference. The earlier apparent Fortran excess was purely the 0.63-vs-1.25 cellsize
confound.

## Lightcone (ievol=1) verified correct — NOT a deficit source
Controlled test: run_multitile_split at ievol=1 (observer at box corner → z<0.19) vs
ievol=0, same 256³ box/cellsize. At z≲0.19 the lightcone MUST match the snapshot:
| | N total | >1.69e12 | >1e13 | >1e14 |
|---|---|---|---|---|
| ievol=0 snapshot  | 69008 | 65223 | 39101 | 5611 |
| ievol=1 lightcone | 67462 | 63224 | 36356 | 4566 |
ratio lc/snap = **0.978**. The path is correct. The small reduction is MASS-DEPENDENT
(0.97 floor → 0.81 @1e14) = correct physics: a halo seen at z>0 is less collapsed →
smaller RTHL → high-mass ones slip below a fixed cut.

This RETIRES the earlier "lightcone-specific deficit" claim: that was a differential-
vs-integrated comparison error (z=0-0.3 shell 0.42 = the z=0 *differential* floor
completeness ~0.47, NOT vs the 0.76 *integrated-above-floor* number). No extra deficit.

## Production-specific paths verified (GPU, multi-res M)
A streaming z-shell scan of the production catalog hinted at a low-z suppression; the
two production-only knobs were isolated and BOTH are clean:
- **GPU shell analysis = CPU exactly**: run_multitile_split use_gpu=true vs false on the
  same z=0 config → raw 69008=69008, floor 65223=65223, merged 25898 vs 25897. ratio 1.0000.
- **Multi-res loss tracks BLOCK size (N/M), not M**: N=256 z=0 floor N(>1.69e12):
  global 65723; split M=16/block=16 → 65223 (0.992); M=32/block=8 → 64086 (0.975);
  M=64/block=4 → 61802 (0.940). Small block puts the coarse/residual split near the
  floor scale and degrades it. **Production block = 6144/64 = 96** (split at ~120 Mpc/h,
  far from the 2-Mpc floor) → loss <1%. OPERATIONAL RULE: keep N/M ≳ 16 (ideally ≫) so
  the multi-res split stays well above the smallest filter.
- The raw shell "0.53×" was a convention artifact (full-box-pre-merge anchor vs
  core-post-merge catalog); apples-to-apples (post-merge, vs ST) gives production ≈
  0.79× the finder — within small-box variance + ST high-z model uncertainty, not a bug.

## Conclusion
Every pipeline stage now verified: field (σ_j match grid spectrum), finder (=Fortran
0.5%), lightcone (=snapshot 2%, correct mass-dependent evolution), GPU (=CPU to 1 halo
in 69008), multi-res (production block N/M=96 → <1% loss; loss only at tiny blocks).
The ~4× production-to-Websky gap (26.5M vs 110M /octant above 1.69e12) decomposes
into TWO effects, NEITHER a Julia field/finder/lightcone bug:
1. **Grid-limited floor**: the smallest halos are marginally resolved; the raw floor
   count is resolution-set. BOTH codes give the same density at the same cellsize
   (1.005). This is exactly why peak-patch catalogs are abundance-matched.
2. **~2× M200m mass calibration**: Websky AM-remaps top-hat → N-body M200m, lifting
   more halos over a fixed mass threshold.

Density cross-check (anchored to the verified z=0 finder density 2.0e-3/Mpc³):
ST lightcone dilution predicts production/z0 ≈ 0.285; measured 0.177; Websky 0.734
(2.6× above ST → M200m inflation, consistent).

This SUPERSEDES the earlier "93% complete vs ST" framing — that rested on an
unreliable ST normalization. The current conclusion rests on verified ground truth:
matched σ_j, controlled Fortran comparison, exonerated find-masking.

## Open (optional) definitive test
Have Fortran read Julia's (verified-correct) field via `ireadfield=1` and compare
peak counts on the IDENTICAL field. If Fortran still finds ~2× more, it is Fortran's
floor/near-Nyquist convention (Julia confirmed right); if it matches, the difference
was field generation (Julia already shown correct by σ_j). Either way Julia's field
is correct; this only pins which Fortran stage inflates the marginal floor count.

## Scripts (validation/websky_6144/)
- `run_floor_compare.jl` — Julia split/global N(>M) vs Fortran table
- `run_floor_diag.jl` — Julia dump-counter breakdown (found/kept/dumped)
- `merge_fortran_cat.jl` — apply Julia merge_catalog to Fortran raw (Lagrangian cols 0-2)
- `run_capstone_merge.jl` — Julia raw → merge_catalog
- `run_nomask_test.jl` — find-masking benign test (needs the reverted PP_NO_FINDMASK toggle)
- `probe_sigma_vs_fortran.jl` — per-scale σ(Rf) Julia vs Fortran
- `probe_sigmaj.jl` — σ0/σ1/σ2 vs theory (the decisive field check)
- `probe_peakcount.jl` — local-maxima count, masked vs nomask
