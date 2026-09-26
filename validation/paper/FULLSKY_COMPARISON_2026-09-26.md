# Full-sky frozen-campaign comparison: first pass (2026-09-26)

Plan: `docs/paper_comparison_plan_2026-09.md`. Error model A3 was fixed before any numbers.

## What was built today

**AMv2 catalogs.** All 8 octants re-abundance-matched with the top-halo fix (jobs
5692696–703). Each file is the same size as its raw catalog. N(>1.7e12) = 5.076–5.078e7 and
N(>1e14) = 6.251–6.255e4.

**Halo maps** (`production/paint_octant.jl` + `halo_profiles.jl`; HEALPix Nside 4096).
- The physics is ported verbatim from the validated cap painters: κ construction B,
  Battaglia-2012 y, and Battaglia gas kSZ with W (plain) and Wc (compensated).
- Every halo is exact per halo: a residual against its analytic integral goes into the
  centre pixel.
- Threads own disjoint ring bands, so there are no pixel races.
- Self-tests: disc sets equal `queryDiscRing` at Nside 64, 1024 and 4096; totals exact to
  1e-9; band-parallel output equals serial.
- Cost: 9.5–10 min and 40 GB per octant. Output is in `/project/.../websky/halomaps_prod/`.

**CIB** (`production/paint_cib.jl`): XGPaint `CIB_Planck2013` at 6 frequencies with the
Websky completeness cut, 48 min and 54 GB per octant. Pixel accumulation is serial because
XGPaint's HEALPix `paint!` has a threaded `+=` race on shared pixels. That finding is noted
here only; nothing was sent upstream. oct000 is done and 7 octants are queued
(jobs 5695072–78).

**Full-sky assembly** (`assemble_fullsky.jl`) is in `/project/.../websky/fullsky_prod/`.

**Spectra tooling** (`spectra.jl`, `compare_auto.jl`): Healpix.jl `map2alm` at Nside 4096,
lmax 8192, takes 2.5 s on 32 threads.

## First full-sky results (`results/auto_*.txt`)

| product | result |
|---|---|
| tSZ (compared at Nside 2048) | mean y ratio **1.009**; C_ℓ 0.76–0.80 at ℓ=300–1000, 0.94 at ℓ≈3500 |
| kSZ total (field + Wc) | 1.1–1.3 at ℓ=50–150; 1.3–1.42 at ℓ=200–1000 (the documented mid-ℓ construction difference); 1.03 at ℓ≈3500; 0.89 at ℓ≈8000 |
| ISW | **fails**: 2× at ℓ≈40 rising to 745× at ℓ≈850; monopole 3.7 μK vs 0.27 |

Per-octant tSZ, C(ℓ=400–800): ours 6.5–7.5e-18 in every octant against Websky's
7.7–9.0e-18. So the deficit is systematic, not realization scatter. The oct000-only
masked check gave 0.93–1.0.

## ⚠️ Geometry defect in the frozen campaign: an unsimulated slab on the octant seams

**Cause.** The production tiling has ntile=16 × nsub=382 = **6112 core cells**. The periodic
box has N=6144 cells, because `src/MultiResolution.jl:554` sets `N = nsub*ntile + 2nbuff`.
Buffers wrap periodically (`mod1(…, N)`), but the 2·nbuff = 32 cells between the last and
first cores are **never in any tile core**. The cores span ±2604.4 Mpc/h while the
observer sits at the box corner, ±2618.1 Mpc/h. Each octant is therefore missing a
**13.6 Mpc/h slab against each of its 3 bounding planes**. Over the full sky that is a
27 Mpc/h slab through the observer along each plane, about 1.6% of the volume.

**Evidence.**
- *Catalog* (oct000, first 4M halos): counts per 5 Mpc/h of distance to the nearest
  plane are 232, 4.6k, 32k, then 70–85k per bin beyond 15 Mpc/h.
- *Maps*, mean signal against angular distance to the nearest plane, relative to the
  value beyond 2°:

  | distance | κ | y |
  |---|---|---|
  | <0.1° | ≈0 | ≈0.01 |
  | 0.2–0.5° | 0.83 | 0.48 |
  | 0.5–1° | — | 0.80 |

  The angular gap is 13.6/χ, so it covers several degrees at z≲0.1.
- *Exactly-zero pixels*: 190k in the κ field and 1.5M in ISW, all within 1° of a plane.
- **Websky's released maps are flat across the planes** (tsz, ksz and cib545 are all
  within noise), so the reference has no gap. Its Fortran tiles cover the full periodic box.

**Consequences.**
- **ISW:** the seams create the catastrophic ℓ>30 excess, because the octant-mean potential
  steps across holes.
- **κ:** holes where the +1.3 mean drops to 0 would dominate C_ℓ at high ℓ, so κ is not
  comparable until this is fixed.
- **tSZ:** a large part of the ~20% C_ℓ deficit at ℓ<1000, since the missing nearby volume
  sits exactly where tSZ comes from.
- **kSZ:** a small effect.
- **Catalog statistics:** octant-volume counts are low by about 1.6%.
- The earlier Tier-A/Tier-B validation used 8° caps or apodized single octants away from
  the planes, which is why this was never seen.

**Fix** (applied 2026-09-26, opt-in `[grid] periodic_cores = true`; see the Decision
section): let the cores tile the whole periodic box, `N = nsub*ntile`, with buffers wrapped
periodically. The index machinery already supported this. For the Websky grid that means
nsub=384, n=416 (nbuff 16), N=6144, and cf32 gives M = ntile·cf = 512 and block 12. The
observer stays at the box corner (±2618), which is now also the core corner.

## Other findings from today's parallel work (details in their notes)

**Performance** (`validation/performance/PERFORMANCE_2026-09-26.md`).
- Two-thirds of each ~7 h catalog job is the serial merge, which uses a fixed 256³ grid of
  20.9 Mpc/h cells.
- With 5 Mpc/h cells the merge is 26× faster with identical survivors, bringing an
  octant to about 2.6 h.
- Scaling from 1 to 4 L40S is 3.15×; one H100 is 1.69× one L40S.

**Tiling** (`validation/tiling/TILING_INVARIANCE_2026-09-26.md`).
- Noise and the coarse field are bit-identical across tilings; the fine fields are not
  (ψ differs by 7–13% of σ).
- Against an exact global FFT there is a ±10% P(k) feature near the coarse Nyquist.
- Merge tie-break: about 22% of raw peaks tie exactly, because of the discrete shell
  radii in `RadialShell.jl`. As a result, merged catalogs are not reproducible run to run
  or across GPU counts. A total-order sort fixes this.

**Theory** (`validation/paper_theory/THEORY_ANCHORS_2026-09-26.md`).
- Field κ is 1.02–1.09 × linear Limber at ℓ=60–2700.
- Field kSZ is 1.14–1.19 × linear theory at ℓ=73–245 in every octant. This is
  unexplained, and the old z_max explanation is refuted (+5–7% only). The deciding test is
  the 3D velocity P(k) against linear theory.

## Decision (user, 2026-09-26): rerun as production-v2

- **Seam:** fixed with `[grid] periodic_cores = true`, using `grid_layout()` in the GPU
  multires and FieldMap paths. The CPU MultiTile and MPI paths refuse the option.
  - Tests: `test_fieldmap.jl` checks that a box-corner observer sees every cell, and that
    the legacy layout shows the deficit and a warning.
  - Pre-flight on a GPU octant: no depletion at the planes.
- **Splice bump:** diagnosed and fixed (`coarse_compensation`, user-approved;
  `../tiling/SPLICE_COMPENSATION_2026-09-26.md`).
- **Bundled:** merge tie-break and `auto_hash_nc`; φ_ij Nyquist zeroing on GPU and CPU.
- **Test status:** `Pkg.test` and the GPU tests pass (job 5696014). That job ran before
  compensation was added; the default path is unchanged.
- **Cleanup:** 258 GB of superseded catalogs deleted (list in
  `/home/yguan/scratch/websky_6144/cleanup_2026-09-26.txt`).
- **Campaign:** `../websky_6144/production_v2/README.md`, submitted 2026-09-26.
- **Old-campaign CIB preview** (seam-affected, for tooling only):
  - mean intensity 1.24 × Websky at every ν, consistent with the Tier-B 1.254
  - C_ℓ about 6× at ℓ≈100, falling to about 1.2 at ℓ≈3000
  - Poisson tail 0.71 at 545 GHz and 0.49 at 857 GHz

  Interpret only after v2 is in, and apply Websky's 400 mJy flux cut to both maps for
  the Poisson comparison.
