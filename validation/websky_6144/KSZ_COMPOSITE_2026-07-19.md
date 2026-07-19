# Composite kSZ vs ksz.fits (2026-07-19)

First absolute-units test of the pipeline (κ was all dimensionless ratios; ksz.fits
is μK). Composite built exactly like the settled κ construction — scheme B =
full-matter field kSZ map (job 4319537, Nside 4096) + halo ΔT = −T_CMB·(v_r/c)·τ(θ)
with truncated-NFW electrons compensated by the Δ=3 sphere (zero net), τ constants
verbatim from FieldMap.jl eq-3.22 kernel; halo v_r reconstructed from the stored
2LPT displacements (finalize_eulerian convention). ksz.fits is z<4.5 websky only
(ksz_patchy.fits is a separate uncorrelated reionization model), so apples-to-apples
with our z<4.6 octant with no tail correction.

## Result (job 4321331): we OVERSHOOT by ~2× in power

| ℓ    | A/ref | B/ref ± scat | fldAll/ref | halo0/ref |
|------|-------|--------------|------------|-----------|
| 165  | 2.47  | 2.28 ± 0.26  | 2.16       | 0.21      |
| 452  | 1.86  | 1.58 ± 0.20  | 1.17       | 0.57      |
| 884  | 1.69  | 1.46 ± 0.12  | 0.88       | 0.70      |
| 1730 | 1.67  | 1.56 ± 0.07  | 0.71       | 0.88      |
| 3383 | 1.92  | 1.97 ± 0.09  | 0.59       | 1.27      |

Band mean 150–2500: A = 2.13, B = 1.92. Cap estimator validated against full-sky
anafast of ksz.fits (caps/full = 0.9–1.35); map std: ours ~4.9 μK vs ksz.fits 3.38 μK.
Whole-octant pseudo-C_ℓ (apodized mask, f_sky=1/8) confirms map-wide, not cap
variance: field/ksz.fits = 2.1, 2.4, 3.2, 1.8, 1.06, 0.93, 0.82 at
ℓ = 118, 163, 226, 320, 449, 626, 873.

## Why — TWO distinct causes, both understood from Stein+2020 §3.2.1/§4.4.2 + release

1. **High-ℓ (halo term): construction difference, not a bug.** Websky's kSZ halo map
   uses the **Battaglia et al. AGN-feedback GAS profiles** (generalized NFW fit: F0,
   xc, β functions of M,z; α=1, γ=−0.3; UPDATES 01-APR-2019 "Battaglia (2017) gas
   density profiles"), painted ONLY for halos with **M200c > 1e13 M☉ and r200c > 0.5′**
   — everything smaller goes into their 2LPT field. Their own release comparison
   (kszcomp.pdf) shows their halo term is tiny: ~5% of C_ℓ at ℓ~1000, ~25% at ℓ=4000.
   Our composite paints ALL 197.5M halos (≥1.2e12 M☉/h) with cuspy c=7 truncated-NFW
   electrons at full f_b — our halo term is ~10× theirs at ℓ~900. Physically theirs
   is the better electron model (feedback pushes gas out of groups); ours was the κ
   construction transplanted. To reproduce ksz.fits the halo painting must switch to
   Battaglia profiles + their mass/size cut (XGPaint has Battaglia machinery).

2. **Low-ℓ (field term): REAL anomaly in our field velocities.** Their low-ℓ kSZ is
   the "Doppler" term from broken LOS-velocity cancellation at the sharp z=4.5
   lightcone cutoff (their eq 4.3, dominant at large scales; their fig-6 dashed curve
   matches their field map). Our field map shows the same character but 2–3× stronger
   at ℓ=118–320, converging to ~1.0 by ℓ≈450 (whole-octant measurement above; z_max
   4.6 vs 4.5 explains only ~10%). Map rendered visually clean — no tile seams/blobs.

   **Suspected mechanism**: the multi-resolution coarse grid is M = ntile×coarse_factor
   = 64³ → coarse Nyquist k ≈ 0.038 h/Mpc. The k ~ 0.03–0.1 band (EXACTLY the
   Doppler-regime modes: ℓ≈150–350 at the τ-weighted χ~3–5 Gpc/h) is the coarse↔fine
   transition, carried by per-tile fine noise (tile 327 Mpc/h, HR residual). Velocity
   ∝ ψ is long-range and directly painted — errors/decoherence in that band inflate
   the Doppler power. Density projections (κ, ξ, mass function) are insensitive
   (coherent displacement errors don't change local density), which is why every
   prior validation passed. Halo v-statistics validated only at r<50 Mpc/h (v12) or
   1-pt (σ_vr) — also blind to this band.

   **Test in flight**: job 4321358 reruns the full-matter octant with
   coarse_factor=16 (256³, Nyquist 0.154) at Nside 2048; same seed (coarse noise is
   block-averaged from the same fine noise, so same realization). Analysis:
   `compare_cf16_ksz.jl` — if ksz_cf16/ksz.fits → ~1 at ℓ<350 while kappa_cf16/cf4
   stays ~1 (density control), mechanism confirmed and the production fix is a
   larger coarse grid.

## Notes

- Self-tests: ∫τdΩ exact (plain), 0 (compensated); anchor τ(1e15 M☉/h, z=0.5, 1′)
  = 5.2e-3 (real-cluster scale), kSZ@300 km/s = 14.1 μK.
- In kszcomp.pdf their FIELD sits ABOVE their TOTAL at ℓ≲300: the compensated halo
  term anti-correlates with the field at low ℓ.
- If (2) confirms, the coarse_factor choice affects any velocity-sensitive Tier-B
  product (kSZ, ISW/moving-lens); κ/density products are unaffected.
