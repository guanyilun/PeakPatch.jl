# Composite κ vs kap.fits: full-map reproduction at the ~10% level (2026-07-18)

Assembled κ maps from our own pipeline components and compared both Websky
double-counting schemes against the released kap.fits (job 4310996,
`compare_composite_kappa.jl`; 6 gnomonic 8° caps, flat-sky C_ℓ, cap scatter = error bar):

- **A (paper §3.1.3 literal)**: field map with halo Lagrangian spheres excluded
  (job 4304169) + plain truncated-NFW halo κ painted from the finecell AM catalog
  (197.5M halos, exact per-pixel angles, per-halo mass-exact residual deposits).
- **B (paper §3.2.4 literal)**: full-matter field map (job 4297504) + NFW halos
  compensated by a uniform Δ=3 sphere of the same total mass.

| ℓ    | A/kap ± scat | B/kap ± scat | fldE/kap | halo/kap |
|------|--------------|--------------|----------|----------|
| 165  | 1.00 ± 0.24  | 0.72 ± 0.18  | 0.40     | 0.19     |
| 323  | 1.02 ± 0.18  | 0.75 ± 0.13  | 0.43     | 0.20     |
| 632  | 0.90 ± 0.04  | 0.70 ± 0.03  | 0.38     | 0.22     |
| 1237 | 0.87 ± 0.05  | 0.73 ± 0.04  | 0.32     | 0.31     |
| 2419 | 0.84 ± 0.02  | 0.78 ± 0.02  | 0.21     | 0.45     |
| 3383 | 0.81 ± 0.01  | 0.79 ± 0.02  | 0.15     | 0.54     |

**Band mean ℓ=150–2500: A/kap = 0.920, B/kap = 0.729.**

## Conclusions

1. **Scheme A is Websky's construction** — the field component excludes lattice sites
   inside halo Lagrangian radii, and halos are painted as plain (uncompensated)
   truncated NFW. The §3.2.4 Δ=3 compensation language evidently describes the kSZ
   *electron* treatment, not the κ map assembly; applied to κ (scheme B) it removes
   too much power (~27%).
2. With scheme A, our pipeline reproduces the released kap.fits at 0.94–1.06 for
   ℓ=165–450 and within 10–20% to ℓ=3400. Component split is physical: field
   dominates at low ℓ, the 1-halo term takes over by ℓ~2500.
3. Remaining high-ℓ deficit (~15–20% at ℓ≳2000) candidates: uncorrected pixel
   windows (ours Nside 2048 vs kap's 4096), 2LPT field missing deep-nonlinear power,
   NFW profile detail differences, and our z≤4.6 vs kap's z<4.5+Gaussian tail.

## Bugs found on the way (fixed in 74bace1)

- The standalone halo kernel divided by χ instead of multiplying (Born kernel
  W ∝ (1+z)·χ·(1−χ/χ*)): halos suppressed by χ²~10⁷, composites collapsed onto the
  bare field curves. The original self-test used the same wrong W — circular. Now
  anchored independently: κ(10¹⁵ Msun/h, z=0.5, 1′) = 0.334.
- Pixel-center sampling misses the NFW cusp for sub-pixel halos: fixed with per-halo
  residual deposits (painted mass exact per halo; probe on a nearby-halo subsample
  showed 0.44 of analytic before the fix).
- Login-node note: processes are silently killed between 4–6 GB RSS; the script
  streams (maps one at a time, per-cap compact halo lists) and the production run
  uses a minimal 1-GPU-slice SLURM allocation (job 4310996 backfilled in 41 s).
