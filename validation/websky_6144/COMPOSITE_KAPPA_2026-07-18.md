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
3. **The high-ℓ deficit (0.84 at ℓ=2419, 0.81 at ℓ=3383 — i.e. 16–19%, with cap
   scatter only ±1–2% → a real systematic) is NOT yet understood.** Quantified
   so far (2026-07-18 follow-up analysis):
   - **Pixel windows explain only +1–2%**, not the deficit. Per-component Gaussian
     w_ℓ² (our field patches carry Nside-2048's 0.885/0.79 at ℓ=2419/3383; the halo
     patches are painted at 0.33′ so w²≈0.99; kap.fits carries Nside-4096's
     0.970/0.942): correcting the measured decomposition moves A/kap 0.839→0.854
     at ℓ=2419 and 0.811→0.82 at ℓ=3383. Where the window is large (field), the
     component is small (~15% of total); where the component is large (halo), the
     window is negligible.
   - **The halo side looks healthy**: our halo/kap fraction (0.45→0.54 at
     ℓ=2400–3400) matches the split measured from Websky's own catalog, and our
     NFW painting matched their halos at 1.02–1.10 over ℓ=2000–6900
     (KAPPA_PAINTED_2026-07-16.md). The deficit is concentrated in the
     **field + cross terms at small scales**.
   - Live suspects, unquantified: (a) **exclusion radius** — we excise the full
     Lagrangian sphere (AM RTHL); a smaller excision leaves more field mass near
     halos, boosting exactly the missing cross/1-halo-scale power (a real modeling
     freedom); (b) our field ran the CPU paint path with subdiv 3 vs Websky's n≤5
     splitting; (c) nearest-neighbor HEALPix→gnomonic resampling aliasing in the
     comparison itself.

## Follow-up plan (2026-07-18)

1. **Nside-4096 field rerun** (gpu_paint=true, subdiv 5 — Phase B built for this):
   removes suspects (b) and most of the window mismatch; `run_fieldmap_oct000_field4096.slurm`.
2. Composite script: print a window-corrected A/kap column from the measured
   per-band decomposition (Gaussian w_ℓ per map's actual Nside + flat halo grid).
3. If the ℓ≳2000 gap persists: **exclusion-radius scan** (R = 0.7·R_L vs R_L, one
   octant field job each) — the knob with genuine physics freedom.
4. Then: kSZ field validation (τ/kSZ maps exist from job 4304169) and the z>4.5
   Gaussian κ tail.

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
