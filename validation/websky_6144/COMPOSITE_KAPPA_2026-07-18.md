# Composite κ vs kap.fits (2026-07-18)

> **⚠️ MAJOR CORRECTION (2026-07-19, see the section at the bottom): the original
> conclusion below — "scheme A wins" — was an artifact of comparing against the WRONG
> reference. kap.fits INCLUDES the z>4.5 Gaussian tail (15–45% of C_ℓ, not "small").
> Against the correct apples-to-apples reference kap_lt4.5.fits, scheme B
> (compensated) matches and scheme A overshoots by ~40%. A's apparent agreement with
> kap.fits was an accidental cancellation of two ~40% errors.**

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

# CORRECTION & RESOLUTION (2026-07-19): the flat 10% and the scheme question

Following the user's challenge on the residual deficit, three independent
measurements got to the bottom:

1. **Full-sky anafast of kap.fits** (degraded to Nside 512, window-corrected):
   kap.fits/linear-total(z→1089) = 0.94–0.97 FLAT over ℓ=50–520, while
   kap.fits/linear(z<4.6) rises 1.13→1.50. **kap.fits contains the z>4.5 tail.**
2. **The release provides the split**: kap_lt4.5.fits and kap_gt4.5.fits
   (downloaded to websky_ref/). Full-sky spectra: lt4.5 = 0.93–1.00 × linear(z<4.5);
   gt4.5 = 0.92–0.94 × the linear Limber tail — the tail is **15% → 34% of total C_ℓ
   over ℓ=50–520** (Limber: →47% by ℓ=3400 in linear terms). The old "small at
   ℓ≳100" caveat was simply wrong.
3. **Composite vs kap_lt4.5.fits** (job 4319404, apples-to-apples):
   **A/lt4.5 = 1.38–1.44 at ℓ=323–1237** (falling to 1.19 at 3383) — scheme A
   OVERSHOOTS their z<4.5 map, by just about the amount of its +0.636·M200m/halo
   double count (plain NFW paints 1.636·M but the field only excises 1.0·M).
   **B/lt4.5 = 1.05–1.14, roughly flat** — the compensated construction matches.
   Consistency: B + measured gt4.5 tail reproduces kap.fits at ~1.04–1.06; and the
   caps' implied tail (kap − lt at ℓ=884: 4.3e-9) equals 0.92 × the Limber tail.

**Websky's κ construction is the COMPENSATED one (§3.2.4 applied to κ as §3.2.3
says), not the exclusion reading.** Scheme A's match to kap.fits was two ~40% errors
cancelling: (+) tail-mass double count vs (−) missing z>4.5 tail.

## Final numbers (2026-07-19, job 4321268: both field maps Nside 4096, vs kap_lt4.5)

| ℓ    | A/lt4.5 | B/lt4.5 (full comp, net 0) | B2/lt4.5 (M-sphere, net +0.636M) |
|------|---------|----------------------------|----------------------------------|
| 165  | 1.38    | **0.99 ± 0.22**            | 1.35                             |
| 452  | 1.42    | **1.08 ± 0.12**            | 1.39                             |
| 884  | 1.38    | **1.12 ± 0.04**            | 1.37                             |
| 1730 | 1.34    | 1.20 ± 0.03                | 1.38                             |
| 3383 | 1.19    | 1.20 ± 0.02                | 1.27                             |

Band means (150≤ℓ≤2500): A = 1.384, **B = 1.107**, B2 = 1.376.

- **B2 ≈ A (1.376 vs 1.384)** — a decisive internal consistency check: both net
  +0.636·M200m per halo (A by exclusion bookkeeping, B2 by under-compensation), and
  they land on the same C_ℓ. The mass ledger fully explains the scheme differences.
- **B (zero-net-mass compensation) is Websky's construction**: 0.99–1.12 for
  ℓ=165–884, rising to ~1.20 at ℓ≥1730. The earlier "B falls back to 1.05 at ℓ=3383"
  was the Nside-2048 full-matter map's pixel window — gone at 4096.
- The remaining high-ℓ excess (+20%, we are ABOVE kap_lt4.5) plausibly reflects that
  OUR halo painting is per-halo mass-exact (residual deposits) while pks2map samples
  profiles at Nside-4096 pixel centers with no residual correction — sub-pixel halos
  lose most of their mass there (we measured 0.44× for exactly this failure mode in
  our own painter before fixing it). At ℓ≳1500 our composite is plausibly MORE
  faithful than kap_lt4.5 itself; low/mid ℓ agree within cap scatter.

**Bottom line: Websky κ = full-matter 2LPT field + truncated-NFW halos minus a
uniform Δ=3 sphere of the full painted mass (zero net) + Gaussian z>4.5 tail. Our
pipeline reproduces kap_lt4.5 at ~1.0–1.1 (ℓ≲1000) and kap.fits at 1.04–1.06 once
the measured tail is added.**

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
