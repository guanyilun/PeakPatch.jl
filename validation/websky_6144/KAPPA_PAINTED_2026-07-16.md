# Painted κ vs Websky kap.fits (2026-07-16) — first paper-figure-level test

Uses the new `NFWKappaProfile` in the XGPaint.jl fork (branch `lensing-kappa`,
github.com/guanyilun/XGPaint.jl): Websky-recipe NFW c=7 (r200m) + r⁻² tail × Born CMB
lensing kernel. Script: `compare_kappa_painted.jl` (run with `--project=~/work/XGPaint.jl`).

Three maps on the matched footprint (7.08° square inscribed in the shared 5° cap,
0.25′ CAR pixels in a cap-centered rotated frame → uniform pixel area):
1. **ours** — oct000 finecell_AM halos (404,048 painted, M>2e12 Msun/h, Eulerian)
2. **wsky_halo** — Websky 10×10-patch halos (395,090), identical painting
3. **kap.fits** — Websky's released convergence map sampled on the same footprint in the
   Websky frame — the **same sky** as (2), so (3)−(2) is exactly their field component
   (+ sub-cut halos + z>4.5 Gaussian tail)

## Results

| ℓ | C_ours | C_wsky_halo | C_kap.fits | ours/whalo | whalo/kap |
|---|---|---|---|---|---|
| 322 | 7.68e-9 | 6.53e-9 | 5.17e-8 | 1.18 | 0.13 |
| 807 | 2.96e-9 | 2.98e-9 | 1.44e-8 | 1.00 | 0.21 |
| 1487 | 1.61e-9 | 1.64e-9 | 5.35e-9 | 0.99 | 0.31 |
| 2742 | 9.00e-10 | 8.64e-10 | 1.87e-9 | 1.04 | 0.46 |
| 6865 | 2.89e-10 | 2.79e-10 | 3.68e-10 | 1.04 | 0.76 |

- **Catalog test: mean ours/wsky_halo = 0.986 over ℓ=300–2000**, and 1.02–1.10 per bin at
  ℓ=2000–6900. With realistic NFW profiles, our catalog's painted κ power matches the
  Websky catalog's to ~1–2% in the mean — much tighter than the point-mass test, because
  the profile-smoothed 1-halo term is less dominated by the few nearest halos. Map rms:
  ours 0.1035 vs wsky_halo 0.1021 (+1.4%).
- **Halo fraction physics (same-sky, catalog-independent): whalo/kap ≈ 0.13–0.2 at
  ℓ≈300–800 rising to 0.76 at ℓ≈6900.** The low-ℓ level is exactly the expected 2-halo
  suppression (f_coll·b_eff)² ≈ (0.25·1.5)² ≈ 0.14 for halos above 2e12; the rise is the
  1-halo term taking over. The high-ℓ ratio is additionally suppressed by kap.fits'
  Nside=4096 pixel window (~0.8 in C_ℓ at ℓ~6900, not corrected here) while our painted
  maps have a 0.25′ grid.

## Implications

- The **halo side of κ painting is done and validated**: profile + kernel + catalog give
  the same κ power as Websky's own halos.
- Reproducing the actual kap.fits (and the paper's lensed-CMB Fig 13) now hinges entirely
  on the **field component** — matter outside halos from the LPT displacement field —
  which is ~80% of C_ℓ^κκ at ℓ≲1000. That is the planned PeakPatch.jl feature
  (field-particle/mass-shell lightcone output), also needed for the kSZ halo/field split.
