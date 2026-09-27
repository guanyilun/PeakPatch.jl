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

   **✅ CONFIRMED (job 4321559, coarse_factor=16 = 256³ coarse grid, Nyquist 0.154,
   same underlying noise → same realization; 65 min on ONE L40S at Nside 2048).**
   Whole-octant pseudo-C_ℓ (`compare_cf16_ksz.jl`):

   | ℓ    | ksz cf4/ref | ksz cf16/ref | κ cf16/cf4 |
   |------|-------------|--------------|------------|
   | 118  | 2.13        | **1.44**     | 1.10       |
   | 163  | 2.35        | **1.40**     | 1.19       |
   | 226  | 3.16        | **1.29**     | 1.16       |
   | 320  | 1.78        | **1.16**     | 1.11       |
   | 449  | 1.06        | 1.01         | 1.03       |
   | ≥626 | 0.93→0.70   | 0.94→0.77    | 0.95-1.04  |

   The low-ℓ excess collapses (2.1-3.2× → 1.2-1.4×) with the density control moving
   only 10-19% — and in the direction of MORE power in the same transition band,
   i.e. the 64³ coarse grid was also mildly suppressing density there (κ field/CAMB
   was 0.94-0.95; cf16 lands ~1.03-1.13). Residual 1.2-1.4× at ℓ≲250: coarse
   Nyquist 0.154 still marginal for the band (tricubic interp degrades near
   Nyquist), z_max 4.6 vs 4.5 (+10% in u0²), octant realization variance.
   **✅ CONVERGED (job 4321879, coarse_factor=32 = 512³, Nyquist 0.31, 1h12m on
   one L40S): cf32/ref = 1.22/1.20/1.13/1.05/0.97 at ℓ=118/163/226/320/449 (cf16
   was 1.44/1.40/1.29/1.16/1.01), κ cf32/cf16 = 0.97-1.04 at low ℓ.** The residual
   ~10-20% at ℓ≲230 matches z_max 4.6-vs-4.5 (+10% in u0²) + single-octant
   realization variance — field-velocity story CLOSED at coarse_factor≈32.**

   Production implication: coarse_factor=4 is fine for catalogs/density statistics
   but NOT for velocity-sensitive painted products; field-map (and eventually
   catalog-velocity) runs should use coarse_factor ≥16 — cost is negligible
   (65 min single-GPU octant at 2048 including the bigger coarse FFTs).

## Notes

- Self-tests: ∫τdΩ exact (plain), 0 (compensated); anchor τ(1e15 M☉/h, z=0.5, 1′)
  = 5.2e-3 (real-cluster scale), kSZ@300 km/s = 14.1 μK.
- In kszcomp.pdf their FIELD sits ABOVE their TOTAL at ℓ≲300: the compensated halo
  term anti-correlates with the field at low ℓ.
- If (2) confirms, the coarse_factor choice affects any velocity-sensitive Tier-B
  product (kSZ, ISW/moving-lens); κ/density products are unaffected.

## cf32-catalog composite (2026-07-26, job 4418280): velocity attribution REFUTED

The coherent-velocity (cf32) catalog did NOT shrink the mid-ℓ excess — halo/ref
rose slightly (0.69/0.77/0.77 at ℓ=323/452/632 vs 0.63/0.70/0.71); W = 1.83,
Wc = 1.38 band mean. Catalog-velocity decoherence was NOT the mid-ℓ driver.

**Decisive patch test** (websky's OWN 10×10 halos — their pos+masses+velocities —
through our verbatim painter): their-catalog halo kSZ C_ℓ = 0.26-0.43× OURS at
ℓ=632-3383. So the excess halo term splits into:
1. **Catalog abundance (~2×in count → ~2.5-3.8× in C_ℓ)**: our selected
   (mh>1e13 + 0.5′) surface density is 2.04× websky's (1090 vs 535 per deg²);
   the dN/dz diag shows the excess concentrated at z≳2 (2.3-3.8×) — our per-shell
   Tinker AM vs websky's real high-z abundance. NEXT: measure ours-vs-websky
   N(M>1e13, z) directly; consider whether websky's high-z M>1e13 halos are
   incomplete (their completeness file only covers the low-mass floor) or whether
   our AM overproduces at high z.
2. **Compensation (~2-4×)**: even their halos give halo/total ≈ 23% at ℓ~900 vs
   ~5-10% in kszcomp.pdf — their released halo map is evidently COMPENSATED
   (§3.2.1 Δ=3 mean-density sphere of the halo mass), suppressing mid-ℓ halo
   power; our W is uncompensated (verbatim per-halo painting has no compensation
   in the Fortran path we read — production must apply it elsewhere, e.g. via the
   table like κ's gas+DM−1). Our Wc (full-τ sphere, R=4rvir) is the right class:
   Wc = 1.38 vs W = 1.83. NEXT: implement their exact compensation (Δ=3 MEAN
   density sphere, mass = halo mass → Rc = (3M/(4π·3·ρ̄m))^{1/3} ≈ 4-6·rvir).
   Also verified: Fortran vrad = v·r̂/c, same as ours.

Bottom line: painter faithful; residual = catalog high-z abundance (real, ours-vs-
theirs difference) × compensation scheme (their production detail). ksz.fits total
remains matched at 1.09-1.26 for ℓ≥2419 (Wc) and the field at ℓ≤450; mid-ℓ is the
convolution of these two identified factors.

## Careful re-check (2026-07-26, job 4419407 + audits): corrections and final state

**CORRECTION to the patch test**: the 10×10 patch is NOT centered on its mean halo
direction (corner near origin; spans [0,10°]×[−7.5°,2.6°]) — the centered square
caught only 16.6% of halos instead of 32.1% → the 0.26-0.43 ratio was a ×1.93
geometry undercount. RECENTERED: websky's own halos through our painter =
**1.16-1.75× OUR halo term** (single-patch realization variance; consistent).

**Abundance question CLOSED**: cf32_AM catalog matches websky N(mh>1e13, z) at
**1.00 total, 0.93-1.05 per Δz=0.25 bin** (with and without the 0.5′ cut; the
"2×" CIB-era count excess was the old finecell_AM's mass calibration + the same
patch-geometry error). Catalog, painter amplitude, and velocities are now ALL
cross-validated between the two pipelines.

**We variant (Websky-exact literal: Δ=3-mean sphere, gas mass fb·mh)**:
W/ref = 1.84, Wc/ref = 1.38, **We/ref = 1.58** band mean — between W and Wc, does
NOT close mid-ℓ. Audit findings on the literal reading:
- gasr(mh,z) = painted-gas/(fb·mh) < 1 at high z (0.56 at 1e13, z=2): the "same
  mass" sphere OVER-compensates high-z halos (net negative) — either their real
  behavior (would help explain their tiny halo term) or the prescription differs.
- Their maptable rt-axis caps transverse radius at 4 Mpc COMOVING (not 4·rvir):
  drops 10-29% of cluster τ (M≥2e14) — a production-side suppression we don't have.
- Sphere-mass ambiguity (SIS mh vs M_RTH): ≤20% in comp amount.

**Honest final state of kSZ mid-ℓ**: with catalog+painter+velocities verified
equivalent, a faithful Battaglia painting of EITHER catalog produces a halo term
~5-10× larger at ℓ~500-1000 than the released map's own decomposition
(kszcomp.pdf), and no uniform-sphere compensation variant (We/Wc/W) reaches
their ~5-10% halo fraction there. The remaining gap is a websky-production detail
not recoverable from the paper, the released Fortran, or websky_model (which is
only a number-density projector): plausibly a much higher effective production
mmin, a stronger compensation than documented, and/or the fixed-4-Mpc + table/
pixel truncations stacking. Composite scorecard vs ksz.fits (Wc, physically
consistent zero-net): 1.09-1.17 at ℓ≥2400, 1.31-1.54 at ℓ=450-1750, 1.32-1.41
at ℓ≤330 (low-ℓ dominated by known zmax-4.6 + octant-realization effects).
For OUR production maps this is a documented model choice, not an open bug:
painter, catalog, and field are each independently validated.

## THE BOTTOM (2026-07-27): their code run, their bug found, the excess identified

We compiled and ran WEBSKY'S OWN pks2map (repo Fortran + bundled HEALPix + their
bbps table; cosmology set to websky values) on their halos_10x10 catalog:

1. **Repo code ≠ paper ≠ released map**: mmin = 2.5e10 HARDCODED (no 1e13 cut, no
   0.5′ cut anywhere); cosmology hardcoded Planck18 (not websky). Their own binary
   (all halos, uncompensated) produces a halo map **~25× the released halo
   component** at ℓ~900 — the released production applied cuts + compensation not
   present in this source.
2. **REAL BUG in their pks2map**: the redshift-cut compaction (pks2map.f90:182-186)
   copies posxyz/rth but NOT vrad (computed at load, pksc.f90:79) → for every halo
   past the first cut index, positions pair with the WRONG halo's velocity. ksz.fits
   is documented z<4.5 from a z≤4.6 catalog → cut active in production. DEMONSTRATED
   with their binary: zmax=4.5 vs 6.0 on identical input → map RMS 5.90e-8 vs
   9.87e-8 (**40% RMS loss**), C_ℓ suppressed 20-33% at ℓ≤1237.
3. **Checkmate test (job 4419750)**: our composite with velocities deliberately
   shuffled (emulating their bug): W 1.84→1.58 (cross term killed) but
   **Wc 1.383→1.373 — unchanged**. So our Wc mid-ℓ excess is NOT the halo-field
   cross and NOT 2-halo velocity coherence: it is the **1-halo (shot) power of the
   compensated Battaglia halos** — velocity-shuffle-insensitive (⟨v²⟩ preserved).
   Patch shuffle test confirms: cut-sample halo C_ℓ at ℓ≥632 is shot-dominated.

**Physical adjudication**: at ℓ~900 our Wc gives D̃ℓ ≈ 1.3-1.4 μK² — at the level
of the hydro-simulation band (Shaw ~1.2, Battaglia ~1.1-1.3 per their own fig 6)
— while ksz.fits sits at ~0.95. The released map's tiny mid-ℓ halo term requires
suppression (cuts ÷3.4 × strong compensation × velocity bug) that their available
code does not document. The truth for the Battaglia gas model plausibly lies at or
near OUR composite; ksz.fits is LOW at mid-ℓ.

**Production disposition**: our kSZ construction (cf32 field + Battaglia halos,
Wc zero-net compensation) is defensible and likely MORE correct than ksz.fits at
mid-ℓ; the reference itself carries a demonstrated velocity-association bug and
undocumented production filtering. Remaining open PHYSICS question (not a
reproduction question): the exact 1-halo/field double-count treatment at mid-ℓ —
bracketed by field-only (low) and Wc (high); best settled against hydro kSZ
templates, not ksz.fits.

## Battaglia-τ implementation audit across codes (2026-07-27)

Comparing Battaglia-2016 gas-density τ across the three implementations:

| | outer exponent | truncation | amplitude convention |
|---|---|---|---|
| B16 paper (eq A1) | −(β−γ)/α = **−4.58** | — | ρ̄fit·(ref. density; fb convention TBD) |
| Fortran/Websky | plain −β = **−3.83** (B16 form COMMENTED OUT in bbps_profile.f90) | x≤4 sphere | ×fb, ρ̄m-referenced, mean-Δ radius in amplitude |
| XGPaint | (β−γ)/α (neg-β storage) = **−4.13** (γ sign slip; should be (β+γ)/α) | NONE (∫ to ∞) | no fb, ρ_crit-comoving; measured 2.3-2.8× Fortran τ at z≤0.5 center, ×4-8 in outskirts |

Point-wise measurement (xgp_tau_check): XGPaint/ours = 1.0-2.8 at x_b=0.25
(z-dependent, tracks ρcr/ρ̄m), growing to 2.3-8 at x_b=3.5, ∞ beyond x=4.
**Arbitration anchor**: B16's own τ scaling (ln τ0 = −6.23, z=0.3, Θ=1.3′ →
τ̄≈2e-3) ≈ Fortran-convention amplitudes (ours ~2-3e-3 for matching clusters);
XGPaint would sit ~2.5× above. TENTATIVE: B16 shape + Fortran-anchored amplitude
= intended physics; both public codes deviate (differently). Worth reporting
upstream to XGPaint (slope sign + amplitude) alongside the pks2map velocity bug.

**Fork port ✅ DONE (2026-07-27, XGPaint fork commit 2059a6d, branch
lensing-kappa)**: `WebskyTauProfile` with (a) convention flag `:websky` (plain −β,
x≤4, Fortran amplitude) vs `:b16` (−(β−γ)/α paper bracket — production default);
(b) delta_comp zero-net compensation (uniform sphere of the FULL painted gas mass
at Δ×fb·ρ̄m, mirroring NFWKappaProfile); (c) kSZ via the existing velocity-aware
`paint!` with proj_v_over_c = −v_rad/c; (d) 20 anchor tests (test_websky_tau.jl),
all passing.

**AMPLITUDE ARBITRATION SETTLED by the τ0 anchor test**: aperture τ̄(Θ<1.3′,
z=0.3, M200c=3e14) vs B16's own ln τ0 = −6.23 → 1.97e-3:
  :b16 = 2.28e-3 (**1.16× — excellent**); :websky = 3.25e-3 (1.65×, the fat −3.83
  tail); XGPaint BattagliaTauProfile ≈ 2.5× (fails). So Fortran amplitude
  bookkeeping + B16 paper shape reproduces B16's own scaling relation — the
  TENTATIVE verdict above is now CONFIRMED. Gas ledger Mgas/(fb·M200c):
  :b16 = 1.09 (1e13) / 1.29 (1e15) (sane); :websky = 1.25 / 2.84.
  Central τ(3e14, z=0.55): 6.9e-3 (:websky) / 6.0e-3 (:b16);
  kSZ(1e14, z=0.5, 300 km/s) = 10.0 μK.
Implementation gotcha for reuse: quadgk on the compensated (zero-net) profile
needs an atol — pure rtol on a ≈0 integral never converges (caused a silent hang).
