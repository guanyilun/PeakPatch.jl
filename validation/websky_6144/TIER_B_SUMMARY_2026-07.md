# Tier-B (map painting) validation summary — final state 2026-07-26

Definitive reference for the map-level Websky reproduction campaign. Detailed
narratives (with the dead ends and corrections) live in the per-product MDs;
this file records the settled constructions, final numbers, and open items.
Catalog-level (Tier-A) validation: see `websky_comprehensive_validation` memory
note — six statistics match Websky at the 1-5% level.

Conventions throughout: octant 000 of the finecell run (box 5236 Mpc/h = 7700 Mpc,
cellsize 0.852 Mpc/h, seed 12345, z<4.6), observer at (−2618)³, AM'd catalogs.
References: the released Websky v0.0 maps (Nside 4096) at mocks.cita.utoronto.ca.
Comparisons: 6 gnomonic 8° caps (flat-sky C_ℓ; cap scatter = error bar) and/or
whole-octant pseudo-C_ℓ (apodized mask, f_sky = 1/8). Different realizations —
ratios ~1 with realization-level scatter are the pass criterion.

## Scorecard

| Product | Construction | Result vs released map | Status |
|---|---|---|---|
| κ (lensing) | full-matter 2LPT field + truncated-NFW halos − Δ=3 sphere of FULL painted mass (zero net) + Gaussian z>4.5 tail | **0.99–1.12 (ℓ=165–884)** vs kap_lt4.5; ~1.20 at ℓ≥1730 (plausibly their sub-pixel mass loss); +tail → 1.04–1.06 vs kap.fits | ✅ closed |
| kSZ field (low ℓ) | full-matter :ksz kernel, coarse_factor=32 | 1.05–1.22 at ℓ≤320, ~1.0 at ℓ=450 (residual = z_max 4.6-vs-4.5 +10% + octant realization) | ✅ closed |
| kSZ composite (high ℓ) | + Battaglia AGN gas halos (pks2map-verbatim), M200c>1e13 + 0.5′ cuts | 1.09–1.17 at ℓ≥2400 (Wc) | ✅ closed |
| kSZ composite (mid ℓ) | same | 1.31–1.54 at ℓ=450–1750 (Wc) | ⚠️ documented model choice (see below) |
| tSZ (Compton-y) | halo-only Battaglia-2012 pressure (pks2map-verbatim), all halos | mean-y ratio **1.096**; C_ℓ 1.18–1.28 in well-measured bands (ref pixwin uncorrected) | ✅ closed (model-faithful) |
| CIB (545 GHz) | XGPaint CIB_Planck2013 (= Websky model) + Websky completeness cut | mean intensity **1.254** vs cib_nu0545; websky-halos-through-XGPaint baseline = 1.18 → our catalog contributes **1.06** | ✅ closed |
| ISW | new :isw FieldMap kernel — per-tile ∇⁻²δ₀ painted as Lagrangian cell value | same-mask **0.66–1.2 over ℓ=11–257**, tracking 5 decades of C_ℓ | ✅ closed |

## Key constructions decoded (things the paper does not state plainly)

- **κ = compensated, not excluded**: the §3.2.4 Δ=3 compensation applies to κ, with
  the sphere carrying the FULL painted NFW mass (1.636·M200m), zero net per halo.
  The §3.1.3 "exclusion" reading, and the paper-literal "same mass" sphere (B2),
  both overshoot by exactly the +0.636·M/halo mass ledger. Decided empirically
  against kap_lt4.5.fits (the release conveniently splits the z>4.5 tail — which
  is 15–34% of kap.fits C_ℓ at ℓ=50–520, NOT small).
- **kSZ halo model** = Battaglia+ AGN gas-density gNFW (bbps_rhotilde: P0=4e3·
  m14^0.29·zp^−0.66, αρ=0.88·m14^−0.03·zp^0.19, β=3.83·m14^0.04·zp^−0.025,
  xc=0.5, γ=0.2), x=r/R200c, spherical truncation x≤4, amplitude
  −(v_r/c)·tau0persigma·(Ωm h²)^(2/3)·(1+z)²·(mh/200)^(1/3)·fb·Σ̃, with
  mh = √(Δcrit(z)/200)·M_RTH (SIS proxy, Bryan–Norman). All verbatim from
  ~/work/peakpatch/src/pks2map/. XGPaint's BattagliaTauProfile has identical
  parameters (independent cross-check).
- **ksz.fits is z<4.5 late-time only** (ksz_patchy is a separate uncorrelated
  reionization model) — apples-to-apples with our octant, no tail correction.
- **Websky's low-ℓ kSZ is the eq-4.3 "Doppler" artifact** of the sharp z=4.5
  cutoff (broken LOS-velocity cancellation), dominant at ℓ≲300.
- **CIB**: XGPaint CIB_Planck2013 IS the Websky model (Shang HOD, z-dependent
  Td). XGPaint's HEALPix paint! outputs MJy/sr directly (fluxes are MJy — do NOT
  divide by Ω_pix again; that division is a no-op at Nside 1024 by coincidence).
- **ISW**: linear-potential LOS integral; must be painted at LAGRANGIAN positions
  (Eulerian deposit imprints spurious δ×φ power at ℓ≳50, up to 900×).

## Pipeline defects found & fixed during this campaign (all at source)

1. **Coarse-grid velocity decoherence** (the big one): coarse grid M=ntile×cf;
   at cf=4 (64³, Nyquist 0.038 h/Mpc) the k~0.03–0.1 band is carried by per-tile
   noise → 2–3× kSZ Doppler excess at ℓ≲350. Density/κ blind to it (which is why
   all prior validation passed). **Production rule: coarse_factor≈32 for any
   velocity/potential-weighted product; cf=4 fine for catalogs/density.**
2. ISW Eulerian-deposit δ×φ contamination → Lagrangian deposit for pw kernels.
3. 8× multi-worker fieldmap overcount (Julia closure scoping trap; see
   julia_closure_scoping_trap memory).
4. Halo-κ kernel /χ-vs-×χ + circular self-test → independent amplitude anchors
   now required in all painters.
5. Sub-pixel profile undersampling → per-halo mass-exactness residual deposits
   (pks2map lacks this — plausibly why kap_lt4.5 is LOW at ℓ≥1730 vs us).
6. CIB double division by Ω_pix (units audit via sum(map)·Ω_pix == sum(fluxes)).

## The one open item: kSZ mid-ℓ (ℓ=450–1750, Wc/ref = 1.31–1.54)

Every input is now cross-validated between the pipelines: catalog N(>1e13, z)
ratio 1.00 (0.93–1.05 per Δz=0.25); websky's own halos through our painter give
1.16–1.75× our halo term (patch realization variance); velocities validated
(σ_vr, v12, their vrad = v·r̂/c convention identical). Yet a faithful Battaglia
painting of EITHER catalog produces a halo term 5–10× larger at ℓ~500–1000 than
the released map's own decomposition (kszcomp.pdf: halo ≈ 5–10% of total there).
Compensation variants tested: W (none) = 1.84, We (paper-literal Δ=3-mean sphere,
gas mass fb·mh) = 1.58, Wc (zero-net full-τ sphere) = 1.38 band mean. None reach
their tiny halo fraction. Audit of the literal prescription: gasr = painted-gas/
(fb·mh) < 1 at z≳1.5 (over-compensation, net negative halos); their maptable also
truncates profiles at 4 Mpc transverse (−10–29% cluster τ). CONCLUSION: the
remaining gap is an unreleased production detail (actual mmin, stronger
compensation, and/or stacked truncations) — not recoverable from the paper, the
Fortran source, or websky_model. For our production maps, Wc (physically
consistent, zero-net electrons) is the documented choice.

**Lessons that generalized** (worth keeping front-of-mind):
- End-to-end ratios ≈ 1 can hide compensating errors (scheme A vs kap.fits =
  two ~40% errors cancelling); always test against component-split references.
- Attribution ≠ diagnosis: two plausible attributions (halo velocities; catalog
  abundance) were REFUTED by direct tests. Run the discriminating experiment.
- Check comparison-patch geometry before trusting densities (the 10×10 patch is
  corner-anchored, not centered — caused a clean ×1.93 phantom).
- Login node kills silently at 4–6 GB; 1-GPU SLURM slices backfill in ~1 min
  vs hours for 4-GPU requests; job time limits: cf32 catalog needs >3h (7h13m).

## Product recipes (for the production campaign)

- Catalogs: run_gpu_octant + apply_abundance_match (cf=4 OK; cf=32 if halo
  velocities feed velocity-sensitive statistics). New convention: Eulerian + km/s.
- Field maps: run_fieldmap_octant, kernels [:kappa,:mass,:tau,:ksz,:isw],
  coarse_factor=32, gpu_paint, Nside 4096 (~1.3 h/octant on one L40S).
- κ map: field κ + XGPaint NFWKappaProfile(delta_comp) halos + kap_gt4.5-style
  Gaussian tail (or the released kap_gt4.5.fits itself).
- kSZ: field ksz + compare_composite_ksz_battaglia.jl-style Battaglia halos (Wc).
- tSZ: compare_tsz.jl painter (all halos). CIB: run_cib.jl via XGPaint (wcut
  optional — only for websky comparison; omit for Tinker-complete skies).
- ISW: :isw kernel map directly (Nside 1024 suffices).

Files: COMPOSITE_KAPPA_2026-07-18.md, KSZ_COMPOSITE_2026-07-19.md,
FIELDMAP_KAPPA_2026-07-17.md, KAPPA_PAINTED_2026-07-16.md, compare_*.jl,
run_*.slurm (this directory); docs/simulation_landscape_2026-07.md (HalfDome/
Backlight context).
