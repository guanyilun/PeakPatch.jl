# Websky paper: comparisons needed before we call it done (2026-09-26)

Review of `paper/` (§1–6 drafted, §7–10 stubs) against `docs/paper_plan.md`,
`validation/websky_6144/TIER_B_SUMMARY_2026-07.md` and the production README.

## Where we stand

**Solid (keep, re-render):** Fortran finder equivalence at matched cellsize (0.5%);
Tier-A six statistics vs Websky (MF, dN/dz, σ_vr, ξ/bias, b(M), v12); Tier-B
scorecard (κ, kSZ field/composite, tSZ mean-y, CIB, ISW); cf convergence for kSZ;
RNG large-angle null test (PASSED); code-level tests (1061 + new finalize/AM/fieldmap).

**The central problem:** every number in §7–8 comes from the *validation octant*
(oct000 finecell, z<4.6, pre-fix AM, cf4/cf16/cf32 mix), not the frozen campaign
(8 octants, z<4.5, cf32, AMv2). The frozen campaign has raw catalogs, AMv2
(running) and 40 field maps, but **no halo-painted maps yet** (κ composite, kSZ
composite, tSZ, CIB). The paper cannot quote its dataset until they exist.

**Data risk:** Websky reference files live on purged scratch
(`scratch/websky_6144/websky_ref/`); `kap.fits`, `kap_lt4.5.fits`, `kap_gt4.5.fits`
are already gone; `/project/.../websky_ref/` is empty. Only a 10°×10° Websky
halo patch exists locally (no full catalog).

## What the original Websky paper validated (Stein+2020, arXiv:2001.08787)

Source: `~/work/peakpatch/2001.08787.txt`, §4. **Every Websky validation is
map-level.** The paper has no halo-catalog plot at all: no mass function, dN/dz or
clustering. For those it defers to Stein+2019 and the Euclid comparison project.
Most of its references are observational data or hydrodynamical simulations, not
theory.

| Websky item | Compared against | Their result | Our status → action |
|---|---|---|---|
| Fig 5 tSZ C_ℓ^yy, ℓ=10–10⁴ | Planck 2015 y, Bolliet+18 re-analysis, ACT D₃₀₀₀(148)=3.4±1.4, SPT D₃₀₀₀(143)=4.08 | "good agreement"; slight excess at ℓ=3000 | have map-vs-map only → **add the same data points** |
| Fig 6 kSZ total/halo/field | 6 hydro sims at ℓ=3000 (Battaglia10, Trac11, Shaw12, Dolag16, Roncarelli17, Park18, CSF-scaled via Park+18 Tab.1) + analytic Doppler (eq 4.3) + OV (eq 4.4) | field dominates at all ℓ; halo ≈5/6 of total at ℓ=5000; "rough agreement" | hydro overlay was optional B5 → **must-do**; **compute eq 4.3/4.4 for our cosmology and z_max=4.5** (a theory anchor for the field kSZ) |
| Fig 7 CIB C_ℓ at 143/217/353/545/857 | Planck 2013 CIB (radio shot noise subtracted), HerMES (Viero+13, 600/857), Lenz+19 | good at all ν; slight Poisson excess at 857; L₀ tuned to Planck 545 at ℓ=500; **400 mJy pixel flux cut** | we have 545 mean intensity only → **5 frequencies, same flux cut, C_ℓ** |
| Table 1 CIB decoherence ⟨C^{νν'}/√(C^{νν}C^{ν'ν'})⟩, 150<ℓ<1000 | Planck 2013, Lenz+19 | similar | **new, cheap once the CIB maps exist** |
| Fig 8 κ C_ℓ + halo/field/z>4.5 components | CAMB Halofit (Takahashi) Limber, cosmic-variance errors | few % at ℓ<1000; **~20% suppression at the smallest ℓ** | **key**: our κ is ~1.20× theirs at ℓ≥1730, matching the suppression they state themselves → compare *both* to Halofit (neutral resolution of item B4) |
| Fig 9 lensed CMB TT/EE/BB | CAMB lensed spectra | excellent at ℓ<5000 | not done → optional (pixell lensing of our κ) |
| Fig 10 CIB×φ, 100–857 GHz | Planck 2013, Lenz+19 × Planck lensing | good; low-ℓ turnover matches | **must-do**: it tests correlations between maps |
| Fig 11 y×T_ν (T=CIB+y), 143–857 | Planck 2015 CIB-cleaned (Planck XXIII 2016) | broad agreement | **must-do** |
| §4.1 cost | 3.84 h on 1128 Skylake cores = 4336 core-h; peak memory 7.67 TB; 5.9 TB of stored ICs | — | direct comparison point for §9 (ours: ~8×7 h×4 L40S at cf32 on 512 GB nodes; no IC storage because the fields are regenerated from the counter RNG) |
| §4.1 geometry | eight corners of one periodic 6144³ box, identical to ours | *claims* periodic effects are only at ℓ<10, but **never tested** | our A6 replication test checks this claim for both codes |
| §4.2 AM | Tinker M200m, Δz=0.1 bins, bilinear interpolation, keep halos with >10 particles before AM; M_min table released | 1.2e12 (z<4) → 4e12 (z=4.6) | same scheme → compare our completeness to their `halo_mass_completion.txt` |

**What this implies for our paper:**
1. Structure §8 as "Websky Figs 5–11 reproduced, with our curve added". Our map against
   theirs is one line in each panel, alongside the same data, hydro and theory references.
   This makes the comparison paper-to-paper and needs no claims about their code.
2. §7 (catalog statistics), ISW, the Fortran finder equivalence, tiling invariance and the
   coarse-factor criterion are **new relative to Websky**. Present them as added validation.
3. Revise section C: overlaying the *same* observational points Websky used is cheap and
   expected, so it is no longer out of scope. Fitting to data still is.
4. Our z_max=4.5 maps match their maps (halos to z=4.6, maps to z=4.5). The halo
   selection is the same too: kSZ halos M200c>1e13 with r200c>0.5′; κ halos subtending >1 pixel.

Reference data to collect (in `validation/websky_6144/refdata/`, with provenance):
- Planck 2015 y C_ℓ; Bolliet+18; ACT/SPT D₃₀₀₀ values.
- Park+18 Table 1 hydro kSZ points, with the Shaw+12 CSF scaling.
- Planck 2013 CIB auto (XXX) and CIB×φ; HerMES Viero+13; Lenz+19 auto and ×κ.
- Planck 2015 y×T (XXIII).

## A. Must-do (blocking)

| # | Item | Why | Cost |
|---|---|---|---|
| A0 | Re-download Websky refs (kap*, ksz, tsz, cib at several ν, isw; ideally full `halos.pksc`) to `/project/.../websky_ref/`, with checksums | needed by every map comparison | download only |
| A1 | Build frozen-campaign map suite from AMv2: κ composite, kSZ composite (Wc), tSZ, CIB (545 + 2–3 more ν) for all 8 octants, then assemble full-sky maps | the paper's dataset; F13 gallery | ~8×(few GPU-h) per product, b1 backfill |
| A2 | Rerun all Tier-A statistics on AMv2 (per octant plus full sky) and all Tier-B spectra on **full-sky vs full-sky** (replaces caps and octant pseudo-C_ℓ) → regenerated scorecard T1 | final numbers | CPU, mostly scripts that already exist |
| A3 | **Error model and pass criteria**, written down before A2: Knox/Gaussian errors plus octant-to-octant scatter for C_ℓ; jackknife for catalog statistics; two realizations (ours vs Websky) → the difference variance is 2× | reviewers will ask "within what?" | small |
| A4 | **Theory anchors independent of Websky**: κ C_ℓ vs CAMB+HMcode Limber (z<4.5); raw (pre-AM) MF vs Tinker per z-bin (finder performance; AM hides it); b(M) vs Tinker10; IC field P(k) vs input; ψ1/ψ2/velocity power vs linear/2LPT theory | shows agreement with Websky is not a shared systematic | small–medium |
| A5 | **Algorithm figures (§3, the core contribution)**: F2 coarse/fine splice cross-section and P(k) continuity across tile boundaries; **tiling independence at catalog level** (same seed, ntile 8 vs 16 on a mid box: matched-halo fraction, mass and position scatter) | the paper's main claim is shown only indirectly now | 2–3 small GPU jobs |
| A6 | **Full-sky assembly checks**: octant seams (map continuity and low-ℓ power at boundaries); **box-replication/periodic-image** correlations (shell cross-correlation between octants beyond χ≈2618; Websky shares the geometry, so compare both) | we release full skies; open item since the RNG test | medium; needs design (see Q2) |
| A7 | **Performance** (F14/T3): per-stage timing of the frozen campaign from logs; strong scaling 1→4 L40S; one H100 node; peak memory vs tile size; Fortran/CPU reference cost (CITA numbers or a stated estimate) | §9 is empty; weakest material | ~1 day queue |
| A8 | Convergence appendix: nbuff (16/24, the oct000 runs exist), filter-bank count/spacing, cellsize (reframe the factor-h runs as a controlled test), coarse_factor for catalogs as well as kSZ | referee staples | mostly existing runs + 1–2 new |

### A3 adopted (2026-09-26, fixed before any production-campaign numbers were computed)

**Estimator.** Full-sky C_ℓ from the 8-octant assembled maps against the full-sky
released maps. No mask: both skies are complete. Pixel windows are divided out for
both, and the Websky tSZ map is Nside 2048. Bins are Δℓ/ℓ ≈ 0.1 with unit weight per mode.

**Errors (two independent realizations).**
- Gaussian: σ²(Ĉ_ℓ) = 2C_ℓ² / [(2ℓ+1)Δℓ] for each map. The ratio
  r = Ĉ^ours/Ĉ^ref then has σ_r ≈ r·√(4/[(2ℓ+1)Δℓ]). Cross-spectra use the Knox form
  (C_AB² + C_AA·C_BB)/[(2ℓ+1)Δℓ].
- Non-Gaussian (tSZ, kSZ halo term, CIB shot noise, κ at high ℓ): an **octant jackknife**,
  i.e. the scatter of the eight per-octant C_ℓ (f_sky = 1/8 each, apodized mask), divided
  by √8. This is applied to *each* map and the larger of the Gaussian and jackknife
  errors is used.
- Caveat: our octants are views of one box. Beyond χ≈2618 Mpc/h they share structure,
  so the jackknife underestimates the variance at ℓ ≲ 30. There we quote the Gaussian
  error only and flag it.

**Pass criteria.**
- *Reproduction claims* (same model painted on both catalogs, e.g. κ, tSZ, CIB,
  kSZ field): band-averaged |r−1| < 3σ_r over the pre-declared ranges below, and no
  single band beyond 4σ.
- *Theory claims* (κ vs Halofit, kSZ field vs Doppler+OV, P(k) vs input): the same
  test, using our realization's error only.
- *Construction-dependent differences* (kSZ mid-ℓ halo term, κ ℓ≳1700 vs the released
  map): not pass/fail. Report them with the discriminating test.
- Pre-declared ranges: κ ℓ=30–3000; tSZ ℓ=100–5000; kSZ field ℓ=30–500; kSZ total
  ℓ=500–8000; CIB ℓ=100–3000 (clustered) and 3000–8000 (Poisson); cross-spectra
  ℓ=100–2000 (the Planck range in Websky Figs 10–11); ISW ℓ=10–250.

**Catalog statistics.** Mass function and dN/dz get Poisson errors plus octant-jackknife
sample variance; ξ(r), b(M) and v12 get the octant jackknife. Pass: 3σ in each
pre-declared bin (as in the Tier-A notes).

## B. Should-do (strengthens; cheap relative to value)

- B1 **Cross-spectra vs Websky**: κ×tSZ, κ×CIB, tSZ×CIB, κ×halo-density. These test
  whether the products are correlated correctly, which auto-spectra cannot.
- B2 **One-point statistics**: κ and y PDFs, tSZ/κ peak counts (tSZ is strongly non-Gaussian).
- B3 tSZ C_ℓ with the pixel window corrected. The current 1.18–1.28 is quoted uncorrected.
- B4 κ high-ℓ excess (~1.20 at ℓ≥1730): resolve with a sub-pixel/pixel-window test or state it as a resolution effect.
- B5 kSZ mid-ℓ hydro-template overlay (Shaw+12, Battaglia+10), as planned for F10.
- B6 Full-sky halo-catalog comparison, if we get the full Websky `halos.pksc`
  (MF/dN/dz over 4π instead of a 10°×10° patch).

B1 (cross-spectra) and B5 (hydro overlay) are now must-do: Websky Figs 10, 11 and 6
(see above). B4 is resolved by comparing to Halofit (Websky Fig 8).

## C. Out of scope (state it in the paper)

*Fitting* to observations (retuning L₀ or pressure profiles); patchy kSZ; radio sources;
post-Born lensing; z>4.5 beyond the Gaussian κ tail. We overlay the observational points
Websky used, for context, but make no new claims about data.

## Order

1. A0 and A3 now, while AMv2 finishes. Also move the remaining scratch refs to `/project` today.
2. A1 (queue-bound) runs in parallel with A4, A5 and A8 (small jobs, independent of AMv2).
3. A2 and B1–B3 once A1 lands. This regenerates T1 and F5–F13.
4. A6 and A7 (A6 needs a design decision). Then §7–9 prose and the neutral-framing pass.

## Open questions for the user

- Q1: Pass criteria: ratio within N·σ of the combined two-realization error, or keep the
  current "≈1 within realization scatter" wording?
- Q2: Replication test design: full-sky shell cross-correlation between octants
  vs. the same statistic on Websky (null = Websky's own level)?
- Q3: Include B-items in v1 of the paper, or defer to a follow-up?
