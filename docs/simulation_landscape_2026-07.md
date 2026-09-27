# Simulation landscape: HalfDome, Backlight vs Websky/PeakPatch.jl (2026-07-25)

Notes from reviewing the algorithms of the neighboring mock-sky suites, and where our
pipeline sits. Sources at bottom.

## HalfDome (arXiv:2407.17462, JCAP05(2025)016; data.cmb-s4.org/halfdome.html)

Purpose: joint Stage-IV analyses (Rubin LSST × SO/CMB-S4/LiteBIRD × DESI/Euclid/
SPHEREx/Roman/PFS). 11 fixed-cosmology realizations + one f_NL=20 run; ~300 TB
released (full-sky lightcones + halo catalogs, z=0–4). IllustrisTNG cosmology.

Algorithm:
- **Gravity: FastPM** (quasi-PM): 6144³ particles, 3.75 Gpc/h box (m_p = 1.95e10
  Msun/h), force mesh 2× particle grid (12288³), **60 steps linear in a** from z=9.
  FastPM's modified kick/drift forces large-scale growth to track the exact
  linear/Zel'dovich solution regardless of step count; halo interiors unresolved.
- **ICs**: CLASS linear P(k); start z=9.
- **Lightcone: box REPLICATION** — box tiled ~2.6× per dimension to cover z=0–4;
  particles placed on the lightcone on the fly, interpolated between time steps;
  redshift-dependent downsampling targeting ℓ_max ≈ 10⁴.
- **Halos**: relaxed-FoF (RFoF) on the fly; **abundance-matched per z-slice to
  Tinker08 M200m**. Reliable above ~320 particles → floor ~6e12 Msun/h (numerical
  artifacts below). ROCKSTAR planned for later releases.
- **Validation**: P(k) within 4% of Aemulus-ν to k~1 h/Mpc; C_ℓ ~10% vs Halofit to
  ℓ~5000; halo power vs TNG300-1-Dark; MF vs Watson/Tinker.

## Backlight (no release paper as of 2026-07; cited in ACT arXiv:2401.13033)

- "Upcoming suite of non-Gaussian sky simulations", CITA/Websky lineage.
- Defining feature: **ensemble size — >1000 realizations planned** (40 existed in
  early 2024), each with CMB-lensing maps + halo catalogs. ACT used it for
  patchy-screening stacking-bias tests precisely because Websky/Agora provide only
  ONE sky each (can't separate bias from noise fluctuations).
- Presumably peak-patch-based (personnel + the ~1000× speed economics; a full
  N-body ensemble of that size is not feasible).

## Comparison

| | HalfDome | Backlight | Websky / PeakPatch.jl |
|---|---|---|---|
| Method | FastPM (60-step PM) | peak-patch (presumed) | peak-patch + 2LPT field |
| Mass floor | 6e12 Msun/h | ? | **~1.2e12 Msun** |
| Lightcone | box replication (2.6³) | ? | single box, 8 corner octants, no repetition |
| Realizations | 11 | **>1000** | 1 (Websky); ours cheap to scale |
| Halo masses | RFoF + Tinker AM | ? | peak-patch RTHL + Tinker AM |
| Foregrounds | future paper | non-Gaussian sky maps | Battaglia/Shang painting (ours: κ/kSZ/tSZ/CIB/ISW validated 2026-07) |

Takeaways for us:
1. **Everyone abundance-matches to Tinker** — HalfDome does per-slice AM on RFoF
   masses exactly as we do on RTHL masses. AM is the field-standard calibration,
   not a peak-patch crutch: no fast method gets absolute halo masses natively.
2. **Lightcone geometries trade different artifacts.** Replication repeats
   structure every 3.75 Gpc/h along the LOS (spurious radial correlations at the
   box scale — directly relevant to kSZ-Doppler LOS-cancellation physics, cf. our
   coarse-grid velocity-coherence saga). The octant trick avoids repetition but
   caps z at the box diagonal and yields one realization per box.
3. **Mass reach**: HalfDome's floor is ~5× shallower than ours — they cannot
   directly paint the CIB (Shang Mpeak ≈ 2e12 Msun); we can (validated vs
   cib_nu0545 at 1.06 relative to the websky-halos-through-XGPaint baseline).
4. **Our GPU pipeline sits in Backlight's niche**: ~3h/octant catalog + ~1.3h for
   the 5-kernel map set on ONE L40S → a ~100-realization full-sky ensemble is a
   modest allocation. Strongest scientific argument for the production hardening
   (frozen config, coarse_factor=32, acceptance suite, 8-octant campaign).

Sources: arXiv:2407.17462 (HalfDome), data.cmb-s4.org/halfdome.html,
halfdomesims.github.io, arXiv:2401.13033 (ACT screening; Backlight description),
arXiv:2001.08787 (Websky), arXiv:1810.07727 (mass-peak-patch algorithm).
