# Same-IC N-body calibration of Peak Patch — plan

Drafted 2026-10-02. Status: PROPOSAL (nothing run yet).

## 0. Motivation

Peak Patch is a hand-built, few-parameter map from the linear field to halos to observables:

    δ_L ──[filter bank]──> F(R), strain ──[ellipsoidal collapse]──> candidate peaks
        ──[exclusion / volume reduction]──> halos (q, R_TH)
        ──[2LPT]──> Eulerian x, v ──[AM]──> masses ──[profiles]──> maps

Several stages contain hand-set choices whose effect we have already measured at the several-%
level in the Websky campaign: filter bank (~1/3 of the full-sky ξ gap), R_TH cap / nbuff, exclusion
rules, get_homel details, and a blanket Tinker AM that absorbs everything else
(see memory notes: clustering-excess attribution, nbuff mass cap, MATCHED_FORTRAN_2026-10.md).

Until now everything has been tuned with **summary statistics** (HMF, ξ, b(M)). An N-body run from
the **same initial conditions** lets us compare **halo by halo**. That separates selection errors
from mass errors and displacement errors, which summary statistics cannot do.

### Science questions

- **Q1 (fidelity):** At fixed ICs, what are Peak Patch's completeness, purity, mass scatter and
  position/velocity errors vs a full N-body run, as a function of mass and z?
- **Q2 (attribution):** Which stage dominates each error (selection, mass, displacement)?
  In particular, what drives the ~5% ξ excess we see against Websky?
- **Q3 (calibration):** How much does a few-parameter calibration of the existing stages remove,
  and do the calibrated parameters transfer across seeds, redshift, resolution and **cosmology**?
- **Q4 (learned residual, stretch):** What does a small learned correction on per-peak features
  add beyond Q3, and which features does it use? That points at missing physics.

## 1. Data: Quijote

Primary set: **`fiducial_HR`**: 100 realizations, 1024³ particles, L = 1 Gpc/h, 2LPTic ICs at z=127,
Gadget-III, snapshots at z = 0, 0.5, 1, 2, 3.
Cosmology: Ωm=0.3175, Ωb=0.049, h=0.6711, ns=0.9624, σ8=0.834.

Transfer set: **`latin_hypercube_HR`**: 2000 cosmologies at 1024³ (seed per sim in `ICs/2LPT.param`).
Resolution-dependence check: fiducial 512³ (same seeds as HR? **verify** in the 2LPT.param files).

Per simulation:
- `ICs/`: ics.* (positions, velocities, IDs at z=127), `2LPT.params`, `inputspec_ics.txt`, generator log
- FoF catalogs; **`FoF_id/`** stores member particle IDs
- Rockstar catalogs plus `out_X_pid.list` (member IDs)
- Particle IDs are consistent across snapshots, so each halo member can be traced back to its IC (Lagrangian) position.

| quantity | Quijote HR | Websky production |
|---|---|---|
| cell / interparticle spacing | 0.977 Mpc/h | 0.852 Mpc/h |
| Rf,min = 2 a_latt | 1.95 Mpc/h | 1.70 Mpc/h |
| PP mass floor (≈ M(Rf,min)) | ~2.7e12 M☉/h | ~1.2e12 M☉/h |
| N-body particle mass | ~8.2e10 M☉/h | — |
| 100-particle halo | ~8e12 M☉/h | — |

So the HR grid nearly matches our production resolution. That's a nice coincidence: a
calibration done on HR transfers to production with only a ~15% change in resolution.
Halos with ≳50–100 particles (≳4–8e12) are robust, which covers the tSZ/kSZ mass range. The
CIB range (~1e12) is **not** testable at HR. That would need AbacusSummit or our own run (Phase 5).

**Access:** run on **rusty** (decided 2026-10-02); Quijote lives there. Killarney has no Globus
collection (Alliance page lists it as TBA) and no Globus tools. Size: an IC snapshot at 1024³ is ~1.07e9 × (12+12+8 B) ≈ 34 GB
before compression. Start with **3 realizations** (~150 GB including halo catalogs and IDs).

## 2. Phase 0: Reconstruct δ_L on the 1024³ grid

Two independent routes, which cross-check each other:

- **Route A (exact phases):** build 2LPTic, run it with each simulation's `2LPT.params`, and patch it to
  dump δ(k) (or the 1LPT displacement field) before particle generation. This gives bit-identical phases.
  Requires reading 2LPTic's seed table / k-plane RNG ordering. No Julia port is needed; just dump to HDF5.
- **Route B (inversion):** from ics.* get lattice q (from the ID ordering, **verify**; or round(x/a),
  since the z=127 displacement ≪ a), compute ψ = x − q (minimum image), then δ₁ = −∇·ψ₁ in Fourier space.
  Remove the 2LPT piece iteratively (ψ₂/ψ₁ ~ (3/7)D ~ 3e-3 at z=127; one iteration suffices).
  Positions are lossy-compressed (6 bits nulled at 1024³); the error is ≪ a, so it's harmless.

Scale to z=0 with D(0)/D(127). Quijote builds its ICs by rescaling CAMB's z=0 P(k) with linear
growth, so this is self-consistent by construction. Use the same D(z) they used and check against
`inputspec_ics.txt`.

**Gates (must pass before Phase 2):**
1. P(k) of the reconstructed δ_L matches inputspec × D² within cosmic variance up to k_Ny/2.
2. A vs B cross-correlation r(k) > 0.999 below k_Ny/2. Handle Nyquist the way
   NOTES_2LPT_NYQUIST_2026-09-24.md says.
3. Our `displacements_1lpt/2lpt(δ_L)` reproduces the measured IC ψ to < 1% RMS (closes the sign/normalization loop).

## 3. Phase 1: Plumbing in PeakPatch.jl

Hook point: `run_tile` in `src/Pipeline.jl`. Phase 1 there generates `delta` with
`generate_grf[_lcg]`. Everything downstream (rfft, LPT, filters, shells, merge, write) is
field-agnostic.

1. **External field input:** `[files] delta_in = "..."` (HDF5, Float32, Mpc/h, z=0 linear). If set,
   skip the GRF and load it instead. Assert grid size and boxsize against the config.
2. **Periodic single box:** right now `find_peaks` skips an `nbuff` border (PeakFind.jl:27-31) and shell
   analysis assumes buffered tiles. Quijote is periodic, so add `periodic=true` and use wrapped (mod1)
   indexing in the peak finder and radial shells. Alternative: pad the smoothed fields periodically.
   Without this we lose ~(1 − (1−2·nbuff/1024)³) ≈ 13% of the volume at nbuff=22, and edge halos
   become un-matchable.
3. **Configs:** a Quijote `pk` file from `inputspec_ics.txt`; a regenerated filter bank with
   Rf,min = 1.95 Mpc/h; keep our production spacing as the *baseline*. Check whether HomelTab needs
   regeneration (it should be cosmology-independent apart from the Ωm entering the collapse-time mapping, **verify**).
   ievol=0 at z = 0, 0.5, 1, 2 (matching the snapshots).
4. **Output:** the extended catalog (ioutshear=1) with Lagrangian q, R_TH, F, e, p, Rf, plus Eulerian
   x, v from `finalize_eulerian`. Also keep the **pre-exclusion candidate list** with flags, which
   Phase 3 needs to re-run exclusion cheaply.
5. **Fortran baseline:** run hpkvd on the same δ_L using the MATCHED_FORTRAN infrastructure, so every
   result has a Julia / Fortran / N-body triad.

Cost: 1024³ Float32 = 4.3 GB per field, ~20–30 live fields → fits one 512 GB CPU node. Expect
minutes per run on a GPU, tens of minutes on CPU. Cheap enough for hundreds of calibration evaluations.

## 4. Phase 2: Halo-by-halo comparison (the core measurement)

### N-body side
For each host halo (FoF b=0.2 from `FoF_id`; Rockstar hosts as the second definition), map member
IDs to IC lattice positions and compute:
- Lagrangian center of mass q_nb and R_L = (3M/4πρ̄)^{1/3}
- Lagrangian inertia tensor (shape e, p, orientation), compared with PP's strain e, p
- the Eulerian center and bulk velocity at the snapshot

Mass definitions: FoF M_fof, plus Rockstar M200m (matches the Tinker/AM convention). Report both;
the definitional spread is part of the answer.

### Matching
- Primary: **Lagrangian overlap**. f_halo = fraction of halo particles with |q − q_pk| < R_TH;
  f_peak = fraction of the PP sphere covered by halo particles. Match on max f_halo·f_peak with a
  greedy, mass-ordered one-to-one assignment.
- Secondary (robustness): |q_pk − q_nb| < R_L.
- Classify each object: matched, PP-only (false positive), N-body-only (missed), fragmented, or merged.

### Metrics (per mass bin, per z, mean ± scatter over realizations)
- completeness C(M_nb), purity P(M_pp)
- mass ratio ln(M_pp/M_nb): bias and scatter, before and after AM
- Lagrangian offset |Δq|/R_L; Eulerian offset |Δx| (2LPT vs N-body); LOS velocity residual
- shape: PP (e, p) vs the Lagrangian inertia of the N-body patch
- population level: HMF; halo auto-P(k) and **cross-correlation coefficient r(k)** for
  number-density-matched samples; b(k) ratio; ξ(r) on 1–20 Mpc/h (compare the 1–3 Mpc/h residual seen vs Websky)

### Attribution by hybrid catalogs (answers Q2)
Build catalogs that swap one ingredient at a time and measure the change in ξ / b:

| catalog | selection | mass | position |
|---|---|---|---|
| H0 | N-body | N-body | N-body |
| H1 | PP-matched → N-body | N-body | N-body |
| H2 | N-body | PP | N-body |
| H3 | N-body | N-body | PP 2LPT (from q_nb) |
| H4 | PP | PP (AM) | PP |

The steps from H0 through H4 split the clustering error into selection, mass-ranking and
displacement terms. This tests directly whether our Websky ξ excess is selection, as we
concluded on 2026-10-02.

## 5. Phase 3: Calibration (few parameters first)

Free parameters, roughly in order of expected leverage:

1. **Filter bank:** Rf,min multiplier, spacing ratio, Rf,max / R_TH cap (discrete → grid search)
2. **Collapse threshold:** a smooth deformation δ_c → δ_c · (1 + α₁ e² + α₂ p + α₃ ν⁻¹), or a
   multiplicative scaling of F_ev. Allow z dependence only if the data demand it.
3. **Exclusion / volume reduction:** the overlap-acceptance threshold in `lagrangian_exclusion!` /
   `volume_reduction!` (src/Merger/Exclusion.jl)
4. **Mass map (replaces blanket AM):** ln M = ln M_pp + g(ν, e, p, F) with a low-order g. The goal is
   that AM becomes a small correction instead of the main calibration.
5. **Displacement:** the smoothing scale applied to ψ at the halo location (multiplier on R_TH);
   optionally a 2LPT growth-coefficient tweak.

**Loss:** field-level terms, Σ_M w_M[(1−C) + (1−P) + σ²_lnM + bias²_lnM], plus summary terms
(HMF, b(M)) as regularizers. Each evaluation re-runs only the stages downstream of the changed
parameter (cache smoothed fields and candidate lists).

**Optimizer:** gradient-free (the pipeline has discrete steps). Grid search for the discrete
choices, Bayesian optimization / Nelder–Mead for the continuous ones.

**Validation protocol:**
- train on fiducial_HR realizations 0–4; test on 5–9 (seed transfer)
- check z transfer: train at z=0, test at z=0.5–2
- check resolution transfer: fiducial 512³ vs HR 1024³ (do the parameters scale with Rf/a_latt?)
- check cosmology transfer: ~20 LH_HR cosmologies spanning σ8 and Ωm. Ideally the calibrated
  parameters are functions of ν = δ_c/σ(R) and are otherwise cosmology-independent.
  If not, that's an important negative result for the "universal" Peak Patch claim.

### Phase 3b (stretch, Q4): learned residual
A small model (GBDT first, for interpretability) on per-peak features (F, e, p, Rf, ν,
neighbour density / tidal environment, distance to the nearest larger peak) that predicts
P(real halo) and a mass correction. Keep the physics skeleton; feature importances indicate
what is missing. Only consider a CNN on the Lagrangian patch if the GBDT residual is large and structured.

## 6. Phase 4: Feed back into production and the paper

- Apply the calibrated parameters to the Websky 6144 config (single octant first) and re-run the
  full-sky Tier-A comparison. Does the 0.948 ξ ratio close? Does AM shrink?
- Paper: Phase 2 metrics give a **methods-paper validation section** independent of Websky
  (which has no N-body counterpart). Phases 3–3b are probably a separate follow-up paper.

## 7. Phase 5 (later): reaching the CIB mass range, and painting

- Low masses (~1e12): AbacusSummit (m_p ≈ 2e9; ICs reproducible, but PLT corrections complicate
  "the linear field"), or our **own run**: write our δ_L, generate particles, run Gadget-4 /
  PKDGRAV3 at ~500 Mpc/h with 2048³. That gives zero convention risk.
- Painting stage: CAMELS (hydro/DMO twins, same ICs; small boxes but a full feedback-parameter grid)
  and FLAMINGO (Panphasia ICs, large volumes; check data availability). This is a separate plan.

## 8. Milestones and decision gates

| # | deliverable | gate | rough effort |
|---|---|---|---|
| M0 | data access confirmed (rusty/Globus), 1 realization staged | — | 1 day |
| M1 | δ_L reconstructed (Route B first, A later) | Phase 0 gates 1–3 | 2–4 days |
| M2 | `delta_in` + periodic mode in run_tile; tests | PP on Quijote runs; HMF sane vs Tinker | 2–3 days |
| M3 | matching + metrics pipeline on 1 realization, z=0 | C, P, mass scatter plots | 3–5 days |
| M4 | 3 realizations × 4 z; hybrid-catalog attribution; Fortran triad | **decision: is calibration worth it?** | 1 week |
| M5 | few-parameter calibration + transfer tests | parameters transfer across seeds | 2–3 weeks |
| M6 | Websky re-run with calibrated parameters | Tier-A ξ gap | 1 week |

Stop/continue at M4: if PP already gets C, P > 90% and σ_lnM < 0.2 above ~1e13, with clustering
errors under 2–3%, calibration is a polish step and Phase 2 is mostly a paper section. If the errors
are large and structured, Phases 3/3b become the main project.

## 9. Risks

- **δ_L reconstruction:** Nyquist handling and the D(z) convention. Mitigated by having two routes plus gates.
- **Halo-definition ambiguity:** FoF vs SO changes "truth" by ~10–20% in mass. Report both;
  calibrate to the definition the downstream painters assume (M200m → Tinker/Battaglia conventions).
- **Statistics at the high-mass end:** a 1 (Gpc/h)³ box holds only a few M > 1e15 halos. Stack realizations.
- **Fixed-z snapshots, not a lightcone:** fine for calibration. The lightcone machinery
  (per-tile threshold) is already validated separately against Websky/Fortran.
- **Overfitting / loss of cosmology dependence:** guarded by the LH transfer test. Prefer
  ν-parameterized corrections.

## 10. Questions for Francisco Villaescusa-Navarro (Quijote lead; shares an office with Yilun)

These replace most of the **verify** items above:
1. Paths on rusty for `fiducial_HR` ICs, FoF_id and Rockstar pid lists. Is a Globus endpoint to
   Killarney preferred instead? (Quijote is NOT mounted on Killarney, checked 2026-10-02.)
2. Does the particle ID order follow the IC lattice order (ID → q directly)?
3. Seeds: does fiducial_HR sim i share its seed with fiducial (512³) sim i? Which 2LPTic version and
   RNG/seed-table convention was used? Is there a known way to dump δ(k) (Route A)?
4. The exact D(z) / P(k) rescaling convention behind `inputspec_ics.txt` (z=127 vs z=0 normalization).
5. Recommended halo definition for this use: FoF b=0.2 vs Rockstar M200m; known problems at low particle number.
6. Interest in collaborating or co-authoring the calibration paper? Can he share the LH_HR subset and future
   higher-resolution runs?

## 11. Open decisions for Yilun

1. ~~Data access~~ DECIDED 2026-10-02: **run on rusty** (Flatiron), next to the Quijote data, so no transfer.
   CPU path (`run_tile`), 1024³ fits on one node. Clone the repo there and use the root project env
   (CUDA is a weak dependency, not needed). Killarney GPUs only if we scale to LH_HR / AbacusSummit.
2. Truth definition: FoF (b=0.2) or Rockstar M200m as the primary target?
3. Scope: Phase 2 as a section of the current methods paper, with calibration as a follow-up? (recommended)
