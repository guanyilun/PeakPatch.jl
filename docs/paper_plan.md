# PeakPatch.jl methods-paper plan (drafted 2026-07-27)

Decisions (with user): **full methods paper** (MNRAS/arXiv style, successor to the
Stein+ 2019 mass-peak-patch and Stein+ 2020 Websky papers); the **Websky
reproduction campaign is the core validation section**; reference-code
discrepancies get **neutral framing** (report our constructions and measured
ratios factually; describe convention differences between implementations without
characterizing other codes as buggy — pending the user's own validation of those
findings, which stay in internal notes).

Working title (strawman): *"PeakPatch.jl: GPU-accelerated peak-patch simulations
and CMB secondary-anisotropy maps in Julia"*.

---

## 1. Paper skeleton (section → existing material)

1. **Introduction** — need for many-realization mock skies (SO/CMB-S4/CCAT
   foregrounds, covariance ensembles); peak-patch heritage (BM96, Stein+ 2019);
   why a rewrite: GPU throughput, single-language stack, per-node ensembles.
   Landscape positioning: HalfDome, Backlight, Agora, Websky2
   (→ `docs/simulation_landscape_2026-07.md`).
2. **The peak-patch algorithm** (recap, cite Stein+ 2019 for details) — Gaussian
   ICs, top-hat filter bank, homogeneous-ellipsoid collapse (strain-based),
   binary exclusion+merging, 2LPT displacements. State our conventions:
   Mpc/h units (→ `CONVENTIONS.md`), Rf,min = 2 a_latt, filter spacing.
3. **Multi-resolution tile decomposition** — the core algorithmic contribution:
   MUSIC-style coarse-field + independent fine tiles, no global FFT; buffer
   geometry (nbuff), coarse–fine splicing, Threefry counter-based RNG
   (determinism = same seed → same fields at any tiling); **coarse-grid
   velocity-coherence requirement** (coarse_factor ≈ 32 for velocity/potential
   products vs 4 for density/catalogs — the cf4/cf16/cf32 kSZ convergence study
   is a *figure*, not a bug report: it derives a resolution criterion
   k_Nyq,coarse ≳ 0.3 h/Mpc from first principles + measurement).
4. **Implementation** — Julia/CUDA.jl design: kernel inventory (FFT convolutions,
   collapse table walks, peak find/merge, 2LPT), per-tile GPU residency, multi-GPU
   task parallelism, MPI extension (distributed FFT path), memory budget vs tile
   size; catalog output conventions (Eulerian positions + km/s velocities);
   abundance matching step (Tinker, Om_total).
5. **Lightcone construction** — single-box observer-at-corner octants, ievol=1,
   per-shell collapse thresholds, z_max; extended (shear) outputs.
6. **Map making** —
   a. Field maps (`src/FieldMap.jl`): eq-3.11-style LOS kernels :kappa/:mass/
      :tau/:ksz/:isw; GPU HEALPix painting (own ang2pix + atomics); sub-cell
      splitting; **Lagrangian vs Eulerian deposit** for potential-weighted
      kernels (ISW) — a genuine methods finding (Eulerian deposit imprints δ×φ).
   b. Halo painting via XGPaint fork: NFWKappaProfile (c=7 truncated NFW +
      r⁻² tail, Δ=3 full-mass compensation), WebskyTauProfile (Battaglia gas
      density; **convention flag :b16/:websky presented neutrally** as "two
      conventions exist in the literature/codes; we implement both and anchor
      to B16's published τ₀ scaling relation"), Battaglia-2012 tSZ pressure,
      CIB via XGPaint CIB_Planck2013.
   c. Composite recipes (field + compensated halos + Gaussian z>4.5 tail) —
      from TIER_B_SUMMARY product recipes.
7. **Validation I: catalogs (Tier-A)** — vs Fortran peak-patch at matched
   config (finder identical to 0.5%) and vs released Websky catalog: mass
   function (~1% with AM), dN/dz, σ_vr (1–2%), ξ(r)/bias (1.06; matched-geometry
   Landy-Szalay), b(M) (≤5% in 4 bins), pairwise v12 (0.997 in kSZ window).
8. **Validation II: maps (Tier-B)** — the scorecard vs released Websky v0.0
   maps: κ 0.99–1.12 (ℓ=165–884, vs kap_lt4.5); kSZ field 1.05–1.22 (ℓ≤320,
   cf32); kSZ high-ℓ 1.09–1.17 (Wc); tSZ mean-y 1.096; CIB 1.06 (catalog
   contribution, after the 1.18 painter baseline); ISW 0.66–1.2 (ℓ=11–257).
   **Neutral treatment of residuals**: kSZ mid-ℓ presented as "the halo-term
   amplitude at ℓ~500–1700 depends on compensation and mass-cut choices not
   fully specified by the reference release; our zero-net construction sits in
   the hydro-simulation band (Shaw/Battaglia templates)" — with the
   catalog-level cross-checks (N(>1e13,z)=1.00, painter cross-validation) shown
   as evidence the *inputs* agree. No claims about reference-code defects.
9. **Performance** — time-to-solution per octant / per full sky, single-GPU and
   multi-GPU scaling, memory footprint, comparison to Fortran-on-CPU cost
   (needs new benchmark runs — see §3 gaps).
10. **Summary & code availability** — GitHub, Zenodo DOI, docs, example configs;
    XGPaint fork status.
- **Appendices**: A) collapse-table construction; B) filter-bank/AM details;
  C) map-kernel derivations (incl. ISW weight, Doppler term of the sharp z-cut);
  D) reproduction configs for every figure.

## 2. Figure & table plan (~14 figures, 3 tables)

| # | Figure | Data status |
|---|---|---|
| F1 | Tiling/buffer schematic + pipeline flowchart | to draw (TikZ/inkscape) |
| F2 | Coarse–fine field splice cross-section + power-spectrum continuity | quick run needed |
| F3 | Filter bank σ(R) + collapse-threshold table | exists (data files) |
| F4 | Finder equivalence: Julia vs Fortran MF at matched cellsize (0.5%) | exists (2026-06-14) |
| F5 | Mass function vs Websky + Tinker (with/without AM) | exists |
| F6 | dN/dz octant vs Websky | exists |
| F7 | ξ(r) + b(M) + v12 panels | exists (compare_xi.jl etc.) |
| F8 | Velocity-coherence study: kSZ C_ℓ at cf=4/16/32 + κ control | exists (jobs 4321559/4321879) |
| F9 | κ composite C_ℓ vs kap_lt4.5 (+components) | exists (job 4321268) |
| F10 | kSZ composite vs ksz.fits (+hydro band overlay) | exists; hydro-template overlay to add |
| F11 | tSZ + CIB panels | exists |
| F12 | ISW same-mask C_ℓ (5 decades) | exists (job 4387665) |
| F13 | Map gallery (κ/kSZ/tSZ/CIB/ISW cutouts) | from production campaign |
| F14 | Performance/scaling curves | **new benchmark runs needed** |
| T1 | Validation scorecard (the TIER_B table, prettified) | exists |
| T2 | Websky-6144 configuration | exists |
| T3 | Timing/resource table per product | needs benchmarks |

## 3. Gaps to close before/while writing (work items)

1. **Frozen production campaign** (the paper's dataset): 8 octants, N=6144,
   box=5236 Mpc/h, cf=32, z_max=4.5 (not 4.6 — removes the known +10% low-ℓ kSZ
   residual), new-convention catalogs + full map suite at Nside 4096, one frozen
   config committed to the repo. ~1 day of Killarney queue time.
2. **Benchmarks** (F14/T3): strong scaling 1→4 L40S per octant; H100 single-node
   number; catalog vs fieldmap timing; peak memory; a Fortran reference timing
   (CITA numbers from NOTES.md or a fresh Sunnyvale run) — currently the weakest
   material.
3. **Convergence appendices**: nbuff sensitivity, filter-count sensitivity,
   cellsize study (partially exists via the factor-h saga — re-frame as a
   controlled convergence test, not the bug story).
4. **Hydro-template overlay for F10** (Shaw+, Battaglia+ kSZ templates digitized)
   — also settles the remaining physics question (1-halo/field double-count) and
   strengthens the neutral mid-ℓ narrative.
5. **Neutral-framing rewrite** of everything imported from the forensic MDs —
   scrub "bug/defect/wrong" language regarding reference codes; keep our own
   fixed defects (they're honest methods content: closure trap, D/a, 2LPT sign,
   Lagrangian deposit) framed as implementation validation.
6. **Code-release hygiene**: README/docs pass, example small-box config that runs
   on one GPU in minutes, CI status, Zenodo DOI, LICENSE check; decide whether
   the XGPaint painters live in the fork long-term or get a small standalone pkg.
7. **Decide**: does `fsc_of_z`/lightcone-threshold treatment match what the
   production campaign uses? (Old 2026-04 finding predates current validated
   runs — verify current code path once, cite in §5.)

## 4. Writing infrastructure

- `paper/` directory in PeakPatch.jl repo (or private overleaf-synced repo if
  collaborators join): MNRAS template, `figures/` generated ONLY by scripts in
  `paper/figscripts/` reading from the frozen campaign outputs — every figure
  reproducible from one `make figures`.
- Reuse the validation scripts (`validation/websky_6144/compare_*.jl`) as
  figscript backbones; port cap/pseudo-C_ℓ estimator into a small shared module.
- BibTeX: BM96, Stein+ 2019/2020, Websky2, MUSIC (Hahn & Abel 2011), CUDA.jl,
  Battaglia 2012/2016, XGPaint, Tinker 2008/2010, HalfDome, Agora, 2LPT refs.

## 5. Suggested phasing

- **P1 (now)**: skeleton + §2–6 prose (implementation is fully known, no new
  data needed); F1/F3 schematics; verify item 3.7.
- **P2**: frozen production campaign + benchmarks (queue-bound; launches early,
  writes §9 later).
- **P3**: validation sections §7–8 — mostly transcription of TIER_B_SUMMARY +
  existing figures re-rendered publication-grade; neutral-framing pass.
- **P4**: hydro overlay, convergence appendices, gallery, polish, internal
  review; decide author list & acknowledgements (Killarney/SciNet, CITA).

Dependencies: P3 figures can use existing octant-000 data even before P2
finishes; only F13/T3 and final numbers wait on the frozen campaign.
