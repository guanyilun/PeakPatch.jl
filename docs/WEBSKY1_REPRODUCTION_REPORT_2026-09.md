# PeakPatch.jl reproduction of WebSky 1.0: fidelity and performance

*Status as of 2026-09-27. Branch `websky2-resolution`. Audience: internal (PI) first; after PI review, Bond, Carlson, Stein et al.; input to the methods paper (`paper/`).*

> **Distribution gate: internal draft.** §6.2 and the reference-side parts of §4 (the tSZ-painter, kSZ-construction and CIB Table-1 checks) describe observed differences with the Websky reference code and maps. They must **not** go to the Websky authors, or anyone outside the group, until the PI has personally validated them and explicitly approved sharing. An external version should drop §6.2 and the reference-side remarks in §4.

All lengths are in Mpc/h, and masses in M☉/h unless stated. The Fortran peak-patch code and the public Websky products use Mpc. Conversion: `boxsize[Mpc/h] = L[Mpc]·h`. For Websky, 7700 Mpc × 0.68 = 5236 Mpc/h (`CONVENTIONS.md`).

**Result tiers used below.**
- **Current:** the production-v2 campaign (8 octants, nbuff=16) and the nbuff=25 "v3test" (oct000 only).
- **Superseded:** the frozen 2026-07 campaign, which had seam gaps; runs from before the 2026-06 fixes; and the Tier-A single-octant validation catalog (`oct000_finecell_AM`), which is used because nothing newer exists for those statistics.
- **Do not quote:** the "19 min/octant" and "~0.9M halos/octant, 120× fewer than Websky" figures. They come from the 2026-04 configuration, before the bug fixes (`validation/performance/PERFORMANCE_2026-09-26.md` §1).

**Framing.** Where we list differences from the reference Fortran code or the released maps, these are *observed differences with evidence*. The PI has not yet personally validated them, and none has been communicated upstream. Sending this draft to Websky authors would count as communicating them, so the distribution gate above applies.

**Glossary.**
- **Tier-A:** catalog-level statistics (MF, dN/dz, σ_vr, ξ, b(M), v12), measured in 2026-06/07 on the single validation octant `oct000_finecell_AM` against Websky's 10°×10° halo patch.
- **Tier-B:** map-level (painted) statistics from the same 2026-07 period, on oct000 caps or octants.
- **Plan codes** (A0–A8, B1–B6, Q1–Q3): items in `docs/paper_comparison_plan_2026-09.md`. For example, A2 = rerun Tier-A/B on current catalogs; A3 = error model and pass criteria; A4 = theory anchors independent of Websky; A6 = full-sky assembly and replication checks; A8 = convergence appendix; B5 = kSZ hydro-template overlay.
- **kSZ halo constructions** (`KSZ_COMPOSITE_2026-07-19.md`):
  - **W:** uncompensated per-halo Battaglia τ painting;
  - **We:** the Websky-literal variant, with a Δ=3-mean compensation sphere and gas mass f_b·m_h;
  - **Wc:** a zero-net compensated variant, a full-τ sphere with R=4 r_vir.
- **cf:** the coarse_factor of the multi-resolution split (coarse grid M = N/cf).
- **octZYX:** the octant label; the three bits give the observer sign per axis.

**Sources.** Facts sourced only to `memory/*.md` notes are not in the repo and are marked *(memory only; unverified in repo)*.

**Ranges.** Ratio ranges in §1, §3.2 and §4 were re-read from `validation/paper/results/*.txt`. Where they differ from the `V2_RESULTS_2026-09-27.md` scorecard, the results files take precedence. The known differences are:

| Quantity | V2_RESULTS | results files (used here) |
|---|---|---|
| tSZ, ℓ 270–1400 | 0.82–0.88 | 0.804–0.883 |
| ISW, ℓ 27–300 | 0.85–1.05 | 0.67–1.05 (0.671 at ℓ=30) |
| κ×CIB | 1.25–1.33 | 1.21–1.39 |
| κ×y | 1.03–1.10 | 0.98–1.15 |
| CIB Poisson | per-ν values differ | results-file values |

V2_RESULTS itself has not been edited.

---

## 1 Summary

**Verdict.**
- **Halo catalog.**
  - The statistics that AM does not enforce (clustering ξ, bias b(M), σ_vr, v12) agree with Websky at 1–6% over the stated ranges (ξ 3–15 Mpc/h, v12 5–30 Mpc/h, b(M) 7.9e12–5e13).
  - Scope: this was measured on **one superseded validation octant** (`oct000_finecell_AM`) against Websky's 10°×10° patch, about 1/51 of an octant, and has not been repeated on v2/v3.
  - Known out-of-band residuals:
    - ξ at 1–3 Mpc/h is 1.2–1.4× (unexplained);
    - v12 falls to 0.90 at 31 Mpc/h and 0.25 at 58 Mpc/h;
    - b(M) is 0.952 at 1.4e14;
    - both catalogs sit 5–13% below Tinker10 in b(M);
    - the top N(>M) bin (3e14) is 0.93.
  - Post-AM N(>M) agreement (1.00–1.02) is **largely enforced by AM**, because both catalogs are abundance-matched to the same Tinker08 M200m function. It tests the AM implementation and volume bookkeeping, not the finder.
  - The independent finder evidence is the Fortran/Julia kept-peak test (limited scope, below). The raw pre-AM MF against Tinker per z-bin (plan A4) has not been done (§7).
- **Seams:** fixed in v2. The signal profile across the octant planes tracks Websky's. This is our own fix, not a reproduced Websky statistic.
- **CIB decoherence:** agrees with the released maps (post-2022 CIB update) to ≤0.025 for all 10 pairs, with no σ attached (143×857 differs by the full 0.025). This uses a same-family CIB model (XGPaint `CIB_Planck2013`), so it is not an independent check. The printed Stein+ Table 1 differs from both (§6.2).
- **ISW at low ℓ:** consistent within large cosmic variance, using an unofficial σ computed for this report (§4).
- **Not yet reproduced to the pre-declared criteria:**
  - κ (6–13% high against `kap_lt4.5`);
  - full-sky tSZ, which is low at nbuff=16 for an identified cause; the fix has been demonstrated on one octant only;
  - CIB, which differs by a mostly painter-driven amplitude: mean 1.23, and 1.18 on Websky's own halos through the same model. There is also a frequency-dependent residual at 857 GHz and in the Poisson regime, where the flux cut has not yet been applied.
- **Differs by construction:** kSZ halo power at mid-ℓ, where the result depends on how the halos are painted.
- **Cost:**
  - The full sky takes about 115 L40S-GPU-h, with no stored initial conditions (ICs).
  - Peak memory is ~30× lower **per octant run** (0.24 TB, against 7.67 TB for Websky's full-box run). For the full sky with 8 concurrent octants it is ≈1.9 TB, ~4× lower in aggregate.
  - The v2 cost is for a mass-capped (nbuff=16) catalog. This is **not** a like-for-like comparison (§5).

**Headline results.**
- **Phase space: the AM-insensitive statistics** (Tier-A octant, superseded catalog; see the §3 caveat):
  - σ_vr agrees to about 1–2% (*memory only*: `memory/catalog_velocity_units.md`). The repo measurement in `BIAS_PAIRWISE_V12_2026-07-16.md:32` is 267.5 vs 248.0 km/s, a ratio of 1.079, in a small 5° cap;
  - ξ bias ratio is 1.06 (3–15 Mpc/h);
  - b(M) ratio is 0.95–1.02 in 4 bins;
  - v12 ratio averages 0.997 over 5–30 Mpc/h.

  Main source: `validation/websky_6144/BIAS_PAIRWISE_V12_2026-07-16.md`. None of these has yet been recomputed on v2.
- **Finder equivalence (limited scope).**
  - Kept-peak densities are 3.925e-3 (Julia) and 3.946e-3 (Fortran) per (Mpc/h)³, where both codes use the same numeric cellsize, 1.2533. Fortran/Julia = 1.005 (`validation/websky_6144/FORTRAN_COMPARISON_2026-06-14.md`).
  - Scope: a 320.8 Mpc/h box; a z=0 snapshot (no lightcone); the old 1.25 Mpc/h cell; pre-v2 code. It is statistical (different seeds) and compares kept-peak density only, not the mass function or matched halos.
- **Mass function after AM** (Tier-A octant; largely enforced by AM, see the Verdict). Ours/Websky N(>M) is 1.00–1.02 from 1.23e12 to 1e14, and N(>1.69e12) = 5.08e7 per octant against Websky's ≈5.09e7 (`validation/websky_6144/LOW_MASS_COMPLETENESS_2026-06-15.md`).
- **Seams fixed in v2.** κ-field and ISW maps have 0 exactly-zero pixels, against 190k and 1.5M in the frozen campaign. The signal profile across the octant planes tracks Websky (`validation/paper/V2_RESULTS_2026-09-27.md`).
- **tSZ.**
  - v2 full sky gives 0.80–0.88× Websky at ℓ=270–1400.
  - Cause: the nbuff=16 shell-search cap (R_TH ≤ 12.8 Mpc/h, raw M ≤ 7.4e14). This leaves N(>1e15, z<1) = 0 per octant, against 51.5 for Tinker.
  - With nbuff=25 on **oct000 only**: N(>1e15) = 52 and C_ℓ ratio 1.06–1.21 at ℓ=300–1400.
- **CIB.**
  - Mean intensity is 1.23× the released maps at every frequency. Websky's own halos through the same XGPaint model give 1.18, so most of the offset comes from the painter/model and not from the catalog.
  - Clustered C_ℓ is 1.44–1.46 at 100–353 GHz, falling to 1.40 (545 GHz) and 1.26 (857 GHz). Poisson-regime ratios are 1.22→0.70 across frequency. Both are frequency-dependent rather than a single (1.23)² factor.
  - Decoherence agrees with the released maps to ≤0.025 for all 10 frequency pairs (no σ).
- **Performance** (v2, one 4×L40S node per octant; nbuff=16, mass-capped catalog):
  - catalog+AM takes 2.97–3.21 h per octant;
  - host memory is about 135 GiB. GPU memory is about 22 GiB per GPU in the n=414 benchmark; it was not measured at n=416 (v2);
  - peak memory is ~30× lower per octant run (0.24 TB vs 7.67 TB for the full-box Websky run). The full sky at 8 concurrent octants is ≈1.9 TB, ~4× lower in aggregate;
  - IC storage is 0 B, against 5.9 TB for Websky.

---

## 2 What we implemented, and how it differs from the Fortran/Websky production

### 2.1 Production-v2 configuration

Sources: `validation/websky_6144/production_v2/README.md` and `config_v2_oct000.toml`.

| Item | Value |
|---|---|
| Seed, N, box, cell | 12345; 6144; 5236 Mpc/h (= 7700 Mpc); 0.852213 Mpc/h |
| Tiling | ntile=16, tile mesh n=416, nbuff=16, nsub=384, `periodic_cores=true` (N = nsub·ntile). `[grid] boxsize` is the **per-tile** box, 354.52 Mpc/h |
| Coarse grid | coarse_factor=32 (M=512, block 12), `coarse_compensation=true` |
| Lightcone | ievol=1, z_out=0, z_max=4.5; observer at ±2618 Mpc/h per axis (the core corner, set by the octZYX bits) |
| Physics | ilpt=2, ioutshear=1, wsmooth=1, rmax2rs=0.0 |
| Cosmology | Ωm=0.31, Ωb=0.049, ΩΛ=0.69, h=0.68; P(k) from `pk_websky.dat` (CAMB, σ8=0.8100) |
| Filters | `filters_websky_finecell.dat`: 23 top-hats, Rf = 1.406–30.44 Mpc/h, ratio 1.15 |
| Collapse table | `HomelTab_websky.dat`, 50×20×20 |

### 2.2 Differences from the Fortran code and the Websky production

| Aspect | Fortran / Websky (Stein+2020) | PeakPatch.jl v2 | Source |
|---|---|---|---|
| Units | Mpc | Mpc/h everywhere, because chi(z) uses H0=100h | `CONVENTIONS.md` |
| White noise | 48-bit LCG, filled sequentially per MPI slab; the realization depends on the task count | Counter-based Threefry2x (20 rounds) with key (seed,0), plus Box-Muller. Each cell depends only on (seed, global index). A bit-exact LCG port is also available | `src/InitialConditions/RandomField.jl:36-103`; `LCG.jl` |
| Initial conditions | Global slab-decomposed FFTW3, 7 values per site, 5.9 TB stored | No global fine FFT. The coarse M³ periodic FFT is spliced to a per-tile isolated residual convolution with Catmull-Rom interpolation (MUSIC/Hoffman-Ribak style); D/T compensation in v2; nothing stored | `docs/multi_resolution_fft.md:60-104`; `src/MultiResolution.jl:584-735` |
| Buffer | ≈ the largest Lagrangian halo radius (~40 Mpc) | nbuff=16 cells = 13.6 Mpc/h. The shell search is capped at nhunt = min(nbuff−1, ⌊1.75 Rf_max/a⌋), the same formula as `hpkvd.f90:207-209` | `src/MultiResolution.jl:630` |
| Filter bank | Paper: Rf,min = 2 a_latt; Rf,max = 36 Mpc (24.5 Mpc/h); count and spacing not stated. `filter_gen.py` default: 1.65 a_latt, spacing 1.15 (the same as ours). Production bank not located; the evidence conflicts (paper 2.0; the dense-bank residual suggests ≈1.2) | Rf,min = 1.65 a_latt, Rf,max = 30.44 Mpc/h, 23 filters, spacing 1.15 | `data/filters_websky_finecell.dat`; `2001.08787.txt:300-303`; `LOW_MASS_COMPLETENESS_2026-06-15.md:66-80` |
| Peak threshold on the lightcone | Constant fcrit = fsc_of_z(0) | Same. Shell analysis uses the per-peak z_pk (1+z_pk, fsc_of_z(z_pk)) | `src/MultiResolution.jl:934-937,1017-1045` |
| Ω_m(a) in the 2LPT coefficient | Standard (Ωm/a³)/(Ωm/a³+ΩΛ) (`hpkvd.f90:762,875`) | Used a³ where the standard form has a⁻³; **fixed 2026-09-27** via `Cosmology.omega_m_a` (§6.1). All v2/v3test catalogs predate the fix | `src/Cosmology/Cosmology.jl` |
| Merge | Exclusion with a fixed NC=256 hash; volume reduction commented out | Same survivors, with a data-sized hash (`clamp(round(cbrt(n)),16,1024)`), a total order (−R_TH,x,y,z), and volume reduction disabled | `src/Merger/Merger.jl:32-80`; `Exclusion.jl:22,57-64` |
| Eulerian positions and velocities | merge_pkvd | Port: x = q+ψ1+ψ2, v = a·100·E·f·(ψ1+2ψ2). Velocity epoch from one chi→z evaluation at the Eulerian distance, with no iteration | `src/Merger/Merger.jl:106-145` |
| AM | Tinker M200m, Δz=0.1 bins, bilinear interpolation; only halos with >10 particles before AM | Tinker 2008 Δ=200m, rank inversion; 46 z-bins (Δz≈0.098), 10k mass bins over 5e11–1e16. **No particle cut.** Top-halo fix applied | `src/AbundanceMatch/AbundanceMatch.jl`; `2001.08787.txt:1044-1052` |
| Catalog format | 10 floats, ~9e8 halos, 33 GB | 33 Float32 (extended pksc), 198.4M halos/octant, 26.2 GB/octant | `src/Merger/Merger.jl:84-162` |
| Field maps (κ, τ, kSZ, ISW) | pks2map / field painting | Separate pass that regenerates the tile fields and paints every 2LPT-displaced lattice cell on the GPU at Nside 4096; ISW is deposited at the Lagrangian position | `src/FieldMap.jl`; `ext/CUDAExt.jl:3441-3515` |
| Halo painting | pks2map (Fortran) | Our CPU painter `paint_octant.jl` / `halo_profiles.jl`: κ is a truncated NFW with a Δ=3 compensation sphere; tSZ is Battaglia12 AGN pressure; kSZ is Battaglia gas with M200c>1e13 and r200c>0.5′, in W and Wc (zero-net) forms. XGPaint is not used for these | `validation/websky_6144/production/paint_octant.jl`, `halo_profiles.jl` |
| CIB | Websky CIB model | XGPaint `CIB_Planck2013`, Websky completeness cut (wcut), serial pixel accumulation | `validation/websky_6144/production/paint_cib.jl` |

### 2.3 Accuracy of the multi-resolution split (field level, N=480, block 12)

Source: `validation/tiling/TILING_INVARIANCE_2026-09-26.md` and `validation/tiling/SPLICE_COMPENSATION_2026-09-26.md`.

| Quantity | Uncompensated (frozen campaign) | Compensated (v2) |
|---|---|---|
| δ rms error vs exact global FFT | 5.10–5.13% of σ (r = 0.99875) | 4.50% (r = 1.0000 below k_Nc) |
| ψx rms error | 15.2–18.0% | 15.2% (near-Nyquist image decoherence) |
| P_split/P_exact, k = 0.05–0.19 h/Mpc | 1.021–1.115 measured (1.020–1.111 predicted by the D(T−D) model) | 0.999–1.003 (k = 0.05–0.29) |
| σ(R) stitched/global | up to 1.035 at R≈11 Mpc/h | 0.993–1.001 |
| δ error vs distance to tile face | flat 4.8–5.6% (no seam) | — |

Further tiling results:
- **Bitwise identical across ntile = 1/2/4/8:** the white noise, the coarse fields and the residual noise.
- **Not bitwise identical:** the fine δ, which differs by 0.28–0.62% of σ_δ.
- **Catalog tiling test** (job 5695201): merged N ratio 0.9993–0.9994; N(>M) within ±0.2% below 1e13 and ±0.5% up to 1e14; median ln(M ratio) 0 ± 0.002; per-halo ln M rms 2.6–9.0%.

**Velocity-coherence rule** (`paper/sections/multires.tex:58-93`):
- The coarse k_Nyq = πM/L must be ≳ 0.3 h/Mpc for products weighted by velocity or potential, i.e. cf≈32.
- cf=4 is enough for catalogs.
- cf32 adds only about 0.5 h to a catalog run (6h39m vs 7h13m, jobs 3968163/4387351; `validation/performance/PERFORMANCE_2026-09-26.md` §1, consistent with sacct 06:38:45 and 07:13:07). Both are frozen-era job walls dominated by the NC=256 merge, so they say nothing about the v2 cost.

---

## 3 Catalog validation

**Caveat that applies to this whole section.**
- All Tier-A statistics below come from one octant, `oct000_finecell_AM`, compared with Websky's public 10°×10° patch (`halos_10x10.pksc`).
- That catalog has several superseded features:
  - z_max=4.6;
  - AM from before the top-halo fix;
  - the order-dependent merge;
  - the seam-gapped nsub=382 layout;
  - the nbuff=16 mass cap;
  - cf=4.
- The raw catalog (job 3968163, 2026-06-16) also pre-dates the Eulerian/velocity and 2LPT-sign fixes (95f04a1, 2026-06-29; 4ae7837, 2026-07-16, GPU path). Its Eulerian positions and velocities were reconstructed at read time from the old-convention fields (`BIAS_PAIRWISE_V12_2026-07-16.md:7-9`). This matters directly for the σ_vr, ξ, b(M) and v12 results.
- Post-AM N(>M) agreement is largely enforced by AM, because both catalogs are matched to Tinker08 M200m, so it is not independent evidence. The AM-insensitive statistics are ξ, b(M), σ_vr and v12, and the Fortran finder test. The raw pre-AM MF against Tinker per z-bin (plan A4) has not been done. The only raw-vs-Tinker numbers are the z<1 tail counts in §3.2.
- **No Tier-A statistic has been rerun on v2 or v3test** (plan item A2, `docs/paper_comparison_plan_2026-09.md`).
- The full Websky `halos.pksc` is not available locally, so a full-sky catalog comparison has not been done.
- Stein+2020 itself shows no halo mass function, dN/dz or clustering plot, so these statistics are additional validation.

### 3.1 Scorecard (ours / Websky unless stated)

| Statistic | Result | Source |
|---|---|---|
| Finder, matched cell 1.2533 Mpc/h, z=0 | Kept-peak density 3.946e-3 (Fortran) vs 3.925e-3 (Julia) per (Mpc/h)³, with the same numeric cellsize in both codes: Fortran/Julia 1.005. Scope: 320.8 Mpc/h box, z=0 snapshot (no lightcone), pre-v2 code, different seeds (statistical), kept-peak density only (no MF or matched halos) | `FORTRAN_COMPARISON_2026-06-14.md:50-60` (labelled "/Mpc³" there) |
| Lightcone vs snapshot (z<0.19) | 0.978; GPU = CPU (raw 69008 = 69008; merged 25898 vs 25897) | same |
| N(>M) after AM, 1.23e12 → 3e14 (largely enforced by AM) | 1.01, 1.01, 1.00, 1.00, 1.00, 1.02, 1.01, 0.93 (last bin: rare-cluster noise) | `LOW_MASS_COMPLETENESS_2026-06-15.md:3-10` (job 4033130) |
| N(>1.69e12) per octant (largely enforced by AM) | 5.08e7 vs ≈5.09e7 | same |
| dN/dz | **Not quantitatively validated** for finecell_AM. "Both peak at z≈1.4" is recorded only in a memory note. The only recorded per-z ratios (M>3e12) are for the superseded coarse-cell AM catalog: ≈0.85→1.0 at z<2 and 1.2–1.36 at z>2.5 | *memory only*: `memory/websky_highz_lightcone_deficit.md`; `LOW_MASS_COMPLETENESS_2026-06-15.md:41-43` (coarse catalog) |
| σ_vr | ~1–2%, flat in mass (*memory only*). In the repo, the small matched cap gives 267.5 vs 248.0 km/s, a ratio of 1.079 | `memory/catalog_velocity_units.md` (`apply_eulerian_velocity.jl` output); `BIAS_PAIRWISE_V12_2026-07-16.md:32` |
| ξ(r), Landy-Szalay with matched cap | bias ratio 1.06 (3–15 Mpc/h); 1–3 Mpc/h ~1.2–1.4× (unexplained; *memory only*) | `BIAS_PAIRWISE_V12_2026-07-16.md:63`; `memory/websky_comprehensive_validation.md` |
| b(M) at 7.9e12 / 2.0e13 / 5.0e13 / 1.4e14 | 1.024 / 1.006 / 0.980 / 0.952. Both catalogs sit 5–13% below Tinker10 | `BIAS_PAIRWISE_V12_2026-07-16.md` §1 |
| v12(r), M>1e13 | mean 0.997 at 5–30 Mpc/h; ours decays faster beyond 30 (0.90 at 31, 0.25 at 58 Mpc/h), not investigated | same §2 |

### 3.2 Current production catalogs (v2, and v3test oct000)

Source: `validation/paper/V2_RESULTS_2026-09-27.md`.

| Quantity | frozen (superseded) | v2 (8 octants) | v3test oct000 (nbuff=25) | Tinker |
|---|---|---|---|---|
| Halos / octant (post-merge) | 195.5M | 198.4M (+1.5% = recovered seam volume) | 198.4M | — |
| Raw N(>1.7e12) | 4.03e7 | 4.09e7 | 4.091e7 | — |
| AM N(>1.7e12) | 5.076–5.078e7 | 5.077e7 | — | — |
| AM N(>1e14) | 6.251–6.255e4 | 6.25e4 | 6.248e4 | — |
| z<1: AM N(>3e14) | — | 3615 | 3610 | 3657 |
| z<1: AM N(>5e14) | — | 750 | 746 | 735 |
| z<1: AM N(>1e15) | — | 0 | 52 | 51.5 |
| Raw / AM max M | — | 7.4e14 / 8.5e14 | 2.6e15 / 1.7e15 | — |

Notes on v3test:
- Raw R_TH tapers smoothly, with no pile-up: the maximum is 19.27 Mpc/h, below the 20.4 cap.
- AM Σ(M/1e14)^(10/3) at z<0.5 is ×1.58 over v2.
- nbuff=25 gives nhunt=24, the largest the GPU shell kernel allows (482 shells against `_MAX_SHELLS_GPU` = 512).
- The v3test section of V2_RESULTS is uncommitted working-tree text.

AM top-halo fix (commit 9f92778), measured on frozen oct000:
- 319 of 195.5M halos change, with median ΔM/M −8.4%;
- ΣM^(10/3) changes by −4.6% overall and −7.9% at z<0.5;
- ΣM^(5/3) changes by −0.2%.

Source: `validation/websky_6144/AM_TOPHALO_FIX_2026-09-25.md`.

### 3.3 RNG and geometry

- **Threefry vs Xoshiro: no detectable generator difference** at ℓ=2–80 with 28 seeds at full N=6144.
  - Exact-permutation family-wise p = 0.70.
  - The one borderline statistic did not replicate out of sample (p = 0.59).
  - Source: `validation/rng/RNG_LARGE_SCALE_TEST_2026-09-23.md`.
- **Box replication is untested.** The full sky is 8 views of one 5236 Mpc/h periodic box, so structure repeats beyond χ≈2618 Mpc/h (z≳1.2). Websky has the same geometry (plan A6).

---

## 4 Map validation

**Setup.**
- Estimator A3: full-sky C_ℓ of the 8-octant v2 maps against the released Websky v0.0 maps, with no mask, pixel windows divided out, and Δℓ/ℓ≈0.1.
- Analysis jobs: 5703683 (ours vs Websky) and 5703908 (vs theory).
- Reference files are in `/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/`.
- Pass criterion, declared in advance (`docs/paper_comparison_plan_2026-09.md:77-105`): the band-averaged |r−1| < 3σ over the declared ℓ range, and no single band above 4σ.
- Construction-dependent differences (kSZ halo term at mid-ℓ, κ at ℓ≳1700) are not scored pass/fail.

> **Caveats on the "status" column.**
> - The σ values were computed for this report as an inverse-variance mean of the `sigma_ratio` column in `validation/paper/results/auto_v2_*.txt`. No repo document states formal verdicts; they should be re-derived with the paper tooling.
> - `compare_auto.jl:28-34` uses max(Gaussian, jackknife) at all ℓ, which is not the plan's Gaussian-only rule at ℓ≲30.
> - The theory and cross-spectrum files carry no σ, so those tests are not scored.

| Product (declared ℓ range) | ours / Websky | ours / theory | Status vs pre-declared criterion |
|---|---|---|---|
| κ, z<4.5 vs `kap_lt4.5` (30–3000) | 1.06–1.08 (ℓ 270–430); 1.10–1.12 (850–1400); 1.13 (2600–3000); 1.02 (8000) | Halofit: 1.00–1.056 (ℓ 150–1650), 1.136 (2894). Field κ vs linear Limber: 0.99–1.04 (152–1000). Websky vs Halofit: 0.935–0.994 (150–1650), so the offset is shared between the two maps | **Fail** at ℓ 30–1700: 1.103 ± 0.0014 (19/43 bands >4σ); ℓ≳1700 construction-dependent |
| tSZ, v2 full sky (100–5000) | mean y 1.021; C_ℓ 0.96–1.29 (40–150), 0.80–0.88 (270–1400), 0.91–0.98 (2600–4100) | — | **Fail**: 0.924 ± 0.0013. Cause identified (nbuff cap, §3.2) |
| tSZ, v3test **oct000 only**, masked | mean y 1.035; C_ℓ 1.06–1.21 (300–1400), 1.03–1.05 (2000–3000). v2 on the same mask: 0.90–0.94 | — | Single octant; not scored. Full-sky v3 pending |
| kSZ total, field + Wc halos (500–8000) | 1.17–1.40 (86–400); ~1.2 (1000); ~1.0 (≥3500) | — | Construction-dependent: 1.122 ± 0.003 |
| kSZ field (30–500) | — | vs exact-LOS Doppler + linear OV: 1.13–1.21 (71–203), 1.05–1.11 (245–360), 1.00–1.02 (434–933) | Not scored (no σ). ℓ 70–200 excess **open** |
| ISW (10–250) | 0.67–1.05 (27–300); ≥1.6 at ℓ≥434 (painting noise floor, outside the declared range) | — | **Pass (weak)**: 0.907 ± 0.075 (1.2σ, max band 1.6σ), with the unofficial σ computed for this report; still passes with Gaussian-only errors at ℓ<30. The range is dominated by ℓ<30, with large two-realization cosmic variance and shared octant structure: individual bands at ℓ<15 run 0.49–1.34. Ours vs linear-theory ISW would be a stronger anchor (not done) |
| CIB mean intensity | 1.23 at all ν. Baseline: Websky's halos through the same XGPaint model give 1.18, so the offset is mostly from the painter/model | — | Normalization **open** |
| CIB clustered (100–3000) | 1.46/1.46/1.45/1.44/1.40/1.26 at 100/143/217/353/545/857 GHz | — | **Fail** on amplitude: ≈1.23² at 100–353 GHz, plus a frequency-dependent fall to 1.26 at 857 GHz |
| CIB Poisson (3000–8000) | 1.22/1.20/1.15/1.04/0.89/0.70 | — | No band >4σ, but the 400 mJy flux cut is applied to neither map. **Open** |
| CIB decoherence (Websky Table 1, 150<ℓ<1000) | All 10 pairs within ≤0.025 of the released maps (post-2022 CIB update), e.g. 545×857 0.956 vs 0.951; 143×857 0.815 vs 0.790 sits at the limit | — | Agrees with the released maps (no σ). Same-family CIB model, so not an independent check. Stein+ printed Table 1 differs (§6.2) |
| κ × CIB (100–2000) | 1.21–1.39 (≈√(1.07·1.5)) | — | Not scored; tracks the CIB normalization |
| y × CIB | 1.19–1.37 (ℓ 190–2940) | — | Not scored; tracks the CIB normalization |
| κ × y | 0.98–1.15 (ℓ 108–2940) | — | Not scored |
| Seams | κ-field and ISW: 0 exactly-zero pixels (frozen: 190k / 1.5M); CIB 120k = empty-pixel floor (frozen: 1.9M). y across the planes 1.22/1.16/1.14 vs Websky 1.18/1.16/1.11 | — | Fixed (our own bug fix); the octant-plane profile tracks Websky's |

Sources: `validation/paper/V2_RESULTS_2026-09-27.md`, `validation/paper/results/{auto,theory,cross}_v2_*.txt`, and `validation/paper_theory/THEORY_ANCHORS_2026-09-26.md`.

**Supporting checks.**
- **tSZ painter on Websky's own halos.** We painted `halos_10x10.pksc` with our painter and compared it with the released `tsz_2048` on a 4.5° disc.
  - C_ℓ ratio: 1.21/1.16/1.11/1.06/1.04/1.06/1.13/1.25 at ℓ=204–3666.
  - Cross-correlation ≥0.992; mean y 1.006.
  - Caveat: this is a single 4.5° disc. tSZ there is dominated by a handful of clusters, so the variance is large and non-Gaussian, and no σ is attached.
  - Within that caveat, the v3test excess of 1.06–1.21 is of the same size as the difference between our painter and the reference on identical input.
  - V2 offers a hypothesis for the painter difference: the reference table is truncated at 4 Mpc transverse. This has not been demonstrated.
- **Reference map versions.**
  - The release `UPDATES` file (mocks.cita.utoronto.ca/data/websky/v0.0/UPDATES, quoted in `V2_RESULTS_2026-09-27.md:87-95`) records:
    - 17-FEB-2022: a ~7% tSZ profile-normalization fix;
    - 18-FEB-2022: a CIB satellite update;
    - 08-APR-2022: a fix for the high-mass profile centre.
  - V2_RESULTS states that our `tsz_2048` reference is the post-fix version.
  - Checksums of the downloaded files are in `websky_ref/md5sums.txt`. The versions are otherwise established only from download dates.
  - This context applies to every tSZ and CIB comparison with the reference.
- **Inference for full-sky v3 (untested).** oct000 is the highest octant: its masked v2 ratio (0.90–0.94) sits above the full-sky value (0.80–0.88). Applying the oct000 v3/v2 factor to the full-sky v2 ratio suggests roughly 0.95–1.1. This has not been tested.
- **kSZ construction.**
  - The catalog input is cross-validated: N(M_h>1e13, z) ratio 1.00.
  - A faithful Battaglia painting of *either* catalog gives a halo term ~5–10× the released decomposition at ℓ~500–1000.
  - Tier-B band means vs `ksz.fits`: W 1.84, We 1.58, Wc 1.38.
  - Shuffling velocities leaves Wc almost unchanged (1.383→1.373), so the mid-ℓ excess is 1-halo shot power.
  - The mid-ℓ halo power differs by construction. The hydro points in Stein+ Fig 6 are ℓ=3000 values (Park+18 Table 1, CSF-scaled; `2001.08787.txt:1258-1267`), so they do not constrain ℓ≈900 directly. At ℓ≈2700–3500 our total is 1.02–1.07× `ksz.fits`. The hydro-template overlay (B5, Fig 6) is pending and is the intended arbiter.
  - Source: `validation/websky_6144/KSZ_COMPOSITE_2026-07-19.md`.
- **Field-kSZ coarse-factor convergence** (Tier-B oct000, vs `ksz.fits`, ℓ=118–449): cf4 1.06–3.16; cf16 1.01–1.44; cf32 0.97–1.22.
- **Estimator transfer** (8 Gaussian skies, masked per-octant pseudo-C_ℓ): 1.00 ± 0.01–0.03 at ℓ≥70 (`THEORY_ANCHORS_2026-09-26.md:57-72`).
- **Not done yet:**
  - the observational and hydro-template overlays for Websky Figs 5–7 and 10–11;
  - the lensed-CMB check (Fig 9).

---

## 5 Performance and resources

> **The comparison is not like-for-like.**
> - Hardware: Websky ran on 2019 Skylake CPU nodes; we run on 2024-era L40S GPU nodes.
> - Catalog size: we keep all ~198.4M halos/octant (≈1.59e9 full sky), against Websky's ~9e8 after its cut.
> - Duplicated volume: each octant processes 2413 of 4096 tiles, so the full sky processes ~4.71 box volumes of fine tiles. 2413 is confirmed for v2: the log `v2_cat_oct000_5696115.err` reports "Tile k/2413".
> - Scope: Websky's 3.84 h figure is "the run time of the simulation"; the paper gives no cost for field maps or painting.
> - Node counts assume 40 cores per Niagara node.
> - **The v2 cost is for a mass-capped catalog.** The nbuff=16 shell cap removes all M>7.5e14 halos, and Websky's ~40 Mpc buffer had no such cap. The like-for-like configuration is v3 (nbuff=25), whose cost comes from one octant and one run so far.
> - **The v3test/v2 timing mixes three changes:** nbuff 16→25, the shell early exit (off in v2, on in v3test), and a different node (kn020 vs kn110). The v2 cost with the early exit is unknown.
> - **Early-exit exactness** was checked bit-for-bit only on the small example config (~610k records, job 5708923) and in a unit test (job 5711281), not at production scale (`SHELL_EARLY_EXIT_2026-09-27.md:26-33`).
>
> The robust claims are zero IC storage, the time to solution on a single node, and the memory footprint *per octant run* (~30×). In aggregate over the full sky the memory saving is ~4×.

| Resource | Websky (Stein+2020 §4.1) | PeakPatch.jl v2 | Source |
|---|---|---|---|
| Hardware per run | 1128 Skylake cores (~28 nodes, assumed) | 1 node, 4×L40S, per octant | `2001.08787.txt:1033-1037` |
| Catalog wall | 3.84 h (full box) | 2.97–3.21 h per octant for catalog+AM (mean 3.04); 3.21 h full sky with 8 nodes in parallel; ~24.3 h sequentially on one node | sacct 5696115–5696143 |
| Catalog compute | 4336 core-h (≈108 node-h, assumed) | 24.32 node-h = 97.3 L40S-GPU-h allocated (58.7 GPU-h of summed stage time) | same |
| Full product set | not stated | 115.0 GPU-h ≈ 28.7 node-h, of which 97.3 catalog+AM, 10.44 field maps, 1.36 halo paint, 5.89 CIB | sacct |
| Peak memory | 7.67 TB (whole box) | ≈135 GiB host + 4×~22 GiB GPU (GPU figure from the n=414 benchmark; not measured at n=416) ≈ 0.24 TB **per node-job/octant**, ~30× less than Websky's 7.67 TB per full-box run. Aggregate for 8 concurrent octants ≈1.9 TB, ~4× less | `PERFORMANCE_2026-09-26.md` §4; sacct MaxRSS |
| Stored ICs | 5.9 TB | 0 B (regenerated from the seed) | same |
| Catalog on disk | 33 GB (10 floats, ~9e8 halos) | 26.2 GB/octant (33 floats); raw+AM × 8 = 390.3 GiB. An 11-float format would be ~8.6 GB/octant | `ls` on /project |

**v2 oct000 stage profile** (max per worker, 4 workers; log `v2_cat_oct000_5696115.err`):

| Stage | Time | Share |
|---|---|---|
| Shell analysis | 5308 s | 80.4% |
| Peak find | 427 s | 6.5% |
| Residual generation | 321 s | 4.9% |
| δ interpolation | 249 s | 3.8% |
| Stage total | 6607 s | — |

- The pipeline itself takes 121.5–126.3 min per octant.
- The block after the pipeline (merge, finalize, 26 GB write) took about 42–52 min, measured from file mtimes. It was not timed per stage.
- The earlier ~10–15 min merge estimate was a projection from a microbenchmark: on 9.81M halos, 346 s at NC=256 against 13 s with 5 Mpc/h cells.
- AM takes about 15 min.

**Scaling** (N=1560, 64 tiles; pre-v2 layout and merge; `PERFORMANCE_2026-09-26.md` §2):

| Hardware | Pipeline time | Speed-up | Efficiency |
|---|---|---|---|
| 1 L40S | 709.7 s | 1.00× | — |
| 2 L40S | 374.6 s | 1.89× | 95% |
| 4 L40S | 225.0 s | 3.15× | 79% (tail imbalance) |
| 1 H100 | 419.1 s | 1.69× one L40S | — |

- Production load balance is 99%.
- Mean GPU utilization was 19–34%, which leaves headroom.

**GPU memory vs tile size.**
- 21.2–23.7 GiB at n=414 and 44.4 GiB at n=576 (the ceiling).
- n=416 (v2) and n=434 (v3test) were not measured.

**nbuff=25 cost** (v3test oct000, job 5710456, same host RSS of 135.8 GiB):
- **Without** the shell early exit: about 20 s/tile, projecting to about 13 h per octant (job 5704175, cancelled).
- **With** the exact early exit (commit a887192):

| | v3test | v2 oct000 |
|---|---|---|
| Job wall | 2:38:03 | 2:58:43 |
| Halo-finding step | 2:22:41 | 2:43:20 |
| Pipeline | 100.0 min | 121.6 min |

  V2_RESULTS reads this as 2h22m vs 2h02m, but that compares the full step with the pipeline alone.
- This is one octant on a different node (kn020 vs kn110). The comparison mixes nbuff 16→25, early exit off→on, and the node change (see the caveats above).
- The v3test log is stamped with commit 1c11cb3, but the early-exit code was uncommitted when it ran. The evidence: the `v3test_cat_oct000_5710456.out` header shows commit 1c11cb3 and a start at 08:11 on 2026-09-27, while `git log` shows a887192 authored at 11:56 that day.

**Superseded (frozen campaign).**
- About 7.05 h per octant; catalog walls sum to 59.4 node-h, including one 10.0 h TIMEOUT.
- Merge+finalize+write took 4.61–4.85 h (7.66 h for oct100, on the slow node kn103) with the GPUs idle, because of the fixed NC=256 hash.
- Most of this cost was that hash, not cf32.

---

## 6 Issues found along the way

### 6.1 In our code (fixed)

| Issue | Effect | Fix |
|---|---|---|
| `fsc_of_z` returned 1.686·D(z) | Broke all lightcone runs | Collapse-table bisection (commit 91c4145) |
| chi(z) a factor h too small | Lightcone stopped at z≈1.96 (~31% of the volume) | 70e0beb |
| `pk_websky.dat` missing /(2π)³ | Field 16× too strong; fake giants (pile-up at RTHL 18.67) | `INVESTIGATION_2026-06-14.md` |
| Box 7700 Mpc/h instead of 7700 Mpc | Cells 1.47× too coarse; mass floor ~3e12 instead of ~1.2e12 | box 5236 Mpc/h plus finecell bank (`LOW_MASS_COMPLETENESS_2026-06-15.md`) |
| AM used Ωm = 0.359 (baryons counted twice) | AM/Websky 0.14 at z=4.6 | Ωm = 0.31 |
| No Eulerian/velocity stage; D/a instead of D; 2LPT +3/7 instead of −3/7 | σ_vr off by 40×; wrong positions | 95f04a1 (CPU); 4ae7837 (GPU path) |
| MPI path: D/a, 2LPT sign, trace-identity Nyquist | — | 8aa4e9d |
| AM top-halo identity default; CLI finalize | Top 1–2 halos per z-bin kept near raw mass | 9f92778 |
| Legacy tile layout N = nsub·ntile + 2nbuff | 13.6 Mpc/h unsimulated slab at each octant plane (1.6% of volume); broke full-sky ISW, κ and tSZ | `periodic_cores` (`FULLSKY_COMPARISON_2026-09-26.md`) |
| Splice P(k) bump (Catmull-Rom vs piecewise-constant residual) | Up to +11% P(k) at 0.5 k_Nc; +3.5% σ(R) | `coarse_compensation` (D/T) |
| Merge order depended on input (22% of peaks tie in R_TH) | ~3% of merged records order-dependent; GPU run-to-run differences of 109 records | Total order (−R,x,y,z) |
| Fixed NC=256 merge hash | 4.6–4.9 h merge | Hash sized to the data |
| GPU tile-local 2LPT Nyquist mismatch | ≲1e-3 of displacement | k=0-only φ_ij zeroing |
| nbuff=16 shell cap | No M>7.5e14; tSZ deficit | nbuff=25 demonstrated on oct000; full campaign pending |
| Stale conclusion that the nhunt clamp was harmless | It held only while the field was wrong | Superseded by V2_RESULTS |

**Open issue in our code (not fixed).**

| Issue | Effect | Status |
|---|---|---|
| Ω_m(a) in the 2LPT coefficient: Julia uses `Omnr*a^3/(Omnr*a^3+OL)` (`src/MultiResolution.jl:1117,1184`; the same form in `Pipeline.jl:270`, `MultiTile.jl:367,719`, `FieldMap.jl:336`, `MPIExt.jl:606`). The standard Ω_m(a), used by Fortran (`hpkvd.f90:762,875`), has a⁻³, so Julia's value → 0 at high z where it should → 1 | Through Ω_m(a)^(−1/143), the ψ2 coefficient is ×1.012/1.019/1.029/1.042 at z=0.5/1/2/4.5. The impact on positions and velocities is expected to be small, because ψ2 is second order (not yet quantified) | **Fixed 2026-09-27**: one helper, `Cosmology.omega_m_a`, is used at all 7 sites, with a regression test in `test/test_cosmology.jl`. v2 and v3test catalogs were made before the fix; the v3 campaign will include it |

**Stale text in the paper draft and comments** (to correct):
- `paper/sections/algorithm.tex:46-47` says partial overlaps reduce R_TH, but volume reduction is disabled (`Merger.jl:68`).
- `paper/sections/maps.tex:54` and `intro.tex:58` say halo profiles are painted with XGPaint. Only CIB uses XGPaint.
- `paper/sections/performance.tex:8-9` still quotes the stale "~19 min on 4×L40S".
- `src/FieldMap.jl` comments are stale on three points: the +3/7 sign, ISW painted at displaced positions, and D/a.
- `MultiResolution.jl:591-593` says shell analysis stays on the CPU; it runs on the GPU.
- `run_gpu_octant.jl:92` logs "volume reduction", which is disabled.
- `TILING_INVARIANCE_2026-09-26.md:112` still calls the canonical merge order "not applied".
- `PERFORMANCE_2026-09-26.md` gives "20 filters"; production uses 23.

### 6.2 Observed differences with the reference (neutral; not yet validated by the PI)

| Observation | Evidence | Source |
|---|---|---|
| Our tSZ painter on Websky's own halos exceeds the released `tsz_2048` by 4–25% in C_ℓ, with cross-correlation ≥0.992 and mean y 1.006 | Hypothesis: the reference table is truncated at 4 Mpc transverse (not demonstrated) | `V2_RESULTS_2026-09-27.md` |
| In the local build, pks2map defaults to mmin=2.5e10 with no 1e13 or 0.5′ cut. On halos_10x10 it gives ~25× the released kSZ halo component at ℓ~900 | The settings used for the released map are not documented in the repo code. Observed in the local (hand-patched) clone; whether the released maps were produced with this code path or these settings is unknown; not validated by the PI | `KSZ_COMPOSITE_2026-07-19.md:167-203` |
| In the local build, the output depends on the z-cut at fixed input. As read, the compaction at `pks2map.f90:182-186` copies posxyz and rth but not vrad | In the local build, zmax 4.5 vs 6.0 on identical input gives 40% lower map RMS and 20–33% lower C_ℓ at ℓ≤1237. Observed in the local (hand-patched) clone; whether the released maps were produced with this code path, binary or zmax setting is unknown; not validated by the PI | same |
| κ: `kap_lt4.5`/Halofit(z<4.5) is 0.935–0.994 at ℓ=150–1650 and 1.00 at ℓ≈2900; ours/Halofit is 1.00–1.056 | The ~20% smallest-scale suppression described in Stein+ (`2001.08787.txt:1604-1606`, Fig 8) lies at ℓ beyond our theory file (which ends at ℓ=2894). It also refers to the total κ, including the z>4.5 Gaussian component, not `kap_lt4.5` alone. Not tested here. Tier-B suggested sub-pixel mass loss (untested) | `results/theory_v2_kappa.txt` |
| CIB decoherence printed in Stein+ Table 1 is lower than we measure on the released maps (e.g. 857×545 0.933 printed vs 0.951 measured) | Candidate: the 18-FEB-2022 CIB satellite update (untested). The garbled text extraction needs checking against the PDF | `2001.08787.txt:1464-1520`; `results/cross_v2_cib_cib.txt` |
---

## 7 Open items and next steps

| # | Item | Status / next action |
|---|---|---|
| 1 | **v3 full-sky campaign at nbuff=25** (8 octants) | Needs user approval. Estimated ~2.6 h job per octant (b3). /project is at 4517/5000 GiB; v3 needs ~370 GB more (*memory only; unverified in repo*: `memory/nbuff_mass_cap.md`). Then redo the tSZ, κ, kSZ-halo, CIB and cross-spectra scorecard |
| 2 | Tier-A statistics on current catalogs | Rerun MF, dN/dz, σ_vr, ξ, b(M) and v12 on v2/v3 with the octant jackknife (A2/A3). None exists yet |
| 3 | Field kSZ excess, 1.13–1.21× linear at ℓ=70–200 | Candidates: ψ2 in the velocities, shared box modes, the z=4.5 cutoff, residual velocity decoherence at cf32. Deciding test: 3D velocity P(k) of the ICs vs linear theory (not run) |
| 4 | CIB normalization (1.23× in mean) | Decide whether to rescale L0 to the released mean or quote the ratio |
| 5 | CIB 400 mJy flux cut | Apply to both maps before scoring the Poisson regime |
| 6 | κ offset (1.06–1.14× `kap_lt4.5`) | Unexplained; shared between the two maps relative to Halofit |
| 7 | gradrf_x/y/z NaN (~0.005% of halos; fields 31–33) | Masses and positions unaffected. Root cause untraced (zero-radius mrf shell suspected) |
| 8 | Ω_m(a) in the 2LPT coefficient: a³ in place of a⁻³ in our code (§6.1) | Code fixed 2026-09-27. Still to do: quantify the ψ2 impact (coefficient ×1.01–1.04 at z=0.5–4.5) by comparing v3 with v3test oct000 |
| 9 | Theory-test and cross-spectrum σ | Add σ to `theory_v2_*`/`cross_v2_*`; implement the Gaussian-only rule at ℓ≲30 in `compare_auto.jl`; re-derive verdicts with the paper tooling |
| 10 | 10-particle pre-AM cut | Undecided. The 0.6% sub-cell (M<1e11) fraction was measured only on the example config |
| 11 | AM below Mmin=5e11 (flat extrapolation) | Effect on the low-mass tail not checked |
| 12 | Observational and hydro-template overlays; lensed CMB | Not produced (Figs 5–7, 9–11) |
| 13 | Box replication (A6); v12 tail beyond 30 Mpc/h; 1–3 Mpc/h ξ excess; shared b(M) shortfall vs Tinker10 | Untested or uninvestigated |
| 14 | Measurements still missing | Timing of merge, finalize and write in v2; GPU memory at n=416/434; v2 cost with the shell early exit; early-exit exactness at production scale; worker-count reproducibility with the canonical merge; MPI global-FFT vs multires agreement; whether the chi(z) handling in finalize matches merge_pkvd exactly; 2LPT sign convention of raw ψ2 vs Fortran |
| 15 | Provenance | The v3test tSZ map is in purgeable scratch (`/home/yguan/scratch/websky_6144/v3test_paint/`); the v3test V2_RESULTS section is uncommitted; the reference map versions are established only from download dates and `websky_ref/md5sums.txt` (see the UPDATES entries in §4) |
| 16 | Paper-draft corrections | See §6.1 |
| 17 | Raw (pre-AM) MF vs Tinker per z-bin (A4) | Not done. It is the finder-level anchor that AM hides; so far only the z<1 tail counts exist (§3.2) |
| 18 | Other A4 theory anchors: IC field P(k) vs input; ψ1/ψ2/velocity power vs linear/2LPT theory | Not run (the velocity P(k) test is also the deciding test for item 3) |
| 19 | Completeness vs Stein+ §4.2 `halo_mass_completion.txt` | The file is on disk in `websky_ref/`; comparison not done |
| 20 | A8 convergence appendix: nbuff (16/25 oct000 runs exist), filter-bank count/spacing, cellsize (the factor-h runs as a controlled test), catalog coarse_factor | Not assembled |
| 21 | Lightcone, v2-code Fortran/Julia finder comparison (the existing test is z=0, a 320 Mpc/h box, pre-v2 code) | Not done |
| 22 | B1 κ × halo-density cross-spectrum | Not done (the other B1 crosses are in §4) |
| 23 | B2 one-point statistics: κ and y PDFs, tSZ/κ peak counts | Not done |
| 24 | Plan questions Q1 (pass-criteria wording), Q2 (replication test design; blocks item 13), Q3 (B-items in paper v1 or deferred) | Undecided; need user decisions |

---

## 8 WebSky 2.0 outlook (brief)

Source: `docs/deep_lightcone_roadmap.md`, revised 2026-09-27.

**What WebSky2.0 is.** WebSky2.0 CMB/LSS runs use the same corner-observer octant geometry as WebSky1.0, at finer cells. They are not radial shells.
- Cells about 0.33–0.39 Mpc/h on 12288–16384 grids.
- nbuff 64–69, Planck18.
- About 3e9 halos per octant.

**Work items, in order:**
1. **Lift the GPU shell limit.** `_MAX_SHELLS_GPU` = 512 caps nhunt at 24. The plan is to stream shells from global memory; the exact early exit keeps the mean cost near the typical halo radius. The same limit sets the current tSZ mass cap.
2. **Slim catalog format.** At 33 floats, ~3e9 halos is about 450 GB per octant. Offer a 10-float option.
3. **Finer single-box octant** at N=12288 (~8× cost). This needs `periodic_cores` support in the MultiTile/MPI paths for multi-node runs.
4. **Planck18 option.**
5. **Optional offset-box deep lightcones** for line-intensity mapping (LIM).
6. **4-field 2LPT**, blocked until the modern Fortran on gw has been read.

---

## Appendix: source index

**Configs and drivers.**
- `validation/websky_6144/production_v2/` (README, `config_v2_oct*.toml`, `config_v3test_oct000.toml`, `run_catalog.slurm`, `run_fieldmap.slurm`)
- `validation/websky_6144/production/` (`paint_octant.jl`, `halo_profiles.jl`, `paint_cib.jl`, `run_paint.slurm`, `run_cib_prod.slurm`)
- `validation/websky_6144/run_gpu_octant.jl`, `run_fieldmap_octant.jl`, `apply_abundance_match.jl`

**Notes, current.**
- Units and results: `CONVENTIONS.md`; `validation/paper/V2_RESULTS_2026-09-27.md`; `validation/paper/FULLSKY_COMPARISON_2026-09-26.md`
- Tiling and splice: `validation/tiling/TILING_INVARIANCE_2026-09-26.md`; `validation/tiling/SPLICE_COMPENSATION_2026-09-26.md`
- Performance: `validation/performance/PERFORMANCE_2026-09-26.md`; `validation/performance/SHELL_EARLY_EXIT_2026-09-27.md`
- Theory, 2LPT and RNG: `validation/paper_theory/THEORY_ANCHORS_2026-09-26.md`; `validation/NOTES_2LPT_NYQUIST_2026-09-24.md`; `validation/rng/RNG_LARGE_SCALE_TEST_2026-09-23.md`
- AM fix: `validation/websky_6144/AM_TOPHALO_FIX_2026-09-25.md`
- Plans: `docs/paper_comparison_plan_2026-09.md`; `docs/multi_resolution_fft.md`; `docs/deep_lightcone_roadmap.md`

**Notes, Tier-A/B validation** (partly superseded):
- `validation/websky_6144/FORTRAN_COMPARISON_2026-06-14.md`
- `LOW_MASS_COMPLETENESS_2026-06-15.md`
- `BIAS_PAIRWISE_V12_2026-07-16.md`
- `TIER_B_SUMMARY_2026-07.md`
- `KSZ_COMPOSITE_2026-07-19.md`
- `INVESTIGATION_2026-06-14.md`
- `INVESTIGATION_2026-04-16.md` (superseded)

**Results.** `validation/paper/results/auto_v2_*.txt`, `theory_v2_*.txt`, `cross_v2_*.txt`

**Reference.**
- Stein+2020: `~/work/peakpatch/2001.08787.txt`
- Released maps: `/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/`
- Local Fortran clone: `~/work/peakpatch` (hand-patched near the cosmology block)

**Job IDs.**

| Purpose | Jobs |
|---|---|
| v2 catalog+AM | 5696115, 119, 123, 127, 131, 135, 139, 143 |
| v2 field maps | 5696740–5696747 |
| v2 halo paint | 5696116, 120, …, 144 |
| v2 CIB | 5696117, 121, …, 145 |
| v2 analysis / theory | 5703683 / 5703908 |
| v2 tail check | 5704125 |
| Seam pre-flight | 5696038, 5696088, 5696647 |
| v3test catalog | 5710456 (no early exit: 5704175, cancelled) |
| v3test tail check | 5715228 |
| v3test tSZ paint | 5715273 |
| Early-exit validation | 5708923, 5711281 |
| Tiling / merge determinism | 5695201, 5695214 |
| Benchmarks | 5692900, 5692991–5692994, 5693005 |
| RNG | 5635142, 5637544, 5638894, 5641162 |
| Tier-A catalog history | 3968163 (finecell raw), 4033130 (finecell AM), 4387351 (cf32) |
| Frozen campaign (superseded) | 4438182–4438197, 5627673; AMv2 5692696–5692703 |

**Commits.** 91c4145 (fsc_of_z), 70e0beb (chi), 95f04a1 and 4ae7837 (D/a, 2LPT sign), 8aa4e9d (MPI), 9f92778 (CLI finalize, AM top halo), 285c72c (ISW at Lagrangian positions), 1c11cb3 (v2: periodic cores, compensation, deterministic merge), a887192 (shell early exit, nbuff diagnosis).
