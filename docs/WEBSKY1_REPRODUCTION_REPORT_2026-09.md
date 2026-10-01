# PeakPatch.jl reproduction of WebSky 1.0: fidelity and performance

*Status as of 2026-10-01 (updated to the production-v4 / v4fs campaign). Branch `websky2-resolution`. Audience: internal (PI) first; after PI review, Bond, Carlson, Stein et al.; input to the methods paper (`paper/`).*

> **Distribution gate: internal draft.** §6.2 and the reference-side parts of §4 (the tSZ-painter, kSZ-construction and CIB Table-1 checks) describe observed differences with the Websky reference code and maps. They must **not** go to the Websky authors, or anyone outside the group, until the PI has personally validated them and explicitly approved sharing. An external version should drop §6.2 and the reference-side remarks in §4.

All lengths are in Mpc/h, and masses in M☉/h unless stated. The Fortran peak-patch code and the public Websky products use Mpc. Conversion: `boxsize[Mpc/h] = L[Mpc]·h`. For Websky, 7700 Mpc × 0.68 = 5236 Mpc/h (`CONVENTIONS.md`).

**Result tiers used below.**
- **Current: "v4fs".** The production-v4 catalogs (8 octants; `validation/websky_6144/production_v4/`) are v3 (nbuff=25, shell early exit, Ω_m(a) fix) plus the two Fortran-faithful catalog steps found missing on 2026-09-28: the `merge_pkvd` **volume reduction** and the `hpkvd` **per-tile lightcone peak threshold** (§2.2, `validation/paper/CLUSTERING_EXCESS_2026-09-28.md`). They are abundance-matched with one full-sky table and `tail_N`=10. The field maps are the v3 field maps, which depend on neither change.
- **Superseded, quoted for comparison:**
  - **v3fs** (exclusion-only merge, constant threshold): ~7% more biased than Websky; includes ~86M extra sub-1e12 halos per octant.
  - **v3** (per-octant AM): v3fs catalogs, but each z-bin's top halo is pinned at M(N=1), so N(>2e15) = 0.
  - **v2** (nbuff=16): the shell-search cap removes every M>7.5e14 halo.
  - the frozen 2026-07 campaign, which had seam gaps, and runs from before the 2026-06 fixes.
- The 2026-06/07 Tier-A numbers on the single validation octant `oct000_finecell_AM` are superseded (§3), except the Fortran finder test.
- **Do not quote:** the "19 min/octant" and "~0.9M halos/octant, 120× fewer than Websky" figures. They come from the 2026-04 configuration, before the bug fixes (`validation/performance/PERFORMANCE_2026-09-26.md` §1).

**Framing.** Where we list differences from the reference Fortran code or the released maps, these are *observed differences with evidence*. The PI has not yet personally validated them, and none has been communicated upstream. Sending this draft to Websky authors would count as communicating them, so the distribution gate above applies.

**Glossary.**
- **Tier-A:** catalog-level statistics (MF, dN/dz, σ_vr, ξ, b(M), v12) against Websky's 10°×10° halo patch; first measured in 2026-06/07 on one validation octant, rerun with 8 caps on v3fs (2026-09-28) and v4fs (2026-10-01).
- **Tier-B:** map-level (painted) statistics from the same 2026-07 period, on oct000 caps or octants.
- **Plan codes** (A0–A8, B1–B6, Q1–Q3): items in `docs/paper_comparison_plan_2026-09.md`. For example, A2 = rerun Tier-A/B on current catalogs; A3 = error model and pass criteria; A4 = theory anchors independent of Websky; A6 = full-sky assembly and replication checks; A8 = convergence appendix; B5 = kSZ hydro-template overlay.
- **kSZ halo constructions** (`KSZ_COMPOSITE_2026-07-19.md`):
  - **W:** uncompensated per-halo Battaglia τ painting;
  - **We:** the Websky-literal variant, with a Δ=3-mean compensation sphere and gas mass f_b·m_h;
  - **Wc:** a zero-net compensated variant, a full-τ sphere with R=4 r_vir.
- **cf:** the coarse_factor of the multi-resolution split (coarse grid M = N/cf).
- **octZYX:** the octant label; the three bits give the observer sign per axis.

**Sources.** Facts sourced only to `memory/*.md` notes are not in the repo and are marked *(memory only; unverified in repo)*.

**Ranges.** Ratio ranges in §1, §3 and §4 are re-read from `validation/paper/results/*_v4fs_*.txt` and `tierA_v4fs.txt` (analysis job 5752162, Tier-A job 5752164). v3fs and v2 values are from `*_v3fs_*` and `*_v2_*`. Where these differ from the prose notes (`V2_RESULTS_2026-09-27.md`, `V3_RESULTS_2026-09-28.md`, `V4_RESULTS_2026-10-01.md`), the results files take precedence.

---

## 1 Summary

**Verdict.**
- **Halo catalog: reproduced.** All Tier-A statistics but two (plus two marginal bins, below) agree with Websky within 2σ of our cap-to-cap scatter (Websky patch against the same cap placed in each of our 8 octants; §3.1):
  - N(>M|z) after AM 0.97–1.05 (M≥3e12, z<3; largely enforced by AM), and 0.987 ± 0.007 at the 1.2e12 floor for z 3–4.5;
  - dN/dz (M>3e12) 0.95–1.07 at z<4.25;
  - σ_vr 0.96–0.99;
  - ξ(3–15 Mpc/h) 0.958 ± 0.026 (1.6σ); ξ at 1–3 Mpc/h 0.98–1.01;
  - b(M) 0.97 / 1.00 / 1.01 / 1.07 in 4 bins (all ≤1.1σ);
  - v12 0.963 ± 0.040.

  The two exceptions are unchanged from v3fs: the last dN/dz bin at our z_max edge (z 4.25–4.5, 2.5σ), and M>1e13 dN/dz at z 2.25–3 (W/ours 0.90–0.93, i.e. ours 8–12% high; 2.8σ, 1.8σ, 2.6σ; §7). Marginal: ξ at r = 8.6 Mpc/h 0.936 (−2.2σ) and N(>1e13) at z 3–4.5 0.913 (−2.0σ).
- **What made it match:**
  - **Volume reduction.** v3fs was ~7% more biased (ξ 0.877 ± 0.029, 4σ). Our merge skipped the Fortran volume-reduction pass; a same-raw-catalog A/B reproduced the small-scale excess exactly (ξ(1–3) 0.796 vs 0.79) and ~70% of the ξ(3–15) offset (0.916). The remaining ~4% (0.96) persists in v4fs, at 1.6σ (§3.1).
  - **Per-tile threshold.** It makes no difference above 1e12, but it removes most sub-1e12 halos that Websky's catalog does not contain.
- **Raw (pre-AM) mass function vs Tinker (full sky):**
  - 0.90–0.92 at 1e13 (z<1.5) and 0.99–1.02 at 1e14 (z<1);
  - 0.62–0.71 at 1e12 (z<3), falling to 0.54 at z 3–4.5.
  - AM matches Tinker to ≤3% above 3e12 except in sparse bins (0.94–0.97 for 1e14 at z 2–3, 5e14 at z 1–1.5, 1e15 at z<0.25), and to 0.92–1.07 at the 1e12 floor.
- **High-mass tail (full sky, z<1):** N(>1e15) = 411 and N(>2e15) = 10, against Tinker's 412.3 and 12.5. v2 had 0 above 1e15 because of the nbuff=16 shell cap.
- **Seams:** fixed since v2.
- **CIB decoherence:** agrees with the released maps (post-2022 CIB update) to ≤0.022 for all 10 pairs, with no σ attached. This uses a same-family CIB model (XGPaint `CIB_Planck2013`), so it is not an independent check. The printed Stein+ Table 1 differs from both (§6.2).
- **ISW at low ℓ:** consistent within large cosmic variance, using an unofficial σ computed for this report (§4).
- **Field kSZ vs linear theory:** 1.01–1.10 at ℓ = 65–200 and 0.96–1.02 up to ℓ≈930 (v3 field maps). The v2 excess of 1.13–1.21 at ℓ = 70–200 is gone, but **which v2→v3 change removed it has not been established** (§7).
- **tSZ, close but not formally passing:**
  - mean y **0.990** (v3fs 1.035, which included the sub-1e12 halos Websky does not have);
  - band average over ℓ = 100–5000 is 1.024 ± 0.002 (3/40 bands >4σ);
  - C_ℓ is 1.01–1.06 at ℓ = 1000–4100 and 1.07–1.16 at ℓ = 500–1000, rising to 1.17–1.33 at ℓ = 270–500 and above 1.8 at ℓ < 150, where a few clusters dominate and σ is large;
  - our painter alone is 1.04–1.16 above the reference on identical halos (§4), which accounts for all of the excess at ℓ≳800, most at ℓ≈500, and about half at ℓ≈300.
- **Not yet reproduced to the pre-declared criteria:**
  - κ (6–13% high against `kap_lt4.5`; ours is 1.00–1.05 × Halofit, the reference 0.94–0.99);
  - CIB, which differs by a mostly painter-driven amplitude: mean 1.23, and 1.18 on Websky's own halos through the same model. There is also a frequency-dependent residual at 857 GHz and in the Poisson regime, where the flux cut has not yet been applied.
- **Differs by construction:** kSZ halo power at mid-ℓ, where the result depends on how the halos are painted.
- **Cost (v4):**
  - The catalog takes 1.62–1.95 h per octant on one 4×L40S node: the full sky in about 1.95 h on 8 concurrent nodes (inferred), or 13.4 h on one node. Websky took 3.84 h on 1128 cores.
  - The full product set is about 74 L40S-GPU-h, with no stored initial conditions (ICs).
  - Peak memory is 89.7 GiB host + 4 × ~22 GiB GPU ≈ 0.19 TB **per octant run**, ~40× below the 7.67 TB of Websky's full-box run; for 8 concurrent octants ≈1.5 TB.
  - This is **not** a like-for-like hardware comparison (§5).

**Headline results.**
- **Phase space: the AM-insensitive statistics** (v4fs; §3.1): σ_vr 0.96–0.99; ξ(3–15) 0.958 ± 0.026; b(M) 0.97–1.07; v12 0.963 ± 0.040.
- **Finder equivalence (limited scope).**
  - Kept-peak densities are 3.925e-3 (Julia) and 3.946e-3 (Fortran) per (Mpc/h)³ at the same numeric cellsize, 1.2533: Fortran/Julia = 1.005 (`validation/websky_6144/FORTRAN_COMPARISON_2026-06-14.md`).
  - Scope: a 320.8 Mpc/h box, a z=0 snapshot, pre-v2 code, different seeds; kept-peak density only.
- **Seams fixed since v2.** κ-field and ISW maps have 0 exactly-zero pixels, against 190k and 1.5M in the frozen campaign (`validation/paper/V2_RESULTS_2026-09-27.md`).
- **tSZ: the v2 deficit had two causes in our pipeline, both fixed.**
  - v2 full sky gave 0.80–0.88× Websky at ℓ=270–1400, because the nbuff=16 shell-search cap left N(>1e15, z<1) = 0 per octant against 51.5 for Tinker.
  - v3 (nbuff=25) restored the raw tail. Its per-octant AM then capped each z-bin at about 1.7e15; the full-sky AM + `tail_N` removed that cap.
  - v4fs tSZ C_ℓ / Websky is 1.33 (ℓ≈300), 1.17 (480), 1.09 (770), 1.06 (1030), 1.02 (2000).
- **CIB** (band averages 1–3% below v3fs: 0.989 at 100 GHz to 0.973 at 857 GHz).
  - Mean intensity is 1.23× the released maps at every frequency. Websky's own halos through the same XGPaint model give 1.18, so most of the offset comes from the painter/model and not from the catalog.
  - Clustered C_ℓ band averages are 1.44 (100 GHz) … 1.41 (353), 1.37 (545), 1.23 (857). Poisson-regime ratios run 1.24→0.66 across frequency.
  - Decoherence agrees with the released maps to ≤0.022 for all 10 pairs (no σ).
- **Performance** (v4, one 4×L40S node per octant):
  - catalog 1.62–1.95 h per octant, with the pipeline at 52–54 min (v3: 100 min). The per-tile threshold cuts the pre-merge candidates from 338M to 201M;
  - IC storage is 0 B, against 5.9 TB for Websky.

---

## 2 What we implemented, and how it differs from the Fortran/Websky production

### 2.1 Production-v4 configuration

Sources: `validation/websky_6144/production_v4/README.md` and `config_v4_oct000.toml`. v4 is v3 (v2 with nbuff 16→25, the GPU shell early exit and the Ω_m(a) fix) plus `volume_reduction` and `peak_threshold_per_tile`.

| Item | Value |
|---|---|
| Seed, N, box, cell | 12345; 6144; 5236 Mpc/h (= 7700 Mpc); 0.852213 Mpc/h |
| Tiling | ntile=16, tile mesh n=434, nbuff=25, nsub=384, `periodic_cores=true` (N = nsub·ntile). `[grid] boxsize` is the **per-tile** box, 369.86 Mpc/h. nhunt = 24 (R_TH cap 20.4 Mpc/h) |
| Coarse grid | coarse_factor=32 (M=512, block 12), `coarse_compensation=true` |
| Lightcone | ievol=1, z_out=0, z_max=4.5; observer at ±2618 Mpc/h per axis (the core corner, set by the octZYX bits) |
| Physics | ilpt=2, ioutshear=1, wsmooth=1, rmax2rs=0.0 |
| Catalog steps (v4) | `peak_threshold_per_tile = true` (Fortran per-tile lightcone threshold); `volume_reduction = true` (Fortran merge_pkvd reduction after exclusion) |
| Cosmology | Ωm=0.31, Ωb=0.049, ΩΛ=0.69, h=0.68; P(k) from `pk_websky.dat` (CAMB, σ8=0.8100) |
| Filters | `filters_websky_finecell.dat`: 23 top-hats, Rf = 1.406–30.44 Mpc/h, ratio 1.15 |
| Collapse table | `HomelTab_websky.dat`, 50×20×20 |
| AM (v4fs; as v3fs) | one Tinker08 M200m table from all 8 raw octants (fsky=1), 10⁴ mass × 46 z bins over z<4.5, `tail_N`=10 (`apply_abundance_match_fullsky.jl`) |

### 2.2 Differences from the Fortran code and the Websky production

| Aspect | Fortran / Websky (Stein+2020) | PeakPatch.jl v4fs | Source |
|---|---|---|---|
| Units | Mpc | Mpc/h everywhere, because chi(z) uses H0=100h | `CONVENTIONS.md` |
| White noise | 48-bit LCG, filled sequentially per MPI slab; the realization depends on the task count | Counter-based Threefry2x (20 rounds) with key (seed,0), plus Box-Muller. Each cell depends only on (seed, global index). A bit-exact LCG port is also available | `src/InitialConditions/RandomField.jl:36-103`; `LCG.jl` |
| Initial conditions | Global slab-decomposed FFTW3, 7 values per site, 5.9 TB stored | No global fine FFT. The coarse M³ periodic FFT is spliced to a per-tile isolated residual convolution with Catmull-Rom interpolation (MUSIC/Hoffman-Ribak style); D/T compensation in v2; nothing stored | `docs/multi_resolution_fft.md:60-104`; `src/MultiResolution.jl:584-735` |
| Buffer | ≈ the largest Lagrangian halo radius (~40 Mpc) | nbuff=25 cells = 21.3 Mpc/h (v2: 16). The shell search is capped at nhunt = min(nbuff−1, ⌊1.75 Rf_max/a⌋), the same formula as `hpkvd.f90:207-209`; nhunt=24 is the GPU kernel's shared-memory maximum. Largest raw R_TH over the 8 octants is 17.2–20.4 Mpc/h (one halo at the cap, in oct101) | `src/MultiResolution.jl:630`; `V3_RESULTS_2026-09-28.md` |
| Filter bank | Paper: Rf,min = 2 a_latt; Rf,max = 36 Mpc (24.5 Mpc/h); count and spacing not stated. `filter_gen.py` default: 1.65 a_latt, spacing 1.15 (the same as ours). Production bank not located; the evidence conflicts (paper 2.0; the dense-bank residual suggests ≈1.2) | Rf,min = 1.65 a_latt, Rf,max = 30.44 Mpc/h, 23 filters, spacing 1.15 | `data/filters_websky_finecell.dat`; `2001.08787.txt:300-303`; `LOW_MASS_COMPLETENESS_2026-06-15.md:66-80` |
| Peak threshold on the lightcone | Per tile: fcrit = fsc_of_z(z_tile), z_tile from the tile centre's distance (`hpkvd.f90:510-514` → `get_pks`, `peakvoidsubs.f90:82`) | Same since v4 (`peak_threshold_per_tile`). Up to v3fs we used a constant fsc_of_z(0), which a 2026-04 note had wrongly attributed to Fortran. Shell analysis uses the per-peak z_pk either way. Effect: removes most sub-1e12 candidates (pre-merge 338M → 201M per octant), ≤0.5% above 1e12 (`results/threshold_ab.txt`) | `src/MultiResolution.jl` (tile loop) |
| Ω_m(a) in the 2LPT coefficient | Standard (Ωm/a³)/(Ωm/a³+ΩΛ) (`hpkvd.f90:762,875`) | Same since 058207a (v2 used a³ in place of a⁻³; §6.1). v3 includes the fix | `src/Cosmology/Cosmology.jl` |
| Merge | Lagrangian exclusion (centre inside a larger halo), then **volume reduction**: each halo loses its caps beyond the mid-plane with every overlapping neighbour (`merge_pkvd.f90:139-140`, `exclusion.f90:120-204`; Stein+2020 §2.3). Fixed NC=256 hash | Same since v4 (`volume_reduction`), with a data-sized hash and a total order (−R_TH,x,y,z); each pair visited once (as read, the local Fortran clone visits equal-radius pairs twice; observed, not validated by the PI). Up to v3fs the reduction pass was skipped (a comment wrongly said Fortran had it commented out); that accounts for the small-scale excess and ~70% of the ~7% bias offset (`CLUSTERING_EXCESS_2026-09-28.md`) | `src/Merger/Merger.jl`; `Exclusion.jl` |
| Eulerian positions and velocities | merge_pkvd | Port: x = q+ψ1+ψ2, v = a·100·E·f·(ψ1+2ψ2). Velocity epoch from one chi→z evaluation at the Eulerian distance, with no iteration | `src/Merger/Merger.jl:106-145` |
| AM | Tinker M200m, Δz=0.1 bins, 10⁴ mass bins, bilinear interpolation; **one full-sky table**; above a slice's top halo the table is left at identity, so the bilinear lookup blends the top halos with their raw masses; only halos with >10 particles before AM | Same table method (46 z-bins, Δz≈0.098, over 5e11–1e16), **one full-sky table from all 8 octants**. Above the last mass edge with ≥10 halos per z-bin, the fractional correction is frozen (`tail_N`), which is monotone and keeps ranks (the reference's identity default can invert them). **No particle cut** | `src/AbundanceMatch/AbundanceMatch.jl`; `~/work/peakpatch/python/catalogue_tools/abundance_match/make_abundancematch_table.py`; `2001.08787.txt:1044-1052` |
| Catalog format | 10 floats, ~9e8 halos, 33 GB | 33 Float32 (extended pksc), 112.44–112.53M halos/octant (v4; v3: 198.4M, the difference almost all below 1e12), ~14.8 GB/octant | `src/Merger/Merger.jl:84-162` |
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

**Method and caveats.**
- The Tier-A statistics are measured on the **v4fs** catalogs (`validation/paper/tierA_v3fs.jl` with TIERA_CAMP=v4, job 5752164; `results/tierA_v4fs.txt`; note `V4_RESULTS_2026-10-01.md`). v3fs (job 5735978, `TIERA_V3FS_2026-09-28.md`) is shown for comparison. The comparison target is Websky's public 10°×10° patch (`halos_10x10.pksc`), because the full Websky `halos.pksc` is not available locally.
- The cap is the one inscribed in that patch (half-angle 5.00°, Ω = 0.02395 sr), with the 2026-07 shell (1600–2000 Mpc/h, z≈0.6–0.8), mass bins and windows. The identical cap is placed at the centre of each of our 8 octants; "ours" is the mean ± std over the 8 caps, and z-scores use std·√(1+1/8) (the printed values in the results file are 1.06× larger).
- Post-AM N(>M) agreement is largely enforced by AM, because both catalogs are matched to Tinker08 M200m, so it is not independent evidence. The AM-insensitive statistics are ξ, b(M), σ_vr and v12, the raw MF against Tinker, and the Fortran finder test.
- Our realization is not Websky's, and the Websky patch is one cap, so each comparison is "is Websky's patch a plausible draw from our cap distribution".
- Stein+2020 itself shows no halo mass function, dN/dz or clustering plot, so these statistics are additional validation.

### 3.1 Scorecard (W/ours = Websky patch / mean of our 8 caps)

| Statistic | v3fs W/ours (z) | **v4fs W/ours (z)** | Verdict (v4fs) |
|---|---|---|---|
| N(>M) after AM, M ≥ 3e12, z < 3 | 0.97–1.05 (within ±1.7) | 0.97–1.05; 1.08 for >1e13 at z<0.5 (within ±1.8) | agree (largely enforced by AM) |
| N(>1e14), z = 0–0.5 / 0.5–1 / 1–1.5 | 1.12 / 0.96 / 1.03 | 1.12 / 0.97 / 1.04 (≤1.8) | agree |
| N(>1.2e12), z 3–4.5 | 0.919 ± 0.007 (−10) | **0.987 ± 0.007 (−1.7)** | agree |
| dN/dz, M>3e12, Δz=0.25, z 0–4.25 | 0.95–1.08 | 0.95–1.07 (within ±1.8) | agree |
| dN/dz, M>3e12, z 4.25–4.5 | 1.27 ± 0.10 | 1.27 ± 0.10 (+2.5) | ours low in the last bin, at our z_max=4.5 edge (Websky reaches 4.6) |
| dN/dz, M>1e13, z 2.25–3 | 0.90–0.93 | 0.90–0.93 (−2.8, −1.8, −2.6) | ours 8–12% high; **open**, unchanged by v4 |
| σ_vr, 5 bins 1e12–3.2e14 | 0.954–0.989 | 0.962–0.991 (≤1.2) | agree |
| **ξ(r), M>5e12, mean 3–15 Mpc/h** | 0.877 ± 0.029 (−4.0) | **0.958 ± 0.026 (−1.6)** | agree |
| ξ(r) at r = 1.22 / 1.80 / 2.66 Mpc/h | 0.750 / 0.808 / 0.836 (−4.5 to −6.7) | **0.977 / 1.010 / 0.983** (≤0.6) | agree |
| b(M) ∝ √ξ̄(6–18 Mpc/h), 5e12–1.3e13 / 1.3–3.2e13 / 3.2–7.9e13 / 0.8–2.5e14 | 0.936 / 0.958 / 0.958 / 0.996 | **0.967 / 1.003 / 1.008 / 1.067** (≤1.1) | agree |
| v12(r), M>1e13, mean 5–30 Mpc/h | 0.927 ± 0.041 (+1.7) | 0.963 ± 0.040 (+0.9) | agree |
| Finder, matched cell 1.2533 Mpc/h, z=0 | Kept-peak density 3.946e-3 (Fortran) vs 3.925e-3 (Julia) per (Mpc/h)³, Fortran/Julia 1.005. Scope: 320.8 Mpc/h box, z=0 snapshot, pre-v2 code, different seeds, density only | — | `FORTRAN_COMPARISON_2026-06-14.md:50-60` |
| Lightcone vs snapshot (z<0.19) | 0.978; GPU = CPU (raw 69008 = 69008; merged 25898 vs 25897) | — | same |

**Raw (pre-AM) mass function vs Tinker08, full sky** (plan A4; v4 raw, v4fs AM; v3 raw in brackets):

| M threshold | raw / Tinker at z = 0–0.25, 0.25–0.5, 0.5–1, 1–1.5, 1.5–2, 2–3, 3–4.5 | AM / Tinker |
|---|---|---|
| > 1e12 | 0.62, 0.64, 0.67, 0.70, 0.71, 0.69, 0.54 [0.57–0.73] | 0.92, 0.97, 1.04, 1.07, 1.00, 1.02, 0.92 |
| > 3e12 | 0.80, 0.82, 0.84, 0.85, 0.84, 0.78, 0.62 [0.66–0.92] | 0.97–1.00 |
| > 1e13 | 0.90, 0.91, 0.92, 0.91, 0.87, 0.77, 0.53 [0.56–1.02] | 0.98–1.00 |
| > 1e14 | 1.01, 1.02, 0.99, 0.89, 0.69, 0.47, — [1.20, 1.19, 1.15, 1.03, 0.80, 0.55] | 0.97–1.00 (z<3) |
| > 5e14 | 1.08, 1.04, 0.89, 0.56 (z<1.5) [1.35, 1.26, 1.07, 0.72] | 0.95–1.01 |
| > 1e15 | 0.83, 0.81, 0.62 (z<1) [1.00, 1.05, 0.92] | 0.94, 1.05, 0.96 |

With volume reduction the raw MF sits closer to Tinker at 1e14–5e14 and lower at 1e13 and above 1e15. AM restores Tinker in every case.

**Cause of the v3fs clustering excess** (`validation/paper/CLUSTERING_EXCESS_2026-09-28.md`):
- Our merge skipped the Fortran volume-reduction pass.
- Diagnostics: abundance matching is not the cause (raw-rank and AM selections are 99.9% the same). v3fs absolute bias was 0.97–1.01 × Tinker10 below 3e13, Websky's 0.93–0.99. v3fs had 16–17% more close pairs at u = d/(R_i+R_j) < 0.5 (6σ).
- A/B of our own pipeline (same raw halos, z = 0.7 snapshot, reduction on/off): ξ(1–3) 0.796, ξ(3–15) 0.916, and the same close-pair profile as Websky/ours. With reduction, our low-mass bias is 0.98–1.00 × Tinker10.
- The second Fortran difference, the per-tile lightcone threshold, has no effect above 1e12 (A/B, `results/threshold_ab.txt`: ≤0.5% in N(>M|z), ξ ratio 1.000).
- v4fs includes both. It brings every clustering statistic within 1.6σ (table above); the ξ(3–15) residual 0.958 ± 0.026 equals the ~0.96 the reduction A/B left unexplained, so a ~4% offset may remain, at the edge of what one Websky patch can resolve.

**Comparison with the 2026-07 single-cap Tier-A** (superseded `oct000_finecell_AM`, one cap, no error bar; `BIAS_PAIRWISE_V12_2026-07-16.md`): its "ξ bias ratio 1.06 (3–15 Mpc/h)" was real for exclusion-only catalogs (v3fs: 1.068 ± 0.018). Its b(M) ratios and v12 (0.997) were one draw within the scatter.

### 3.2 Current production catalogs (v4fs, with v2, v3 and v3fs for comparison)

Sources: `validation/paper/V4_RESULTS_2026-10-01.md` (tail log `results/tails_v4fs.log`), `V3_RESULTS_2026-09-28.md` (tail jobs 5728241 for v3, 5735812 for v3fs) and `V2_RESULTS_2026-09-27.md`.

**Per octant.**

| Quantity | v2 (nbuff=16) | v3 / v3fs (nbuff=25) | **v4 / v4fs** |
|---|---|---|---|
| Pre-merge candidates | — | ~337.8M | 200.95–201.17M |
| Halos / octant (post-merge) | 198.4M | 198.39–198.52M | 112.44–112.53M (mean R_TH 1.6525 → 1.6323 Mpc/h after reduction; 25–38 halos removed by it) |
| Raw N(>1.7e12) | 4.09e7 | 4.09–4.10e7 | 3.887–3.893e7 |
| AM N(>1e14) | 6.25e4 (forced equal) | v3: 6.25e4 (forced equal); v3fs: 6.22–6.27e4 | 6.22–6.28e4 |
| Raw / AM max M | 7.4e14 / 8.5e14 | raw 1.8–3.1e15; AM v3fs 1.84–2.91e15 | raw 1.78–2.75e15; AM 1.80–2.72e15; max raw R_TH 17.0–19.7 Mpc/h (below the 20.4 cap) |

**Full sky, z<1 (sum over 8 octants).**

| | raw v3 | v2 AM | v3 AM (per octant) | v3fs AM | raw v4 | **v4fs AM** | Tinker |
|---|---|---|---|---|---|---|---|
| N(>3e14) | 34910 | ~28920 | 28874 | 28861 | 28971 | 28851 | 29255 |
| N(>5e14) | 7118 | ~6000 | ~5937 | 5942 | 5829 | 5929 | 5878 |
| N(>1e15) | 415 | **0** | 408 | 416 | 321 | **411** | 412.3 |
| N(>2e15) | 11 | 0 | **0** | 9 | 9 | **10** | 12.5 |

(v2 AM sums are 8× the oct000 value.)

By z-bin, v4fs N(>1e15) is 120 / 210 / 81 at z = 0–0.25 / 0.25–0.5 / 0.5–1, against Tinker 127.7 / 200.3 / 84.2; N(>2e15) is 5 / 4 / 1 against 6.0 / 5.5 / 1.0.

Notes:
- nbuff=25 gives nhunt=24, the largest the GPU shell kernel allows (482 shells against `_MAX_SHELLS_GPU` = 512).
- **Sub-1e12 halos.** v4 has ~86M fewer halos per octant than v3, almost all below 1e12, from the per-tile threshold (pre-merge 338M → 201M; the reduction itself removes only 25–38 halos per octant). On oct000, *abundance-matched* v4fs/v3fs counts are 0.07–0.12 at 1e11–1e12, 0.85–1.03 at 1–3e12 and 1.000 above 3e12, so ΣM^{5/3} (∝ mean y) falls by 3.8% (0.998 above 1e12; `validation/paper/ymoments_v3fs_v4fs.jl`, `results/ymoments_oct000_v3fs_v4fs.txt`). The measured full-sky mean y falls by 4.3%. This brings our mass range closer to Websky's, which keeps only halos above 10 particles before AM.
- **Ω_m(a) fix, size on the catalog** (oct000, all 198,432,491 halos matched one-to-one between pre-fix v3test and v3 on fix-invariant keys; `check_omfix_oct000.jl`, job 5720824): masses unchanged; rms position shift 0.007 Mpc/h (max 0.06); rms velocity shift 1.2 km/s, i.e. 0.15% at χ<1300 rising to 0.34% at χ>3500 Mpc/h.
- Catalog record order follows tile completion order and is **not** canonical between runs; any positional comparison must match halos first.

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
- Estimator A3: full-sky C_ℓ of the 8-octant v4fs maps against the released Websky v0.0 maps, with no mask, pixel windows divided out, and Δℓ/ℓ≈0.1.
- Analysis job: 5752162 (v4fs; ours vs Websky, cross-spectra and theory). v3fs: 5731508; v3: 5728240; v2: 5703683 / 5703908.
- Reference files are in `/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/`.
- Pass criterion, declared in advance (`docs/paper_comparison_plan_2026-09.md:77-105`): the band-averaged |r−1| < 3σ over the declared ℓ range, and no single band above 4σ.
- Construction-dependent differences (kSZ halo term at mid-ℓ, κ at ℓ≳1700) are not scored pass/fail.

> **Caveats on the "status" column.**
> - The σ values were computed for this report as an inverse-variance mean of the `sigma_ratio` column in `validation/paper/results/auto_v4fs_*.txt` (the same computation reproduces the v2 values quoted previously). No repo document states formal verdicts; they should be re-derived with the paper tooling.
> - `compare_auto.jl:28-34` uses max(Gaussian, jackknife) at all ℓ, which is not the plan's Gaussian-only rule at ℓ≲30.
> - The theory and cross-spectrum files carry no σ, so those tests are not scored. The analysis log does print a fractional Gaussian σ for the κ and field-kSZ theory tables (e.g. field kSZ at ℓ=78: 1.103 ± 0.043, 2.4σ), so the field-kSZ test could be scored.

| Product (declared ℓ range) | ours / Websky (v4fs) | ours / theory | Status vs pre-declared criterion |
|---|---|---|---|
| κ, z<4.5 vs `kap_lt4.5` (30–3000) | 1.06–1.08 (ℓ 270–430); 1.10–1.11 (850–1400); 1.13 (2600–3000); 1.00–1.01 (7600–8100). v4fs/v3fs 0.974–1.002 (the halo term changes slightly) | Halofit: 0.995–1.053 (ℓ 150–1650), 1.13 (2894). Field κ vs linear Limber: 0.98–1.04 (150–1650). Websky vs Halofit: 0.935–0.994 (150–1650), so the offset is shared between the two maps | **Fail** at ℓ 30–1700: 1.099 ± 0.0014 (19/43 bands >4σ); ℓ≳1700 construction-dependent |
| tSZ, full sky (100–5000) | mean y 0.990 (v3fs 1.035); C_ℓ 1.79–3.23 (40–150), 1.31–1.79 (150–270), 1.17–1.33 (270–500), 1.07–1.16 (500–1000), 1.03–1.06 (1000–1400), 1.01–1.03 (1400–4100). v2: 0.80–0.88 (270–1400) | — | **Fail** formally: 1.024 ± 0.0017 (3/40 bands >4σ, max 6.6σ); v3fs 1.025 (5/40), v2 0.924 (25/40). Dividing out the painter difference on identical halos (below) leaves ≈1.12 / 1.04 / 1.02 / 1.00 at ℓ ≈ 330 / 530 / 850 / 1240 |
| kSZ total, field + Wc halos (500–8000) | 0.87–1.27 (20–50); 1.06–1.30 (86–400); 1.19–1.26 (450–1100); 1.07–1.19 (1100–2500); 0.95–1.00 (≥3500) | — | Construction-dependent: 1.105 ± 0.0035 (v3fs 1.112, v2 1.122) |
| kSZ field (30–500) | — | vs exact-LOS Doppler + linear OV (v3 field maps, reused in v4fs): 1.01–1.10 (65–200), 0.98–1.02 (200–400), 0.96–1.00 (430–930). v2: 1.13–1.21 (71–203) | Not scored (no σ). The v2 ℓ = 70–200 excess is gone; **cause not attributed** (§7) |
| ISW (10–250) | 0.68–1.03 (27–300; 0.82–0.98 at 100–300); ≥1.09 at ℓ≥434 rising to 214 at ℓ≈930 (painting noise floor, outside the declared range; 2.7× lower than v2 at ℓ≈800) | — | **Pass (weak)**: 0.851 ± 0.073 (2.0σ, max band 1.6σ; v3 field maps), with the unofficial σ computed for this report; v2 was 0.907 ± 0.075. The range is dominated by ℓ<30, with large two-realization cosmic variance and shared octant structure: individual bands at ℓ<15 run 0.56–1.37. Ours vs linear-theory ISW would be a stronger anchor (not done) |
| CIB mean intensity | 1.226/1.227/1.227/1.228/1.230/1.232 at 100–857 GHz. Baseline: Websky's halos through the same XGPaint model give 1.18, so the offset is mostly from the painter/model | — | Normalization **open** |
| CIB clustered (100–3000) | band averages 1.444 / 1.440 / 1.433 / 1.412 / 1.366 / 1.228 at 100–857 GHz (v4fs/v3fs 0.95–1.00 per band) | — | **Fail** on amplitude: 1.41–1.44 at 100–353 GHz (cf. the mean-intensity ratio squared, 1.23² = 1.51), falling with frequency to 1.23 at 857 GHz |
| CIB Poisson (3000–8000) | 1.18–1.24 (100), 1.15–1.21 (143), 1.10–1.17 (217), 0.99–1.08 (353), 0.84–0.94 (545), 0.66–0.75 (857 GHz) | — | No band >4σ, but the 400 mJy flux cut is applied to neither map. **Open** |
| CIB decoherence (Websky Table 1, 150<ℓ<1000) | All 10 pairs within ≤0.022 of the released maps (post-2022 CIB update), e.g. 545×857 0.955 vs 0.951; 143×857 0.812 vs 0.790 sits at the limit | — | Agrees with the released maps (no σ). Same-family CIB model, so not an independent check. Stein+ printed Table 1 differs (§6.2) |
| κ × CIB (100–2000) | 1.19–1.38 (≈√(1.07·1.5)) | — | Not scored; tracks the CIB normalization |
| y × CIB | 1.12–1.37 (ℓ 190–2940) | — | Not scored; tracks the CIB normalization and the tSZ tail |
| κ × y | 0.89–1.24 (ℓ 108–2940) | — | Not scored |
| Seams | κ-field and ISW: 0 exactly-zero pixels (frozen: 190k / 1.5M); CIB 130k = empty-pixel floor (v3fs 120k; frozen: 1.9M). near-plane enhancement (0–0.05° vs 2–5°): y 1.283 (Websky 1.18), CIB545 1.104 (Websky 1.078), i.e. the same shape, somewhat stronger in the innermost bin | — | Fixed (our own bug fix) |

Sources: `validation/paper/V4_RESULTS_2026-10-01.md`, `validation/paper/results/{auto,theory,cross}_v4fs_*.txt` (comparisons: `*_v3fs_*`, `*_v2_*`; `V3_RESULTS_2026-09-28.md`, `V2_RESULTS_2026-09-27.md`), and `validation/paper_theory/THEORY_ANCHORS_2026-09-26.md`.

**v2 → v3 → v3fs → v4fs, what changed.**
- nbuff 16→25 restores the tail; this raised tSZ from 0.80–0.82 to 1.03–1.09 at ℓ = 480–1000.
- Full-sky AM + `tail_N` then raised tSZ further by ×1.09 (ℓ≈480), ×1.04 (1000), ×1.01 (3000) and ×1.14–1.35 at ℓ = 100–300.
- κ, kSZ and CIB band averages change by ≤0.5% between v3 and v3fs (single CIB 857 GHz bands by up to 1.6%).
- The field maps (κ_field, kSZ field, ISW) changed between v2 and v3; two changes are confounded there: the Ω_m(a) coefficient and the tile layout (n 416→434, nbuff 16→25).
- v3fs → v4fs (volume reduction + per-tile threshold): mean y −4.3% (sub-1e12 halos removed); tSZ C_ℓ ×0.97 at ℓ≈300, ×0.99 at 800, ×1.00 at 3000; κ ×0.974–1.002; kSZ ×0.99–1.00; CIB ×0.95–1.00 per band. Field maps unchanged.

**Supporting checks.**
- **tSZ painter on Websky's own halos.** We painted `halos_10x10.pksc` with our painter and compared it with the released `tsz_2048` on a 4.5° disc.
  - C_ℓ ratio: 1.21/1.16/1.11/1.06/1.04/1.06/1.13/1.25 at ℓ=204–3666.
  - Cross-correlation ≥0.992; mean y 1.006.
  - Caveat: this is a single 4.5° disc. tSZ there is dominated by a handful of clusters, so the variance is large and non-Gaussian, and no σ is attached.
  - Within that caveat, dividing the v4fs full-sky tSZ ratio at the nearest matching bands (ℓ = 326 / 526 / 848 / 1241: 1.294 / 1.157 / 1.076 / 1.039) by these painter factors leaves a catalog-side ratio of about 1.12 / 1.04 / 1.02 / 1.00 (v3fs: 1.14 / 1.06 / 1.02 / 1.00). So the painter accounts for all of the excess at ℓ≳800, most of it at ℓ≈500, and about half at ℓ≈300. Our sky is a different random realization from Websky's, so the few dominant clusters do not correspond one-to-one.
  - V2 offers a hypothesis for the painter difference: the reference table is truncated at 4 Mpc transverse. This has not been demonstrated.
- **Reference map versions.**
  - The release `UPDATES` file (mocks.cita.utoronto.ca/data/websky/v0.0/UPDATES, quoted in `V2_RESULTS_2026-09-27.md:87-95`) records:
    - 17-FEB-2022: a ~7% tSZ profile-normalization fix;
    - 18-FEB-2022: a CIB satellite update;
    - 08-APR-2022: a fix for the high-mass profile centre.
  - V2_RESULTS states that our `tsz_2048` reference is the post-fix version.
  - Checksums of the downloaded files are in `websky_ref/md5sums.txt`. The versions are otherwise established only from download dates.
  - This context applies to every tSZ and CIB comparison with the reference.
- **Per-octant tSZ** (v4fs, C(400–800)): ours 7.9–10.7e-18 against Websky 7.7–9.0e-18; v2 had 6.5–7.5e-18. At ℓ = 100–200 the octant-to-octant ratio spans 0.69–4.8 (v3fs 0.76–5.1), i.e. a handful of clusters.
- **kSZ construction.**
  - The catalog input is cross-validated: N(M_h>1e13, z) ratio 1.00.
  - A faithful Battaglia painting of *either* catalog gives a halo term ~5–10× the released decomposition at ℓ~500–1000.
  - Tier-B band means vs `ksz.fits`: W 1.84, We 1.58, Wc 1.38.
  - Shuffling velocities leaves Wc almost unchanged (1.383→1.373), so the mid-ℓ excess is 1-halo shot power.
  - The mid-ℓ halo power differs by construction. The hydro points in Stein+ Fig 6 are ℓ=3000 values (Park+18 Table 1, CSF-scaled; `2001.08787.txt:1258-1267`), so they do not constrain ℓ≈900 directly. At ℓ≈2700–3600 our total is 1.00–1.04× `ksz.fits` (v4fs). The hydro-template overlay (B5, Fig 6) is pending and is the intended arbiter.
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
> - Catalog size: we keep all ~112.5M halos/octant (≈9.0e8 full sky, v4), against Websky's ~9e8 after its cut (v3: 198.4M/octant).
> - Duplicated volume: each octant processes 2413 of 4096 tiles, so the full sky processes ~4.71 box volumes of fine tiles. 2413 is confirmed for v2–v4 (logs report "Tile k/2413").
> - Scope: Websky's 3.84 h figure is "the run time of the simulation"; the paper gives no cost for field maps or painting.
> - Node counts assume 40 cores per Niagara node.
> - **The costs below are for v4** (nbuff=25, volume reduction, per-tile threshold), the configuration that reproduces both the tail and the clustering. v3 and v2 are kept for comparison.
> - **v4 vs v3:** the per-tile threshold removes the z>0 sub-1e12 candidate flood before shell analysis (pre-merge 201M vs 338M per octant), which roughly halves the pipeline. AM moved out of the catalog job into the full-sky step.
> - **v3 vs v2 timing mixes two changes:** nbuff 16→25 and the shell early exit (off in v2, on in v3). The v2 cost with the early exit is unknown.
> - **Early-exit exactness** was checked bit-for-bit on the small example config (~610k records, job 5708923) and in a unit test (job 5711281). No production-scale early-exit-off run exists to compare against (both v3test and v3 ran with it on; `SHELL_EARLY_EXIT_2026-09-27.md`).
>
> The robust claims are zero IC storage, the time to solution on a single node, and the memory footprint *per octant run* (~40×). In aggregate over the full sky (8 concurrent octants) the memory saving is ~5×.

| Resource | Websky (Stein+2020 §4.1) | PeakPatch.jl v4 | Source |
|---|---|---|---|
| Hardware per run | 1128 Skylake cores (~28 nodes, assumed) | 1 node, 4×L40S, per octant | `2001.08787.txt:1033-1037` |
| Catalog wall | 3.84 h (full box) | 1.62–1.95 h per octant (mean 1.68; halo finding, merge, finalize, write); ≈1.95 h full sky if the 8 octants run concurrently (inferred from the longest job; in this campaign they did not, oct000 started ~12 h earlier); 13.4 h sequentially on one node. v3 (incl. per-octant AM): 2.64–2.83 h; v2: 2.97–3.21 h | sacct 5752129–5752136; v3 5716328–5716356 (step 4); v2 5696115–5696143 |
| Catalog compute | 4336 core-h (≈108 node-h, assumed) | 13.4 node-h = 53.7 L40S-GPU-h allocated (v3: 87.1, v2: 97.3) | same |
| Full product set | not stated | ≈73.5 GPU-h ≈ 18.4 node-h: 53.7 catalogs, 1.45 full-sky AM (table 0.42 + 8 applies 1.03), 1.00 halo paint, 5.40 CIB, 11.9 field maps (v3 field maps, reused). v3: 106.4 GPU-h | sacct 5752129–5752164; field maps 5716331 + 4k |
| Peak memory | 7.67 TB (whole box) | 89.7 GiB host (MaxRSS, oct000; v3 134.5) + 4×~22 GiB GPU (GPU figure from the n=414 benchmark; not measured at n=434) ≈ 0.19 TB **per node-job/octant**, ~40× less than Websky's 7.67 TB per full-box run. Aggregate for 8 concurrent octants ≈1.5 TB, ~5× less | `PERFORMANCE_2026-09-26.md` §4; sacct MaxRSS |
| Stored ICs | 5.9 TB | 0 B (regenerated from the seed) | same |
| Catalog on disk | 33 GB (10 floats, ~9e8 halos) | 14.8 GB/octant (33 floats, 112.5M halos); raw+AM × 8 ≈ 221 GiB. An 11-float format would be ~5 GB/octant | `ls` |

**v4 per-octant split** (logs `v4_cat_oct*_5752*.{out,err}`):
- pipeline (halo finding) 52.4–54.3 min, against 99.6–101.5 in v3: the per-tile threshold removes ~40% of the candidates, almost all at z>0 and below 1e12, before the shell analysis;
- merge (exclusion + volume reduction), finalize and the 15 GB write: the rest of the 1.62–1.95 h job, ≈45–63 min including startup; not timed per stage;
- full-sky AM afterwards: table 0.42 h, then 0.12–0.14 h per octant.

**v3 per-octant split** (logs `v3_cat_oct*_57163*.{out,err}`): pipeline 99.6–101.5 min (v2 121.5–126.3: the early exit more than pays for the larger nbuff); merge, finalize and write 42–53 min; per-octant AM 15–16 min.

**v2 oct000 stage profile** (max per worker, 4 workers; log `v2_cat_oct000_5696115.err`):

| Stage | Time | Share |
|---|---|---|
| Shell analysis | 5308 s | 80.4% |
| Peak find | 427 s | 6.5% |
| Residual generation | 321 s | 4.9% |
| δ interpolation | 249 s | 3.8% |
| Stage total | 6607 s | — |

- The earlier ~10–15 min merge estimate was a projection from a microbenchmark: on 9.81M halos, 346 s at NC=256 against 13 s with 5 Mpc/h cells.
- A v3/v4 stage profile has not been extracted.

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
- n=416 (v2) and n=434 (v3, v4) were not measured.

**nbuff=25 without the early exit** (job 5704175, cancelled): about 20 s/tile, projecting to about 13 h per octant. With the exact early exit (commit a887192) the whole v3 campaign ran at 2.64–2.83 h per octant.

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
| nbuff=16 shell cap | No M>7.5e14; tSZ deficit | nbuff=25 (v3 campaign); full-sky tail matches Tinker |
| Per-octant AM table (fsky=1/8) + cap above each z-bin's top halo | Each z-bin's rank-1 halo pinned at M(N=1) ≈ 1.7e15; N(>2e15) = 0 against 12.5 | One full-sky table from all octants + `tail_N` (942930a); v3fs: 9 vs 12.5 |
| Ω_m(a) in the 2LPT coefficient used a³ where the standard form has a⁻³ (7 sites) | ψ2 coefficient ×1.01–1.04 at z = 0.5–4.5; measured catalog effect 0.007 Mpc/h rms in position and 0.15–0.34% in velocity | `Cosmology.omega_m_a` + regression test (058207a) |
| Stale conclusion that the nhunt clamp was harmless | It held only while the field was wrong | Superseded by V2_RESULTS |
| Merge skipped the Fortran volume-reduction pass (comment wrongly said Fortran had it commented out) | ~7% higher halo bias than Websky (ξ(3–15) W/ours 0.877 ± 0.029, 4σ); 16–17% excess close pairs | `merge_catalog(...; volume_reduction=true)` (4ba9ef8), on in v4: ξ 0.958 ± 0.026 |
| Constant lightcone peak threshold fsc_of_z(0), where Fortran uses fsc_of_z(z_tile) (2026-04 note misread the Fortran) | Extra sub-1e12 candidates at z>0 (338M vs 201M pre-merge per octant), ~2× pipeline cost; ≤0.5% above 1e12 | `[run] peak_threshold_per_tile` (4ba9ef8), on in v4 |

**Stale text in the paper draft and comments** (to correct):
- `paper/sections/algorithm.tex:46-47` says partial overlaps reduce R_TH. That is the Fortran algorithm and is what v4 does; keep it, and add the per-tile lightcone threshold.
- `paper/sections/maps.tex:54` and `intro.tex:58` say halo profiles are painted with XGPaint. Only CIB uses XGPaint.
- `paper/sections/performance.tex:8-9` still quotes the stale "~19 min on 4×L40S".
- `src/FieldMap.jl` comments are stale on three points: the +3/7 sign, ISW painted at displaced positions, and D/a.
- `MultiResolution.jl:591-593` says shell analysis stays on the CPU; it runs on the GPU.
- `TILING_INVARIANCE_2026-09-26.md:112` still calls the canonical merge order "not applied".
- `PERFORMANCE_2026-09-26.md` gives "20 filters"; production uses 23.

### 6.2 Observed differences with the reference (neutral; not yet validated by the PI)

| Observation | Evidence | Source |
|---|---|---|
| Our tSZ painter on Websky's own halos exceeds the released `tsz_2048` by 4–25% in C_ℓ, with cross-correlation ≥0.992 and mean y 1.006 | Hypothesis: the reference table is truncated at 4 Mpc transverse (not demonstrated) | `V2_RESULTS_2026-09-27.md` |
| In the local build, pks2map defaults to mmin=2.5e10 with no 1e13 or 0.5′ cut. On halos_10x10 it gives ~25× the released kSZ halo component at ℓ~900 | The settings used for the released map are not documented in the repo code. Observed in the local (hand-patched) clone; whether the released maps were produced with this code path or these settings is unknown; not validated by the PI | `KSZ_COMPOSITE_2026-07-19.md:167-203` |
| In the local build, the output depends on the z-cut at fixed input. As read, the compaction at `pks2map.f90:182-186` copies posxyz and rth but not vrad | In the local build, zmax 4.5 vs 6.0 on identical input gives 40% lower map RMS and 20–33% lower C_ℓ at ℓ≤1237. Observed in the local (hand-patched) clone; whether the released maps were produced with this code path, binary or zmax setting is unknown; not validated by the PI | same |
| κ: `kap_lt4.5`/Halofit(z<4.5) is 0.935–0.994 at ℓ=150–1650 and 1.00 at ℓ≈2900; ours/Halofit is 0.995–1.053 (v4fs) | The ~20% smallest-scale suppression described in Stein+ (`2001.08787.txt:1604-1606`, Fig 8) lies at ℓ beyond our theory file (which ends at ℓ=2894). It also refers to the total κ, including the z>4.5 Gaussian component, not `kap_lt4.5` alone. Not tested here. Tier-B suggested sub-pixel mass loss (untested) | `results/theory_v4fs_kappa.txt` |
| CIB decoherence printed in Stein+ Table 1 is lower than we measure on the released maps (e.g. 857×545 0.933 printed vs 0.951 measured on the released maps (ours 0.955)) | Candidate: the 18-FEB-2022 CIB satellite update (untested). The garbled text extraction needs checking against the PDF | `2001.08787.txt:1464-1520`; `results/cross_v4fs_cib_cib.txt` |
---

## 7 Open items and next steps

| # | Item | Status / next action |
|---|---|---|
| 1 | **Attribute the field-map change v2 → v3** | The field-kSZ excess at ℓ = 70–200 (1.13–1.21 → 1.01–1.10) and the ISW changes coincide with two changes: the Ω_m(a) coefficient and the tile layout (nbuff 16→25). Deciding test: one oct000 field map with one change and not the other (~1.5 h on 1 GPU). Until then, do not claim which change resolved it |
| 2 | **Clustering excess: resolved in v4** | Cause: the missing Fortran volume-reduction pass (`CLUSTERING_EXCESS_2026-09-28.md`). v4fs Tier-A: ξ(3–15) 0.958 ± 0.026, ξ(1–3) 0.98–1.01, b(M) 0.97–1.07, v12 0.963 ± 0.040, all ≤1.6σ (§3.1); the ~4% ξ(3–15) residual left by the reduction A/B persists at 1.6σ. Remaining: whether the reduction rule is more *physical* than exclusion-only needs a matched N-body comparison (bias ±4%, ξ(1–3) ±20% depend on it) |
| 3 | `tail_N` choice | 10 is a judgement call (below ~10 halos per slice the relative Poisson scatter is ≳30%); not derived in any note, and the sensitivity to tail_N (e.g. 3 / 10 / 30) has not been measured |
| 4 | CIB normalization (1.23× in mean) | Decide whether to rescale L0 to the released mean or quote the ratio |
| 5 | CIB 400 mJy flux cut | Apply to both maps before scoring the Poisson regime |
| 6 | κ offset (1.06–1.13× `kap_lt4.5`) | Unexplained; shared between the two maps relative to Halofit |
| 7 | gradrf_x/y/z NaN (~0.005% of halos; fields 31–33) | Masses and positions unaffected. Root cause untraced (zero-radius mrf shell suspected) |
| 8 | tSZ low-ℓ excess (1.31–3.2 at ℓ = 40–270; single bands larger at ℓ < 40) | Dominated by a handful of clusters (octant ratios 0.69–4.8 at ℓ = 100–200); σ is large. A matched-mass comparison (drop or cap the top clusters in both maps) would separate realization variance from a systematic |
| 9 | Theory-test and cross-spectrum σ | Add σ to `theory_*`/`cross_*`; implement the Gaussian-only rule at ℓ≲30 in `compare_auto.jl`; re-derive verdicts with the paper tooling |
| 10 | 10-particle pre-AM cut | Undecided. v4 already removes most sub-1e12 halos (counts 0.07–0.12 of v3 at 1e11–1e12), which moved mean y 1.035 → 0.990; an explicit Websky-style cut would make the mass range exactly like for like |
| 11 | AM near the floor (flat extrapolation below Mmin=5e11) | Measured: AM/Tinker at M>1e12 is 0.92–1.07 in v4fs (v3fs 0.94–1.17), non-monotonic in z (§3.1); ≤3% above 3e12 |
| 12 | Observational and hydro-template overlays; lensed CMB | Not produced (Figs 5–7, 9–11) |
| 13 | Box replication (A6); v12 tail beyond 30 Mpc/h; dN/dz M>1e13 at z 2.25–3 (ours 8–12% high, up to 2.8σ, unchanged v3fs→v4fs) | Untested or uninvestigated |
| 14 | Measurements still missing | Per-stage timing of merge, finalize and write; a v3/v4 stage profile; GPU memory at n=434; v2 cost with the shell early exit; an early-exit-off production run; worker-count reproducibility with the canonical merge; MPI global-FFT vs multires agreement; whether the chi(z) handling in finalize matches merge_pkvd exactly; 2LPT sign convention of raw ψ2 vs Fortran |
| 15 | **Provenance / storage** | All current catalogs and maps are on **purgeable scratch**: v4 raw + v4fs AM catalogs (`scratch/websky_6144/catalogs_v4/`, 16 × 13.8 GiB) and v4fs maps (`scratch/websky_6144/v4fs/`, via /project symlinks `{halomaps,cibmaps,fullsky}_v4fs` and the symlink dir `fieldmaps_v4fs`, which points at `fieldmaps_v3` on /project). Superseded v3/v3fs catalogs and maps (~470 GiB) are also on scratch. Copy v4fs to /project (4.8/5.0 TiB used) after deleting superseded data there; the deletion is the PI's decision. Reference map versions are established only from download dates and `websky_ref/md5sums.txt` |
| 16 | Paper-draft corrections | See §6.1 |
| 17 | Raw (pre-AM) MF vs Tinker per z-bin (A4) | Done (§3.1). v4: 0.90–0.92 at 1e13 (z<1.5), 0.99–1.02 at 1e14 (z<1), 0.62–0.71 at 1e12 |
| 18 | Other A4 theory anchors: IC field P(k) vs input; ψ1/ψ2/velocity power vs linear/2LPT theory | Not run (the velocity P(k) test would also help with item 1) |
| 19 | Completeness vs Stein+ §4.2 `halo_mass_completion.txt` | The file is on disk in `websky_ref/`; comparison not done |
| 20 | A8 convergence appendix: nbuff (full 16/25 campaigns exist), AM variant (per-octant vs full-sky vs tail_N), merge variant (exclusion-only v3fs vs reduction v4fs), filter-bank count/spacing, cellsize (the factor-h runs as a controlled test), catalog coarse_factor | Not assembled |
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
- `validation/websky_6144/production_v4/` (README, `config_v4_oct*.toml`, `run_catalog.slurm`, `submit_all.sh`)
- `validation/websky_6144/production_v3/` (README, `config_v3_oct*.toml`, `run_catalog.slurm`, `run_fieldmap.slurm`, `submit_all.sh`, `submit_am_fullsky.sh`)
- `validation/websky_6144/production_v2/` (README, `config_v2_oct*.toml`, `config_v3test_oct000.toml`)
- `validation/websky_6144/production/` (`paint_octant.jl`, `halo_profiles.jl`, `paint_cib.jl`, `run_paint.slurm`, `run_cib_prod.slurm`)
- `validation/websky_6144/run_gpu_octant.jl`, `run_fieldmap_octant.jl`, `apply_abundance_match.jl`, `apply_abundance_match_fullsky.jl`

**Notes, current.**
- Units and results: `CONVENTIONS.md`; `validation/paper/V4_RESULTS_2026-10-01.md`; `validation/paper/CLUSTERING_EXCESS_2026-09-28.md`; `validation/paper/TIERA_V3FS_2026-09-28.md`; `validation/paper/V3_RESULTS_2026-09-28.md`; `validation/paper/V2_RESULTS_2026-09-27.md`; `validation/paper/FULLSKY_COMPARISON_2026-09-26.md`
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

**Results.** `validation/paper/results/{auto,theory,cross}_v4fs_*.txt`, `tierA_v4fs.txt`, `tails_v4fs.log` (current); `*_v3fs_*`, `*_v3_*`, `*_v2_*`, `tierA_v3fs.txt` for comparison; `merge_ab.txt`, `threshold_ab.txt`, `tierA_diag_v3fs.txt` (the clustering investigation)

**Reference.**
- Stein+2020: `~/work/peakpatch/2001.08787.txt`
- Released maps: `/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/`
- Local Fortran clone: `~/work/peakpatch` (hand-patched near the cosmology block)

**Job IDs.**

| Purpose | Jobs |
|---|---|
| v4 catalogs | 5752129–5752136 |
| v4fs AM table / applies / paint / CIB | 5752137 / 5752138 + 3k / 5752139 + 3k / 5752140 + 3k |
| v4fs analysis / tails / Tier-A | 5752162 / 5752163 / 5752164 |
| Clustering investigation: Tier-A diagnostics / reduction A/B / threshold A/B | 5736375 / 5736480 / 5750925 |
| v3 catalog+AM | 5716328 + 4k (k = 0…7) |
| v3 halo paint / CIB / field maps | 5716329 / 5716330 / 5716331 + 4k |
| v3 analysis / tail check | 5728240 / 5728241 |
| v3fs AM table / applies | 5731483 / 5731484 + 3k |
| v3fs paint / CIB | 5731485 + 3k / 5731486 + 3k |
| v3fs analysis / tail check | 5731508 / 5735812 |
| Ω_m(a) fix size (oct000) | 5720824 |
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
| Tier-A v3fs rerun (8 caps + raw MF vs Tinker) | 5735978 |
| Tier-A catalog history (superseded) | 3968163 (finecell raw), 4033130 (finecell AM), 4387351 (cf32) |
| Frozen campaign (superseded) | 4438182–4438197, 5627673; AMv2 5692696–5692703 |

**Commits.** 91c4145 (fsc_of_z), 70e0beb (chi), 95f04a1 and 4ae7837 (D/a, 2LPT sign), 8aa4e9d (MPI), 9f92778 (CLI finalize, AM top halo), 285c72c (ISW at Lagrangian positions), 1c11cb3 (v2: periodic cores, compensation, deterministic merge), a887192 (shell early exit, nbuff diagnosis), 058207a (Ω_m(a) in the 2LPT coefficient), b6f8064 (production-v3), 942930a (full-sky AM + tail_N), 4ba9ef8 (volume reduction and per-tile threshold options), 756fa4b (production-v4).
