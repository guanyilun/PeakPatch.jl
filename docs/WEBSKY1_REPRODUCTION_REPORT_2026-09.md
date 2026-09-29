# PeakPatch.jl reproduction of WebSky 1.0: fidelity and performance

*Status as of 2026-09-28 (updated to the production-v3 / v3fs campaign). Branch `websky2-resolution`. Audience: internal (PI) first; after PI review, Bond, Carlson, Stein et al.; input to the methods paper (`paper/`).*

> **Distribution gate: internal draft.** §6.2 and the reference-side parts of §4 (the tSZ-painter, kSZ-construction and CIB Table-1 checks) describe observed differences with the Websky reference code and maps. They must **not** go to the Websky authors, or anyone outside the group, until the PI has personally validated them and explicitly approved sharing. An external version should drop §6.2 and the reference-side remarks in §4.

All lengths are in Mpc/h, and masses in M☉/h unless stated. The Fortran peak-patch code and the public Websky products use Mpc. Conversion: `boxsize[Mpc/h] = L[Mpc]·h`. For Websky, 7700 Mpc × 0.68 = 5236 Mpc/h (`CONVENTIONS.md`).

**Result tiers used below.**
- **Current: "v3fs".** The production-v3 catalogs (8 octants, nbuff=25, shell early exit, Ω_m(a) fix; `validation/websky_6144/production_v3/`), abundance-matched with **one full-sky table and `tail_N`=10**, as in the reference procedure (§2.2). The field maps are the v3 field maps, which don't depend on AM.
- **Superseded, quoted for comparison:**
  - **v3** (per-octant AM): identical catalogs, but each z-bin's top halo is pinned at M(N=1), so N(>2e15) = 0.
  - **v2** (nbuff=16): the shell-search cap removes every M>7.5e14 halo.
  - the frozen 2026-07 campaign, which had seam gaps, and runs from before the 2026-06 fixes.
- The 2026-06/07 Tier-A numbers on the single validation octant `oct000_finecell_AM` are superseded by the v3fs rerun (§3), except the Fortran finder test.
- **Do not quote:** the "19 min/octant" and "~0.9M halos/octant, 120× fewer than Websky" figures. They come from the 2026-04 configuration, before the bug fixes (`validation/performance/PERFORMANCE_2026-09-26.md` §1).

**Framing.** Where we list differences from the reference Fortran code or the released maps, these are *observed differences with evidence*. The PI has not yet personally validated them, and none has been communicated upstream. Sending this draft to Websky authors would count as communicating them, so the distribution gate above applies.

**Glossary.**
- **Tier-A:** catalog-level statistics (MF, dN/dz, σ_vr, ξ, b(M), v12) against Websky's 10°×10° halo patch; first measured in 2026-06/07 on one validation octant, rerun on v3fs with 8 caps on 2026-09-28.
- **Tier-B:** map-level (painted) statistics from the same 2026-07 period, on oct000 caps or octants.
- **Plan codes** (A0–A8, B1–B6, Q1–Q3): items in `docs/paper_comparison_plan_2026-09.md`. For example, A2 = rerun Tier-A/B on current catalogs; A3 = error model and pass criteria; A4 = theory anchors independent of Websky; A6 = full-sky assembly and replication checks; A8 = convergence appendix; B5 = kSZ hydro-template overlay.
- **kSZ halo constructions** (`KSZ_COMPOSITE_2026-07-19.md`):
  - **W:** uncompensated per-halo Battaglia τ painting;
  - **We:** the Websky-literal variant, with a Δ=3-mean compensation sphere and gas mass f_b·m_h;
  - **Wc:** a zero-net compensated variant, a full-τ sphere with R=4 r_vir.
- **cf:** the coarse_factor of the multi-resolution split (coarse grid M = N/cf).
- **octZYX:** the octant label; the three bits give the observer sign per axis.

**Sources.** Facts sourced only to `memory/*.md` notes are not in the repo and are marked *(memory only; unverified in repo)*.

**Ranges.** Ratio ranges in §1, §3.2 and §4 are re-read from `validation/paper/results/*_v3fs_*.txt` (analysis job 5731508). v2 values are from `*_v2_*.txt`. Where these differ from the prose notes (`V2_RESULTS_2026-09-27.md`, `V3_RESULTS_2026-09-28.md`), the results files take precedence.

---

## 1 Summary

**Verdict.**
- **Halo catalog** (v3fs, Tier-A rerun 2026-09-28 with errors: the Websky patch against the same cap placed in each of our 8 octants; `validation/paper/TIERA_V3FS_2026-09-28.md`).
  - **Counts and kinematics agree:** N(>M|z) after AM 0.97–1.05 (z<3, M≥3e12; largely enforced by AM), dN/dz 0.95–1.08 (z<4.25), σ_vr 0.95–0.99, all within ≲2σ of our cap-to-cap scatter.
  - **Ours is more clustered than Websky:** ξ(3–15 Mpc/h) W/ours = 0.877 ± 0.029 (4.0σ), i.e. bias ≈ 1.07× Websky's; b(M) 4–7% higher below 3e13 (2.5σ, 1.7σ) and equal above; |v12| ~7% higher (1.7σ); ξ at 1–3 Mpc/h 1.2–1.33× (4.5–6.7σ). The 2026-07 single-cap "bias ratio 1.06" is confirmed, now with an error bar; the cause is not identified (§7).
  - **Raw mass function vs Tinker (plan A4, full sky):** within 2% at 1e13 for z<1.5; +15–35% above 1e14 at z<1; 0.55–0.8 at z>2 and 0.6–0.9 at 1–3e12. AM matches Tinker to ≤3% above 3e12, but is off by up to 17% at the 1e12 floor (shared with the Websky catalog).
  - The independent finder evidence is the Fortran/Julia kept-peak test (limited scope, below).
- **High-mass tail (v3fs, full sky, z<1):** N(>1e15) = 416 against Tinker's 412.3, and N(>2e15) = 9 against 12.5. v2 had 0 halos above 1e15 because of the nbuff=16 shell cap.
- **Seams:** fixed since v2. The signal profile across the octant planes tracks Websky's. This is our own fix, not a reproduced Websky statistic.
- **CIB decoherence:** agrees with the released maps (post-2022 CIB update) to ≤0.025 for all 10 pairs, with no σ attached (143×857 differs by 0.024). This uses a same-family CIB model (XGPaint `CIB_Planck2013`), so it is not an independent check. The printed Stein+ Table 1 differs from both (§6.2).
- **ISW at low ℓ:** consistent within large cosmic variance, using an unofficial σ computed for this report (§4).
- **Field kSZ vs linear theory:** 1.01–1.10 at ℓ = 65–200 and 0.96–1.02 up to ℓ≈930 in v3. The v2 excess of 1.13–1.21 at ℓ = 70–200 is gone, but **which v2→v3 change removed it has not been established** (§7).
- **tSZ, now close but not formally passing:**
  - band average over ℓ = 100–5000 is 1.025 ± 0.002 (5/40 bands >4σ), against 0.924 in v2;
  - C_ℓ is 1.01–1.06 at ℓ = 1000–4100 and 1.08–1.17 at ℓ = 500–1000, rising to 1.2–1.4 at ℓ = 270–500 and above 1.8 at ℓ < 150, where a few clusters dominate and σ is large;
  - our painter alone is 1.04–1.16 above the reference on identical halos (§4): that accounts for all of the excess at ℓ≳800, most at ℓ≈500, and about half at ℓ≈300.
- **Not yet reproduced to the pre-declared criteria:**
  - κ (6–14% high against `kap_lt4.5`; ours is 1.00–1.06 × Halofit, the reference 0.94–0.99);
  - CIB, which differs by a mostly painter-driven amplitude: mean 1.23–1.24, and 1.18 on Websky's own halos through the same model. There is also a frequency-dependent residual at 857 GHz and in the Poisson regime, where the flux cut has not yet been applied.
- **Differs by construction:** kSZ halo power at mid-ℓ, where the result depends on how the halos are painted.
- **Cost (v3, the like-for-like configuration):**
  - Catalog+AM takes 2.64–2.83 h per octant on one 4×L40S node; the full sky takes about 2.8 h on 8 nodes, or 21.8 h on one node.
  - The full product set is about 106 L40S-GPU-h, with no stored initial conditions (ICs).
  - Peak memory is ~30× lower **per octant run** (0.24 TB, against 7.67 TB for Websky's full-box run). For the full sky with 8 concurrent octants it is ≈1.9 TB, ~4× lower in aggregate.
  - This is **not** a like-for-like hardware comparison (§5).

**Headline results.**
- **Phase space: the AM-insensitive statistics** (v3fs, Websky patch vs our 8 caps; §3.1):
  - σ_vr 0.95–0.99 (≤1.4σ);
  - ξ(3–15 Mpc/h) 0.877 ± 0.029 (ours more clustered, 4σ);
  - b(M) W/ours 0.94 / 0.96 / 0.96 / 1.00 in 4 bins;
  - v12 W/ours 0.927 ± 0.041 over 5–30 Mpc/h.
- **Finder equivalence (limited scope).**
  - Kept-peak densities are 3.925e-3 (Julia) and 3.946e-3 (Fortran) per (Mpc/h)³, where both codes use the same numeric cellsize, 1.2533. Fortran/Julia = 1.005 (`validation/websky_6144/FORTRAN_COMPARISON_2026-06-14.md`).
  - Scope: a 320.8 Mpc/h box; a z=0 snapshot (no lightcone); the old 1.25 Mpc/h cell; pre-v2 code. It is statistical (different seeds) and compares kept-peak density only, not the mass function or matched halos.
- **Mass function after AM** (largely enforced by AM, see the Verdict). v3fs ours/Websky N(>M|z) is 0.97–1.05 for M≥3e12 at z<3 (§3.1); at the 1.2e12 floor for z 3–4.5 ours has 9% more halos.
- **Seams fixed since v2.** κ-field and ISW maps have 0 exactly-zero pixels, against 190k and 1.5M in the frozen campaign. The signal profile across the octant planes tracks Websky (`validation/paper/V2_RESULTS_2026-09-27.md`).
- **tSZ: the v2 deficit had two causes in our pipeline, both now fixed.**
  - v2 full sky gave 0.80–0.88× Websky at ℓ=270–1400, because the nbuff=16 shell-search cap (R_TH ≤ 12.8 Mpc/h, raw M ≤ 7.4e14) left N(>1e15, z<1) = 0 per octant against 51.5 for Tinker.
  - v3 (nbuff=25) restores the raw tail: 43–60 per octant above 1e15 at z<1. Its per-octant AM then capped each z-bin at about 1.7e15, so N(>2e15) = 0 in every octant.
  - v3fs (full-sky AM + `tail_N`) gives 416 above 1e15 and 9 above 2e15 over the full sky, against 412 and 12.5 for Tinker. tSZ C_ℓ / Websky is 1.37 (ℓ≈300, σ 0.09), 1.19 (480, σ 0.05), 1.10 (770), 1.06 (1030), 1.02 (2000); mean y 1.035.
- **CIB** (band averages unchanged from v2 to ≤0.5%; single bands move by up to ~2.5%).
  - Mean intensity is 1.23–1.24× the released maps at every frequency. Websky's own halos through the same XGPaint model give 1.18, so most of the offset comes from the painter/model and not from the catalog.
  - Clustered C_ℓ is 1.43–1.46 at 100–353 GHz, falling to 1.40 (545 GHz) and 1.26 (857 GHz). Poisson-regime ratios are 1.22→0.70 across frequency. Both are frequency-dependent rather than a single (1.23)² factor.
  - Decoherence agrees with the released maps to ≤0.025 for all 10 frequency pairs (no σ).
- **Performance** (v3, one 4×L40S node per octant; nbuff=25 with the exact shell early exit):
  - catalog+AM takes 2.64–2.83 h per octant (v2 at nbuff=16 without the early exit: 2.97–3.21 h);
  - host memory is about 134.5 GiB. GPU memory is about 22 GiB per GPU in the n=414 benchmark; it was not measured at n=434 (v3);
  - peak memory is ~30× lower per octant run (0.24 TB vs 7.67 TB for the full-box Websky run). The full sky at 8 concurrent octants is ≈1.9 TB, ~4× lower in aggregate;
  - IC storage is 0 B, against 5.9 TB for Websky.

---

## 2 What we implemented, and how it differs from the Fortran/Websky production

### 2.1 Production-v3 configuration

Sources: `validation/websky_6144/production_v3/README.md` and `config_v3_oct000.toml`. v3 is v2 with nbuff 16→25, the GPU shell early exit and the Ω_m(a) fix.

| Item | Value |
|---|---|
| Seed, N, box, cell | 12345; 6144; 5236 Mpc/h (= 7700 Mpc); 0.852213 Mpc/h |
| Tiling | ntile=16, tile mesh n=434, nbuff=25, nsub=384, `periodic_cores=true` (N = nsub·ntile). `[grid] boxsize` is the **per-tile** box, 369.86 Mpc/h. nhunt = 24 (R_TH cap 20.4 Mpc/h) |
| Coarse grid | coarse_factor=32 (M=512, block 12), `coarse_compensation=true` |
| Lightcone | ievol=1, z_out=0, z_max=4.5; observer at ±2618 Mpc/h per axis (the core corner, set by the octZYX bits) |
| Physics | ilpt=2, ioutshear=1, wsmooth=1, rmax2rs=0.0 |
| Cosmology | Ωm=0.31, Ωb=0.049, ΩΛ=0.69, h=0.68; P(k) from `pk_websky.dat` (CAMB, σ8=0.8100) |
| Filters | `filters_websky_finecell.dat`: 23 top-hats, Rf = 1.406–30.44 Mpc/h, ratio 1.15 |
| Collapse table | `HomelTab_websky.dat`, 50×20×20 |
| AM (v3fs) | one Tinker08 M200m table from all 8 raw octants (fsky=1), 10⁴ mass × 46 z bins over z<4.5, `tail_N`=10 (`apply_abundance_match_fullsky.jl`) |

### 2.2 Differences from the Fortran code and the Websky production

| Aspect | Fortran / Websky (Stein+2020) | PeakPatch.jl v3fs | Source |
|---|---|---|---|
| Units | Mpc | Mpc/h everywhere, because chi(z) uses H0=100h | `CONVENTIONS.md` |
| White noise | 48-bit LCG, filled sequentially per MPI slab; the realization depends on the task count | Counter-based Threefry2x (20 rounds) with key (seed,0), plus Box-Muller. Each cell depends only on (seed, global index). A bit-exact LCG port is also available | `src/InitialConditions/RandomField.jl:36-103`; `LCG.jl` |
| Initial conditions | Global slab-decomposed FFTW3, 7 values per site, 5.9 TB stored | No global fine FFT. The coarse M³ periodic FFT is spliced to a per-tile isolated residual convolution with Catmull-Rom interpolation (MUSIC/Hoffman-Ribak style); D/T compensation in v2; nothing stored | `docs/multi_resolution_fft.md:60-104`; `src/MultiResolution.jl:584-735` |
| Buffer | ≈ the largest Lagrangian halo radius (~40 Mpc) | nbuff=25 cells = 21.3 Mpc/h (v2: 16). The shell search is capped at nhunt = min(nbuff−1, ⌊1.75 Rf_max/a⌋), the same formula as `hpkvd.f90:207-209`; nhunt=24 is the GPU kernel's shared-memory maximum. Largest raw R_TH over the 8 octants is 17.2–20.4 Mpc/h (one halo at the cap, in oct101) | `src/MultiResolution.jl:630`; `V3_RESULTS_2026-09-28.md` |
| Filter bank | Paper: Rf,min = 2 a_latt; Rf,max = 36 Mpc (24.5 Mpc/h); count and spacing not stated. `filter_gen.py` default: 1.65 a_latt, spacing 1.15 (the same as ours). Production bank not located; the evidence conflicts (paper 2.0; the dense-bank residual suggests ≈1.2) | Rf,min = 1.65 a_latt, Rf,max = 30.44 Mpc/h, 23 filters, spacing 1.15 | `data/filters_websky_finecell.dat`; `2001.08787.txt:300-303`; `LOW_MASS_COMPLETENESS_2026-06-15.md:66-80` |
| Peak threshold on the lightcone | Constant fcrit = fsc_of_z(0) | Same. Shell analysis uses the per-peak z_pk (1+z_pk, fsc_of_z(z_pk)) | `src/MultiResolution.jl:934-937,1017-1045` |
| Ω_m(a) in the 2LPT coefficient | Standard (Ωm/a³)/(Ωm/a³+ΩΛ) (`hpkvd.f90:762,875`) | Same since 058207a (v2 used a³ in place of a⁻³; §6.1). v3 includes the fix | `src/Cosmology/Cosmology.jl` |
| Merge | Exclusion with a fixed NC=256 hash; volume reduction commented out | Same survivors, with a data-sized hash (`clamp(round(cbrt(n)),16,1024)`), a total order (−R_TH,x,y,z), and volume reduction disabled | `src/Merger/Merger.jl:32-80`; `Exclusion.jl:22,57-64` |
| Eulerian positions and velocities | merge_pkvd | Port: x = q+ψ1+ψ2, v = a·100·E·f·(ψ1+2ψ2). Velocity epoch from one chi→z evaluation at the Eulerian distance, with no iteration | `src/Merger/Merger.jl:106-145` |
| AM | Tinker M200m, Δz=0.1 bins, 10⁴ mass bins, bilinear interpolation; **one full-sky table**; above a slice's top halo the table is left at identity, so the bilinear lookup blends the top halos with their raw masses; only halos with >10 particles before AM | Same table method (46 z-bins, Δz≈0.098, over 5e11–1e16), **one full-sky table from all 8 octants**. Above the last mass edge with ≥10 halos per z-bin, the fractional correction is frozen (`tail_N`), which is monotone and keeps ranks (the reference's identity default can invert them). **No particle cut** | `src/AbundanceMatch/AbundanceMatch.jl`; `~/work/peakpatch/python/catalogue_tools/abundance_match/make_abundancematch_table.py`; `2001.08787.txt:1044-1052` |
| Catalog format | 10 floats, ~9e8 halos, 33 GB | 33 Float32 (extended pksc), 198.4M halos/octant (928.4M full sky in 5e11–1e16 at z<4.5), 26.2 GB/octant | `src/Merger/Merger.jl:84-162` |
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
- The Tier-A statistics were rerun on the **v3fs** catalogs on 2026-09-28 (`validation/paper/tierA_v3fs.jl`, job 5735978; note `TIERA_V3FS_2026-09-28.md`). The comparison target is still Websky's public 10°×10° patch (`halos_10x10.pksc`), because the full Websky `halos.pksc` is not available locally.
- The cap is the one inscribed in that patch (half-angle 5.00°, Ω = 0.02395 sr), with the 2026-07 shell (1600–2000 Mpc/h, z≈0.6–0.8), mass bins and windows. The identical cap is placed at the centre of each of our 8 octants; "ours" is the mean ± std over the 8 caps, and z-scores use std·√(1+1/8).
- Post-AM N(>M) agreement is largely enforced by AM, because both catalogs are matched to Tinker08 M200m, so it is not independent evidence. The AM-insensitive statistics are ξ, b(M), σ_vr and v12, the raw MF against Tinker, and the Fortran finder test.
- Our realization is not Websky's, and the Websky patch is one cap, so each comparison is "is Websky's patch a plausible draw from our cap distribution".
- Stein+2020 itself shows no halo mass function, dN/dz or clustering plot, so these statistics are additional validation.

### 3.1 Scorecard (W/ours = Websky patch / mean of our 8 caps, v3fs)

| Statistic | W/ours | z | Verdict |
|---|---|---|---|
| N(>M) after AM, M ≥ 3e12, z < 3 | 0.97–1.05 (1.08 for >1e13 at z<0.5) | within ±1.7 | agree (largely enforced by AM) |
| N(>1e14), z < 1.5 | 0.96 / 1.12 / 1.03 | ≤1.6 | agree |
| N(>1.2e12), z 3–4.5 | 0.919 ± 0.007 | −10 | ours 9% more at the floor at high z (near-floor AM) |
| dN/dz, M>3e12, Δz=0.25, z 0–4.25 | 0.95–1.08 | within ±1.9 | agree |
| dN/dz, M>3e12, z 4.25–4.5 | 1.27 ± 0.10 | +2.5 | ours low in the last bin, at our z_max=4.5 edge (Websky reaches 4.6) |
| σ_vr, 1e12–1e13 / 1e13–3.2e14 | 0.989, 0.988 / 0.980, 0.969, 0.954 (± 0.02–0.03) | −0.5 / −0.8 to −1.4 | agree |
| ξ(r), M>5e12, mean 3–15 Mpc/h | **0.877 ± 0.029** | **−4.0** | **ours more clustered** (bias ≈ 1.068 ± 0.018) |
| ξ(r) at 1–3 Mpc/h | 0.75–0.84 | −4.5 to −6.7 | ours 1.2–1.33× |
| ξ(r) at 19–41 Mpc/h | 0.72–0.93 (± 0.1–0.3) | ≤0.9 | consistent (noisy) |
| b(M) ∝ √ξ̄(6–18 Mpc/h), 5e12–1.3e13 / 1.3–3.2e13 / 3.2–7.9e13 / 0.8–2.5e14 | 0.936 / 0.958 / 0.958 / 0.996 | −2.5 / −1.7 / −0.9 / −0.1 | ours 4–7% higher below 3e13 |
| v12(r), M>1e13, mean 5–30 Mpc/h | 0.927 ± 0.041 | +1.7 | ours ~7% stronger infall (≤2σ) |
| Finder, matched cell 1.2533 Mpc/h, z=0 | Kept-peak density 3.946e-3 (Fortran) vs 3.925e-3 (Julia) per (Mpc/h)³, Fortran/Julia 1.005. Scope: 320.8 Mpc/h box, z=0 snapshot, pre-v2 code, different seeds, density only | — | `FORTRAN_COMPARISON_2026-06-14.md:50-60` |
| Lightcone vs snapshot (z<0.19) | 0.978; GPU = CPU (raw 69008 = 69008; merged 25898 vs 25897) | — | same |

**Raw (pre-AM) mass function vs Tinker08, full sky** (plan A4; v3 raw, v3fs AM):

| M threshold | raw / Tinker at z = 0–0.25, 0.25–0.5, 0.5–1, 1–1.5, 1.5–2, 2–3, 3–4.5 | AM / Tinker |
|---|---|---|
| > 1e12 | 0.64, 0.66, 0.69, 0.72, 0.73, 0.71, 0.57 | 0.94, 1.03, 1.17, 1.00, 1.01, 1.00, 1.07 |
| > 3e12 | 0.87, 0.89, 0.91, 0.92, 0.90, 0.84, 0.66 | 0.97–1.00 |
| > 1e13 | 1.01, 1.02, 1.02, 1.00, 0.95, 0.83, 0.56 | 0.98–1.00 |
| > 1e14 | 1.20, 1.19, 1.15, 1.03, 0.80, 0.55, — | 0.96–1.00 (z<3) |
| > 5e14 | 1.35, 1.26, 1.07, 0.72 (z<1.5) | 0.97–1.01 |
| > 1e15 | 1.00, 1.05, 0.92 (z<1) | 0.96–1.05 |

**Comparison with the 2026-07 single-cap Tier-A** (superseded `oct000_finecell_AM`, one cap, no error bar; `BIAS_PAIRWISE_V12_2026-07-16.md`): its "ξ bias ratio 1.06 (3–15 Mpc/h)" is confirmed (now 1.068 ± 0.018). Its b(M) ratios (ours/W 1.024 / 1.006 / 0.980 / 0.952) and v12 (0.997) are **not** reproduced by the 8-cap mean: that single cap was one draw within the scatter.

### 3.2 Current production catalogs (v3fs, with v2 and v3 for comparison)

Sources: `validation/paper/V3_RESULTS_2026-09-28.md` (tail jobs 5728241 for v3, 5735812 for v3fs) and `V2_RESULTS_2026-09-27.md`.

**Per octant.**

| Quantity | v2 (nbuff=16) | v3 / v3fs (nbuff=25) |
|---|---|---|
| Halos / octant (post-merge) | 198.4M | 198.39–198.52M (same catalogs for v3 and v3fs) |
| Raw N(>1.7e12) | 4.09e7 | 4.09–4.10e7 |
| AM N(>1e14) | 6.25e4 (forced equal in every octant) | v3: 6.25e4 (forced equal); v3fs: 6.22–6.27e4 (varies between octants, as sample variance should) |
| Raw / AM max M | 7.4e14 / 8.5e14 | raw 1.8–3.1e15; AM v3 1.66–1.71e15; AM v3fs 1.84–2.91e15 |

**Full sky, z<1 (sum over 8 octants).**

| | raw v3 | v2 AM | v3 AM (per octant) | **v3fs AM** | Tinker |
|---|---|---|---|---|---|
| N(>3e14) | 34910 | ~28920 | 28874 | 28861 | 29255 |
| N(>5e14) | 7118 | ~6000 | ~5937 | 5942 | 5878 |
| N(>1e15) | 415 | **0** | 408 | **416** | 412.3 |
| N(>2e15) | 11 | 0 | **0** | **9** | 12.5 |

(v2 AM sums are 8× the oct000 value.)

By z-bin, v3fs N(>1e15) is 122 / 206 / 88 at z = 0–0.25 / 0.25–0.5 / 0.5–1, against Tinker 127.7 / 200.3 / 84.2; N(>2e15) is 5 / 3 / 1 against 6.0 / 5.5 / 1.0.

Notes:
- The raw v3 tail is itself close to Tinker above 1e15 (415 vs 412), so AM moves the top only a little. It moves the 3e14–5e14 range down by ~17%.
- Raw R_TH tapers smoothly below the 20.4 Mpc/h cap in all 8 octants; one halo in oct101 sits at the cap.
- nbuff=25 gives nhunt=24, the largest the GPU shell kernel allows (482 shells against `_MAX_SHELLS_GPU` = 512).
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
- Estimator A3: full-sky C_ℓ of the 8-octant v3fs maps against the released Websky v0.0 maps, with no mask, pixel windows divided out, and Δℓ/ℓ≈0.1.
- Analysis job: 5731508 (v3fs; ours vs Websky, cross-spectra and theory). v2: 5703683 / 5703908; v3: 5728240.
- Reference files are in `/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/`.
- Pass criterion, declared in advance (`docs/paper_comparison_plan_2026-09.md:77-105`): the band-averaged |r−1| < 3σ over the declared ℓ range, and no single band above 4σ.
- Construction-dependent differences (kSZ halo term at mid-ℓ, κ at ℓ≳1700) are not scored pass/fail.

> **Caveats on the "status" column.**
> - The σ values were computed for this report as an inverse-variance mean of the `sigma_ratio` column in `validation/paper/results/auto_v3fs_*.txt` (the same computation reproduces the v2 values quoted previously). No repo document states formal verdicts; they should be re-derived with the paper tooling.
> - `compare_auto.jl:28-34` uses max(Gaussian, jackknife) at all ℓ, which is not the plan's Gaussian-only rule at ℓ≲30.
> - The theory and cross-spectrum files carry no σ, so those tests are not scored. The analysis log does print a fractional Gaussian σ for the κ and field-kSZ theory tables (e.g. field kSZ at ℓ=78: 1.103 ± 0.043, 2.4σ), so the field-kSZ test could be scored.

| Product (declared ℓ range) | ours / Websky (v3fs) | ours / theory | Status vs pre-declared criterion |
|---|---|---|---|
| κ, z<4.5 vs `kap_lt4.5` (30–3000) | 1.06–1.08 (ℓ 270–430); 1.10–1.12 (850–1400); 1.13–1.14 (2600–3000); 1.02–1.03 (7600–8100). Unchanged from v2 (v3fs/v2 = 0.994–1.003 at ℓ ≥ 100; 0.987–1.016 including ℓ < 100) | Halofit: 1.00–1.06 (ℓ 150–1650), 1.14 (2894). Field κ vs linear Limber: 0.98–1.04 (150–1650). Websky vs Halofit: 0.935–0.994 (150–1650), so the offset is shared between the two maps | **Fail** at ℓ 30–1700: 1.103 ± 0.0014 (19/43 bands >4σ); ℓ≳1700 construction-dependent |
| tSZ, full sky (100–5000) | mean y 1.035; C_ℓ 1.86–3.38 (40–150), 1.35–1.83 (150–270), 1.19–1.37 (270–500), 1.08–1.17 (500–1000), 1.04–1.06 (1000–1400), 1.01–1.03 (1400–4100). v2: 0.80–0.88 (270–1400) | — | **Fail** formally: 1.025 ± 0.0017 (5/40 bands >4σ, max 7.1σ); v2 was 0.924 ± 0.0013 (25/40). Dividing out the painter difference on identical halos (below) leaves ≈1.14 / 1.06 / 1.02 / 1.00 at ℓ ≈ 330 / 530 / 850 / 1240 |
| kSZ total, field + Wc halos (500–8000) | 0.87–1.27 (20–50); 1.06–1.31 (86–400); 1.20–1.27 (450–1100); 1.08–1.20 (1100–2500); 0.96–1.01 (≥3500) | — | Construction-dependent: 1.112 ± 0.0035 (v2 1.122) |
| kSZ field (30–500) | — | vs exact-LOS Doppler + linear OV: 1.01–1.10 (65–200), 0.98–1.02 (200–400), 0.96–1.00 (430–930). v2: 1.13–1.21 (71–203) | Not scored (no σ). The v2 ℓ = 70–200 excess is gone; **cause not attributed** (§7) |
| ISW (10–250) | 0.68–1.03 (27–300; 0.82–0.98 at 100–300); ≥1.09 at ℓ≥434 rising to 214 at ℓ≈930 (painting noise floor, outside the declared range; 2.7× lower than v2 at ℓ≈800) | — | **Pass (weak)**: 0.851 ± 0.073 (2.0σ, max band 1.6σ), with the unofficial σ computed for this report; v2 was 0.907 ± 0.075. The range is dominated by ℓ<30, with large two-realization cosmic variance and shared octant structure: individual bands at ℓ<15 run 0.56–1.37. Ours vs linear-theory ISW would be a stronger anchor (not done) |
| CIB mean intensity | 1.229/1.229/1.230/1.233/1.237/1.245 at 100–857 GHz. Baseline: Websky's halos through the same XGPaint model give 1.18, so the offset is mostly from the painter/model | — | Normalization **open** |
| CIB clustered (100–3000) | band averages 1.46 (100 GHz) … 1.26 (857 GHz), as in v2 (v3fs/v2 within 0.5%) | — | **Fail** on amplitude: ≈1.23² at 100–353 GHz, plus a frequency-dependent fall to 1.26 at 857 GHz |
| CIB Poisson (3000–8000) | 1.18–1.26 (100), 1.15–1.24 (143), 1.10–1.19 (217), 0.99–1.10 (353), 0.84–0.95 (545), 0.66–0.75 (857 GHz) | — | No band >4σ, but the 400 mJy flux cut is applied to neither map. **Open** |
| CIB decoherence (Websky Table 1, 150<ℓ<1000) | All 10 pairs within ≤0.025 of the released maps (post-2022 CIB update), e.g. 545×857 0.956 vs 0.951; 143×857 0.814 vs 0.790 sits at the limit | — | Agrees with the released maps (no σ). Same-family CIB model, so not an independent check. Stein+ printed Table 1 differs (§6.2) |
| κ × CIB (100–2000) | 1.20–1.41 (≈√(1.07·1.5)) | — | Not scored; tracks the CIB normalization |
| y × CIB | 1.18–1.45 (ℓ 190–2940) | — | Not scored; tracks the CIB normalization and the tSZ tail |
| κ × y | 0.93–1.29 (ℓ 108–2940) | — | Not scored |
| Seams | κ-field and ISW: 0 exactly-zero pixels (frozen: 190k / 1.5M); CIB 119k = empty-pixel floor (frozen: 1.9M). near-plane enhancement (0–0.05° vs 2–5°): y 1.245 (Websky 1.18), CIB545 1.109 (Websky 1.078), i.e. the same shape, somewhat stronger in the innermost bin | — | Fixed (our own bug fix) |

Sources: `validation/paper/V3_RESULTS_2026-09-28.md`, `validation/paper/results/{auto,theory,cross}_v3fs_*.txt` (v2 comparison: `*_v2_*.txt`, `V2_RESULTS_2026-09-27.md`), and `validation/paper_theory/THEORY_ANCHORS_2026-09-26.md`.

**v2 → v3 → v3fs, what changed.**
- nbuff 16→25 restores the tail; this raised tSZ from 0.80–0.82 to 1.03–1.09 at ℓ = 480–1000.
- Full-sky AM + `tail_N` then raised tSZ further by ×1.09 (ℓ≈480), ×1.04 (1000), ×1.01 (3000) and ×1.14–1.35 at ℓ = 100–300.
- κ, kSZ and CIB band averages change by ≤0.5% between v3 and v3fs (single CIB 857 GHz bands by up to 1.6%).
- The field maps (κ_field, kSZ field, ISW) changed between v2 and v3; two changes are confounded there: the Ω_m(a) coefficient and the tile layout (n 416→434, nbuff 16→25).

**Supporting checks.**
- **tSZ painter on Websky's own halos.** We painted `halos_10x10.pksc` with our painter and compared it with the released `tsz_2048` on a 4.5° disc.
  - C_ℓ ratio: 1.21/1.16/1.11/1.06/1.04/1.06/1.13/1.25 at ℓ=204–3666.
  - Cross-correlation ≥0.992; mean y 1.006.
  - Caveat: this is a single 4.5° disc. tSZ there is dominated by a handful of clusters, so the variance is large and non-Gaussian, and no σ is attached.
  - Within that caveat, dividing the v3fs full-sky tSZ ratio at the nearest matching bands (ℓ = 326 / 526 / 848 / 1241: 1.326 / 1.173 / 1.083 / 1.043) by these painter factors leaves a catalog-side ratio of about 1.14 / 1.06 / 1.02 / 1.00. So the painter accounts for all of the excess at ℓ≳800, most of it at ℓ≈500, and about half at ℓ≈300. Our sky is a different random realization from Websky's, so the few dominant clusters do not correspond one-to-one.
  - V2 offers a hypothesis for the painter difference: the reference table is truncated at 4 Mpc transverse. This has not been demonstrated.
- **Reference map versions.**
  - The release `UPDATES` file (mocks.cita.utoronto.ca/data/websky/v0.0/UPDATES, quoted in `V2_RESULTS_2026-09-27.md:87-95`) records:
    - 17-FEB-2022: a ~7% tSZ profile-normalization fix;
    - 18-FEB-2022: a CIB satellite update;
    - 08-APR-2022: a fix for the high-mass profile centre.
  - V2_RESULTS states that our `tsz_2048` reference is the post-fix version.
  - Checksums of the downloaded files are in `websky_ref/md5sums.txt`. The versions are otherwise established only from download dates.
  - This context applies to every tSZ and CIB comparison with the reference.
- **Per-octant tSZ** (v3fs, C(400–800)): ours 8.0–10.8e-18 against Websky 7.7–9.0e-18; v2 had 6.5–7.5e-18. At ℓ = 100–200 the octant-to-octant ratio spans 0.76–5.1, i.e. a handful of clusters.
- **kSZ construction.**
  - The catalog input is cross-validated: N(M_h>1e13, z) ratio 1.00.
  - A faithful Battaglia painting of *either* catalog gives a halo term ~5–10× the released decomposition at ℓ~500–1000.
  - Tier-B band means vs `ksz.fits`: W 1.84, We 1.58, Wc 1.38.
  - Shuffling velocities leaves Wc almost unchanged (1.383→1.373), so the mid-ℓ excess is 1-halo shot power.
  - The mid-ℓ halo power differs by construction. The hydro points in Stein+ Fig 6 are ℓ=3000 values (Park+18 Table 1, CSF-scaled; `2001.08787.txt:1258-1267`), so they do not constrain ℓ≈900 directly. At ℓ≈2700–3500 our total is 1.01–1.06× `ksz.fits` (v3fs). The hydro-template overlay (B5, Fig 6) is pending and is the intended arbiter.
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
> - Duplicated volume: each octant processes 2413 of 4096 tiles, so the full sky processes ~4.71 box volumes of fine tiles. 2413 is confirmed for v2 and v3 (logs report "Tile k/2413").
> - Scope: Websky's 3.84 h figure is "the run time of the simulation"; the paper gives no cost for field maps or painting.
> - Node counts assume 40 cores per Niagara node.
> - **The costs below are for v3** (nbuff=25, no mass cap below ~3e15), the configuration that reproduces the tail. v2 (nbuff=16, mass-capped) is kept for comparison.
> - **v3 vs v2 timing mixes two changes:** nbuff 16→25 and the shell early exit (off in v2, on in v3). The v2 cost with the early exit is unknown.
> - **Early-exit exactness** was checked bit-for-bit on the small example config (~610k records, job 5708923) and in a unit test (job 5711281). No production-scale early-exit-off run exists to compare against (both v3test and v3 ran with it on; `SHELL_EARLY_EXIT_2026-09-27.md`).
>
> The robust claims are zero IC storage, the time to solution on a single node, and the memory footprint *per octant run* (~30×). In aggregate over the full sky the memory saving is ~4×.

| Resource | Websky (Stein+2020 §4.1) | PeakPatch.jl v3 | Source |
|---|---|---|---|
| Hardware per run | 1128 Skylake cores (~28 nodes, assumed) | 1 node, 4×L40S, per octant | `2001.08787.txt:1033-1037` |
| Catalog wall | 3.84 h (full box) | 2.64–2.83 h per octant for catalog+AM (mean 2.72); 2.83 h full sky with 8 nodes in parallel; 21.8 h sequentially on one node. v2: 2.97–3.21 h (mean 3.04) | sacct 5716328–5716356 (step 4); v2 5696115–5696143 |
| Catalog compute | 4336 core-h (≈108 node-h, assumed) | 21.8 node-h = 87.1 L40S-GPU-h allocated. v2: 24.3 node-h = 97.3 GPU-h | same |
| Full product set | not stated | 106.4 GPU-h ≈ 26.6 node-h: 87.1 catalog+AM, 11.9 field maps, 1.36 halo paint, 5.95 CIB. Full-sky AM (v3fs) adds 0.64 (table) + 1.67 (8 applies) GPU-h and a repaint (1.45 + 6.0 GPU-h) | sacct |
| Peak memory | 7.67 TB (whole box) | ≈134.5 GiB host (MaxRSS, oct000) + 4×~22 GiB GPU (GPU figure from the n=414 benchmark; not measured at n=434) ≈ 0.24 TB **per node-job/octant**, ~30× less than Websky's 7.67 TB per full-box run. Aggregate for 8 concurrent octants ≈1.9 TB, ~4× less | `PERFORMANCE_2026-09-26.md` §4; sacct MaxRSS |
| Stored ICs | 5.9 TB | 0 B (regenerated from the seed) | same |
| Catalog on disk | 33 GB (10 floats, ~9e8 halos) | 26.2 GB/octant (33 floats); raw+AM × 8 = 390 GiB. An 11-float format would be ~8.6 GB/octant | `ls` |

**v3 per-octant split** (logs `v3_cat_oct*_57163*.{out,err}`):
- pipeline (halo finding) 99.6–101.5 min, against 121.5–126.3 min in v2: the early exit more than pays for the larger nbuff;
- merge, finalize and the 26 GB write: 42–53 min (halo-finding step 142–153 min minus the pipeline), as in v2; not timed per stage;
- AM (per octant, inside the job) 15–16 min.

**v2 oct000 stage profile** (max per worker, 4 workers; log `v2_cat_oct000_5696115.err`):

| Stage | Time | Share |
|---|---|---|
| Shell analysis | 5308 s | 80.4% |
| Peak find | 427 s | 6.5% |
| Residual generation | 321 s | 4.9% |
| δ interpolation | 249 s | 3.8% |
| Stage total | 6607 s | — |

- The earlier ~10–15 min merge estimate was a projection from a microbenchmark: on 9.81M halos, 346 s at NC=256 against 13 s with 5 Mpc/h cells.
- A v3 stage profile has not been extracted.

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
- n=416 (v2) and n=434 (v3) were not measured.

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
| κ: `kap_lt4.5`/Halofit(z<4.5) is 0.935–0.994 at ℓ=150–1650 and 1.00 at ℓ≈2900; ours/Halofit is 1.00–1.06 (v3fs 0.996–1.058) | The ~20% smallest-scale suppression described in Stein+ (`2001.08787.txt:1604-1606`, Fig 8) lies at ℓ beyond our theory file (which ends at ℓ=2894). It also refers to the total κ, including the z>4.5 Gaussian component, not `kap_lt4.5` alone. Not tested here. Tier-B suggested sub-pixel mass loss (untested) | `results/theory_v3fs_kappa.txt` |
| CIB decoherence printed in Stein+ Table 1 is lower than we measure on the released maps (e.g. 857×545 0.933 printed vs 0.951 measured) | Candidate: the 18-FEB-2022 CIB satellite update (untested). The garbled text extraction needs checking against the PDF | `2001.08787.txt:1464-1520`; `results/cross_v3fs_cib_cib.txt` |
---

## 7 Open items and next steps

| # | Item | Status / next action |
|---|---|---|
| 1 | **Attribute the field-map change v2 → v3** | The field-kSZ excess at ℓ = 70–200 (1.13–1.21 → 1.01–1.10) and the ISW changes coincide with two changes: the Ω_m(a) coefficient and the tile layout (nbuff 16→25). Deciding test: one oct000 field map with one change and not the other (~1.5 h on 1 GPU). Until then, do not claim which change resolved it |
| 2 | **Ours ~7% more biased than Websky** (ξ 3–15 Mpc/h 1.14×, 4σ; b(M) +4–7% below 3e13; 1–3 Mpc/h 1.2–1.33×) | Tier-A rerun done (§3.1). Cause not identified; candidates: mass-rank scatter at fixed AM mass, merge/exclusion, Eulerian displacement convention, tile layout. A test: b(M) at fixed *raw* mass rank, and ξ of the raw catalog |
| 3 | `tail_N` choice | 10 is a judgement call (below ~10 halos per slice the relative Poisson scatter is ≳30%); not derived in any note, and the sensitivity to tail_N (e.g. 3 / 10 / 30) has not been measured |
| 4 | CIB normalization (1.23× in mean) | Decide whether to rescale L0 to the released mean or quote the ratio |
| 5 | CIB 400 mJy flux cut | Apply to both maps before scoring the Poisson regime |
| 6 | κ offset (1.06–1.14× `kap_lt4.5`) | Unexplained; shared between the two maps relative to Halofit |
| 7 | gradrf_x/y/z NaN (~0.005% of halos; fields 31–33) | Masses and positions unaffected. Root cause untraced (zero-radius mrf shell suspected) |
| 8 | tSZ low-ℓ excess (1.35–3.4 at ℓ = 40–270; single bands up to 8 at ℓ < 40) | Dominated by a handful of clusters (octant ratios 0.76–5.1 at ℓ = 100–200); σ is large. A matched-mass comparison (drop or cap the top clusters in both maps) would separate realization variance from a systematic |
| 9 | Theory-test and cross-spectrum σ | Add σ to `theory_*`/`cross_*`; implement the Gaussian-only rule at ℓ≲30 in `compare_auto.jl`; re-derive verdicts with the paper tooling |
| 10 | 10-particle pre-AM cut | Undecided. The 0.6% sub-cell (M<1e11) fraction was measured only on the example config |
| 11 | AM near the floor (flat extrapolation below Mmin=5e11) | Measured: AM/Tinker at M>1e12 is 0.94–1.17 and non-monotonic in z (§3.1), shared with the Websky catalog; ≤3% above 3e12 |
| 12 | Observational and hydro-template overlays; lensed CMB | Not produced (Figs 5–7, 9–11) |
| 13 | Box replication (A6); v12 tail beyond 30 Mpc/h; 1–3 Mpc/h ξ excess; shared b(M) shortfall vs Tinker10 | Untested or uninvestigated |
| 14 | Measurements still missing | Per-stage timing of merge, finalize and write; a v3 stage profile; GPU memory at n=434; v2 cost with the shell early exit; an early-exit-off production run; worker-count reproducibility with the canonical merge; MPI global-FFT vs multires agreement; whether the chi(z) handling in finalize matches merge_pkvd exactly; 2LPT sign convention of raw ψ2 vs Fortran |
| 15 | **Provenance / storage** | The v3 raw catalogs, the v3fs AM catalogs (`scratch/websky_6144/catalogs_v3/`, 16 × 24.4 GiB) and all v3fs maps (`scratch/websky_6144/v3fs/`, 79 GiB, reached via /project symlinks `{halomaps,cibmaps,fullsky}_v3fs`) are on **purgeable scratch**; copy to /project (4.8/5.0 TiB used) or delete superseded catalogs before a purge. The reference map versions are established only from download dates and `websky_ref/md5sums.txt` (see the UPDATES entries in §4) |
| 16 | Paper-draft corrections | See §6.1 |
| 17 | Raw (pre-AM) MF vs Tinker per z-bin (A4) | Done (§3.1): within 2% at 1e13 (z<1.5); +15–35% above 1e14 at z<1; 0.55–0.8 at z>2 |
| 18 | Other A4 theory anchors: IC field P(k) vs input; ψ1/ψ2/velocity power vs linear/2LPT theory | Not run (the velocity P(k) test would also help with item 1) |
| 19 | Completeness vs Stein+ §4.2 `halo_mass_completion.txt` | The file is on disk in `websky_ref/`; comparison not done |
| 20 | A8 convergence appendix: nbuff (full 16/25 campaigns exist), AM variant (per-octant vs full-sky vs tail_N), filter-bank count/spacing, cellsize (the factor-h runs as a controlled test), catalog coarse_factor | Not assembled |
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
- `validation/websky_6144/production_v3/` (README, `config_v3_oct*.toml`, `run_catalog.slurm`, `run_fieldmap.slurm`, `submit_all.sh`, `submit_am_fullsky.sh`)
- `validation/websky_6144/production_v2/` (README, `config_v2_oct*.toml`, `config_v3test_oct000.toml`)
- `validation/websky_6144/production/` (`paint_octant.jl`, `halo_profiles.jl`, `paint_cib.jl`, `run_paint.slurm`, `run_cib_prod.slurm`)
- `validation/websky_6144/run_gpu_octant.jl`, `run_fieldmap_octant.jl`, `apply_abundance_match.jl`, `apply_abundance_match_fullsky.jl`

**Notes, current.**
- Units and results: `CONVENTIONS.md`; `validation/paper/V3_RESULTS_2026-09-28.md`; `validation/paper/V2_RESULTS_2026-09-27.md`; `validation/paper/FULLSKY_COMPARISON_2026-09-26.md`
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

**Results.** `validation/paper/results/{auto,theory,cross}_v3fs_*.txt` (current); `*_v3_*`, `*_v2_*` for comparison

**Reference.**
- Stein+2020: `~/work/peakpatch/2001.08787.txt`
- Released maps: `/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/`
- Local Fortran clone: `~/work/peakpatch` (hand-patched near the cosmology block)

**Job IDs.**

| Purpose | Jobs |
|---|---|
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

**Commits.** 91c4145 (fsc_of_z), 70e0beb (chi), 95f04a1 and 4ae7837 (D/a, 2LPT sign), 8aa4e9d (MPI), 9f92778 (CLI finalize, AM top halo), 285c72c (ISW at Lagrangian positions), 1c11cb3 (v2: periodic cores, compensation, deterministic merge), a887192 (shell early exit, nbuff diagnosis), 058207a (Ω_m(a) in the 2LPT coefficient), b6f8064 (production-v3), 942930a (full-sky AM + tail_N).
