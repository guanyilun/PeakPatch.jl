# Theory anchors for the paper (2026-09-26)

Plan item A4 (theory part) of `docs/paper_comparison_plan_2026-09.md`. These are curves
independent of Websky, used to show that agreement with Websky is not a shared systematic.
They also reproduce the theory lines of Websky Figs 6 and 8.

## What and how

`theory_anchors.py` (CAMB 1.5.0, local CVMFS wheel, venv `/home/yguan/scratch/theory_env`;
the old `camb_env` was gutted by the scratch purge). Conventions match `src/FieldMap.jl`:
- lengths in Mpc/h;
- flat ΛCDM Ωm=0.31, ΩΛ=0.69, no radiation, for χ(z), D and f;
- P0 = `pk_websky.dat`×(2π)³, i.e. the exact IC spectrum;
- dτ/dχ = σ_T n_e0 x_e(z)(1+z)², with f_e=0.9 inside n_e0, Y=0.245, He²⁺ at z<3 and He⁺ above.

CAMB with the `generate_pk_camb.py` setup reproduces the IC file to within 0.14% (`pk_check.txt`).

| file | content |
|---|---|
| `results/kappa_limber.txt` | C_ℓ^κκ Limber (Websky eq 3.21), ℓ=2–10⁴, for linear / Halofit-Takahashi / HMcode2020 × {0→z*, 0→4.5, 4.5→z*}; χ*=CAMB χ(1089)=9390.7; plus linear-with-our-P0 and a no-radiation χ*=9451 variant |
| `results/ksz_doppler.txt` | linear Doppler kSZ D_ℓ [μK²]: Websky eq 4.3 boundary term (z_max 4.5, 4.6); **exact linear LOS integral** with a sharp z=4.5 cut, ℓ≤1000; boundary + Limber bulk term (approximation) |
| `results/ksz_ov.txt` | linear Ostriker–Vishniac (eq 4.4), Limber, z<4.5 |
| `measure_fieldmap_cl.jl` → `results/fieldmap_cl_prod_octXXX.txt` | absolute pseudo-C_ℓ of the frozen-campaign field maps (apodized octant, Nside 4096→2048) |
| `compare_fieldmaps_theory.py` → `results/fieldmaps_vs_theory_prod.txt` | field maps vs theory, with χ*=9656 as painted |
| `estimator_transfer.jl` → `results/estimator_transfer.txt` | Gaussian-sky test of the pseudo-C_ℓ estimator |

**Derivation notes.**
- Exact linear Doppler: Δ_ℓ(k) = ∫dχ G(χ) j_ℓ′(kχ)/k with G = (dτ/dχ)·aHfD/c, and
  C_ℓ = (2/π)∫dk k² P0 Δ_ℓ². It is computed by parts; checked against the direct j′ form at ℓ=10 (0.1%).
  Integrating by parts, the boundary term is exactly eq 4.3. The bulk term includes the He step at z=3.
- OV: P_q⊥ = (aHf)²D⁴ ∫d³k′/(2π)³ P0(|k−k′|)P0(k′)(1−μ²)k(k−2k′μ)/(k′²|k−k′|²). This is Websky's eq 4.4
  integrand, which I re-derived (both contractions). In Limber, C_ℓ = ½∫dχ/χ² (dτ/dχ)² P_q⊥(ℓ/χ)/c².

## Numbers

**κ (full 0→z*, Halofit):**
- C_100 = 1.58e-7; HMcode/Halofit = 1.00 at ℓ=100 and 0.96 at ℓ=3000–5000.
- Halofit/linear = 1.07 at ℓ=500, 1.20 at 1000, 1.89 at 3000.
- The z>4.5 fraction is 19% at ℓ=100 and 32–33% at ℓ=500–1000. This matches the 15–34% measured from the
  released kap_gt4.5 split (`COMPOSITE_KAPPA_2026-07-18.md`), an independent check of the Gaussian-tail split.
- The kernel's χ* choice (9391 / 9451 / 9656) matters at ≲1%.

**Doppler kSZ (z_max=4.5):**
- The eq 4.3 boundary term gives D_ℓ = 4.3 / 3.6 / 2.2 / 1.3 μK² at ℓ = 2 / 13 / 44 / 82.
  This matches Websky's "large-scale power of a few μK²".
- z_max 4.6 vs 4.5 raises it by 5–7%. The memo figure of "+10% in u0²" was an overestimate.
- The **exact LOS** result is *lower* than the boundary term: 1.26 / 2.34 / 1.84 / 1.14 μK² at the same ℓ.
  The ratio full/eq 4.3 is 0.29 at ℓ=2, 0.65 at ℓ=13, 0.83 at ℓ=44, 0.89 at ℓ=82, and 0.93–0.99 for ℓ≳150.
  So eq 4.3 is a good approximation only at ℓ≳100. Boundary plus bulk-Limber overshoots
  (1.34× at ℓ=100, 1.11× at 300), because the cross term it drops is negative.

**Linear OV (z<4.5):** D_ℓ = 0.09 / 0.35 / 0.80 / 0.87 / 0.78 / 0.60 μK² at
ℓ = 100 / 300 / 1000 / 3000 / 5000 / 10⁴. It overtakes the Doppler term at ℓ≈230.

## Sanity: frozen-campaign field maps vs theory (all 8 octants)

`measure_fieldmap_cl.jl` computes an apodized **geometric** octant mask from the pixel directions
(octZYX bit 1 → that axis looks toward −, the same convention as `validation/paper/spectra.jl`).
Power outside the geometric octant is 0.35–2.5e-9 of the total for every octant.
A first pass used a `map != 0` mask instead: near-observer cells displaced across the octant planes
paint a few pixels outside the octant, and that mask included them. This gave κ/linear up to 38× at
ℓ≈20 for octants 010/011/100/111. That pass is discarded.

**Estimator transfer** (`estimator_transfer.jl` → `results/estimator_transfer.txt`, 8 Gaussian skies drawn
from the theory spectra, same mask and banding):

| ℓ range | transfer |
|---|---|
| ℓ≥70 | 1.00±0.01 to 1.00±0.03 for both κ and kSZ |
| ℓ=20–60 | 0.95–1.12, with 10–38% single-sky scatter |

The estimator is unbiased where it matters.

Mean over 8 octants, with sd = octant-to-octant scatter (`results/fieldmaps_vs_theory_prod.txt`).
The octants are eight views of **one** box, so the sd understates realization variance for modes the
octants share.

| ℓ_eff | κ/linear (sd) | κ/Halofit | kSZ D_ℓ [μK²] | Doppler_lin | OV_lin | kSZ/(D+OV) (sd) |
|---|---|---|---|---|---|---|
| 40 | 1.15 (0.15) | 1.15 | 2.08 | 1.93 | 0.02 | 1.06 (0.13) |
| 73 | 1.08 (0.06) | 1.08 | 1.55 | 1.28 | 0.05 | 1.16 (0.06) |
| 109 | 1.06 (0.03) | 1.05 | 1.04 | 0.80 | 0.10 | 1.16 (0.07) |
| 163 | 1.04 (0.05) | 1.02 | 0.69 | 0.42 | 0.17 | 1.19 (0.05) |
| 245 | 1.03 (0.03) | 1.00 | 0.54 | 0.20 | 0.28 | 1.14 (0.04) |
| 365 | 1.06 (0.02) | 1.00 | 0.58 | 0.09 | 0.43 | 1.12 (0.02) |
| 546 | 1.07 (0.02) | 0.95 | 0.69 | 0.03 | 0.59 | 1.11 (0.02) |
| 815 | 1.08 (0.01) | 0.87 | 0.82 | 0.01 | 0.73 | 1.10 (0.02) |
| 1217 | 1.03 (0.01) | 0.72 | 0.91 | 0.00 | 0.84 | 1.08 (0.01) |
| 1817 | 1.02 (0.01) | 0.57 | 1.00 | 0.00 | 0.89 | 1.12 (0.01) |
| 2713 | 1.13 (0.00) | 0.48 | 1.07 | 0.00 | 0.88 | 1.21 (0.01) |

**κ field.** At ℓ=60–2700 the map is **1.02–1.09 × linear Limber**; at ℓ≤450 it sits between linear and
Halofit, as a 2LPT lattice should. Its 1-halo deficit relative to Halofit (0.48 at ℓ=2700) is expected:
the composite κ supplies that power through the painted halos. The field-matter synthesis and the
Born κ kernel are anchored to theory at the few-% level, independently of Websky.

**kSZ field.**
- At high ℓ (≳400, OV-dominated) the map is 1.08–1.21 × *linear* OV. That excess is expected,
  because the 2LPT density carries more small-scale power than linear theory (κ field/linear reaches
  1.13 at ℓ=2700). Linear OV is a lower bound.
- At ℓ=73–245 (Doppler-dominated) the map is **1.14–1.19 × linear theory, consistently in all 8 octants**
  (sd 0.04–0.07). An excess also appears in the table at ℓ=163 (1.19), where OV is only 30% of the total.
- The old explanation for the ~1.2 residual vs ksz.fits (z_max 4.6-vs-4.5 "+10% in u0²") is
  superseded on two counts. These maps are already z<4.5, and the z_max effect is only 5–7%
  (see Numbers above).
- The remaining ~15% Doppler-regime excess is **unexplained**. Candidates:
  - (a) the ψ2 part of the velocities (v ∝ D ψ1 + 2 D2 ψ2), which linear theory omits;
  - (b) shared low-k modes of the single 5236 Mpc/h box (the octant sd cannot probe them; k~ℓ/χ_max≈0.01–0.03,
    so ~10–25 k_f);
  - (c) residual coarse/fine velocity decoherence at cf32 (the July cf16→cf32 study was still moving
    at ℓ≲230: 1.29→1.13 at ℓ=226).
- The discriminating test is the 3D velocity power spectrum of the ICs at k=0.005–0.05 (ψ1 and ψ1+ψ2)
  against linear theory, at cf32 and a larger coarse_factor. This is plan item A4, "velocity power vs
  linear/2LPT theory"; it is not done here.

## Caveats

- Linear OV is a lower bound for a 2LPT field. No nonlinear-OV variant was computed.
- The theory is continuous; the periodic 5236 Mpc/h box has k_f = 1.2e-3 h/Mpc and discrete modes.
  This matters for the lowest-ℓ Doppler (k ~ ℓ/χ_max).
- Pseudo-C_ℓ uses a Gaussian pixel-window approximation (the transfer test covers mask and mode coupling,
  not the pixel window). For the paper's final numbers, use the full-sky tooling in `validation/paper/spectra.jl`.
- Websky eq 4.3 is printed with f_e to the first power. Here f_e sits inside n_e0 and enters squared.
  Read literally, eq 4.3 would be 1/0.9 = 1.11× higher.
