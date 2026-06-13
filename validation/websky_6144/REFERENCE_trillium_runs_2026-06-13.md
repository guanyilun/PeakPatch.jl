# Reference: Nate's Trillium / WebSky2 Fortran run log (recovered 2026-06-13)

Source: Google spreadsheet "Trillium Runs" (Nate's canonical Peak Patch Fortran
runs). Recovered while hunting for the authoritative parameters behind the 120x
halo deficit in our Julia Websky 6144 GPU run.

Trillium node: 192 cores (2x 96), 755 GiB RAM, max 128 nodes/run.
FFTW constraint: n_ext must be divisible by N_tasks (n_slab = n_ext/N_tasks integer).
**All distances in Mpc (NOT Mpc/h).** iLPT=2 (2LPT) throughout.

## Canonical Websky-style production runs (the rows that matter)

| Run path | box [Mpc] | n_mesh | n_buff | n_tile | n_ext | n_phys | cellsize [Mpc] | max smooth [Mpc] | buffer [Mpc] | cenz | M_min [Msun] | N_halos | elapsed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 26.04.16_WebSky2/octant_+++_test1 | 6000 | 588 | 69 | 27 | 12288 | 12150 | 0.494 | 40 | 34.1 | 3000 | 2.32e11 | **3.03e9** | 5h39m |
| 26.04.16_WebSkyCO/Doga_zoom1 | 3000 | 588 | 69 | 27 | 12288 | 12150 | 0.247 | 20 | 17.0 | 5000 | 2.89e10 | **4.01e9** | 5h10m |
| 25.11.25/ng0_nm400_nt32_NN86 | 10000 | 400 | 40 | 32 | 10320 | 10240 | 0.977 | 40 | 39.1 | 0 | 1.79e12 | **1.876e9** | 1h16m |
| 25.11.25/ng0_nm636_nt32_NN128 | 8000 | 636 | 64 | 32 | 16384 | 16256 | 0.492 | 40 | 31.5 | 0 | 2.29e11 | 1.876e9 (OOM) | — |
| 25.11.25/example_trillium_20node_z0 | 1100 | 320 | 27 | 21 | 5640 | 5586 | 0.197 | 5 | 5.32 | 0 | 1.47e10 | 2.84e8 | 27m |

Later WebSky2 octant geometries being tuned (no halo counts yet, in-progress rows):
- box 8000, n_mesh 604, n_buff 64, n_tile 32, n_ext 15360, cellsize 0.525, cenz 4000, M_min 2.79e11
- box 8645, n_mesh 604, n_buff 64, n_tile 32, n_ext 15360, cellsize 0.568, cenz 4322, M_min 3.51e11

## The two rules this table makes explicit

1. **Buffer ≈ max smoothing scale.** In EVERY production row, `n_buff * cellsize`
   (the "buffer [Mpc]" column) is set to roughly the max smoothing scale:
   34 vs 40, 17 vs 20, 39 vs 40, 31 vs 40, 5.3 vs 5. The buffer is sized so the
   largest filter fits inside the tile halo.

2. **cellsize sets M_min sets N_halos.** M_min(top-hat) ∝ cellsize³ (since
   R_f,min = 2·cellsize, Stein rule). Halo counts are in the **billions** for
   these cellsizes (0.25–1.0 Mpc).

## Why our Julia run got 120x too few halos (quantified)

### CORRECTION (2026-06-13): our resolution ALREADY matches original Websky.

An earlier draft of this section claimed the 120x was driven by grid resolution
(cellsize too coarse, "need N≈10240"). **That was wrong** — it anchored on Nate's
*newer, finer WebSky2* runs (box 6000 Mpc, cellsize ~0.49 Mpc) instead of the
*original* Stein+2020 Websky. The Websky paper (2001.08787 §4.1, lines 1035-1037)
is explicit:

> "creating periodic initial conditions for a **(7.7 Gpc)³, 6,144³ particle
>  simulation**, placing an observer at each of the eight corners ... the octant
>  method has the advantage of doubling the one-dimensional resolution ... raising
>  the volume resolution by a factor of 8."

So each octant simulation IS 6144³ at 7700 Mpc/h — **identical to our config**. The
effective 12288³ / 15.4 Gpc lightcone comes from the octant (corner-observer) trick,
not from a finer grid.

Mass floor check (h-units, ρ_m = 2.775e11·Ωm):
- particle mass m_p = ρ_m·(7700/6144)³ = **1.69e11 Msun/h**
- Websky retains halos > **10 particles** = **1.69e12 Msun/h** (paper §4.2),
  → post-abundance-match completeness ~1.2e12 Msun. **We match this floor.**

My earlier "M_min = 5.68e12" was the top-hat mass at R_f,min=2.507 Mpc/h — the
WRONG proxy. Peak-patch halo masses go below the smallest-filter top-hat mass
(ellipsoidal collapse gives R_halo < R_f), so the real floor is ~10 particles, not
the top-hat at R_f,min. Resolution and R_f,min (= 2·cellsize) both match Websky.

## Bottom line (corrected)

The 120x is **NOT a resolution problem** — our N=6144 / 7700 Mpc/h grid reproduces
Websky's resolution and ~1.7e12 Msun/h mass floor exactly. The 7.35M result was
produced **before** the fixes (fsc_of_z, nbuff, rmax2rs, volume reduction in commits
91c4145 / c9d3df9) and has **not been re-measured since**. The remaining suspects are
algorithmic/config, per INVESTIGATION_2026-04-16:
  - **filter bank**: Websky uses optimal σ(R)-based spacing (paper §2, line 300);
    ours is 20 fixed 1.15x filters over the same range (R_f,min≈2·cellsize → R_f,max=36).
  - **fsc_of_z lightcone bug** (fixed, not yet re-run at ievol=1).
  - **nbuff / merging / volume exclusion / peak threshold** differences.

NEXT STEP: re-run one octant with the FIXED code at the current (Websky-matched)
resolution and measure the new halo count vs Websky's ~1.1e8 halos/octant
(9e8 total / 8 octants). The "120x" is a stale number from the buggy run.

Nate's WebSky2 / Trillium runs above remain a useful *cross-check* (and show realistic
billion-halo counts), but they are a finer-grid, different-box LIM project — NOT the
resolution our original-Websky comparison should match.

## Other notes from the sheet
- WebSky2 (26.04.16) is a line-intensity-mapping target (CO/[CII]); bottom-right block
  maps observer shells a/b/c/X/Y/Z/W/V/U/T to comoving distance ranges for CCAT EoR-Spec
  [CII] (3.5<z<8) and CO surveys. Different science target than original Websky tSZ/kSZ.
- All runs iLPT=2, NG model 0 (Gaussian) for the baselines; fNL variants at ±1e2, ±1e3.
- Memory failures scale with peak count (cenz=0 finds more peaks → OOM vs cenz=7500).
