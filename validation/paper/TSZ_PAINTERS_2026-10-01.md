# Which painter reproduces the released Websky tSZ map? (2026-10-01)

**Question.** Is XGPaint, rather than our port of the Fortran `pks2map` physics, what made the
released maps? If so, should we paint tSZ with XGPaint?

**Test.** Websky's own halos (`halos_10x10.pksc`, 2.08M halos, z 0.003–4.61) were painted three
ways at Nside 2048. Each map was compared with the released post-2022 `tsz_2048.fits` on the same
apodized 4.5° and 3° discs used by `check_tsz_painter_wsky.jl`.
- Scripts: `paint_tsz_xgpaint_wsky.jl` (XGPaint env) and `compare_tsz_three_painters.jl`.
- Job 5844292; `results/tsz_three_painters_wskypatch.txt`, log `results/tsz_three_painters_job5844292.log`.

The three painters:
- **ours:** `production/paint_octant.jl`, our port of `pks2map` (Battaglia+12 AGN pressure, Y0 from
  `maptable.f90`, the Websky mass used as the "M200c" proxy, spherical truncation x ≤ 4).
- **xg_asis:** upstream XGPaint `Battaglia16ThermalSZProfile`, the same Battaglia+12 AGN parameters.
  This code was not modified by us; our fork's additions are the κ and τ profiles only. Websky mass
  passed as M200c, i.e. the same convention as ours. Cosmology Ωc = 0.261, Ωb = 0.049, h = 0.68.
- **xg_nfw7:** the same, with M200c from M200m for an NFW with c = 7 (Stein+2020 §3.1.1). The median
  M200c/M200m is 0.966.

## Result (4.5° disc; the 3° disc agrees)

| | ours | xg_asis | xg_nfw7 |
|---|---|---|---|
| mean y / released | **1.006** | 1.501 | 1.103 |
| C_ℓ / released, ℓ = 204 / 330 / 534 / 837 | 1.21 / 1.16 / 1.11 / 1.06 | 3.00 / 2.58 / 2.28 / 2.04 | 1.18 / 1.14 / 1.15 / 1.16 |
| C_ℓ / released, ℓ = 1266 / 1872 / 2730 / 3666 | 1.04 / 1.06 / 1.13 / 1.25 | 1.92 / 1.89 / 2.08 / 2.39 | 1.18 / 1.23 / 1.44 / 1.75 |
| cross-correlation r_ℓ, ℓ = 204 → 3666 | **1.000 → 0.992** | 0.998 → 0.915 | 0.999 → 0.900 |

## Conclusion

**The released tSZ map is reproduced much better by our port of the Fortran physics than by
XGPaint**, under either mass convention:
- **Mean y:** ours 1.006, against 1.10–1.50 for XGPaint.
- **Map-level correlation:** ≥0.992 at all ℓ for ours, against 0.90–0.92 at ℓ ≈ 3700 for XGPaint.
  The XGPaint profile *shape* differs from the released map on small scales.
- **C_ℓ:** ours is closer at every ℓ ≥ 534, and comparable at ℓ = 200–330.

So the released v0.0 tSZ map comes from the Fortran path, or something very close to it, not from
the XGPaint implementation. We keep our painter for tSZ.

Our remaining 4–21% C_ℓ difference from the released map is a smaller, shape-preserving detail of
that path. The candidate, not demonstrated, is the reference table's 4 Mpc transverse truncation.

The XGPaint mass convention matters a lot: passing M200m as M200c gives ×1.5 in mean y, against
×1.0 for the Fortran normalization at the same mass. So XGPaint's amplitude convention differs
from `maptable.f90`'s Y0 at fixed mass. This is an observed difference between two codes, not yet
validated, and not to be reported upstream without the PI's go-ahead.

**Not tested:** κ, because the only XGPaint κ profile is our own addition to the fork, so it is not
independent. CIB is already painted with XGPaint.

## Follow-up: emulating the local Fortran painter (job 5846125; `results/tsz_painter_variants_wskypatch.txt`)

We added switches to `production/paint_octant.jl`, default off:
- `TSZ_TRUNC4MPC=1`: zero beyond 4 Mpc transverse, the `maptable.f90` rmaxt;
- `TSZ_NO_NORM=1`: pixel-centre sampling only, with no exact-integral deposit, as in
  `haloproject.f90`.

The same Websky halos were painted with each, on the 4.5° disc, ratio to the released map:

| ℓ | ours (rerun) | trunc4 + norm | no-norm | Fortran (trunc4 + no-norm) |
|---|---|---|---|---|
| mean y | 1.006 | 1.005 | 1.007 | **0.992** |
| 204 / 330 / 534 | 1.21 / 1.16 / 1.11 | 1.43 / 1.53 / 1.70 | 1.21 / 1.16 / 1.13 | 1.20 / 1.14 / 1.11 |
| 837 / 1266 / 1872 | 1.06 / 1.04 / 1.06 | 1.89 / 2.06 / 2.44 | 1.08 / 1.07 / 1.11 | 1.08 / 1.10 / 1.12 |
| 2730 / 3666 | **1.13 / 1.25** | 3.25 / 5.28 | 1.29 / 1.56 | 1.28 / 1.57 |
| r_ℓ at 3666 | **0.992** | 0.583 | 0.902 | 0.902 |

- **The rerun is identical to the September map**, so the baseline is reproducible.
- **trunc4 + norm is an artifact.** The y cut beyond 4 Mpc is deposited into the centre pixel, which
  gives spikes. It is not a candidate.
- **Full Fortran emulation** gives the closest mean y, 0.992, but does **not** remove the ℓ = 200–330
  excess: 1.20 / 1.14, against 1.21 / 1.16 for ours. It is worse at ℓ ≳ 1800 (r falls to 0.90).
- **So the 4 Mpc truncation hypothesis is refuted** as the cause of the large-scale excess. The
  released map is also not the local Fortran code's pixel-centre sampling. Stein+2020 eq. 3.7
  describes an integral normalization like ours.

**What remains:** a 15–20% excess at ℓ ≈ 200–330, which falls to ~4% by ℓ ≈ 1300, with r ≥ 0.999.
It is carried by the extended outskirts of the most massive clusters: the released map has less y
there than a 4 r200c spherical cut. The candidate is a smaller production cut radius or the 2022
high-mass-centre fix (`UPDATES`, 08-APR-2022). The next test is a cut radius of 2 and 3 r200c, which
needs the Battaglia projection tables rebuilt for each x_max.
