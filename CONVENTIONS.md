# Unit & coordinate conventions (READ BEFORE WRITING A CONFIG)

## TL;DR
**The Julia pipeline works in `Mpc/h` for ALL lengths.** The legacy Fortran peak-patch
(`hpkvd`) and the public Websky/WebSky2 catalogs work in **`Mpc`**. The two differ by a
factor of `h`. Mixing them caused a real production bug (see
`memory/websky_cellsize_factor_h.md`): a Websky "(7.7 Gpc)³" box was entered as
`boxsize = 7700` (interpreted as Mpc/h) when Websky means **7700 Mpc**, making our cells
`1/h ≈ 1.47×` too coarse and raising the halo mass floor from ~1e12 to ~3e12.

## The rule
To match a Fortran/Websky box specified as `L` **Mpc**, set:

```
boxsize_[Mpc/h]  =  L_[Mpc] * h
```

Example (Websky, h=0.68):  `7700 Mpc * 0.68 = 5236 Mpc/h`  →  `boxsize = 5236.0`.
Resulting cellsize = `5236 / 6144 = 0.852 Mpc/h` (= `1.2533 Mpc`), matching Websky.

## What is in which unit
| quantity | Julia pipeline | Fortran / Websky |
|---|---|---|
| `boxsize`, `cellsize`, observer `cenx/y/z` | **Mpc/h** | Mpc |
| halo positions `x,y,z` in output `.pksc` | **Mpc/h** | Mpc |
| `RTHL` (Lagrangian radius) | **Mpc/h** | Mpc |
| filter scales `Rf` in `filter*.dat` | **Mpc/h** | Mpc |
| power spectrum `k` / `P(k)` | `h/Mpc` / `(Mpc/h)³` | same |
| comoving distance `chi(z)` | **Mpc/h** (uses H0=100h) | Mpc |

**Why Mpc/h is forced:** `Cosmology.chi(z)` returns `(2.998e5/100)∫dz/E` = Mpc/h. The
lightcone assigns each halo a redshift from its distance via this `chi`, so positions/box
MUST be Mpc/h for the lightcone to be self-consistent. Do not "fix" this by feeding Mpc.

## Cross-comparing catalogs (the other place this bites)
When comparing OUR catalog (Mpc/h) to a Websky catalog (Mpc), convert Websky lengths first:
- distance: `r_[Mpc/h] = sqrt(x²+y²+z²)_[Mpc] * h`
- mass: identical formula `M = (4π/3)·(2.775e11·Ωm)·R³`; Websky `R` is Mpc so use `R*h`
  to get Mpc/h before cubing (yields Msun/h). Websky's own `readhalos.py` uses
  `rho = 2.775e11·Ωm·h²` with `R` in Mpc → Msun, i.e. `M[Msun/h] = M[Msun]·h`.

## Runtime guard
`run_multitile*` logs (verbose) print cellsize and box in **both** Mpc/h and Mpc at Phase 0,
e.g. `cellsize=0.8522 Mpc/h (=1.2533 Mpc), full box=5236.0 Mpc/h (=7700.0 Mpc)`. If you
intended a Fortran/Websky box of `L Mpc`, check the "(… Mpc)" value reads `L`.

## Config generators
`validation/websky_6144/generate_octant_configs.py` derives `boxsize` from a physical
`BOX_FULL_MPC` value times `h` — keep that pattern so configs are correct-by-construction.
