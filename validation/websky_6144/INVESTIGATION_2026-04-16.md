# Investigation: 120x Fewer Halos Than Original Websky

Date: 2026-04-16

## Context

Our GPU production run found **7.35 million halos** across 8 octants (full sky).
The original Websky (Stein+ 2020, arXiv:2001.08787, Section 4.4.3) reports
**~9 x 10^8 (900 million) halos**. This is a factor of ~120x fewer.

This document describes a systematic comparison of our Julia code against the
Fortran reference at `../peakpatch/` and the published paper.

---

## Finding 1: `fsc_of_z` Bug — Wrong Critical Overdensity at z > 0

**Severity: HIGH — correctness bug affecting all lightcone (ievol=1) runs**

### What `fsc_of_z` does

This function returns the critical linear overdensity for spherical collapse by
redshift z. It is used in two places:

1. **Peak finding** — cells with delta < fcrit are skipped (PeakFind.jl:40)
2. **Shell analysis** — determines the largest collapsing radius around each
   peak (RadialShell.jl:319, ShellAnalysisGPU.jl:152)

### The Fortran implementation (correct)

File: `../peakpatch/src/hpkvd/peakvoidsubs.f90:25-51`

```fortran
real function fsc_of_z(z)
  use TabInterp
  real, intent(in) :: z
  real fv, dfv, zvir

  fv = TabInterpX2          ! start at max log10(Frho) in collapse table
  dfv = 1e-4
  zvir = 100

  do while(zvir > 1+z)
     zvir = TabInterpInterpolate(fv, 0., 0.)   ! collapse redshift for spherical (e=0, p=0)
     fv = fv - dfv
  enddo

  fsc_of_z = 10**fv         ! return Frho = 10^(log10(Frho))
end function
```

This walks the collapse table from high to low overdensity, finding the
**minimum linear overdensity Frho for which spherical collapse occurs by
redshift z**. For the standard cosmology:

| z     | D(z)  | Fortran fsc_of_z(z) | Physical meaning                    |
|-------|-------|---------------------|-------------------------------------|
| 0.0   | 1.0   | ~1.686              | Standard spherical collapse at z=0  |
| 1.0   | ~0.5  | ~3.4                | Must be denser to collapse earlier  |
| 2.0   | ~0.3  | ~5.6                | Even denser                         |
| 4.0   | ~0.15 | ~11                 | Only extreme peaks collapse by z=4  |

The physical intuition: the linear density field is extrapolated to z=0. For a
region to have collapsed by an earlier epoch z > 0, its z=0 linear overdensity
must be **higher** than 1.686 (since D(z) < 1, the actual overdensity at epoch
z is delta_0 * D(z), which must reach delta_c = 1.686).

### The Julia implementation (WRONG)

File: `src/HaloFinder/RadialShell.jl:191-196`

```julia
function fsc_of_z(z::Float64, tables)
    a = 1.0 / (1.0 + z)
    dlin, = Dlinear_ab(a, tables)
    return Float64(1.686 * dlin)       # <-- BUG: multiplies by D(z) instead of dividing
end
```

This returns `1.686 * D(z)`, which **decreases** with z:

| z     | D(z)  | Julia fsc_of_z(z) | Error vs Fortran |
|-------|-------|-------------------|------------------|
| 0.0   | 1.0   | 1.686             | 0% (agrees)      |
| 1.0   | ~0.5  | ~0.84             | 4x too low       |
| 2.0   | ~0.3  | ~0.51             | 11x too low      |
| 4.0   | ~0.15 | ~0.25             | 44x too low      |

The formula is the INVERSE of the correct one. It should be `1.686 / D(z)`,
not `1.686 * D(z)`.

### Impact on shell analysis (the halo-count-affecting path)

The shell analysis uses fcrit to find the largest radius where the mean
enclosed overdensity Fbar(R) >= fcrit. This determines the halo's Lagrangian
radius (RTHL) and mass.

**At z=2, for a peak with central delta ~ 3.0:**
- **Fortran**: fcrit = 5.6. Since delta_peak = 3.0 < 5.6, the peak is
  REJECTED (Fbar never exceeds fcrit). Physically correct — this peak hasn't
  collapsed by z=2.
- **Julia**: fcrit = 0.51. Since delta_peak = 3.0 > 0.51, the peak is
  ACCEPTED and assigned a large radius (wherever Fbar drops to 0.51).
  Physically wrong — gives an unrealistically large halo.

The net effect is complex:
- At high z, Julia accepts peaks that should be rejected → some extra halos
- At high z, accepted halos get inflated radii → more overlap exclusion
- The inflated radii corrupt mass estimates for all lightcone halos at z > 0

### Additional issue: per-tile vs constant threshold in peak finding

The Fortran peak-finding (get_pks) uses a **single constant threshold** for
the entire lightcone:

```fortran
fv = fsc_of_z(z)    ! z = z_out = 0.0 for lightcone, gives fv = 1.686

! The per-cell lightcone code is COMMENTED OUT in the Fortran:
!  if(ievol==1) then
!     rc = sqrt(xc**2+yc**2+zc**2)
!     ac = afn(rc)
!     rdsc = 1/ac - 1
!     if(rdsc > maximum_redshift) cycle
!     fv = fsc_of_z(rdsc)
!  endif

if(ff.lt.fv) cycle   ! fv = 1.686 for ALL cells
```

Our Julia code applies a **per-tile threshold** in lightcone mode:

```julia
# MultiResolution.jl:866-873
if ievol == 1
    fcrit_tile = Float32(fsc_of_z(z_tile, growth_tables))  # varies per tile
    fill!(fcrits_per_filter, fcrit_tile)
```

This is a double divergence from Fortran:
1. The threshold value is wrong (multiply vs divide by D)
2. The threshold varies per tile instead of being constant

---

## Finding 2: Filter Bank Configuration (Primary Count Discrepancy)

### What the paper says

Section 2.2 (lines 298-302):
> "The filter bank can consist of logarithmically-spaced filter radii R, or
> linear or logarithmic spacing in sigma(R), with optimal filter spacings to
> maximize both accuracy and efficiency presented in [68], from a minimum
> radius of Rf,min = 2 a_latt..."

Reference [68] is Stein, Alvarez & Bond (2019), MNRAS 483, 2236
(arXiv:1810.07727), which presents an optimized filter spacing algorithm.

### What we use vs what the Fortran generates

Both our Julia and the Fortran `python/filter_gen.py` use **fixed 1.15x
logarithmic spacing**, producing ~20-21 filters from Rf_min to Rf_max.

The production Websky almost certainly used a different, denser filter bank
from the optimal spacing algorithm of ref [68]. The exact configuration is not
in the Fortran repo we have — it may have been a separate script or
hand-tuned for the Niagara production run.

### Mass threshold comparison

| Source                     | Rf_min (Mpc/h) | M_min (M_sun/h) | Total halos |
|----------------------------|----------------|------------------|-------------|
| Our run                    | 2.507          | 5.7e12           | 7.35e6      |
| Fortran (rmincell=1.65)    | 2.068          | 3.2e12           | (not run)   |
| Paper (10 particles)       | —              | 1.2e12           | ~9e8        |

The halo mass function goes as dn/dM ~ M^{-2} at the low-mass end. A factor
of ~5 reduction in minimum mass (5.7e12 -> 1.2e12) with the exponential Press-
Schechter/Tinker tail can produce 50-100x more halos — accounting for most of
the 120x discrepancy.

Note: the paper says halos are found at radii down to a lattice spacing
(1.25 Mpc/h) via the shell analysis shrinking step (Section 2.2, line 317:
"measurements are performed at decreasing radii until the region collapses,
or is equal to a radius smaller than a lattice size"). The filter scale sets
where PEAKS are found, but the shell analysis determines the actual halo
radius, which can be smaller.

---

## Finding 3: `rmax2rs = 1.0` vs Fortran's `0.0`

Our config: `rmax2rs = 1.0`
Fortran default: `rmax2rs = 0.0` (param.params line 195)

When rmax2rs > 0, the shell analysis limits the search volume to
(rmax2rs * Rf / a_latt)^3 particles. For the smallest filter:

| Setting        | Rf=2.507, a_latt=1.253 | Rf=35.7, a_latt=1.253 |
|----------------|------------------------|-----------------------|
| rmax2rs = 1.0  | 33 particles           | 67k particles          |
| rmax2rs = 0.0  | npartmax (nbuff^3)     | npartmax (nbuff^3)     |

With rmax2rs=0 and nbuff=22 (Fortran): npartmax ~ 44,600 for ALL filters.
With rmax2rs=1 (our run): only 33 particles for the smallest filter.

This severely limits the shell analysis at small scales, reducing its ability
to accurately determine the collapsing radius and potentially rejecting peaks
that should pass.

---

## Finding 4: `nbuff = 8` Is Below Minimum

Our config: `nbuff = 8` → 8 x 1.253 = 10.0 Mpc/h buffer thickness

The Fortran documentation (`peakpatchtools.py:80`):
> "This buffer has to be at least about 32 Mpc in thickness"

The Fortran example: `nbuff = 22` → ~27.5 Mpc/h at its grid resolution.

With nbuff=8:
- Shell analysis for large-filter peaks near tile boundaries is truncated
- Rf_max = 35.7 Mpc/h >> 10 Mpc/h buffer → large-scale peaks corrupted

However, the multi-resolution approach is somewhat less sensitive to nbuff
than the standard tiled FFT, since the isolated FFT handles the long-range
field correctly. The main concern is shell analysis radius limitation near
tile edges, which affects a minority of halos.

---

## Finding 5: Minor Differences

### rmincell: 2.0 vs 1.65
Our Julia: Rf_min = 2.0 * a_latt = 2.507 Mpc/h (20 filters)
Fortran: rmincell = 1.65, giving Rf_min = 1.65 * a_latt = 2.068 Mpc/h (21 filters)

One extra filter at the small end, slightly extending mass coverage.

### Peak threshold "roughly 1.5"
The paper (line 306) says peaks are found with "density contrast above a
threshold determined from spherical collapse, roughly 1.5". Our code and the
Fortran both use delta_c = 1.686. The paper's "roughly 1.5" is likely loose
language, not a different parameter — but worth verifying.

---

## Why Previous Validation Tests Didn't Catch the fsc_of_z Bug

### What WAS tested

**Julia-vs-Fortran (snapshot mode, z=0):** Three runs, all ievol=0:

| Validation run         | ievol | z_out | fsc_of_z(0) Julia | fsc_of_z(0) Fortran | Match? |
|------------------------|-------|-------|-------------------|---------------------|--------|
| 256^3 single-tile      | 0     | 0.0   | 1.686             | 1.686               | YES    |
| 1024^3 single-tile     | 0     | 0.0   | 1.686             | 1.686               | YES    |
| ntile=2 multi-tile     | 0     | 0.0   | 1.686             | 1.686               | YES    |

Results: sub-1% halo count agreement in all three runs.

At z=0, D(0) = 1, so `1.686 * D(0) = 1.686 / D(0) = 1.686`. The multiply-
vs-divide bug cancels exactly when z=0.

**Lightcone mode (ievol=1) WAS used in production:**

| Run                    | ievol | Notes                                    |
|------------------------|-------|------------------------------------------|
| websky_2048 (1760^3)   | 1     | First lightcone run; results never checked |
| websky_6144 (8 octants)| 1     | Production run; yielded 7.35M halos       |

Lightcone mode was NOT skipped — it was deployed. But it was never validated
against Fortran or against theoretical expectations.

### The missed validation

The websky_2048 validation plan (`validation/websky_2048/NOTES_2026-04-13.md`,
lines 167-175) explicitly included the check that would have caught this bug:

> **Per-peak fcrit**: in lightcone mode, fcrit = fsc_of_z(z_peak) varies per
> peak. At z=0, fcrit ~ 1.686. At z=1, fcrit should be **higher** (stronger
> collapse threshold at earlier times). Verify a sample of halos.

This checkbox was left unchecked. The websky_2048 job was submitted
(job 596530) but the results were never analyzed — the work moved directly
to the GPU pipeline and the Killarney 6144 production before anyone verified
the lightcone-specific behavior.

The irony: the validation plan correctly identified the expected behavior
(fcrit increases with z) and proposed exactly the right check, but it was
never executed.

### Why the bug survived to production

1. **No Fortran lightcone reference run was ever set up.** All Julia-vs-Fortran
   comparisons used z=0 snapshots. A Fortran ievol=1 comparison would require
   setting up observer coordinates, lightcone geometry, and per-peak redshift
   computation — significantly more complex than a snapshot comparison.

2. **The websky_2048 validation was abandoned mid-flight.** The first ievol=1
   run was submitted but its validation checklist (which included the key
   fcrit-vs-z check) was never completed. Development pivoted to the GPU
   pipeline and Killarney deployment.

3. **The function signature is misleading.** The Julia fsc_of_z takes a
   `tables` argument (growth tables), suggesting careful cosmology-awareness.
   The `1.686 * dlin` formula LOOKS reasonable at first glance — it's the
   growth-factor-scaled threshold, just with the wrong direction of scaling.

4. **The production run "worked."** The 6144 lightcone produced 7.35M halos
   with consistent counts across 8 octants, correct spatial distribution,
   and ~19 min runtime. Nothing obviously failed — the halos were just wrong
   in mass/size at z>0, and the total count was plausible enough not to
   trigger immediate alarm (the 120x discrepancy was attributed to filter
   bank differences, which IS a real contributing factor).

5. **rmax2rs and nbuff differences were intentional trade-offs for GPU
   tiling.** The GPU multi-resolution approach uses more, smaller tiles
   (ntile=16 vs Fortran's ntile=2-9), so nbuff=8 was chosen to keep the
   padded FFT size manageable on GPU memory. rmax2rs=1.0 was set to limit
   shell analysis cost per tile. These were treated as acceptable
   approximations, not as potential sources of large error.

### The README documents the buggy behavior as intended

README.md line 134 describes lightcone mode:
> "The collapse threshold `fcrit` is computed per-peak from `fsc_of_z(z_peak)`"

This is the Julia behavior (per-peak fcrit in peak finding), but does NOT
match the Fortran, where peak-finding uses a CONSTANT threshold and only
the shell analysis uses per-peak fcrit. The README's Validation section
(lines 348-431) contains zero lightcone entries despite lightcone being
a headline feature (line 14). The validation gap was hiding in plain sight.

### What tests would have caught this

- **Completing the websky_2048 validation checklist** — specifically the
  "Per-peak fcrit" check at line 168 of that file.
- A Fortran vs Julia comparison at z_out > 0 with ievol=0 (snapshot at
  z=1 or z=2) would immediately show divergent halo counts.
- A lightcone (ievol=1) validation on a small grid with known Fortran
  output would have revealed both the fsc_of_z bug and the per-tile
  vs constant threshold difference.
- A unit test of fsc_of_z against the collapse table for z > 0 would have
  caught the multiply-vs-divide error directly.

---

## Summary of Issues and Expected Impact

| Issue                          | Type       | Expected impact on halo count   |
|--------------------------------|------------|---------------------------------|
| fsc_of_z wrong formula         | Bug        | Corrupts halo masses/sizes at z>0; complex effect on counts |
| Per-tile vs constant threshold | Bug/Design | Changes which peaks are found at z>0 |
| Filter bank (20 fixed vs optimal) | Config  | PRIMARY: 50-100x fewer halos (mass function steepness) |
| rmax2rs = 1.0 vs 0.0          | Config     | Limits shell analysis accuracy at small scales |
| nbuff = 8 vs 22+              | Config     | Affects large-scale peaks near tile boundaries |

The **filter bank** is likely the dominant factor in the 120x discrepancy
(accounting for ~50-100x through minimum mass difference and filter coverage).
The **fsc_of_z bug** corrupts lightcone halo properties but its net effect on
total count is harder to predict without a test run. The **rmax2rs** and
**nbuff** settings are secondary but should be corrected.

---

## Recommended Fix Order

1. Fix `fsc_of_z` to use collapse table (matching Fortran exactly)
2. Fix peak-finding to use constant threshold in lightcone mode (match Fortran)
3. Set `rmax2rs = 0.0` in octant configs
4. Rerun single octant to measure impact of fixes 1-3
5. Investigate optimal filter spacing (ref [68]) and denser filter bank
6. Consider increasing nbuff (requires re-tuning nmesh/ntile for GPU memory)
