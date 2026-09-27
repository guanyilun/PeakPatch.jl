# Roadmap: WebSky2.0-class resolution, and optional deep offset-box lightcones

> **Revised 2026-09-27.** The previous version (2026-06-29) said WebSky2 builds a z~8 lightcone
> from stacked radial shells, with 4-field 2LPT "for free" per shell. **That was a misreading
> of Nate's Trillium run table.**
>
> - The WebSky2 rows have `cenz` = 3000/4000/4322 = **box/2**, i.e. the ordinary
>   corner-observer octant geometry we already use, at finer cells. The rows are even named
>   `octant_+++`.
> - The "U/V/W/…" shell names in the old text have no source in our notes.
> - No WebSky2.0 paper exists yet. CITA lists "The WebSky2.0 Mocks for next-generation CMB and
>   LSS surveys – N. J. Carlson et al., in prep".

## What the evidence actually says

1. **WebSky2.0 CMB/LSS mocks.** Source: Nate's "Trillium Runs" spreadsheet, transcribed in
   `validation/websky_6144/REFERENCE_trillium_runs_2026-06-13.md`.
   - Geometry: octant runs like WebSky1.0.
   - Cells: about 0.49–0.57 Mpc, roughly 0.33–0.39 Mpc/h, on n_ext 12288–16384 grids.
   - Buffers: n_buff 64–69, with buffer ≈ R_smooth,max.
   - Cosmology: Planck18, iLPT=2.
   - Output and cost: about 3e9 halos per octant; 128 Trillium nodes for about 5–6 h per octant.
   - Some later test rows list N_fields=4 with iLPT=2. **How that works is unknown**: the local
     Fortran clone still picks 4 fields for 1LPT and 7 for 2LPT (`hpkvd.f90:223-224`).
2. **Deep line-intensity mocks.** Source: arXiv:2510.18312, the WebSky [CII] forecasts by
   Carlson, Bond et al.; text at `~/work/peakpatch/2510.18312.txt`, l.667–679.
   - Two adjacent (1100 Mpc)³ cubes **centred 7.5 and 8.6 Gpc from the observer**, covering
     z = 3.5–8 in a 4 deg² patch.
   - Grid: 5640³ cells per cube (5586³ excluding buffers), run as 21³ sub-volumes.
   - Resolution: 197 kpc cells, M_min 4.3e9 M☉.
   - This is the existing Fortran `cenx/ceny/cenz` offset-box feature
     (`hpkvd.f90:1325-1331, 1531`). It needs no new shell code inside the halo finder.
3. **Fortran reference.**
   - The local clone `~/work/peakpatch` is the 2020-era GitLab code plus local edits.
   - The modern code (with `hpkvdmodule.f90`) and the WebSky2 run directories are on gw:
     `/fs/lustre/project/act/njcarlson/peakpatch` and `.../WebSky2`. These need interactive
     login and have not been read in the current session.
   - **Read these before designing anything WebSky2-specific.**

## Work items, in order

1. **Lift the GPU shell-count limit.** `_MAX_SHELLS_GPU = 512`, in shared memory, gives
   nhunt ≤ 24 cells (`ext/CUDAExt.jl`).
   - WebSky2-class n_buff ≈ 66 needs nhunt up to about 65. The same limit also sets the
     current tSZ mass cap (`validation/paper/V2_RESULTS_2026-09-27.md`).
   - Approach: stream shells from global memory in chunks. The shell early exit keeps the
     mean cost near the typical halo radius.
   - Acceptance: catalogs identical to the current kernel at nbuff ≤ 25, and CPU/GPU agreement
     at larger nhunt.
2. **Slim catalog format.** At about 3e9 halos per octant, our 33-float records come to about
   450 GB per octant. Offer a Websky-style 10-float format, with shear columns optional.
3. **Finer single-box octant.** N=12288 at the current box is about 2× finer and costs about 8×.
   - Config, filter bank and validation against the N=6144 run on a common mass range.
   - MultiTile/MPI must support `periodic_cores` for multi-node runs (they currently error).
4. **Planck18 cosmology option.** Regenerate P(k), the collapse table and the AM cosmology.
   This is routine.
5. **Optional: offset-box deep lightcones** (for LIM science, not WebSky2.0 CMB mocks).
   - Make an observer outside the box a supported mode. `_warn_if_observer_outside_cores`
     currently only warns.
   - Add a z_min/χ_min cut, which the local Fortran does not have in the finder.
   - Assemble adjacent cubes.
6. **4-field 2LPT.** Blocked until the modern Fortran (gw) shows what N_fields=4 with iLPT=2
   actually does.

## What already carries over
- Multi-resolution tiling (O(nmesh³) per tile, no global FFT), periodic cores, splice
  compensation.
- Multi-tile GPU acceleration, and the exact shell early exit.
- The tile cull beyond χ(z_max) (`MultiResolution.jl`).
