# Field-matter lightcone output — design plan (2026-07-16)

## Why

The halo side of Websky map-making is done and validated (catalogs match in 6 statistics;
painted halo-κ matches Websky's own halos to ~1%, see
`validation/websky_6144/KAPPA_PAINTED_2026-07-16.md`). The measured gap to the released
`kap.fits` is the **field component — ~80% of C_ℓ^κκ at ℓ≲1000** (halo fraction 0.13–0.2
at ℓ=300–800). The same field machinery is required for the kSZ field component (Stein+2020
Fig 6 splits kSZ into halo+field). This cannot be built from catalogs: it needs the LPT
displacement fields, so it is a PeakPatch.jl feature. The Fortran repo does NOT contain
Websky's field painter (hpkvd only dumps ψ fields via `ioutfield`; pks2map paints halos and
merely cites the field equations) — we implement from the paper spec.

## The Websky algorithm (Stein+2020 §3.1.3, §3.2.4)

For each Lagrangian lattice site i (outside halos, for field maps) at q_i:
1. z_i = z(|q_i − obs|) — redshift from the **Lagrangian** distance (displacement's effect
   on emission time declared negligible).
2. Eulerian position x_i = q_i + D(z_i)·s⁽¹⁾ + D⁽²⁾(z_i)·s⁽²⁾ (2LPT).
3. Peculiar velocity v_i = a·[Ḋ(z_i)s⁽¹⁾ + Ḋ⁽²⁾(z_i)s⁽²⁾] (for kSZ).
4. Contribution to the pixel containing n̂(x_i):  δF = (a³_latt/Ω_pix)·W_F(z_i, v·q̂)
   (eq 3.11; the optional Eulerian bias factor b_F=1 for κ and kSZ).
5. **Sub-cell splitting**: if the cell subtends more than a pixel, split into n³ sub-volumes
   (n ≤ 5), each carrying 1/n³ of the contribution, displaced independently — suppresses
   near-observer shot noise.

Field kernels (eq 3.22–3.24), verified consistent with our halo-κ convention:
- W_f^τ  = f_e·ρ_b,0·σ_T·(1+z)²/(μ_e·m_p·χ²), f_e=0.9; He once-ionized z>3, fully z<3
- W_f^ksz = −(v_r/c)·W_f^τ
- W_f^κ  = (3/2)·Ω_M·(H0/c)²·(1+z)·(1−χ/χ*)/χ

Double counting (§3.2.1): the **halo** map is compensated by subtracting a uniform sphere
of the same mass/position with radius at overdensity 3 (slightly inside the Lagrangian
radius); field lattice sites inside halos are excluded. Mass/momentum then conserved on
large scales.

## Architecture decision: a separate painting pass

**`run_multitile_fieldmap(cfg; nside, kernels, ...)` — a standalone pass that regenerates
the tile fields and paints, instead of hooking into the halo pipeline.** Justification:
- The RNG is deterministic (Threefry, counter-based): same seed ⇒ bit-identical fields.
  Nothing needs to be stored from the halo run.
- The pass needs ONLY δ→ψ₁ (+2LPT), i.e. Phase-1 of the tile loop
  (`src/MultiResolution.jl:727–834`). It skips the expensive parts entirely: no filter
  loop, no peak finding, no shell analysis, no merge (those dominate the 6–7 h octant).
  Estimated cost ~1 h/octant on 4×L40S.
- Halo exclusion (field-only maps) needs the FINAL merged catalog anyway — a second pass
  reading the catalog is the natural place (a tile cannot know about neighbours' halos
  whose Lagrangian spheres poke into its core).

## Integration facts (from the current code)

- Per-tile fields: `psi_tile[1:3]`, `psi2_tile[1:3]`, `delta_tile` — nmesh³ Float32,
  **device-resident CuArrays** on the GPU path (`MultiResolution.jl:780–834`). The new pass
  reuses `_generate_extended_residual` + `isolated_convolve_gpu_multi` +
  `interpolate_to_tile_gpu` + `compute_2lpt_gpu` exactly as the halo path does, then stops.
- Geometry: `tile_center(it,jt,kt,ntile,dcore_box)` (centered coords); core = interior
  `nmesh − 2·nbuff`; `alatt` = cellsize; observer `obs` and `chi2z` table in scope; tile
  pruning by `chi(z_max)` already exists (`MultiResolution.jl:679–690`). Paint **core cells
  only** (buffers belong to neighbours).
- Workers: tiles are processed concurrently per worker/GPU — accumulate into **per-worker
  maps**, sum at the end (mirrors `local_halos_basic` pattern, `MultiResolution.jl:695`).
- ⚠️ **2LPT sign convention**: the stored/computed ψ₂ carries the +3/7 convention; Eulerian
  reconstruction uses a MINUS (see `Merger.finalize_eulerian` and the 2026-06-24 velocity
  fix). The cell displacement and velocity formulas MUST be copied from
  `finalize_eulerian`, not re-derived. Growth factors D, f (and the D₂=−3/7D², f₂≈2f
  relations) from `Cosmology.Dlinear_tables`/`Dlinear_ab`.

## Output design

- **Kernel maps** (primary): one HEALPix RING map per requested kernel, accumulated
  per-worker in Float64, written as FITS. Kernels: `:kappa`, `:tau`, `:ksz` (needs v_r),
  plus `:mass` (raw Σ a³latt per pixel — generic). Nside configurable, default 4096
  (0.8 GB Float64/map/worker — fine on host, fits device too).
- **Mass shells** (optional, phase C): nshell×map at Nside ≤ 1024 for generic reuse
  (50 shells × 1.2e7 pix × 8 B ≈ 5 GB). Not needed for the κ/kSZ goals.
- New extension `ext/HealpixExt.jl` (Healpix.jl as a weakdep, like HDF5/CUDA/MPI) for
  pixel indexing + FITS I/O; the core pass lives in src/ and works without Healpix for
  flat-sky test output.
- Config: `[run] ioutfieldmap`, `field_nside`, `field_kernels`, `field_mode`
  (`all` | `exclude_halos`), driver `validation/websky_6144/run_fieldmap_octant.jl`.

## Cost analysis (octant, N=6144, alatt=0.852 Mpc/h, Nside 4096)

- Cells: ~512 tiles × 384³ core ≈ 2.9e10.
- GPU per-cell work (Eulerian + z + kernel weight) is trivial; the bottleneck is
  **pixelization + map accumulation**.
- Sub-cell splitting: cell subtends >1 pixel (0.86′≈2.5e-4 rad) for χ ≲ 3400 Mpc/h; with
  the n≤5 cap the extra work is bounded by ~1.3e11 sub-points concentrated at χ<700 Mpc/h.
- **Phase A (CPU pixelization)**: chunked device→host copies of (x,v,w), threaded
  `vec2pixRing` + accumulate. ~3–4e10 points at Nside 2048 ≈ 30 min on 24 threads —
  acceptable for first validation (ℓ≤2000).
- **Phase B (GPU pixelization)**: implement nested ang2pix in the CUDA kernel + atomic
  adds into a device-resident map (0.8 GB at Nside 4096); host only receives final maps.
  Removes the transfer/CPU bottleneck; enables Nside 4096 + full 5³ splitting.

## Phases

- **A — core pass + κ (validate the concept)**: `run_multitile_fieldmap` with `:all` mode
  (full matter → direct kap.fits analogue), CPU pixelization, `:kappa` + `:mass` kernels,
  FITS out. ~300–400 lines + driver.
- **B — performance**: on-device ang2pix + atomics, 5³ splitting, Nside 4096.
- **C — kSZ + field-only**: velocity kernel maps (v_r), `exclude_halos` mode (mask cells
  within the overdensity-3 sphere of each catalog halo via the Merger spatial hash), the
  matching uniform-sphere compensation for the HALO kSZ/κ maps, optional mass shells.

## Validation plan (phase A exit criteria)

1. **Mass conservation**: Σ map(:mass) = ρ̄ × V(octant ∩ z<z_max) to <0.1%.
2. **Mean κ**: map mean ≈ analytic ∫W_κ dχ (unsubtracted ⟨κ⟩ of matter to z_max).
3. **C_ℓ^κκ vs kap.fits**: full-matter octant κ map vs kap.fits on same-size footprints —
   levels should now agree at ℓ≲1000 (different realization ⇒ compare band means);
   and ours_field vs the already-measured (kap.fits − Websky-halo) difference.
4. **Limber cross-check**: C_ℓ vs ∫dχ W² P_lin/χ² at ℓ≈100–500 (linear scales).

## Open questions (to resolve during implementation)

- Exact Ḋ⁽²⁾: use f₂ = dlnD₂/dlna·H with D₂=−3/7D²Ω^(−1/143) or the `finalize_eulerian`
  approximation — must match the catalog convention for consistency.
- Whether to redshift cells by Lagrangian (paper) or Eulerian distance — follow the paper
  (Lagrangian) for fidelity; trivial to flip.
- kap.fits' z>4.5 Gaussian tail: exclude ℓ-range where it matters from comparisons
  (it is sub-dominant at ℓ≳100) or add the same Limber tail.
- Nside for production: 4096 to match Websky vs 2048 for speed at ℓ≤2000.
