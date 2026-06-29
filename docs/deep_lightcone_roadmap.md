# Roadmap: deep (z~8) WebSky2-class lightcone via radial-shell geometry

Status: **planned / future** (agreed 2026-06-15). Not started.

## Goal
Extend the pipeline from the current single-box lightcone (z≤4.6) to a **deep z~8 lightcone**
built from **radial shells**, matching WebSky2-class catalogs (CCAT/CII/CO line-intensity-mapping
science needs z~8). This is the agreed future direction.

## Why shells (and why we don't use them today)
- Today: ONE periodic box (6144³), observer at a corner, per-peak redshift from distance, whole
  z=0→4.6 range inside the box (χ(4.6)≈5258 Mpc/h fits). Simple; matches *original* Websky.
- Single-box couples depth to resolution: a bigger box at fixed N = coarser cells (this is exactly
  what caused the factor-h cellsize bug, see `memory/websky_cellsize_factor_h.md`). To reach
  χ(8)≈8500 Mpc/h at fine cells in one box → impractically huge grid.
- **Radial shells decouple depth from resolution**: each shell is its own appropriately-sized,
  finely-resolved box at a comoving-distance `cenz`, with a narrow Δz. Stack shells → deep lightcone.

## Required changes (deltas from the current pipeline)
1. **Radial-shell lightcone tiling** — the core new capability:
   - Partition the lightcone into comoving-distance shells; map each z-shell → (cenz, box size,
     resolution). (WebSky2 uses shells like U/V/W/X/Y/Z/a/b/c — see Nate's Trillium table.)
   - Per-shell box placement / observer geometry (observer offset at cenz, not box corner).
   - Cross-shell stitching: boundary handling, periodic replication, phase/seed continuity so
     structure is continuous across shell joins.
   - Bookkeeping for 8 octants × N shells.
2. **4-field 2LPT** — comes ~free with shells: narrow Δz ⇒ D(z)≈const across a shell ⇒ combine
   ψ₁+ψ₂ at the shell's central redshift into 3 displacement fields → 4 total (vs our current 7).
   ~43% field-memory saving. (Coupled to shells; not a standalone change — see
   `memory/websky2_geometry_fields.md`.)
3. **Finer resolution per shell** (~0.33 Mpc/h, grids ~12k–16k linear vs our 6144) → bigger FFTs,
   more peaks. WebSky2 cost: ~128 nodes, ~5–6 h/octant, ~3–4×10⁹ halos/octant.
4. **Planck18 cosmology** (trivial config): Omx=0.2645, OmB=0.0493, Omvac=0.6862, h=0.6735,
   ns=0.9649, σ8=0.8111, pk = planck18.
5. **Larger buffer**: nbuff ≈ R_smooth,max (~34 Mpc) so the `nhunt`/shell-analysis radius isn't
   clamped (WebSky2 nbuff~66).

## Already done / reusable (the foundation)
- Multi-resolution (O(nmesh³)/tile) — works *within* each shell box.
- Multi-tile + GPU acceleration (~17× vs CPU; ~3× over Fortran) — per shell.
- Validated finder (=Fortran 0.5%), field (σ_j match), per-box lightcone, collapse threshold
  (δc/D(z) to 2.5%), abundance matching (Tinker, Om=0.31).
- Units guardrails (factor-h fix, `CONVENTIONS.md`, both-units runtime log).

## Suggested phasing
- **Phase 0 (in progress):** single-box validated vs *original* Websky at z=4.6 — finish the
  factor-h confirmation (corrected cellsize 0.852 Mpc/h, run 3962039).
- **Phase 1 — shell MVP:** implement radial-shell partition + per-shell box runs + stitching, still
  at z=4.6. Validate: shell-stitched catalog ≈ single-box catalog (regression test).
- **Phase 2 — go deep:** extend to z~8, finer per-shell resolution, 4-field 2LPT, Planck18.
- **Phase 3 — scale & produce:** multi-node, 8 octants × N shells, production catalogs + maps.

## Open design questions
- Shell thickness Δz vs count vs accuracy of the narrow-Δz 4-field combine (where does it break?).
- Structure continuity across shell boundaries (shared seed/phases vs independent + replication).
- Uniform cells per shell vs finer at low z (more nonlinear) — resolution schedule.
- Reuse WebSky2's exact shell scheme vs design our own.

## References
- WebSky2: Nate's Trillium runs `/scratch/nate/runs/26.04.16_WebSky2/` (params not on gw); cosmology
  + filter from gw Cambridge `njcarlson/peakpatch/.../24.09.12_Cambridge/.../param.params`.
- Peak-patch lightcone method: Stein+2019/2020.
- Field-count source: `njcarlson/peakpatch/src/hpkvd/generate_random_field2.f90:240-241` (n_fields 4/7).
