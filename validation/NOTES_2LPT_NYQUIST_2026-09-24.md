# MPI-path displacement bugs + the 2LPT Nyquist convention — 2026-09-24

## Found while wiring `finalize_eulerian` into the CLI
`ext/MPIExt.jl` (the `use_mpi=true` path of `bin/peakpatch.jl`) had missed BOTH fixes
that were applied to `MultiTile.jl` (95f04a1) and `MultiResolution.jl` (4ae7837):

1. **Growth factor**: `_, _, D_pk = Dlinear_ab(...)` took D/a (3rd return), not D.
2. **2LPT sign**: `-Sbar2 .* (-3/7 …)` = +3/7; Fortran and the other paths use −3/7.
3. **2LPT source (trace identity)**: the MPI path builds src2 = δ²/2 − Σφᵢᵢ²/2 − Σφᵢⱼ²
   (saves memory), which needs (Σφᵢᵢ)² = δ². Its φᵢⱼ kernel zeroed k=0 AND the Nyquist
   planes while δ² kept them, so the identity failed. Serial `LPT.jl` (direct determinant)
   zeroes only k=0. Fixed by zeroing only k=0 in the MPI φᵢⱼ/trace kernels.

Why no test caught it: `test_mpi_multitile.jl` compared positions/RTHL only, at z=0 with
1LPT — where D/a = D (a=1) and there is no 2LPT term. New testset "MPI np=1
displacements ≈ serial (z=1, 2LPT)" compares all six displacement fields.
Verified: with fixes all 6 fields match serial exactly (ratio 1.0000); against the
original MPIExt.jl 6/8 assertions FAIL (so the test catches all three bugs).

## Production (GPU) carries (3) — quantified, negligible
The GPU tile-local 2LPT (`CUDAExt._periodic_kernel!` + `δ²/2` trace term; mirrors
`MultiResolution._apply_kernel_inplace!`) has the same Nyquist mismatch. NOT changed:
the production campaign is frozen. `quantify_2lpt_nyquist.jl` (production cell 0.852
Mpc/h, Websky P(k), seed 12345):

| grid | rms Δψ₂/ψ₂ (cell) | (R=2 Mpc/h) | z=0 total-displacement error (R=2) |
|---|---|---|---|
| 128³ | 3.3% | 1.6% | 2.2e-3 |
| 256³ | 2.3% | 1.1% | 1.4e-3 |

Scales ≈ 1/√n → production 414³ tiles ≲ 1e-3 of total displacement at halo scales
(less at higher z, ψ₂ ∝ D²); ≪ the 1–2% σ_vr agreement with Websky. Fix (zero only
k=0 in the GPU φᵢⱼ kernel, id 3) whenever the GPU path is next re-frozen; re-run
`test/test_multiresolution_gpu.jl` and the Websky σ_vr check after.
