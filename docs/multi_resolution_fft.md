# Multi-Resolution Field Generation for PeakPatch

## Motivation

PeakPatch's Phase 1 (field generation) requires a global forward FFT on the
full N³ grid to convolve white noise with the transfer function √P(k). At
N=6144, this requires ~1 TB of distributed memory per node with MPI
PencilFFTs. This is the binding constraint on grid size — N=12288 is
infeasible with a monolithic global FFT on any current hardware.

The multi-resolution approach eliminates this bottleneck by splitting the
convolution into long-wavelength (global coarse grid) and short-wavelength
(per-tile local FFT) contributions, reducing per-tile memory to O(nmesh³)
instead of O(N³).

## Background

This approach follows the lineage of multi-scale initial condition generators
for cosmological simulations:

- **GRAFIC2** (Bertschinger 2001, [astro-ph/0103301](https://arxiv.org/abs/astro-ph/0103301)):
  Pioneered the Hoffman-Ribak constrained noise split for multi-resolution
  Gaussian random fields. The key insight: split at the white noise level,
  not the power spectrum, making the decomposition mathematically exact.

- **MUSIC** (Hahn & Abel 2011, [arXiv:1103.6031](https://arxiv.org/abs/1103.6031)):
  Improved GRAFIC2 with a three-component decomposition (δ_self + δ_coarse +
  δ_bnd) and isolated (non-periodic) FFT via 2× zero-padding. Achieves
  ~10⁻⁴ relative error vs GRAFIC2's ~10⁻¹.

- **Panphasia** (Jenkins 2013, [arXiv:1306.5968](https://arxiv.org/abs/1306.5968)):
  Hierarchical octree-based Gaussian white noise with 50 levels, enabling
  deterministic per-cell noise generation. PeakPatch's Threefry RNG serves
  a similar role.

## Algorithm

### Noise-Level Split (Hoffman-Ribak)

The density field is a convolution of white noise with the transfer function:

    δ(x) = T * ξ(x)

where T(k) = √P(k) and ξ is Gaussian white noise with unit variance per cell.

We split the noise into a coarse component and a residual:

    ξ(x) = ξ̃₀(x) + ξ₁(x)

where ξ̃₀ is the block-averaged coarse noise spread to the fine grid via
**nearest-neighbor** assignment (each fine cell gets the mean of its coarse
cell), and ξ₁ = ξ - ξ̃₀ is the residual. By construction, ξ₁ has **exactly**
zero mean in each coarse cell (the Hoffman-Ribak constraint). This means
T * ξ₁ is dominated by short-wavelength modes that decay quickly in real
space.

**Important**: The spreading must use nearest-neighbor (not interpolation) to
preserve the exact zero-mean property. Tri-linear or tri-cubic interpolation
would blend across cell boundaries, breaking the constraint and introducing
spurious long-wavelength modes in the residual.

Then:

    δ(x) = T * ξ̃₀(x) + T * ξ₁(x)
          = δ_coarse(x) + δ_self(x)

**No taper function or k_split parameter needed** — the split is determined
entirely by the coarse grid resolution M.

### Phase 1a: Coarse Grid (Global, Cheap)

1. Generate coarse noise by block-averaging Threefry fine-grid noise:
   `coarse[I,J,K] = Σ fine_cells_in_block / √block³`
2. Periodic FFT on M³ grid (M = ntile × coarse_factor)
3. Apply full T(k) = √(P(k) dk³ M³)
4. IFFT → δ_coarse on M³ grid
5. Similarly compute ψ_coarse (1LPT), ψ²_coarse (2LPT), ∇²δ_coarse

Memory: M³ × 4 bytes. For M=36 (ntile=9, cf=4): 187 KB. Trivial.

### Phase 1b: Per-Tile Local Fields

For each tile (it, jt, kt):

1. **Generate residual noise**: Use Threefry to regenerate ξ on the tile's
   nmesh³ region, then subtract the nearest-neighbor spread of coarse block
   means. The residual ξ₁ has exactly zero mean per coarse cell.

2. **Isolated FFT convolution**: Zero-pad ξ₁ to (2×nmesh)³, apply
   T(k) = √(P(k) dk³ n³), IFFT, trim to nmesh³. The 2× padding ensures
   non-periodic (isolated) boundary conditions, avoiding wrap-around.

3. **Interpolate coarse contribution**: Tri-cubic Catmull-Rom interpolation
   of δ_coarse from M³ grid to tile's nmesh³ region.

4. **Combine**: δ_tile = δ_self + δ_coarse_interp

5. **Displacement fields**: Same split for ψ (1LPT kernel ik/k²).
   For 2LPT, compute src2 from the combined δ_tile on the local grid
   (periodic FFT OK since 2LPT is a second-order correction).

Memory per tile: ~(2×nmesh)³ × 4 bytes for the padded FFT.
At nmesh=768: (1536)³ × 4 = 14.5 GB. Fits on a laptop.

### Phases 2-4: Unchanged

Peak finding and shell analysis operate per-tile using the combined fields.
No changes needed — the tile fields have the same format as the global
approach.

## Accuracy

### Error Decomposition

Analysis on a 120³ test grid (ntile=2, cf=5, M=10) reveals the error budget:

| Component | Correlation | RMS Error | Notes |
|-----------|------------|-----------|-------|
| Self term (isolated FFT of residual) | 99.998% | 0.57% | Near perfect |
| Coarse term (M³ FFT + interpolation) | 96.4% | 26.6% | **The bottleneck** |
| Combined (δ_self + δ_coarse_interp) | 99.9% | 3.5-4.3% | Coarse term is ~13% of total variance |

The self term (isolated FFT of the Hoffman-Ribak residual) is essentially
exact. The error is dominated by the coarse grid: the M³ periodic FFT cannot
represent the high-frequency harmonics of the step-function block averaging.
When the coarse field is interpolated to the fine tile grid, these missing
harmonics show up as interpolation error.

### Coarse Factor Tradeoff

The `coarse_factor` parameter controls the balance between two competing
error sources:

| cf | M (ntile=2) | block | δ RMS | ψ₁ RMS | Halo Δ | Notes |
|----|-------------|-------|-------|--------|--------|-------|
| 3  | 6  | 20 | 2.5% | 37% | 2.7% | Best delta, worst displacement |
| 5  | 10 | 12 | 3.7% | 28% | **1.2%** | **Best halo counts** |
| 6  | 12 | 10 | 3.9% | 26% | 5.9% | |
| 10 | 20 | 6  | 6.0% | 21% | 13% | |
| 20 | 40 | 3  | 10.9% | 17% | — | Best displacement, worst delta |

- **Small cf** → large blocks → residual has very short-range correlations →
  isolated FFT is very accurate. But coarse grid has few cells → poor
  interpolation for displacement fields (1/k² kernel amplifies long modes).

- **Large cf** → small blocks → residual retains long-wavelength modes that
  the tile-sized isolated FFT cannot capture. But coarse grid is fine →
  better displacement field interpolation.

- **cf=5** is the optimal compromise for halo count accuracy, achieving
  **1.2% agreement** on a 120³ test grid.

### Why Boundary Correction Doesn't Help

MUSIC's third term (δ_bnd) corrects for residual noise beyond the tile
boundary contributing to the field inside the tile. We investigated this by
extending the isolated FFT to include a shell of residual noise from beyond
the tile (nshell parameter).

Result: negligible improvement. The self term already captures 99.998% of its
contribution. The dominant error source is the coarse grid interpolation, not
the missing boundary noise. The existing tile buffer zone (nbuff cells)
already extends the noise region sufficiently.

### Validation Results

**Field-level comparison** (N=120, ntile=2, cf=5):
- Delta: 99.9% correlation, 3.5-4.3% RMS error across all 8 tiles
- Displacement ψ₁: 94-98% correlation, 20-38% RMS

**Halo count comparison** (N=120, ntile=2, cf=5):
- Global FFT: 2101 halos
- Multi-resolution: 2127 halos
- Difference: **1.2%** (both 1LPT and 2LPT)

On production grids (N=6144, ntile=4+), accuracy should improve because:
1. Tiles are a smaller fraction of the box (better isolated FFT)
2. More coarse cells per dimension (better interpolation)
3. Buffer is a smaller fraction of the tile (less boundary effect)

## Memory Comparison

| Approach | Peak memory/node | MPI required | Max N feasible |
|----------|-----------------|-------------|----------------|
| Global FFT (out-of-place) | ~1108 GB | Yes (PencilFFTs) | ~5000 |
| Global FFT (in-place RFFT!) | ~955 GB | Yes (PencilFFTs) | ~6144 |
| Multi-resolution | ~15 GB/tile | No | Unlimited |

## Implementation

- `src/MultiResolution.jl`: `run_multitile_split` and `compare_fields_split`
- `test/test_multiresolution.jl`: Field-level and halo count validation
- `docs/multi_resolution_fft.md`: This document

### Key Functions

| Function | Purpose |
|----------|---------|
| `_downsample_noise` | Block-average Threefry noise to coarse grid |
| `_generate_extended_residual` | Generate Hoffman-Ribak residual for tile region |
| `_isolated_convolve` | 2×-padded FFT convolution (non-periodic) |
| `_interpolate_to_tile` | Tri-cubic Catmull-Rom interpolation of coarse field |
| `_spread_coarse_to_tile` | Nearest-neighbor spreading (for residual computation) |
| `_periodic_convolve!` | Standard periodic FFT convolution (coarse grid) |
| `compare_fields_split` | Field-level accuracy diagnostic vs global FFT |

### Threefry as the Key Enabler

The Threefry counter-based RNG makes this approach possible: each cell's
noise value is computed from (seed, global_index) in O(1), so any tile can
independently regenerate its own noise region without communication or
storing the full N³ array. This plays the same role as Panphasia's
hierarchical octree in MUSIC.

## Usage

```julia
using PeakPatch

halos = run_multitile_split(cfg; ntile=4, seed=42,
                            coarse_factor=5,  # optimal for halo counts
                            verbose=true)
```

## Future Work

1. **Production integration**: Add `lowmem_split=true` option to
   `bin/peakpatch.jl` for runs that bypass the global FFT entirely.

2. **Parallel execution**: Since tiles are independent, trivially parallelize
   with `Threads.@threads` or distributed across nodes without MPI FFT.

3. **Benchmarking at scale**: Validate on 6144³ against the global MPI FFT
   approach to confirm accuracy improvement at production scales.

4. **Improved coarse interpolation**: The dominant error source is the M³
   coarse grid interpolation. Spectral methods or higher-order schemes could
   reduce this, though the current 1.2% halo count accuracy may be
   sufficient for most applications.

## References

- Bertschinger, E. (2001). "Multiscale Gaussian Random Fields and Their
  Application to Cosmological Simulations." ApJS, 137, 1.
  [astro-ph/0103301](https://arxiv.org/abs/astro-ph/0103301)

- Hahn, O. & Abel, T. (2011). "Multi-scale initial conditions for
  cosmological simulations." MNRAS, 415, 2101.
  [arXiv:1103.6031](https://arxiv.org/abs/1103.6031)

- Jenkins, A. (2013). "A new way of setting the phases for cosmological
  multiscale Gaussian initial conditions." MNRAS, 434, 2094.
  [arXiv:1306.5968](https://arxiv.org/abs/1306.5968)

- Michaux, M. et al. (2021). "Accurate initial conditions for cosmological
  N-body simulations." MNRAS, 500, 663.
  [arXiv:2008.09588](https://arxiv.org/abs/2008.09588)
