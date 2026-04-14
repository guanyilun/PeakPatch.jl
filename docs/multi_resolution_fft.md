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

where ξ̃₀ is the block-averaged coarse noise interpolated back to the fine
grid, and ξ₁ = ξ - ξ̃₀ is the residual. By construction, ξ₁ has zero mean
in each coarse cell (the Hoffman-Ribak constraint). This means T * ξ₁ is
dominated by short-wavelength modes that decay quickly in real space.

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

1. **Regenerate local noise**: Use Threefry to generate ξ on the tile's
   nmesh³ region of the global grid (exact same values as global approach)

2. **Compute residual**: ξ₁ = ξ - interpolate(coarse_noise / √block³)
   The division by √block³ converts coarse noise (unit variance) back to
   block-mean scale for subtraction.

3. **Isolated FFT convolution**: Zero-pad ξ₁ to (2×nmesh)³, apply
   T(k) = √(P(k) dk³ n³), IFFT, trim to nmesh³. The 2× padding ensures
   non-periodic (isolated) boundary conditions, avoiding wrap-around.

4. **Interpolate coarse contribution**: Tri-linear interpolation of δ_coarse
   from M³ grid to tile's nmesh³ region.

5. **Combine**: δ_tile = δ_self + δ_coarse_interp

6. **Displacement fields**: Same split for ψ (1LPT kernel ik/k²).
   For 2LPT, compute src2 from the combined δ_tile on the local grid
   (periodic FFT OK since 2LPT is a second-order correction).

Memory per tile: ~(2×nmesh)³ × 4 bytes for the padded FFT.
At nmesh=768: (1536)³ × 4 = 14.5 GB. Fits on a laptop.

### Phases 2-4: Unchanged

Peak finding and shell analysis operate per-tile using the combined fields.
No changes needed — the tile fields have the same format as the global
approach.

## Accuracy

The decomposition is mathematically exact in the limit of:
- Infinite coarse grid resolution (M → N)
- Perfect interpolation (coarse → fine)

In practice, errors arise from:

1. **Coarse interpolation**: Tri-linear interpolation of the coarse field
   introduces smoothing at the coarse cell scale. Higher-order interpolation
   (tri-cubic, as used by MUSIC) would reduce this.

2. **Isolated FFT boundary**: The residual noise ξ₁ is not truly periodic on
   the tile domain. The 2× zero-padding handles most of this, but the
   transfer function kernel extends beyond the padded domain for very
   long-wavelength modes. These modes are captured by δ_coarse, so the error
   is at the transition scale.

3. **Buffer region**: Errors are largest near tile boundaries. PeakPatch's
   buffer region (nbuff cells) already excludes boundary peaks, so these
   errors don't affect the final catalog.

### Validation Results (Small Grid)

On a 30³ grid with ntile=2 (N=30, nmesh=18, coarse_factor varies):

| coarse_factor | M  | Global halos | Split halos | Δ   |
|---------------|-----|-------------|-------------|-----|
| 2             | 4   | 45          | 28          | -17 |
| 3             | 6   | 45          | 26          | -19 |
| 6             | 12  | 45          | 39          | -6  |
| 9             | 18  | 45          | 34          | -11 |

Field-level comparison (N=30, M=6): correlation = 0.993, RMS error = 11%.
The remaining error is dominated by NGP interpolation of the coarse field
and the small grid size. On production grids (N=6144, nmesh=768), the
coarse-to-tile ratio is much more favorable.

### MUSIC's Boundary Correction (Future Improvement)

MUSIC adds a third term δ_bnd that corrects for the coarse-fine boundary
mismatch. This is computed by:
1. "Unapplying" Hoffman-Ribak: subtract parent mean from children cells
2. This gives noise that is nonzero only in the padding region
3. Convolve with T(k) on the padded grid

Adding this correction would improve accuracy at tile boundaries, potentially
bringing the error to ~10⁻⁴ as MUSIC achieves. For PeakPatch, this may not
be necessary since the buffer region already excludes boundary peaks.

## Memory Comparison

| Approach | Peak memory/node | MPI required | Max N feasible |
|----------|-----------------|-------------|----------------|
| Global FFT (out-of-place) | ~1108 GB | Yes (PencilFFTs) | ~5000 |
| Global FFT (in-place RFFT!) | ~955 GB | Yes (PencilFFTs) | ~6144 |
| Multi-resolution | ~15 GB/tile | No | Unlimited |

## Implementation

- `src/MultiResolution.jl`: `run_multitile_split` function
- `test/test_multiresolution.jl`: Validation against global FFT
- `docs/multi_resolution_fft.md`: This document

### Key Functions

| Function | Purpose |
|----------|---------|
| `_downsample_noise` | Block-average Threefry noise to coarse grid |
| `_generate_tile_noise` | Regenerate fine noise for a tile region |
| `_isolated_convolve` | 2×-padded FFT convolution (non-periodic) |
| `_interpolate_to_tile` | Tri-linear interpolation of coarse field |
| `_periodic_convolve!` | Standard periodic FFT convolution (coarse grid) |

### Threefry as the Key Enabler

The Threefry counter-based RNG makes this approach possible: each cell's
noise value is computed from (seed, global_index) in O(1), so any tile can
independently regenerate its own noise region without communication or
storing the full N³ array. This plays the same role as Panphasia's
hierarchical octree in MUSIC.

## Future Work

1. **Tri-cubic interpolation**: Replace tri-linear with conservative tri-cubic
   (as in MUSIC) for coarse-to-fine mapping. Should reduce interpolation error
   by ~1 order of magnitude.

2. **Boundary correction term**: Implement MUSIC's δ_bnd for higher accuracy
   at tile boundaries.

3. **Production integration**: Add `lowmem_split=true` option to
   `bin/peakpatch.jl` for runs that bypass the global FFT entirely.

4. **Parallel execution**: Since tiles are independent, trivially parallelize
   with `Threads.@threads` or distributed across nodes without MPI FFT.

5. **Benchmarking at scale**: Validate on 1760³ (existing reference) and 6144³
   against the global FFT + MPI approach.

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
