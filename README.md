# PeakPatch.jl

A Julia reimplementation of the Fortran Peak Patch framework. PeakPatch.jl
reproduces the physics of the original ~5000-line Fortran code in ~1900 lines of
idiomatic Julia, with native MPI support via
[PencilFFTs.jl](https://github.com/jipolanco/PencilFFTs.jl) and HDF5 output for
direct integration with
[XGPaint.jl](https://github.com/WebSky-CITA/XGPaint.jl).

This code is developed by a coding agent, so use with more care!

## Features

- **Lightcone evolution** (`ievol=1`): per-peak redshift from observer distance, with tile pruning and horizon cuts
- **Single-tile and multi-tile pipelines** with 1LPT and 2LPT displacement fields
- **MPI-distributed multi-tile** via PencilFFTs.jl pencil decomposition (scales to N² ranks)
- **Deterministic RNG**: Threefry2x counter-based noise (same seed, same field, any rank count)
- **Lagrangian exclusion and volume reduction** (merger) for multi-tile catalogs
- **Extended catalog output**: 33-field records with strain, gradients, Laplacian, formation redshift
- **Local-type non-Gaussianity** (fNL modes 1-2: correlated and uncorrelated)
- **Binary-compatible `.pksc` I/O** (downstream tools unchanged)
- **HDF5 catalog output** with ra/dec/redshift/mass for XGPaint.jl sky map painting
- **TOML-based configuration** and CLI driver
- **Adaptive ODE solver** (OrdinaryDiffEq.jl Tsit5) for ellipsoidal collapse, with RK4 fallback for Fortran validation

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/guanyilun/PeakPatch.jl")
```

Requires Julia 1.11+.

## Quick Start

### From TOML config (recommended)

```sh
julia --project=. bin/peakpatch.jl config.toml --verbose
```

See [`examples/config.toml`](examples/config.toml) for all available options.

### From Julia

```julia
using PeakPatch

# Read parameters from a Fortran-compatible binary file
cfg = PipelineConfig(read_params_bin("hpkvd_params.bin"))
halos = run_tile(cfg; seed=42, verbose=true)

# Or from a TOML config
using TOML
config = TOML.parsefile("config.toml")
cfg = PipelineConfig(config)
halos = run_tile(cfg; seed=42)
```

### Multi-tile

```julia
halos = run_multitile(cfg; ntile=2, seed=42)

# With merger (exclusion + volume reduction)
merged = merge_catalog(halos; verbose=true)
```

### MPI

```sh
mpiexec -np 4 julia --project=. bin/peakpatch.jl config.toml
```

Or programmatically:

```julia
using MPI; MPI.Init()
using PencilFFTs, PencilArrays, PeakPatch

halos = run_multitile_mpi(cfg; ntile=4, seed=42, comm=MPI.COMM_WORLD)
```

### Lightcone mode

In lightcone mode (`ievol=1`), each peak receives its own redshift based on
comoving distance from the observer, rather than all peaks sharing a single
global redshift. This is the mode used by the Websky simulation.

```toml
[grid]
cenx = 0.0    # observer position (Mpc/h), default: box center
ceny = 0.0
cenz = 0.0

[run]
ievol = 1       # enable lightcone evolution
z_out = 0.0     # minimum redshift (observer)
z_max = 4.6     # maximum redshift (horizon cut)
```

In lightcone mode:
- The collapse threshold `fcrit` is computed per-peak from `fsc_of_z(z_peak)`
- Growth factor `D(z)` and velocity scaling are per-peak
- Tiles whose nearest corner exceeds `z_max` are pruned (multi-tile/MPI)
- Peaks beyond `maximum_redshift` are discarded

For full-sky runs (like Websky), run 8 octants with `cenx=ceny=cenz=-boxsize/2`
and stitch the catalogs.

### HDF5 output (for XGPaint.jl)

```julia
using HDF5, PeakPatch

halos = run_tile(cfg; seed=42)
cosmo = CosmologyParams(0.315, 0.049, 0.685, 0.674, 0.965, 0.808)
write_catalog_hdf5("catalog.h5", halos, cosmo)
```

The HDF5 output includes comoving positions, ra/dec, redshift, and M200m, and is
directly compatible with XGPaint.jl's `read_halo_catalog_hdf5()`.

### Collapse table generation

The ellipsoidal collapse ODE uses OrdinaryDiffEq.jl's Tsit5() solver with
adaptive stepping by default. The original hand-rolled RK4 is available via
`solver=:rk4` for Fortran validation:

```julia
using PeakPatch

cosmo = CosmologyParams(0.315, 0.049, 0.685, 0.674, 0.965, 0.808)

# Default: adaptive Tsit5 (more accurate)
ep = EllipsoidParams(cosmo)
table = make_table_threaded(ep, CollapseTableParams(); verbose=true)

# Fortran-compatible: hand-rolled RK4
ep_rk4 = EllipsoidParams(cosmo; solver=:rk4)
table_rk4 = make_table_threaded(ep_rk4, CollapseTableParams(); verbose=true)
```

Standalone CLI script:

```sh
julia --project=. bin/generate_table.jl -o HomelTab.dat -v
julia --project=. bin/generate_table.jl --solver rk4 -o HomelTab.dat -v
```

Or via TOML config with on-the-fly table generation before the pipeline:

```toml
[run]
generate_table = true
ode_solver = "diffeq"   # or "rk4" for Fortran compatibility
```

## Package Structure

```
PeakPatch.jl/
  src/
    PeakPatch.jl                   # Module entry point
    Cosmology/Cosmology.jl         # H(z), D(z), growth factor, comoving distance, chi-to-z inversion
    InitialConditions/
      PowerSpectrum.jl             # P(k) table I/O and log-log interpolation
      RandomField.jl               # Gaussian random field generation (Threefry RNG)
      LPT.jl                       # 1LPT + 2LPT displacement fields
      NonGaussian.jl               # Local-type fNL (modes 1-2)
      LCG.jl                       # Fortran-compatible 48-bit LCG (validation only)
    HaloFinder/
      Filters.jl                   # Window functions and k-space smoothing
      PeakFind.jl                  # 3x3x3 local maximum detection
      RadialShell.jl               # Radial shell integration around peaks
    EllipsoidalCollapse/
      EllipsoidalCollapse.jl       # Triaxial ellipsoid ODE (Tsit5 default, RK4 fallback)
      CollapseTable.jl             # zvir(F, e, p) interpolation table (Interpolations.jl)
    IO/
      Parameters.jl                # PipelineConfig + FortranParams (binary I/O)
      Catalog.jl                   # .pksc catalog I/O (11-field and 33-field)
    Merger/
      Exclusion.jl                 # Spatial hash, sphere overlap, exclusion
      Merger.jl                    # Lagrangian exclusion + volume reduction
    Pipeline.jl                    # Single-tile driver (5-phase pipeline, lightcone support)
    MultiTile.jl                   # Serial multi-tile driver (lightcone + tile pruning)
  ext/
    MPIExt.jl                      # MPI-distributed driver (PencilFFTs)
    HDF5Ext.jl                     # HDF5 catalog output (for XGPaint.jl)
  bin/
    peakpatch.jl                   # CLI driver (TOML config)
  examples/
    config.toml                    # Example configuration
  test/                            # ~509 tests
```

## Pipeline

The single-tile pipeline (`run_tile`) proceeds in five phases:

1. **Initialization** — cosmological parameters, growth tables, collapse table, filter bank
2. **Field generation** — Gaussian random field, FFT, 1LPT/2LPT displacements, optional fNL
3. **Multi-scale peak finding** — smooth at each filter scale (largest first), find local maxima
4. **Shell analysis** — radial profile integration, virialization, halo properties
5. **Catalog output** — write `.pksc` binary catalog

Multi-tile (`run_multitile`) generates the GRF on the full periodic grid and
extracts overlapping tiles for per-tile halo finding. The MPI extension
(`run_multitile_mpi`) distributes all FFT and k-space operations across ranks
using pencil decomposition, with point-to-point tile redistribution.

When lightcone mode is enabled (`ievol=1`), the pipeline computes per-peak
redshifts from comoving distance to the observer. The collapse threshold,
growth factor, and velocity scaling all become redshift-dependent per peak.
Tiles beyond the maximum redshift horizon are pruned before processing.

## Differences from Fortran

Beyond being a direct port, PeakPatch.jl incorporates several modernizations:

- **StaticArrays** for 3-vectors and 3x3 matrices (`SVector{3}`, `SMatrix{3,3}`)
  in `PeakResult`, `kernel_strain`, and `analyse_peak` -- eliminates heap
  allocations and enables value semantics (no defensive `copy()` calls)
- **Pre-allocated RK4 workspace** in the ellipsoidal collapse integrator,
  eliminating up to 30k allocations per halo
- **Cached collapse kernel table** (`atab4`) at module load time instead of
  recomputing the 1101-element vector per peak
- **Threaded peak analysis** (`Threads.@threads`) in both `Pipeline.jl` and
  `MultiTile.jl` with thread-safe atomic diagnostic counters
- **Parametric precision types** -- `PeakGrid{T}`, `PeakCandidate{T}`, and
  generic `smooth_field`/LPT/NonGaussian functions that propagate field
  precision (`Float32` or `Float64`) instead of hardcoding `Float32`
- **Interpolations.jl** for collapse table lookup (replaces hand-coded
  trilinear interpolation)
- **Adaptive ODE solver** (OrdinaryDiffEq.jl Tsit5) as default for ellipsoidal
  collapse, with the original RK4 retained for Fortran validation (`solver=:rk4`)

## Fortran Correspondence

The table below maps each Fortran source file to its Julia equivalent, for
reference when cross-checking implementations.

### Core Pipeline (fully reimplemented)

| Component | Fortran | Julia |
|---|---|---|
| Cosmology & growth factor | `psubs_Dlinear.f90`, `cosmology.f90` | `Cosmology/Cosmology.jl` |
| Power spectrum I/O | `modules/RandomField/pktable.f90`, `read_power_spectrum.f90` | `InitialConditions/PowerSpectrum.jl` |
| Gaussian random field | `modules/RandomField/gaussian_field.f90`, `random.f90` | `InitialConditions/RandomField.jl` |
| 1LPT/2LPT displacements | `modules/RandomField/chi2zeta.f90`, `make_zeta.f90` | `InitialConditions/LPT.jl` |
| Non-Gaussian ICs (modes 1-2) | `modules/RandomField/RandomField.f90` (NonGauss 1-2) | `InitialConditions/NonGaussian.jl` |
| Peak finding | `hpkvd/peakvoidsubs.f90` (`get_pks`) | `HaloFinder/PeakFind.jl` |
| Filter/smoothing | `hpkvd/hpkvd.f90` (smoothing loops) | `HaloFinder/Filters.jl` |
| Radial shell analysis | `hpkvd/peakvoidsubs.f90` (`get_homel`, `icloud`, `atab4`) | `HaloFinder/RadialShell.jl` |
| Ellipsoidal collapse ODE | `modules/HomogeneousEllipsoid/HomogeneousEllipsoid.f90` | `EllipsoidalCollapse/EllipsoidalCollapse.jl` |
| Collapse table (TabInterp) | `modules/TabInterp/TabInterp.f90` | `EllipsoidalCollapse/CollapseTable.jl` |
| Binary params I/O | `modules/GlobalVariables/input_parameters.f90` | `IO/Parameters.jl` |
| Catalog I/O (.pksc) | `hpkvd/io.f90`, `pks2map/pksc.f90` | `IO/Catalog.jl` |
| Merger/exclusion | `merge_pkvd/merge_pkvd.f90`, `exclusion.f90` | `Merger/Exclusion.jl`, `Merger/Merger.jl` |
| Multi-tile domain decomp | `modules/RandomField/tiles.f90` | `MultiTile.jl` |
| MPI distribution | `modules/SlabToCube/SlabToCube.f90`, FFTW MPI | `ext/MPIExt.jl` (PencilFFTs) |
| Python orchestration | `python/peak-patch.py`, `peakpatchtools.py` | `bin/peakpatch.jl` (TOML config) |
| HDF5 catalog output | `python/catalogue_tools/pksc2hdf5.py` | `ext/HDF5Ext.jl` |
| Lightcone evolution (`ievol`) | `hpkvd/hpkvd.f90` (per-peak z, `afn`, horizon cut) | `Pipeline.jl`, `MultiTile.jl`, `MPIExt.jl` + `Cosmology.jl` (`chi_to_z`) |

### Not ported (replaced by XGPaint.jl)

| Component | Fortran | Replacement |
|---|---|---|
| Halo → sky map projection | `pks2map/` (~1500 lines) | [XGPaint.jl](https://github.com/WebSky-CITA/XGPaint.jl) |
| CMB map projection | `pks2cmb/` | XGPaint.jl |
| Profile integration tables | `make_maptable.f90`, `make_cmbtable.f90` | XGPaint.jl |
| HEALPix/FITS I/O | `modules/External/` (HEALPix, CFITSIO) | Healpix.jl / FITSIO.jl via XGPaint |
| BBPS/line profiles | `bbps_profile.f90`, `line_profile.f90` | XGPaint.jl |

### Not yet ported

| Component | Fortran | Notes |
|---|---|---|
| Non-Gaussian modes 3-8 | `modules/RandomField/RandomField.f90` | Spike, deltaN, intermittent, and other bispectrum shapes |
| Distributed fNL (MPI) | — | Modes 1-2 not yet applied in `MPIExt.jl` |
| Octant stitching | `python/peakpatchtools.py` | Full-sky from 8 shifted-observer catalogs |
| Instability module | `src/instability/` | Inflationary perturbation evolution (separate physics) |

## Validation

PeakPatch.jl has been validated component-by-component against the Fortran
reference implementation:

| Component | Agreement |
|---|---|
| LCG random numbers | Bit-for-bit identical |
| Growth factor D(z) | 12 significant figures |
| P(k) interpolation | Machine epsilon at table nodes |
| 1LPT/2LPT divergence | < 10⁻¹² |
| Collapse table | 99.5% of 20k points within 0.35% |
| Radial shell analysis | All outputs < 1% vs Fortran |
| Full tile (142³, LCG) | 10741 vs 10667 halos (0.7%) |
| Multi-tile (ntile=1) | Bit-identical to single-tile |
| MPI (np=1,2,4) | Exact match with serial |

The ~0.7% halo count difference with identical ICs is due to P(k) interpolation
method (Julia uses Float64 log-log on the raw table; Fortran uses Float32 linear
on a rebinned grid). This is expected and irreducible without degrading numerical
accuracy.

## Testing

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

MPI tests (not part of the standard test suite):

```sh
julia --project=. -e 'using MPI; run(`$(MPI.mpiexec()) -np 2 julia --project=. test/test_mpi_multitile.jl`)'
```

## License

TBD
