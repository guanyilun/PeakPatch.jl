module PeakPatch

# Re-export submodules for convenient access
include("Cosmology/Cosmology.jl")
include("InitialConditions/PowerSpectrum.jl")
include("InitialConditions/LCG.jl")
include("InitialConditions/RandomField.jl")
include("InitialConditions/LPT.jl")
include("InitialConditions/NonGaussian.jl")
include("HaloFinder/Filters.jl")
include("HaloFinder/PeakFind.jl")
include("IO/Parameters.jl")
include("IO/Catalog.jl")
include("EllipsoidalCollapse/EllipsoidalCollapse.jl")
include("EllipsoidalCollapse/CollapseTable.jl")
include("HaloFinder/RadialShell.jl")
include("Merger/Exclusion.jl")
include("Merger/Merger.jl")
include("MassFunction/MassFunction.jl")
include("AbundanceMatch/AbundanceMatch.jl")
include("ShellAnalysisGPU.jl")
include("Pipeline.jl")
include("MultiTile.jl")
include("MultiResolution.jl")

# Re-export public API
using .Cosmology: CosmologyParams, E2, H, chi, growth_factor, growth_rate, delta_c,
    DlinearTables, Dlinear_tables, Dlinear_ab, Dfnofa,
    ChiToZTable, build_chi_to_z, chi_to_z, peak_redshift
using .PowerSpectrum: load_pk, load_pk_nongaussian
using .LCG: LCG
using .NonGaussian: apply_fnl_correlated!, apply_fnl_uncorrelated!
using .RandomField: generate_grf, generate_grf_lcg,
    fill_noise_threefry!, fill_noise_threefry_region!
using .LPT: displacements_1lpt, displacements_2lpt,
    displacement_1lpt_component, displacement_2lpt_component,
    compute_src2_k, compute_laplacian_field
using .Filters: gaussian_window, gaussian_window_fortran, tophat_window, read_filterbank
using .PeakFind: PeakCandidate, find_peaks
using .RadialShell: ShellCell, PeakGrid, PeakResult, no_collapse,
    hRinteg, atab4, precompute_shells, analyse_peak, normalize_strain!, normalize_strain,
    fsc_of_z, get_evals, reset_dump_counters!, get_dump_counts
using .Parameters: PipelineConfig, FortranParams, read_params_bin, write_params_bin
using .Catalog: HaloRecord, ExtHaloRecord, write_pksc, read_pksc
using .EllipsoidalCollapse: EllipsoidParams, evolve_ellipse_full,
    get_b_2, _elliptic_rd
using .CollapseTable: CollapseTableParams, CollapseTableInterp,
    make_table, make_table_threaded,
    write_homeltab, read_homeltab, interpolate
using .Exclusion: SpatialHash, build_hash, sphere_overlap,
    lagrangian_exclusion!, volume_reduction!
using .Merger: merge_catalog
using .MassFunction: rho_mean, R_of_M, M_of_R, sigma_R, sigma_M,
    dlnsigma_dlnM, tinker_dndlnM, sheth_tormen_dndlnM,
    cumulative_ngtm, precompute_sigma
using .AbundanceMatch: AbundanceTable, build_abundance_table,
    abundance_match, save_abundance_table, load_abundance_table
using .Pipeline: run_tile
using .MultiTile: run_multitile, run_multitile_lowmem, extract_tile, tile_center
using .MultiResolution: run_multitile_split, compare_fields_split
using .ShellAnalysisGPU: ShellTables, build_shell_tables, analyse_peak_gpu

# MPI extension stub — method defined in ext/MPIExt.jl when MPI+PencilFFTs are loaded
function run_multitile_mpi end

# HDF5 extension stub — method defined in ext/HDF5Ext.jl when HDF5 is loaded
function write_catalog_hdf5 end

# CUDA extension stubs — methods defined in ext/CUDAExt.jl when CUDA.jl is loaded
"""
    shell_fbar_gather_gpu(delta, peaks_i, peaks_j, peaks_k, stab; threads=128)
        -> (Fshell::Matrix{Float32}, nshell::Matrix{Int32})

Gather per-shell delta sums for each peak using a block-per-peak CUDA kernel.
Returns matrices of shape `(nshells, npeaks)` on the host. Validates the GPU
gather primitive against `analyse_peak_gpu`'s CPU shell accumulation.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function shell_fbar_gather_gpu end

"""
    shell_gather_psi_gpu(delta, etax, etay, etaz,
                         peaks_i, peaks_j, peaks_k, stab; threads=128)
        -> (Fshell, nshell, Sshell, Gshell)

Multi-field block-per-peak gather: for each peak and each shell, computes
- `Fshell[s, p]`   = Σ delta[cell]
- `nshell[s, p]`   = count of in-bounds cells
- `Sshell[L, s, p]` = Σ eta_L[cell]                for L ∈ {1,2,3} = {x,y,z}
- `Gshell[L, s, p]` = Σ delta[cell] * d_L           un-normalised (divide by rad[s] at caller)

Matches the `local_S[L]`, `local_G[L]` accumulators of `analyse_peak_gpu`'s
per-shell cooperative gather (see `src/ShellAnalysisGPU.jl:196-232`).

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function shell_gather_psi_gpu end

"""
    shell_gather_strain_gpu(delta, etax, etay, etaz,
                            peaks_i, peaks_j, peaks_k, stab; threads=128)
        -> (Fshell, nshell, Sshell, Gshell, SRshell)

Extends `shell_gather_psi_gpu` with the full 3×3 strain accumulator:
- `SRshell[L, K, s, p]` = Σ eta_L[cell] * d_K            (un-normalised; divide by rad[s] at caller)

The `SRshell` accumulator is the raw per-shell gradient tensor that
`analyse_peak_gpu`'s `kernel_strain` step integrates radially. Not
symmetric per-shell; symmetry enters through `Symmetric(Ebar)` in
`get_evals` after the radial integration.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function shell_gather_strain_gpu end

"""
    shell_gather_full_gpu(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                          peaks_i, peaks_j, peaks_k, stab; threads=128)
        -> (Fshell, nshell, Sshell, S2shell, Gshell, Gfshell, SRshell)

Full-field block-per-peak gather matching the CPU prototype in
`analyse_peak_gpu` (`src/ShellAnalysisGPU.jl:196-232`):
- `Fshell[s, p]`           = Σ delta
- `nshell[s, p]`           = in-bounds cell count
- `Sshell[L, s, p]`        = Σ eta_L
- `S2shell[L, s, p]`       = Σ eta2_L                 (pass zeros to skip 2LPT)
- `Gshell[L, s, p]`        = Σ delta * d_L            (un-normalised)
- `Gfshell[L, s, p]`       = Σ delta[iv3,iv1,iv2] * d_L  (transposed-index gradient)
- `SRshell[L, K, s, p]`    = Σ eta_L * d_K            (raw 3×3 strain accumulator)

The `[iv3, iv1, iv2]` transposed delta read on `Gfshell` mirrors the
existing Fortran/CPU convention used by `kernel_strain`'s `gradf` output.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function shell_gather_full_gpu end

"""
    eig3_symmetric_batch_gpu(mats::AbstractArray{<:Real,3}) -> Matrix{Float32}

Compute eigenvalues of a batch of 3×3 symmetric matrices on GPU using the
analytic Cardano formula. Input shape `(3, 3, N)`; uses only the upper
triangle. Output shape `(3, N)` with eigenvalues sorted ascending.

Test harness for `_eig3_symmetric_f32` — validates against CPU LAPACK.
Not called from the halo finder directly.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function eig3_symmetric_batch_gpu end

"""
    interp3_trilinear_gpu(table, X1, X2, Y1, Y2, Z1, Z2, xs, ys, zs, out_val)
        -> Vector{Float32}

Batched manual trilinear interpolation of a 3D Float32 `table` on an evenly
spaced grid `range(X1, X2, length=size(table,1))` × etc. Out-of-range
queries return `out_val`. Matches `CollapseTable.interpolate`.

Test harness for `_interp3_trilinear_f32` — validates against the CPU
`Interpolations.jl` path. Not called from the halo finder directly.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function interp3_trilinear_gpu end

"""
    kernel_strain_gpu(rad, Gshell, Gshellf, SRshell, mp, mlow, mupp,
                       wRnor, aRnor, hlatt_1, hlatt_2)
        -> (Ebar::SMatrix{3,3,Float32}, grad::SVector{3,Float32},
            gradf::SVector{3,Float32})

GPU-side radial integration. Launches a single-thread kernel that calls
`_kernel_strain_f32` on the supplied per-shell profiles. Intended as a
test harness for validation against CPU `kernel_strain`; the fused
analyse_peak kernel calls `_kernel_strain_f32` inline on thread 0.

- `rad`     : Vector of per-shell radii (length ≥ max shell index)
- `Gshell`  : Matrix (3, nshells) — component × shell gradient
- `Gshellf` : Matrix (3, nshells) — transposed-access gradient
- `SRshell` : Array  (3, 3, nshells) — raw strain accumulator
- `mp, mlow, mupp` : shell indices (1-based)
- `wRnor, aRnor, hlatt_1, hlatt_2` : scaling constants

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function kernel_strain_gpu end

"""
    analyse_peak_gpu_cuda(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                          peaks_i, peaks_j, peaks_k, stab,
                          collapse_table_vals, X1, X2, Y1, Y2, Z1, Z2,
                          alatt, ir2min, ZZon, Rfclvi;
                          fcrit_override=nothing, growth_tables=nothing,
                          rmax2rs=0.0, lapd=nothing, mask=nothing, nbuff=0,
                          threads=128, ct_out_val=-1.0)
        -> NamedTuple with Matrix{Float32} fields per-peak

Full-pipeline GPU halo analysis for a batch of peaks. Launches two CUDA
kernels:
1. The full-field gather kernel (`_shell_gather_full_kernel!`) produces
   per-shell profiles (F, n, S, S², G, Gf, SR).
2. The post-process kernel (`_post_process_kernel!`) runs phases 4-6 of
   `analyse_peak_gpu` on thread 0 per block: Fbar cumulative, fcrit
   crossing, inward virialization walk, RTHL interpolation. Thread 0
   writes a compact result for each peak.

Returns a `NamedTuple` with the following fields (per-peak vectors/matrices):
- `RTHL`         : Vector{Float32} — threshold-crossing radius (-1.0 if no collapse)
- `Fbarx`        : Vector{Float32} — mean-within-RTHL delta
- `e_v`, `p_v`   : Vector{Float32} — ellipticity / prolaticity of final strain
- `strain_final` : 3×3×npeaks       — final strain tensor after RTHL interp
- `eigs`         : 3×npeaks         — ascending eigenvalues of normalised strain
- `Srb`          : Vector{Float32} — energy factor ∫ Fbar(r) r⁴ dr / (Fbarx·RTHL⁵)
- `Sbar`         : 3×npeaks         — mean 1LPT displacement over collapsed shells
- `Sbar2`        : 3×npeaks         — mean 2LPT displacement over collapsed shells
- `gradpk`       : 3×npeaks         — gradient of delta at RTHL
- `gradpkf`      : 3×npeaks         — transposed-axis gradient at RTHL
- `gradpkrf`     : 3×npeaks         — gradient at filter scale (Rfclvi)
- `d2F`          : Vector{Float32} — mean Laplacian within RTHL (0 when `lapd` is omitted)
- `zvir_half`    : Vector{Float32} — formation redshift (collapse z at RTHL/2),
                   max over 10 radii in [RTHL/2, RTHL]; -1.0 when RTHL < 3 or no-collapse
- `mask`         : Array{Int8,3} or nothing — cells within RTHL flagged to 1 when
                   a mutable `mask` keyword is passed.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function analyse_peak_gpu_cuda end

"""
    analyse_peaks_gpu_cuda_multirf(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                                    batches, stab, ct_table,
                                    X1, X2, Y1, Y2, Z1, Z2, alatt, ZZon;
                                    lapd=nothing, mask=nothing, nbuff=0,
                                    growth_tables=nothing, rmax2rs=0.0,
                                    fcrit_override=nothing, threads=128,
                                    ct_out_val=-1.0)
        -> (results::Vector{NamedTuple}, mask::Union{Array{Int8,3},Nothing})

Batched variant of `analyse_peak_gpu_cuda`: uploads the seven input fields
and the mask grid to the GPU **once**, then runs the gather + post-process
kernels separately for each entry of `batches`. Avoids redundant PCIe traffic
when the same tile is analysed at multiple filter scales.

`batches` is a vector of `(peaks_i, peaks_j, peaks_k, Rf, ir2min)` tuples —
one tuple per filter-scale batch. The returned `results` has the same ordering,
each entry a NamedTuple with the same fields as `analyse_peak_gpu_cuda`'s return
(without `mask`).

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function analyse_peaks_gpu_cuda_multirf end

"""
    isolated_convolve_gpu(noise, pk, boxsize_local, n;
                          kernel_fn_id=0, dim1=0, dim2=0,
                          nshell=0, pk_table_n=8192) -> Array{Float32,3}

GPU port of `MultiResolution._isolated_convolve` using cuFFT. Replaces the
zero-padded isolated FFT convolution `√P(k) [* kernel_fn(k)]` with a GPU
pipeline:
1. Upload noise to GPU + zero-pad to `(2*n_input)³`
2. Forward rfft via cuFFT
3. Apply transfer function (P(k) lookup + optional kernel) in a CUDA kernel
4. Inverse irfft via cuFFT
5. Extract central tile region; return as `Array{Float32,3}` on host

`kernel_fn_id` selects the transfer-function variant:
- `0` (default): `√P(k)` — used for δ
- `1`: `i ki / k² · √P(k)` — 1LPT displacement, `dim1 ∈ {1,2,3}`
- `2`: `-i ki / k² · √P(k)` — 2LPT displacement, `dim1 ∈ {1,2,3}`
- `3`: `-ki kj / k² · √P(k)` — `φ_ij`, `dim1, dim2 ∈ {1,2,3}`
- `4`: `k² · √P(k)` — Laplacian

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function isolated_convolve_gpu end

"""
    isolated_convolve_gpu_multi(noise, pk, boxsize_local, n;
                                 kernels=[(0,0,0), (1,1,0), (1,2,0), (1,3,0)],
                                 nshell=0)
        -> Vector{Array{Float32,3}}

Multi-output variant of `isolated_convolve_gpu`. Shares one (upload + zero-pad
+ forward rFFT) of the noise across N transfer kernels, then loops per kernel
running (transfer + Hermitise + irFFT + extract). Used in the per-tile loop
of `run_multitile_split` to compute (δ, ψ₁_x, ψ₁_y, ψ₁_z) from the same
residual noise without redoing the largest cost (zero-padded forward FFT)
four times.

`kernels` is a vector of `(kernel_fn_id, dim1, dim2)` tuples — IDs match
`isolated_convolve_gpu`. Returns one host array per kernel in the same order.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function isolated_convolve_gpu_multi end

"""
    compute_2lpt_gpu(delta_tile, nmesh, boxsize_local) -> Vector{Array{Float32,3}}

GPU port of the tile-local 2LPT block in `run_multitile_split`
(MultiResolution.jl). Given the combined δ field on the `nmesh³` tile grid
(periodic), computes the three 2LPT displacement fields ψ₂_{x,y,z}. Uses
cuFFT for all 11 transforms (1 forward δ rfft, 6 φ_ij irfft, 1 forward
src² rfft, 3 ψ₂ irfft) plus a custom k-space kernel for the φ_ij and
ψ₂ multipliers — matching `_apply_kernel_inplace!` exactly (zero the DC
origin, Nyquist rows/planes).

Returns three `Array{Float32,3}` on the host, indexed `[x,y,z]`, matching
the CPU output of `irfft(psi2_k, nmesh)` for `dim ∈ {1,2,3}`.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function compute_2lpt_gpu end

"""
    compute_laplacian_gpu(delta_tile, nmesh, boxsize_local) -> Array{Float32,3}

GPU port of the tile-local Laplacian block in `run_multitile_split`:
`irfft(k² · rfft(δ))` with the CPU-matching Nyquist zeroing. Single cuFFT
pair + one kernel launch.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function compute_laplacian_gpu end

"""
    interpolate_to_tile_gpu(coarse, it, jt, kt, nsub, nmesh, N, M) -> Array{Float32,3}

GPU port of `MultiResolution._interpolate_to_tile` — a tri-cubic
(Catmull-Rom) interpolation from an M³ coarse periodic grid to the
tile's nmesh³ grid. One CUDA thread per output cell; the coarse field
is uploaded once per call (trivial at O(M³) ≤ ~64KB for typical M ≤ 40).

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function interpolate_to_tile_gpu end

"""
    peak_find_tile_gpu(delta_tile, filters, fcrits, tile_masks_in,
                       xbx, ybx, zbx, alatt, nbuff, wsmooth, ioutshear)
      -> (peaks::Vector{PeakCandidate}, Rf::Vector{Float64},
           FcRf::Vector{Float32},       d2Rf::Vector{Float32},
           mask_out::Array{Int8,3})

GPU port of the per-tile peak-finding loop in
`MultiResolution.run_multitile_split`. Uploads `delta_tile`, computes
rfft(δ) once, then for each filter scale `Rf`:
  1. smooths δ in k-space with window `wsmooth` (0=fortran-Gauss, 1=tophat,
     3=k²-Gauss) — one `_smooth_window_kernel!` + cuFFT irfft
  2. scans the inner cube `(nbuff+1:n-nbuff)³` for strict local maxima above
     `fcrits[i]` and not masked — `_find_peaks_kernel!` with atomic counter
  3. optionally (when ioutshear ≥ 1) also computes a wsmooth=3 field for
     the d²F/dR² diagnostic and samples it at each found peak

The filter loop runs sequentially on the GPU (mask is read+written per
filter), but within a filter all cells are processed in parallel. Output
matches CPU `smooth_field` + `find_peaks` at Float32 precision (the peak
order differs because GPU threads complete in non-deterministic order;
the *set* of peaks is identical).

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function peak_find_tile_gpu end

"""
    set_cuda_device!(id::Int) -> nothing

Bind the current task to CUDA device `id` (0-based). Thin wrapper around
`CUDA.device!` that avoids having to import CUDA into PeakPatch's main
namespace (it's a weakdep). Used by the multi-GPU dispatcher in
`run_multitile_split(...; devices=[...])`.

Defined in `ext/CUDAExt.jl` — requires `using CUDA`.
"""
function set_cuda_device! end

export
    CosmologyParams, E2, H, chi, growth_factor, growth_rate, delta_c,
    DlinearTables, Dlinear_tables, Dlinear_ab, Dfnofa,
    ChiToZTable, build_chi_to_z, chi_to_z, peak_redshift,
    load_pk, load_pk_nongaussian,
    apply_fnl_correlated!, apply_fnl_uncorrelated!,
    generate_grf,
    fill_noise_threefry!, fill_noise_threefry_region!,
    displacements_1lpt, displacements_2lpt,
    displacement_1lpt_component, displacement_2lpt_component,
    compute_src2_k, compute_laplacian_field,
    gaussian_window, gaussian_window_fortran, tophat_window, read_filterbank,
    PeakCandidate, find_peaks,
    ShellCell, PeakGrid, PeakResult,
    hRinteg, atab4, precompute_shells, analyse_peak, normalize_strain!, normalize_strain,
    fsc_of_z, get_evals, reset_dump_counters!, get_dump_counts,
    PipelineConfig, FortranParams, read_params_bin, write_params_bin,
    HaloRecord, ExtHaloRecord, write_pksc, read_pksc,
    EllipsoidParams, evolve_ellipse_full, get_b_2,
    CollapseTableParams, CollapseTableInterp,
    make_table, make_table_threaded,
    write_homeltab, read_homeltab, interpolate,
    sphere_overlap, merge_catalog,
    rho_mean, R_of_M, M_of_R, sigma_R, sigma_M,
    dlnsigma_dlnM, tinker_dndlnM, sheth_tormen_dndlnM,
    cumulative_ngtm, precompute_sigma,
    AbundanceTable, build_abundance_table,
    abundance_match, save_abundance_table, load_abundance_table,
    run_tile,
    run_multitile, run_multitile_lowmem, extract_tile, tile_center,
    run_multitile_split, compare_fields_split,
    run_multitile_mpi,
    write_catalog_hdf5,
    shell_fbar_gather_gpu,
    shell_gather_psi_gpu,
    shell_gather_strain_gpu,
    shell_gather_full_gpu,
    eig3_symmetric_batch_gpu,
    interp3_trilinear_gpu,
    kernel_strain_gpu,
    analyse_peak_gpu_cuda,
    analyse_peaks_gpu_cuda_multirf,
    isolated_convolve_gpu, isolated_convolve_gpu_multi,
    compute_2lpt_gpu, compute_laplacian_gpu,
    peak_find_tile_gpu, interpolate_to_tile_gpu,
    ShellTables, build_shell_tables, analyse_peak_gpu

end # module PeakPatch
