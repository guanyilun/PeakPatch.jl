module MultiResolution

# Multi-resolution field generation following Hahn & Abel (2011, MUSIC).
#
# Instead of a global N³ FFT, splits the convolution T*ξ into:
#   δ_coarse : contribution from coarse-grid noise (captures long modes)
#   δ_self   : contribution from tile-local residual noise (isolated FFT)
#
# The split is at the WHITE NOISE level (not the power spectrum), following
# the Hoffman-Ribak constrained realization.  No taper function needed.

using FFTW
using Printf: @sprintf

using ..Cosmology: CosmologyParams, Dlinear_tables, Dlinear_ab, chi,
    build_chi_to_z, peak_redshift
using ..PowerSpectrum: load_pk
using ..RandomField: _threefry_gaussian, generate_grf, fill_noise_threefry!
using ..LPT
using ..LPT: displacements_1lpt
using ..Filters: smooth_field, gaussian_window_fortran, tophat_window, read_filterbank
using ..PeakFind: PeakCandidate, find_peaks
using ..RadialShell: PeakGrid, precompute_shells, analyse_peak, fsc_of_z
using ..ShellAnalysisGPU: build_shell_tables
using ..Parameters: PipelineConfig
using ..Catalog: HaloRecord, ExtHaloRecord
using ..CollapseTable: CollapseTableInterp, read_homeltab
using ..MultiTile: tile_center, extract_tile

export run_multitile_split, compare_fields_split

# ============================================================
# Core building blocks
# ============================================================

"""
Generate Threefry noise for a tile's region of the global N³ grid.
Returns nmesh³ array with the exact same noise values as the global grid.
"""
function _generate_tile_noise(it::Int, jt::Int, kt::Int,
                               nsub::Int, nmesh::Int, N::Int, seed::Integer)
    nbuff = (nmesh - nsub) ÷ 2
    noise = Array{Float32,3}(undef, nmesh, nmesh, nmesh)
    s = UInt64(seed)
    i0 = (it - 1) * nsub + 1 - nbuff
    j0 = (jt - 1) * nsub + 1 - nbuff
    k0 = (kt - 1) * nsub + 1 - nbuff

    Threads.@threads for lk in 1:nmesh
        @inbounds for lj in 1:nmesh, li in 1:nmesh
            gi = mod1(i0 + li - 1, N)
            gj = mod1(j0 + lj - 1, N)
            gk = mod1(k0 + lk - 1, N)
            idx = gi + (gj - 1) * N + (gk - 1) * N * N
            pair_idx = UInt64((idx - 1) >> 1)
            g1, g2 = _threefry_gaussian(s, pair_idx)
            noise[li, lj, lk] = Float32(iseven(idx - 1) ? g1 : g2)
        end
    end
    return noise
end

"""
Downsample fine-grid noise to coarse grid by block-averaging.
Each coarse cell (I,J,K) = mean of block³ fine cells.
Normalized so that coarse noise has unit variance per cell.
"""
function _downsample_noise(N::Int, M::Int, seed::Integer)
    block = N ÷ M
    coarse = zeros(Float64, M, M, M)
    s = UInt64(seed)

    Threads.@threads for K in 1:M
        for J in 1:M, I in 1:M
            val = 0.0
            for dk in 0:block-1, dj in 0:block-1, di in 0:block-1
                gi = (I - 1) * block + di + 1
                gj = (J - 1) * block + dj + 1
                gk = (K - 1) * block + dk + 1
                idx = gi + (gj - 1) * N + (gk - 1) * N * N
                pair_idx = UInt64((idx - 1) >> 1)
                g1, g2 = _threefry_gaussian(s, pair_idx)
                val += Float64(iseven(idx - 1) ? g1 : g2)
            end
            coarse[I, J, K] = val / sqrt(Float64(block)^3)
        end
    end
    return Float32.(coarse)
end

"""
Interpolate coarse grid to fine tile grid using tri-linear interpolation.
Maps global fine-grid indices to fractional coarse-grid coordinates.
"""
function _interpolate_to_tile(coarse::Array{Float32,3},
                               it::Int, jt::Int, kt::Int,
                               nsub::Int, nmesh::Int, N::Int, M::Int)
    block = N ÷ M
    nbuff = (nmesh - nsub) ÷ 2
    tile = Array{Float32,3}(undef, nmesh, nmesh, nmesh)

    i0 = (it - 1) * nsub + 1 - nbuff
    j0 = (jt - 1) * nsub + 1 - nbuff
    k0 = (kt - 1) * nsub + 1 - nbuff

    Threads.@threads for lk in 1:nmesh
        @inbounds for lj in 1:nmesh, li in 1:nmesh
            # Global fine-grid index → fractional coarse-grid coordinate
            gi = mod1(i0 + li - 1, N)
            gj = mod1(j0 + lj - 1, N)
            gk = mod1(k0 + lk - 1, N)

            # Fractional position in coarse grid (1-based, cell-centered)
            fx = (gi - 0.5) / block + 0.5
            fy = (gj - 0.5) / block + 0.5
            fz = (gk - 0.5) / block + 0.5

            tile[li, lj, lk] = _tricubic(coarse, fx, fy, fz, M)
        end
    end
    return tile
end

"""
Nearest-neighbor spreading of coarse grid to fine tile grid.
Each fine cell gets the value of the coarse cell it belongs to.
This preserves the Hoffman-Ribak property: residual = fine - spread(coarse)
has exactly zero mean per coarse cell.
"""
function _spread_coarse_to_tile(coarse::Array{Float32,3},
                                it::Int, jt::Int, kt::Int,
                                nsub::Int, nmesh::Int, N::Int, M::Int)
    block = N ÷ M
    nbuff = (nmesh - nsub) ÷ 2
    tile = Array{Float32,3}(undef, nmesh, nmesh, nmesh)

    i0 = (it - 1) * nsub + 1 - nbuff
    j0 = (jt - 1) * nsub + 1 - nbuff
    k0 = (kt - 1) * nsub + 1 - nbuff

    Threads.@threads for lk in 1:nmesh
        @inbounds for lj in 1:nmesh, li in 1:nmesh
            gi = mod1(i0 + li - 1, N)
            gj = mod1(j0 + lj - 1, N)
            gk = mod1(k0 + lk - 1, N)
            # Which coarse cell does this fine cell belong to?
            ci = (gi - 1) ÷ block + 1
            cj = (gj - 1) ÷ block + 1
            ck = (gk - 1) ÷ block + 1
            tile[li, lj, lk] = coarse[ci, cj, ck]
        end
    end
    return tile
end

"""Catmull-Rom cubic interpolation weight for offset t ∈ [0,1)."""
@inline function _catmull_rom(t::Float32)
    # Returns weights for points at -1, 0, 1, 2 relative to the cell
    t2 = t * t; t3 = t2 * t
    w0 = -0.5f0 * t3 + t2 - 0.5f0 * t
    w1 =  1.5f0 * t3 - 2.5f0 * t2 + 1.0f0
    w2 = -1.5f0 * t3 + 2.0f0 * t2 + 0.5f0 * t
    w3 =  0.5f0 * t3 - 0.5f0 * t2
    return (w0, w1, w2, w3)
end

"""Tri-cubic (Catmull-Rom) interpolation with periodic wrapping on M³ grid."""
function _tricubic(field::Array{Float32,3}, fx, fy, fz, M)
    ix = floor(Int, fx); dx = Float32(fx - ix)
    iy = floor(Int, fy); dy = Float32(fy - iy)
    iz = floor(Int, fz); dz = Float32(fz - iz)

    wx = _catmull_rom(dx)
    wy = _catmull_rom(dy)
    wz = _catmull_rom(dz)

    val = 0.0f0
    @inbounds for (dkk, wk) in zip(-1:2, wz)
        kk = mod1(iz + dkk, M)
        for (djj, wj) in zip(-1:2, wy)
            jj = mod1(iy + djj, M)
            for (dii, wi) in zip(-1:2, wx)
                ii = mod1(ix + dii, M)
                val += wi * wj * wk * field[ii, jj, kk]
            end
        end
    end
    return val
end

"""
Generate the Hoffman-Ribak residual noise for an extended region around a tile.
The extended region has size (nmesh + 2*nshell)³, centered on the tile.
The residual = fine_noise - nearest_neighbor(coarse_block_mean) has exactly
zero mean per coarse cell, ensuring T*residual is short-range.
"""
function _generate_extended_residual(it::Int, jt::Int, kt::Int,
                                      nsub::Int, nmesh::Int, N::Int, seed::Integer,
                                      coarse_noise::Array{Float32,3}, M::Int,
                                      nshell::Int)
    block = N ÷ M
    nbuff = (nmesh - nsub) ÷ 2
    next = nmesh + 2 * nshell
    inv_sqrt_block3 = Float32(1.0 / sqrt(Float64(block)^3))
    s = UInt64(seed)

    # Origin of extended region in global coordinates
    i0 = (it - 1) * nsub + 1 - nbuff - nshell
    j0 = (jt - 1) * nsub + 1 - nbuff - nshell
    k0 = (kt - 1) * nsub + 1 - nbuff - nshell

    residual = Array{Float32,3}(undef, next, next, next)

    Threads.@threads for lk in 1:next
        @inbounds for lj in 1:next, li in 1:next
            gi = mod1(i0 + li - 1, N)
            gj = mod1(j0 + lj - 1, N)
            gk = mod1(k0 + lk - 1, N)

            # Fine noise from Threefry
            idx = gi + (gj - 1) * N + (gk - 1) * N * N
            pair_idx = UInt64((idx - 1) >> 1)
            g1, g2 = _threefry_gaussian(s, pair_idx)
            fine_val = Float32(iseven(idx - 1) ? g1 : g2)

            # Coarse block mean (nearest-neighbor)
            ci = (gi - 1) ÷ block + 1
            cj = (gj - 1) ÷ block + 1
            ck = (gk - 1) ÷ block + 1
            block_mean = coarse_noise[ci, cj, ck] * inv_sqrt_block3

            residual[li, lj, lk] = fine_val - block_mean
        end
    end
    return residual
end

"""
Isolated (non-periodic) FFT convolution of noise with transfer function √P(k).
Zero-pads to 2× size to avoid circular wrap-around.

When `nshell > 0`, the input noise is an extended (nmesh+2*nshell)³ array
centered on the tile. The result is trimmed to the tile's nmesh³ region,
capturing boundary contributions from the shell.
"""
function _isolated_convolve(noise::Array{Float32,3}, pk, boxsize_local::Float64,
                             n::Int; kernel_fn=nothing, nshell::Int=0)
    n_input = size(noise, 1)  # n (no shell) or n + 2*nshell
    dx = boxsize_local / n    # cell size (same for tile and extended region)
    boxsize_ext = n_input * dx
    n2 = 2 * n_input

    padded = zeros(Float32, n2, n2, n2)
    padded[1:n_input, 1:n_input, 1:n_input] .= noise

    padded_k = rfft(padded)

    dk = 2π / (n2 * dx)  # dk for the padded box
    kx_arr = FFTW.rfftfreq(n2, n2 * dk)
    ky_arr = FFTW.fftfreq(n2, n2 * dk)
    kz_arr = FFTW.fftfreq(n2, n2 * dk)

    for iz in 1:n2, iy in 1:n2, ix in 1:size(padded_k, 1)
        kx = Float64(kx_arr[ix]); ky = Float64(ky_arr[iy]); kz = Float64(kz_arr[iz])
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0
            padded_k[ix, iy, iz] = 0
            continue
        end
        k = sqrt(k2)
        amp = sqrt(pk(k) * dk^3 * n2^3)
        if kernel_fn !== nothing
            padded_k[ix, iy, iz] *= amp * kernel_fn(kx, ky, kz, k2)
        else
            padded_k[ix, iy, iz] *= amp
        end
    end

    result = irfft(padded_k, n2)
    # Extract tile region (offset by nshell if extended)
    s1 = nshell + 1; s2 = nshell + n
    return result[s1:s2, s1:s2, s1:s2]
end

"""Apply √P(k) convolution on a (small) periodic grid. Standard approach."""
function _periodic_convolve!(noise_k, pk, n::Int, boxsize::Float64;
                              kernel_fn=nothing)
    dk = 2π / boxsize
    kx_arr = FFTW.rfftfreq(n, n * dk)
    ky_arr = FFTW.fftfreq(n, n * dk)
    kz_arr = FFTW.fftfreq(n, n * dk)

    for iz in 1:n, iy in 1:n, ix in 1:size(noise_k, 1)
        kx = Float64(kx_arr[ix]); ky = Float64(ky_arr[iy]); kz = Float64(kz_arr[iz])
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0
            noise_k[ix, iy, iz] = 0
            continue
        end
        k = sqrt(k2)
        amp = sqrt(pk(k) * dk^3 * n^3)
        if kernel_fn !== nothing
            noise_k[ix, iy, iz] *= amp * kernel_fn(kx, ky, kz, k2)
        else
            noise_k[ix, iy, iz] *= amp
        end
    end
end

# Kernel functions for LPT
_kernel_1lpt(dim) = (kx, ky, kz, k2) -> begin
    ki = dim == 1 ? kx : dim == 2 ? ky : kz
    return im * ki / k2
end

_kernel_2lpt(dim) = (kx, ky, kz, k2) -> begin
    ki = dim == 1 ? kx : dim == 2 ? ky : kz
    return -im * ki / k2
end

_kernel_phi_ij(di, dj) = (kx, ky, kz, k2) -> begin
    ki = di == 1 ? kx : di == 2 ? ky : kz
    kj = dj == 1 ? kx : dj == 2 ? ky : kz
    return -ki * kj / k2
end

_kernel_laplacian() = (kx, ky, kz, k2) -> k2

# ============================================================
# GPU dispatch helper
# ============================================================
# `isolated_convolve_gpu` is declared as a stub at the PeakPatch top level
# and filled in by `ext/CUDAExt.jl` when CUDA.jl is loaded. We access it
# lazily (by name) so this file does not require a forward reference.
@inline _pp_parent() = parentmodule(@__MODULE__)

"""
    _isolated_convolve_dispatch(use_gpu, noise, pk, boxsize, n, kernel_id, dim1, dim2; nshell)

Shared entry point used by `run_multitile_split`. Calls the CPU
`_isolated_convolve` or the GPU `isolated_convolve_gpu` based on `use_gpu`.
`kernel_id` matches the CUDA kernel IDs:
  0=δ, 1=1LPT (uses `dim1`), 2=2LPT (uses `dim1`),
  3=φ_ij (uses `dim1`, `dim2`), 4=Laplacian.
"""
function _isolated_convolve_dispatch(use_gpu::Bool, noise::Array{Float32,3}, pk,
                                      boxsize::Float64, n::Int,
                                      kernel_id::Int, dim1::Int, dim2::Int;
                                      nshell::Int=0)
    if use_gpu
        fn = getglobal(_pp_parent(), :isolated_convolve_gpu)
        return fn(noise, pk, boxsize, n;
                  kernel_fn_id=kernel_id, dim1=dim1, dim2=dim2, nshell=nshell)
    end
    kf = kernel_id == 0 ? nothing :
         kernel_id == 1 ? _kernel_1lpt(dim1) :
         kernel_id == 2 ? _kernel_2lpt(dim1) :
         kernel_id == 3 ? _kernel_phi_ij(dim1, dim2) :
         kernel_id == 4 ? _kernel_laplacian() :
         error("bad kernel_id $kernel_id")
    return _isolated_convolve(noise, pk, boxsize, n; kernel_fn=kf, nshell=nshell)
end

"""
    _gpu_analyse_peaks_batch(delta_tile, psi_tile, psi2_tile, lapd_tile, mask,
                             peaks, Rf_per_peak, stab_gpu, ct_table_gpu,
                             ct_params, alatt, ZZon, nbuff, rmax2rs, growth_tables)

Batched GPU shell analysis: groups peaks by filter scale Rf (so `ir2min`,
`Rfclvi` are uniform per batch) and invokes `analyse_peak_gpu_cuda` once per
group. Returns a Vector of NamedTuple-like results indexed the same as `peaks`.
Slots for peaks that no-collapse-on-GPU retain `RTHL = -1.0f0`.
"""
function _gpu_analyse_peaks_batch(delta_tile::Array{Float32,3},
                                    psi_x, psi_y, psi_z,
                                    psi2_x, psi2_y, psi2_z,
                                    lapd_tile, mask::Array{Int8,3},
                                    peaks::Vector{PeakCandidate},
                                    Rf_per_peak::Vector{Float64},
                                    stab_gpu, ct_table_gpu,
                                    ct_params, alatt::Float64,
                                    ZZon::Float64, nbuff::Int,
                                    rmax2rs::Float64, growth_tables)
    npeaks = length(peaks)
    # Per-peak outputs, pre-filled with sentinel RTHL = -1
    RTHL      = fill(-1.0f0, npeaks)
    Fbarx     = zeros(Float32, npeaks)
    e_v       = zeros(Float32, npeaks)
    p_v       = zeros(Float32, npeaks)
    Srb       = zeros(Float32, npeaks)
    d2F       = zeros(Float32, npeaks)
    zvir_half = fill(-1.0f0, npeaks)
    Sbar      = zeros(Float32, 3, npeaks)
    Sbar2     = zeros(Float32, 3, npeaks)
    gradpk    = zeros(Float32, 3, npeaks)
    gradpkf   = zeros(Float32, 3, npeaks)
    gradpkrf  = zeros(Float32, 3, npeaks)
    strain_f  = zeros(Float32, 3, 3, npeaks)

    fn_gpu_multi = getglobal(_pp_parent(), :analyse_peaks_gpu_cuda_multirf)
    # Group peaks by Rf. Filters are already processed in descending-Rf order
    # in the caller, so visiting unique Rfs in any order preserves correctness
    # here (the mask side-effect accumulates idempotently).
    unique_Rfs = unique(Rf_per_peak)
    batch_idxs = Vector{Vector{Int}}()
    batches = Any[]
    for Rf in unique_Rfs
        idxs = findall(==(Rf), Rf_per_peak)
        isempty(idxs) && continue
        ir2min = min(floor(Int, (1.75 * Rf / alatt)^2),
                     floor(Int, (40.0 / alatt - 1)^2))
        pi_b = Int32[peaks[i].i for i in idxs]
        pj_b = Int32[peaks[i].j for i in idxs]
        pk_b = Int32[peaks[i].k for i in idxs]
        push!(batches, (pi_b, pj_b, pk_b, Rf, ir2min))
        push!(batch_idxs, idxs)
    end

    if isempty(batches)
        return (RTHL=RTHL, Fbarx=Fbarx, e_v=e_v, p_v=p_v, Srb=Srb, d2F=d2F,
                 zvir_half=zvir_half, Sbar=Sbar, Sbar2=Sbar2,
                 gradpk=gradpk, gradpkf=gradpkf, gradpkrf=gradpkrf,
                 strain_final=strain_f)
    end

    # psi2_* may be nothing (ilpt=1). Substitute zero arrays so the GPU
    # kernel sees the expected input shapes.
    z_arr  = psi2_x === nothing ? zeros(Float32, size(delta_tile)) : psi2_x
    z_arr2 = psi2_y === nothing ? zeros(Float32, size(delta_tile)) : psi2_y
    z_arr3 = psi2_z === nothing ? zeros(Float32, size(delta_tile)) : psi2_z

    # One GPU call, fields uploaded once, kernels launched per-batch.
    batch_res = fn_gpu_multi(
        delta_tile, psi_x, psi_y, psi_z, z_arr, z_arr2, z_arr3,
        batches, stab_gpu, ct_table_gpu,
        ct_params.X1, ct_params.X2,
        ct_params.Y1, ct_params.Y2,
        ct_params.Z1, ct_params.Z2,
        alatt, ZZon;
        growth_tables=growth_tables, rmax2rs=rmax2rs,
        lapd=lapd_tile, mask=mask, nbuff=nbuff)

    # Scatter results back to per-peak slots
    for (bi, idxs) in enumerate(batch_idxs)
        res = batch_res.results[bi]
        for (b, idx) in enumerate(idxs)
            RTHL[idx]      = res.RTHL[b]
            Fbarx[idx]     = res.Fbarx[b]
            e_v[idx]       = res.e_v[b]
            p_v[idx]       = res.p_v[b]
            Srb[idx]       = res.Srb[b]
            d2F[idx]       = res.d2F[b]
            zvir_half[idx] = res.zvir_half[b]
            for L in 1:3
                Sbar[L, idx]     = res.Sbar[L, b]
                Sbar2[L, idx]    = res.Sbar2[L, b]
                gradpk[L, idx]   = res.gradpk[L, b]
                gradpkf[L, idx]  = res.gradpkf[L, b]
                gradpkrf[L, idx] = res.gradpkrf[L, b]
                for K in 1:3
                    strain_f[L, K, idx] = res.strain_final[L, K, b]
                end
            end
        end
    end
    # Copy the mask back into the in-place mask array provided by caller
    if batch_res.mask !== nothing
        copyto!(mask, batch_res.mask)
    end

    return (RTHL=RTHL, Fbarx=Fbarx, e_v=e_v, p_v=p_v, Srb=Srb, d2F=d2F,
             zvir_half=zvir_half, Sbar=Sbar, Sbar2=Sbar2,
             gradpk=gradpk, gradpkf=gradpkf, gradpkrf=gradpkrf,
             strain_final=strain_f)
end

# ============================================================
# Main pipeline
# ============================================================

"""
    run_multitile_split(cfg; ntile, seed, verbose, coarse_factor=4)

Multi-resolution variant of `run_multitile` following the MUSIC approach
(Hahn & Abel 2011). Splits the convolution T*ξ into:

- **δ_coarse**: Long-wavelength contribution computed on a coarse M³ grid
  (M = ntile × coarse_factor), using block-averaged fine-grid noise.
  Periodic FFT — captures modes from the entire box.

- **δ_self**: Short-wavelength contribution computed per-tile using the
  residual noise ξ₁ = ξ - interpolate(ξ₀). Isolated FFT (2× zero-padded)
  — non-periodic, avoids wrap-around artifacts.

The residual noise ξ₁ has zero mean in each coarse cell (Hoffman-Ribak
constraint), so T*ξ₁ is dominated by short-wavelength modes that decay
quickly in real space. This makes the isolated FFT accurate even on a
tile-sized domain.

Memory per tile: O(nmesh³) instead of O(N³). No MPI or distributed FFT.
"""
function run_multitile_split(cfg::PipelineConfig; ntile::Int, seed::Integer=42,
                              verbose::Bool=false, coarse_factor::Int=0,
                              coarse_grid::Int=0, use_gpu::Bool=false,
                              devices::Union{Nothing,AbstractVector{Int}}=nothing,
                              profile::Bool=false)
    # When `use_gpu=true`, the isolated-FFT convolutions in the per-tile loop
    # are routed through `isolated_convolve_gpu` (cuFFT). Shell analysis
    # remains on the CPU so the full PeakResult (incl. zvir_half) is available.
    # When `devices=[d0, d1, ...]` is provided with >=2 entries, the per-tile
    # loop is distributed across those CUDA devices via `Threads.@spawn`, one
    # worker task per device. Requires Julia launched with
    # `--threads >= length(devices)` and no `CUDA_VISIBLE_DEVICES` masking.
    if use_gpu
        haskey(ENV, "CUDA_VISIBLE_DEVICES") ||
            @info "run_multitile_split(use_gpu=true): CUDA_VISIBLE_DEVICES not set; using default GPU"
        isdefined(_pp_parent(), :isolated_convolve_gpu) || error("use_gpu=true requires CUDA.jl to be loaded (ext/CUDAExt.jl)")
    end
    if devices !== nothing && !isempty(devices)
        use_gpu || error("devices=$devices requires use_gpu=true")
        length(devices) > Threads.nthreads() &&
            @warn "devices=$devices length > Threads.nthreads()=$(Threads.nthreads()); workers will share threads and may not run concurrently"
    end
    n_workers = (use_gpu && devices !== nothing && length(devices) >= 1) ? length(devices) : 1
    # ---- Geometry ----
    nmesh = cfg.n
    nbuff = cfg.nbuff
    nsub = nmesh - 2 * nbuff
    N = nsub * ntile + 2 * nbuff
    alatt = cfg.boxsize / nmesh
    boxsize_full = N * alatt
    dcore_box = nsub * alatt

    # ---- Phase 0: Initialization ----
    Om_total = cfg.Omx + cfg.OmB
    cosmo = CosmologyParams(Om_total, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)
    growth_tables = Dlinear_tables(cosmo)
    ct_array, ct_params = read_homeltab(cfg.tabfile)
    ct = CollapseTableInterp(ct_array, ct_params)
    filters = read_filterbank(cfg.filterfile)
    sort!(filters; by=f -> -f[3])
    z_out = cfg.z_out
    a_out = 1.0 / (1.0 + z_out)
    ZZon = 1.0 + z_out
    fcrit = Float32(fsc_of_z(z_out, growth_tables))
    _, _, D_out = Dlinear_ab(a_out, growth_tables)
    Rfclmax = filters[1][3]
    nhunt = min(nbuff - 1, floor(Int, Rfclmax * 1.75 / alatt))
    shells = precompute_shells(nhunt)
    Omnr = cosmo.Om
    ioutshear = cfg.ioutshear
    wsmooth = cfg.wsmooth
    ilpt = cfg.ilpt
    pk = load_pk(cfg.pkfile)

    ievol = cfg.ievol
    obs = (cfg.cenx, cfg.ceny, cfg.cenz)
    z_max = cfg.z_max
    chi2z = ievol == 1 ? build_chi_to_z(cosmo; z_max=z_max + 1.0) : nothing

    # GPU shell-analysis tables — built once, uploaded per-call inside the
    # CUDAExt method. Only meaningful when `use_gpu` and `ievol == 0` (per-peak
    # ZZon in ievol=1 is not batched on GPU yet; those tiles fall back to CPU).
    use_gpu_shell = use_gpu && ievol == 0
    stab_gpu = use_gpu_shell ? build_shell_tables(nhunt) : nothing
    ct_table_gpu = use_gpu_shell ? Float32.(ct_array) : nothing

    # ---- Multi-resolution setup ----
    # Determine coarse grid size M. Priority: coarse_grid > coarse_factor > auto.
    # M must divide N. Optimal block (N/M) is ~nmesh/3 to nmesh/2 for best
    # halo count accuracy (residual is short-range, isolated FFT captures it well).
    if coarse_grid > 0
        M = coarse_grid
    elseif coarse_factor > 0
        M = ntile * coarse_factor
    else
        # Auto-select: find M that gives block ≈ nmesh/3
        target_block = nmesh ÷ 3
        best_M = 0; best_dist = N
        for m in 2:N÷2
            N % m == 0 || continue
            dist = abs(N ÷ m - target_block)
            if dist < best_dist
                best_dist = dist; best_M = m
            end
        end
        M = best_M
    end
    @assert N % M == 0 "N=$N must be divisible by M=$M"
    block = N ÷ M

    verbose && @info "Phase 0: N=$N, ntile=$ntile, M_coarse=$M, block=$block"

    # ---- Phase 1a: Coarse grid fields (global, cheap) ----
    coarse_noise = _downsample_noise(N, M, seed)
    coarse_k = rfft(coarse_noise)

    # δ_coarse: convolve coarse noise with full T(k)
    delta_coarse_k = copy(coarse_k)
    _periodic_convolve!(delta_coarse_k, pk, M, boxsize_full)
    delta_coarse = irfft(delta_coarse_k, M)

    # ψ_coarse: 1LPT displacements on coarse grid
    psi_coarse = Vector{Array{Float32,3}}(undef, 3)
    for dim in 1:3
        psi_k = copy(coarse_k)
        _periodic_convolve!(psi_k, pk, M, boxsize_full; kernel_fn=_kernel_1lpt(dim))
        psi_coarse[dim] = irfft(psi_k, M)
    end

    # 2LPT on coarse grid (src2 via trace identity, then displacement)
    psi2_coarse = nothing
    if ilpt >= 2
        src2 = zeros(Float32, M, M, M)
        # delta² / 2
        src2 .+= delta_coarse .^ 2 .* 0.5f0
        for d in 1:3
            phi_k = copy(coarse_k)
            _periodic_convolve!(phi_k, pk, M, boxsize_full; kernel_fn=_kernel_phi_ij(d, d))
            phi_r = irfft(phi_k, M)
            src2 .-= phi_r .^ 2 .* 0.5f0
        end
        for (di, dj) in ((1,2), (1,3), (2,3))
            phi_k = copy(coarse_k)
            _periodic_convolve!(phi_k, pk, M, boxsize_full; kernel_fn=_kernel_phi_ij(di, dj))
            phi_r = irfft(phi_k, M)
            src2 .-= phi_r .^ 2
        end
        src2_k = rfft(src2)
        psi2_coarse = Vector{Array{Float32,3}}(undef, 3)
        for dim in 1:3
            psi2_k = copy(src2_k)
            _periodic_convolve!(psi2_k, pk, M, boxsize_full; kernel_fn=_kernel_2lpt(dim))
            psi2_coarse[dim] = irfft(psi2_k, M)
        end
    end

    # Laplacian on coarse grid
    lapd_coarse = nothing
    if ioutshear >= 1
        lapd_k = copy(coarse_k)
        _periodic_convolve!(lapd_k, pk, M, boxsize_full; kernel_fn=_kernel_laplacian())
        lapd_coarse = irfft(lapd_k, M)
    end

    # Keep coarse noise for residual computation per tile
    coarse_k = nothing

    verbose && @info "Phase 1a: coarse fields done (M=$M)"

    # ---- Tile list ----
    tile_ids = NTuple{3,Int}[]
    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        push!(tile_ids, (it, jt, kt))
    end
    if ievol == 1
        chi_max = chi(z_max, cosmo)
        half = dcore_box / 2.0
        filter!(tile_ids) do tid
            it, jt, kt = tid
            xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)
            dx = max(abs(xbx - obs[1]) - half, 0.0)
            dy = max(abs(ybx - obs[2]) - half, 0.0)
            dz = max(abs(zbx - obs[3]) - half, 0.0)
            sqrt(dx^2 + dy^2 + dz^2) <= chi_max
        end
    end

    # Per-worker accumulators: when n_workers > 1 each dispatched task pushes
    # into its own halo vector to avoid cross-thread contention. For the
    # single-worker path these are length-1 vectors of the originals.
    local_halos_basic = [HaloRecord[] for _ in 1:n_workers]
    local_halos_ext   = [ExtHaloRecord[] for _ in 1:n_workers]

    # Per-worker timing dicts (only populated when profile=true).
    _new_timing_dict() = Dict{String,Float64}(
        "01_residual_gen" => 0.0,
        "02_iso_fft_delta" => 0.0,
        "03_interp_coarse_delta" => 0.0,
        "04_iso_fft_psi1" => 0.0,
        "05_interp_coarse_psi1" => 0.0,
        "06_2lpt_periodic_fft" => 0.0,
        "07_laplacian_periodic_fft" => 0.0,
        "08_peak_find" => 0.0,
        "09_shell_analysis" => 0.0,
        "10_record_packing" => 0.0,
    )
    local_timings = [_new_timing_dict() for _ in 1:n_workers]

    boxsize_local = nmesh * alatt

    # Per-tile work wrapped in a closure. Captures all enclosing read-only
    # state; writes go into `local_halos_basic[wid]` / `local_halos_ext[wid]`
    # / `local_timings[wid]`. Safe for concurrent invocation across workers
    # as long as each worker uses its own `wid`.
    process_tile! = function (wid::Int, ti::Int, tid::NTuple{3,Int})
        it, jt, kt = tid
        halos_basic = local_halos_basic[wid]
        halos_ext   = local_halos_ext[wid]
        timings     = local_timings[wid]
        _tic()  = profile ? time() : 0.0
        _toc!(key::String, t0::Float64) = profile ? (timings[key] += time() - t0) : nothing

        # ---- Phase 1b: Tile-local fields (MUSIC decomposition) ----

        # 1. Generate residual noise for this tile.
        t0 = _tic()
        residual = _generate_extended_residual(it, jt, kt, nsub, nmesh, N, seed,
                                                coarse_noise, M, 0)
        _toc!("01_residual_gen", t0)

        # 2. δ_self via isolated FFT of residual noise
        t0 = _tic()
        delta_self = _isolated_convolve_dispatch(use_gpu, residual, pk,
                                                  boxsize_local, nmesh, 0, 0, 0)
        _toc!("02_iso_fft_delta", t0)

        # 3. δ_coarse interpolated from global coarse field
        t0 = _tic()
        delta_long = if use_gpu
            fn = getglobal(_pp_parent(), :interpolate_to_tile_gpu)
            fn(delta_coarse, it, jt, kt, nsub, nmesh, N, M)
        else
            _interpolate_to_tile(delta_coarse, it, jt, kt, nsub, nmesh, N, M)
        end
        _toc!("03_interp_coarse_delta", t0)

        # 4. Combined tile field
        delta_tile = delta_self .+ delta_long
        delta_self = nothing; delta_long = nothing

        # 1LPT displacements (same residual + coarse split)
        psi_tile = Vector{Array{Float32,3}}(undef, 3)
        for dim in 1:3
            t0 = _tic()
            psi_self = _isolated_convolve_dispatch(use_gpu, residual, pk,
                                                    boxsize_local, nmesh, 1, dim, 0)
            _toc!("04_iso_fft_psi1", t0)
            t0 = _tic()
            psi_long = if use_gpu
                fn = getglobal(_pp_parent(), :interpolate_to_tile_gpu)
                fn(psi_coarse[dim], it, jt, kt, nsub, nmesh, N, M)
            else
                _interpolate_to_tile(psi_coarse[dim], it, jt, kt, nsub, nmesh, N, M)
            end
            _toc!("05_interp_coarse_psi1", t0)
            psi_tile[dim] = psi_self .+ psi_long
        end

        residual = nothing  # free after all isolated convolutions done

        # 2LPT: compute from combined delta on tile grid (periodic FFT OK here
        # since 2LPT is a second-order correction dominated by short modes)
        psi2_tile = nothing
        if ilpt >= 2
            t0 = _tic()
            if use_gpu
                fn = getglobal(_pp_parent(), :compute_2lpt_gpu)
                psi2_tile = fn(delta_tile, nmesh, boxsize_local)
            else
                delta_tile_k = rfft(delta_tile)
                src2_local = zeros(Float32, nmesh, nmesh, nmesh)
                src2_local .+= delta_tile .^ 2 .* 0.5f0
                for d in 1:3
                    phi_k = copy(delta_tile_k)
                    _apply_kernel_inplace!(phi_k, nmesh, boxsize_local, _kernel_phi_ij(d, d))
                    src2_local .-= irfft(phi_k, nmesh) .^ 2 .* 0.5f0
                end
                for (di, dj) in ((1,2), (1,3), (2,3))
                    phi_k = copy(delta_tile_k)
                    _apply_kernel_inplace!(phi_k, nmesh, boxsize_local, _kernel_phi_ij(di, dj))
                    src2_local .-= irfft(phi_k, nmesh) .^ 2
                end
                src2_local_k = rfft(src2_local)
                psi2_tile = Vector{Array{Float32,3}}(undef, 3)
                for dim in 1:3
                    psi2_k = copy(src2_local_k)
                    _apply_kernel_inplace!(psi2_k, nmesh, boxsize_local, _kernel_2lpt(dim))
                    psi2_tile[dim] = irfft(psi2_k, nmesh)
                end
                delta_tile_k = nothing; src2_local_k = nothing
            end
            _toc!("06_2lpt_periodic_fft", t0)
        end

        # Laplacian (compute from combined delta on tile grid, simpler than split)
        lapd_tile = nothing
        if ioutshear >= 1
            t0 = _tic()
            if use_gpu
                fn = getglobal(_pp_parent(), :compute_laplacian_gpu)
                lapd_tile = fn(delta_tile, nmesh, boxsize_local)
            else
                delta_tile_k_tmp = rfft(delta_tile)
                _apply_kernel_inplace!(delta_tile_k_tmp, nmesh, boxsize_local, _kernel_laplacian())
                lapd_tile = irfft(delta_tile_k_tmp, nmesh)
                delta_tile_k_tmp = nothing
            end
            _toc!("07_laplacian_periodic_fft", t0)
        end

        # Force GC on the CPU path where intermediate FFT buffers pile up;
        # skip on the GPU path where `unsafe_free!` has already released the
        # large device allocations.
        use_gpu || GC.gc()

        # ---- Phase 2: Peak finding ----
        t0 = _tic()
        tile_masks  = zeros(Int8, nmesh, nmesh, nmesh)
        tile_peaks  = PeakCandidate[]
        tile_Rf     = Float64[]
        tile_FcRf   = Float32[]
        tile_d2Rf   = Float32[]

        xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)
        # Per-filter fcrit (ievol=1 depends on tile position; ievol=0 constant)
        fcrits_per_filter = Vector{Float32}(undef, length(filters))
        if ievol == 1
            z_tile = peak_redshift(obs[1], obs[2], obs[3], xbx, ybx, zbx, chi2z)
            fcrit_tile = Float32(fsc_of_z(z_tile, growth_tables))
            fill!(fcrits_per_filter, fcrit_tile)
        else
            fill!(fcrits_per_filter, fcrit)
        end

        if use_gpu
            fn = getglobal(_pp_parent(), :peak_find_tile_gpu)
            pf = fn(delta_tile, filters, fcrits_per_filter, tile_masks,
                    xbx, ybx, zbx, alatt, nbuff, wsmooth, ioutshear)
            tile_peaks = pf.peaks
            tile_Rf    = pf.Rf
            tile_FcRf  = pf.FcRf
            tile_d2Rf  = pf.d2Rf
            copyto!(tile_masks, pf.mask)
        else
            delta_tile_k = rfft(delta_tile)
            for ic in 1:length(filters)
                Rf = filters[ic][3]
                delta_s = smooth_field(delta_tile_k, nmesh, boxsize_local, Rf, wsmooth)

                new_peaks = find_peaks(delta_s, tile_masks, xbx, ybx, zbx,
                                        alatt, nbuff, fcrits_per_filter[ic], Rf)
                append!(tile_peaks, new_peaks)
                append!(tile_Rf, fill(Rf, length(new_peaks)))

                if !isempty(new_peaks)
                    lapd_s = ioutshear >= 1 ?
                        smooth_field(delta_tile_k, nmesh, boxsize_local, Rf, 3) : nothing
                    for pk_cand in new_peaks
                        i, j, k = pk_cand.i, pk_cand.j, pk_cand.k
                        push!(tile_FcRf, delta_s[i, j, k])
                        push!(tile_d2Rf, lapd_s !== nothing ? lapd_s[i, j, k] : 0.0f0)
                    end
                end
            end
            delta_tile_k = nothing
        end
        _toc!("08_peak_find", t0)

        # ---- Phase 3-4: Shell analysis ----
        # (Closure body: early-return instead of `continue` when no peaks.)
        isempty(tile_peaks) && return nothing

        fill!(tile_masks, 0)
        pg = PeakGrid(delta_tile,
                        psi_tile[1], psi_tile[2], psi_tile[3],
                        psi2_tile !== nothing ? psi2_tile[1] : nothing,
                        psi2_tile !== nothing ? psi2_tile[2] : nothing,
                        psi2_tile !== nothing ? psi2_tile[3] : nothing,
                        tile_masks, (nmesh, nmesh, nmesh),
                        lapd_tile)

        # If GPU shell analysis is enabled, run the batched kernel once per
        # unique filter Rf up-front; subsequent per-peak work (record packing)
        # reads from the precomputed `gpu_res` slots instead of calling CPU
        # analyse_peak.
        gpu_res = nothing
        if use_gpu_shell
            t0 = _tic()
            gpu_res = _gpu_analyse_peaks_batch(
                delta_tile,
                psi_tile[1], psi_tile[2], psi_tile[3],
                psi2_tile !== nothing ? psi2_tile[1] : nothing,
                psi2_tile !== nothing ? psi2_tile[2] : nothing,
                psi2_tile !== nothing ? psi2_tile[3] : nothing,
                lapd_tile === nothing ? zeros(Float32, nmesh, nmesh, nmesh) : lapd_tile,
                tile_masks,
                tile_peaks, tile_Rf,
                stab_gpu, ct_table_gpu, ct_params,
                alatt, ZZon, nbuff, cfg.rmax2rs, growth_tables)
            _toc!("09_shell_analysis", t0)
        end

        t_loop_start = _tic()
        shell_cpu_time = 0.0
        for idx in 1:length(tile_peaks)
            peak = tile_peaks[idx]
            Rf = tile_Rf[idx]
            ir2min = min(floor(Int, (1.75 * Rf / alatt)^2),
                         floor(Int, (40.0 / alatt - 1)^2))

            ZZon_pk = ZZon
            if ievol == 1
                z_pk = peak_redshift(obs[1], obs[2], obs[3], peak.x, peak.y, peak.z, chi2z)
                ZZon_pk = 1.0 + z_pk
                z_pk > z_max && continue
            end

            # Shell analysis: CPU per-peak, or batched-GPU slot
            result = if use_gpu_shell
                r = gpu_res
                # Package GPU batch outputs into the same interface the CPU
                # `analyse_peak` result exposes (fields are name-compatible
                # with PeakResult).
                (RTHL       = Float64(r.RTHL[idx]),
                 Fbarx      = Float64(r.Fbarx[idx]),
                 e_v        = Float64(r.e_v[idx]),
                 p_v        = Float64(r.p_v[idx]),
                 Srb        = Float64(r.Srb[idx]),
                 d2F        = Float64(r.d2F[idx]),
                 zvir_half  = Float64(r.zvir_half[idx]),
                 Sbar       = Float64.(r.Sbar[:, idx]),
                 Sbar2      = Float64.(r.Sbar2[:, idx]),
                 gradpk     = Float64.(r.gradpk[:, idx]),
                 gradpkf    = Float64.(r.gradpkf[:, idx]),
                 gradpkrf   = Float64.(r.gradpkrf[:, idx]),
                 strain_mat = Float64.(r.strain_final[:, :, idx]))
            else
                if profile
                    ts = time()
                    r = analyse_peak(pg, peak.ipp, alatt, ir2min, ZZon_pk, Rf, ct, shells;
                                      nbuff=nbuff, growth_tables=growth_tables,
                                      rmax2rs=cfg.rmax2rs)
                    shell_cpu_time += time() - ts
                    r
                else
                    analyse_peak(pg, peak.ipp, alatt, ir2min, ZZon_pk, Rf, ct, shells;
                                  nbuff=nbuff, growth_tables=growth_tables,
                                  rmax2rs=cfg.rmax2rs)
                end
            end

            result.RTHL <= 0 && continue

            a_pk = 1.0 / ZZon_pk
            _, _, D_pk = Dlinear_ab(a_pk, growth_tables)
            RTHL_phys = Float32(result.RTHL * alatt)
            Sbar_vel = result.Sbar .* D_pk

            Sbar2_vel = zeros(3)
            if ilpt >= 2
                Om_a = Omnr * a_pk^3 / (Omnr * a_pk^3 + cosmo.OL)
                Sbar2_vel = -result.Sbar2 .* (-3.0/7.0 * Om_a^(-1.0/143) * D_pk^2)
            end

            if ioutshear >= 1
                sm = result.strain_mat
                push!(halos_ext, ExtHaloRecord(
                    Float32(peak.x), Float32(peak.y), Float32(peak.z),
                    Float32(Sbar_vel[1]), Float32(Sbar_vel[2]), Float32(Sbar_vel[3]),
                    RTHL_phys,
                    Float32(Sbar2_vel[1]), Float32(Sbar2_vel[2]), Float32(Sbar2_vel[3]),
                    Float32(result.Fbarx),
                    Float32(result.e_v), Float32(result.p_v),
                    Float32(sm[1,1]), Float32(sm[2,2]), Float32(sm[3,3]),
                    Float32(sm[2,3]), Float32(sm[1,3]), Float32(sm[1,2]),
                    Float32(result.d2F), Float32(result.zvir_half),
                    Float32(result.gradpk[1]), Float32(result.gradpk[2]), Float32(result.gradpk[3]),
                    Float32(result.gradpkf[1]), Float32(result.gradpkf[2]), Float32(result.gradpkf[3]),
                    Float32(Rf), tile_FcRf[idx], tile_d2Rf[idx],
                    Float32(result.gradpkrf[1]), Float32(result.gradpkrf[2]), Float32(result.gradpkrf[3])
                ))
            else
                push!(halos_basic, HaloRecord(
                    Float32(peak.x), Float32(peak.y), Float32(peak.z),
                    Float32(Sbar_vel[1]), Float32(Sbar_vel[2]), Float32(Sbar_vel[3]),
                    RTHL_phys,
                    Float32(Sbar2_vel[1]), Float32(Sbar2_vel[2]), Float32(Sbar2_vel[3]),
                    Float32(result.Fbarx)
                ))
            end
        end
        if profile
            t_loop_total = time() - t_loop_start
            # CPU path: the shell_cpu_time sum is the analyse_peak share;
            #           remainder is record packing.
            # GPU path: shell time is already in "09_shell_analysis" (the batch
            #           call), so the full loop is post-shell work (packing).
            if use_gpu_shell
                timings["10_record_packing"] += t_loop_total
            else
                timings["09_shell_analysis"] += shell_cpu_time
                timings["10_record_packing"] += t_loop_total - shell_cpu_time
            end
        end

        n_halos = ioutshear >= 1 ? length(halos_ext) : length(halos_basic)
        verbose && @info "  Tile $ti/$(length(tile_ids)) ($it,$jt,$kt) [w=$wid]: $(length(tile_peaks)) peaks, $n_halos worker-local halos"
        return nothing
    end  # end of process_tile!

    # ---- Dispatch: sequential (1 worker) or multi-GPU parallel (n_workers > 1) ----
    if n_workers > 1
        # Round-robin partition of tile_ids across workers (load-balance by
        # interleaving; small tiles won't cluster on one worker).
        wid_of = Int[mod1(i, n_workers) for i in 1:length(tile_ids)]
        tasks = Task[]
        for wid in 1:n_workers
            my_indices = [i for i in 1:length(tile_ids) if wid_of[i] == wid]
            my_device  = devices[wid]
            t = Threads.@spawn begin
                # Bind this task to its assigned CUDA device via the
                # CUDAExt-provided stub (avoids importing CUDA into
                # PeakPatch.jl core).
                set_dev = getglobal(_pp_parent(), :set_cuda_device!)
                set_dev(my_device)
                for idx in my_indices
                    process_tile!(wid, idx, tile_ids[idx])
                end
            end
            push!(tasks, t)
        end
        foreach(wait, tasks)
    else
        for (ti, tid) in enumerate(tile_ids)
            process_tile!(1, ti, tid)
        end
    end

    # Merge per-worker halo vectors. Ordering is not halo-canonical (pksc has
    # no ordering requirement); tests sort by coordinates before comparing.
    halos_basic = reduce(vcat, local_halos_basic; init=HaloRecord[])
    halos_ext   = reduce(vcat, local_halos_ext;   init=ExtHaloRecord[])

    if profile
        # Per-stage aggregation: max across workers approximates wall time for
        # that stage (tiles processed in parallel); sum is total GPU-seconds.
        keys_sorted = sort!(collect(keys(local_timings[1])))
        agg_max = Dict(k => maximum(td[k] for td in local_timings) for k in keys_sorted)
        agg_sum = Dict(k => sum(td[k] for td in local_timings)     for k in keys_sorted)
        wall_proxy = sum(values(agg_max))
        total_gpu_s = sum(values(agg_sum))
        @info "run_multitile_split profile (use_gpu=$use_gpu, ntile=$ntile, nmesh=$(cfg.n), n_workers=$n_workers)"
        for key in keys_sorted
            tmax = agg_max[key]; tsum = agg_sum[key]
            pct = 100 * tmax / max(wall_proxy, 1e-9)
            @info @sprintf("  %-28s  max/worker=%7.2f s  (%5.1f%%)   sum=%8.2f s",
                           key, tmax, pct, tsum)
        end
        @info @sprintf("  %-28s  wall~max=%7.2f s  total-GPU-s=%8.2f s", "TOTAL_stages", wall_proxy, total_gpu_s)
    end

    halos = ioutshear >= 1 ? halos_ext : halos_basic
    verbose && @info "Done: $(length(halos)) halos"
    return halos
end

"""Apply a k-space kernel to an existing k-space array (no P(k), just the kernel)."""
function _apply_kernel_inplace!(arr_k, n::Int, boxsize::Float64, kernel_fn)
    dk = 2π / boxsize
    kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
    ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
    kz_arr = Float64.(FFTW.fftfreq(n, n * dk))
    nyq = n ÷ 2; nk = n ÷ 2 + 1

    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = kx_arr[ix]; ky = ky_arr[iy]; kz = kz_arr[iz]
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0 || ix == nk || iy == nyq + 1 || iz == nyq + 1
            arr_k[ix, iy, iz] = 0
            continue
        end
        arr_k[ix, iy, iz] *= kernel_fn(kx, ky, kz, k2)
    end
end

"""
    compare_fields_split(pk, N, boxsize, seed, ntile, coarse_factor; verbose=false)

Generate delta and 1LPT displacement fields using both global FFT and the
multi-resolution split, then return per-tile comparison metrics.

Returns a vector of NamedTuples with fields:
  tile, corr_delta, rms_delta, corr_psi1, rms_psi1, std_global_delta
"""
function compare_fields_split(pk, N::Int, boxsize::Float64, seed::Int,
                               ntile::Int, coarse_factor::Int;
                               nbuff::Int=8, coarse_grid::Int=0, verbose::Bool=false)
    # Infer tile geometry: N = nsub*ntile + 2*nbuff → nsub = (N - 2*nbuff)/ntile
    nsub = (N - 2 * nbuff) ÷ ntile
    nmesh = nsub + 2 * nbuff
    @assert nsub * ntile + 2 * nbuff == N "N=$N not consistent with ntile=$ntile, nbuff=$nbuff"

    M = coarse_grid > 0 ? coarse_grid : ntile * coarse_factor
    @assert N % M == 0 "N=$N must be divisible by M=$M"
    block = N ÷ M
    alatt = boxsize / N

    verbose && @info "compare_fields_split: N=$N, ntile=$ntile, nmesh=$nmesh, nsub=$nsub, M=$M, block=$block"

    # ---- Global FFT reference ----
    delta_global = generate_grf(N, pk, boxsize, seed)
    delta_k_global = rfft(delta_global)
    psi_global = displacements_1lpt(delta_k_global, N, boxsize)

    # ---- Multi-resolution split ----
    coarse_noise = _downsample_noise(N, M, seed)
    coarse_k = rfft(coarse_noise)

    # Coarse delta
    delta_coarse_k = copy(coarse_k)
    _periodic_convolve!(delta_coarse_k, pk, M, boxsize)
    delta_coarse = irfft(delta_coarse_k, M)

    # Coarse 1LPT displacements
    psi_coarse = Vector{Array{Float32,3}}(undef, 3)
    for dim in 1:3
        psi_k = copy(coarse_k)
        _periodic_convolve!(psi_k, pk, M, boxsize; kernel_fn=_kernel_1lpt(dim))
        psi_coarse[dim] = irfft(psi_k, M)
    end

    boxsize_tile = nmesh * alatt

    results = NamedTuple[]

    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        # Extract global reference for the CORE of this tile (no buffer wrapping issues)
        ig0 = (it - 1) * nsub + 1  # global start of core
        jg0 = (jt - 1) * nsub + 1
        kg0 = (kt - 1) * nsub + 1
        ref_delta_core = delta_global[ig0:ig0+nsub-1, jg0:jg0+nsub-1, kg0:kg0+nsub-1]
        ref_psi1_core  = psi_global[1][ig0:ig0+nsub-1, jg0:jg0+nsub-1, kg0:kg0+nsub-1]

        # Residual noise (tile-only, no boundary shell)
        residual = _generate_extended_residual(it, jt, kt, nsub, nmesh, N, seed,
                                                coarse_noise, M, 0)

        # δ_self (isolated) + δ_coarse (interpolated)
        delta_self = _isolated_convolve(residual, pk, boxsize_tile, nmesh)
        delta_long = _interpolate_to_tile(delta_coarse, it, jt, kt, nsub, nmesh, N, M)
        delta_split = delta_self .+ delta_long

        # ψ₁_self (isolated) + ψ₁_coarse (interpolated)
        psi1_self = _isolated_convolve(residual, pk, boxsize_tile, nmesh;
                                        kernel_fn=_kernel_1lpt(1))
        psi1_long = _interpolate_to_tile(psi_coarse[1], it, jt, kt, nsub, nmesh, N, M)
        psi1_split = psi1_self .+ psi1_long

        # Core-only comparison (skip buffer regions)
        core = (nbuff+1):(nbuff+nsub)
        d_ref = Float64.(ref_delta_core)
        d_spl = Float64.(delta_split[core, core, core])
        p_ref = Float64.(ref_psi1_core)
        p_spl = Float64.(psi1_split[core, core, core])

        corr_d = _correlation(d_ref, d_spl)
        rms_d  = _rms_frac(d_ref, d_spl)
        corr_p = _correlation(p_ref, p_spl)
        rms_p  = _rms_frac(p_ref, p_spl)
        std_d  = sqrt(sum(d_ref .^ 2) / length(d_ref))

        push!(results, (tile=(it,jt,kt), corr_delta=corr_d, rms_delta=rms_d,
                         corr_psi1=corr_p, rms_psi1=rms_p, std_global_delta=std_d))

        verbose && @info "  Tile ($it,$jt,$kt): corr_δ=$(round(corr_d; digits=6)), " *
                         "rms_δ=$(round(rms_d*100; digits=2))%, " *
                         "corr_ψ₁=$(round(corr_p; digits=6)), " *
                         "rms_ψ₁=$(round(rms_p*100; digits=2))%"
    end

    return results
end

function _correlation(a, b)
    ma = sum(a) / length(a); mb = sum(b) / length(b)
    da = a .- ma; db = b .- mb
    return sum(da .* db) / sqrt(sum(da .^ 2) * sum(db .^ 2))
end

function _rms_frac(a, b)
    diff = a .- b
    rms_diff = sqrt(sum(diff .^ 2) / length(diff))
    rms_a = sqrt(sum(a .^ 2) / length(a))
    return rms_a > 0 ? rms_diff / rms_a : 0.0
end

end # module MultiResolution
