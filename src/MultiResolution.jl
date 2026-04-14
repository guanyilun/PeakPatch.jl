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

using ..Cosmology: CosmologyParams, Dlinear_tables, Dlinear_ab, chi,
    build_chi_to_z, peak_redshift
using ..PowerSpectrum: load_pk
using ..RandomField: _threefry_gaussian
using ..LPT
using ..Filters: smooth_field, gaussian_window_fortran, tophat_window, read_filterbank
using ..PeakFind: PeakCandidate, find_peaks
using ..RadialShell: PeakGrid, precompute_shells, analyse_peak, fsc_of_z
using ..Parameters: PipelineConfig
using ..Catalog: HaloRecord, ExtHaloRecord
using ..CollapseTable: CollapseTableInterp, read_homeltab
using ..MultiTile: tile_center, extract_tile

export run_multitile_split

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

    for lk in 1:nmesh, lj in 1:nmesh, li in 1:nmesh
        # Global fine-grid index → fractional coarse-grid coordinate
        gi = mod1(i0 + li - 1, N)
        gj = mod1(j0 + lj - 1, N)
        gk = mod1(k0 + lk - 1, N)

        # Fractional position in coarse grid (1-based, cell-centered)
        fx = (gi - 0.5) / block + 0.5
        fy = (gj - 0.5) / block + 0.5
        fz = (gk - 0.5) / block + 0.5

        tile[li, lj, lk] = _trilinear(coarse, fx, fy, fz, M)
    end
    return tile
end

"""Tri-linear interpolation with periodic wrapping on M³ grid."""
function _trilinear(field::Array{Float32,3}, fx, fy, fz, M)
    ix = floor(Int, fx); dx = Float32(fx - ix)
    iy = floor(Int, fy); dy = Float32(fy - iy)
    iz = floor(Int, fz); dz = Float32(fz - iz)

    i0 = mod1(ix, M);   i1 = mod1(ix + 1, M)
    j0 = mod1(iy, M);   j1 = mod1(iy + 1, M)
    k0 = mod1(iz, M);   k1 = mod1(iz + 1, M)

    @inbounds begin
        c000 = field[i0, j0, k0]; c100 = field[i1, j0, k0]
        c010 = field[i0, j1, k0]; c110 = field[i1, j1, k0]
        c001 = field[i0, j0, k1]; c101 = field[i1, j0, k1]
        c011 = field[i0, j1, k1]; c111 = field[i1, j1, k1]
    end

    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c01 = c001 * (1 - dx) + c101 * dx
    c11 = c011 * (1 - dx) + c111 * dx

    c0 = c00 * (1 - dy) + c10 * dy
    c1 = c01 * (1 - dy) + c11 * dy

    return c0 * (1 - dz) + c1 * dz
end

"""
Isolated (non-periodic) FFT convolution of noise with transfer function √P(k).
Zero-pads to 2× size to avoid circular wrap-around.
Returns the convolved field trimmed back to the original size.
"""
function _isolated_convolve(noise::Array{Float32,3}, pk, boxsize_local::Float64,
                             n::Int; kernel_fn=nothing)
    # Zero-pad to 2n for isolated (non-periodic) convolution
    n2 = 2 * n
    padded = zeros(Float32, n2, n2, n2)
    padded[1:n, 1:n, 1:n] .= noise

    padded_k = rfft(padded)

    dk = 2π / (n2 * boxsize_local / n)  # dk for the padded box
    kx_arr = FFTW.rfftfreq(n2, n2 * dk)
    ky_arr = FFTW.fftfreq(n2, n2 * dk)
    kz_arr = FFTW.fftfreq(n2, n2 * dk)

    # The amplitude normalization is the same regardless of padding:
    # amp = √(P(k) · dk³ · n²) where dk and n refer to the padded grid
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
    return result[1:n, 1:n, 1:n]
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
                              verbose::Bool=false, coarse_factor::Int=4)
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

    # ---- Multi-resolution setup ----
    M = ntile * coarse_factor
    @assert N % M == 0 "N=$N must be divisible by M=$M (try coarse_factor=$(N ÷ ntile))"
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
    # (re-downsample is expensive; cache it)
    coarse_noise_cached = _downsample_noise(N, M, seed)
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

    halos_basic = HaloRecord[]
    halos_ext   = ExtHaloRecord[]

    boxsize_local = nmesh * alatt

    for (ti, tid) in enumerate(tile_ids)
        it, jt, kt = tid

        # ---- Phase 1b: Tile-local fields (MUSIC decomposition) ----

        # 1. Generate fine noise for this tile region
        tile_noise = _generate_tile_noise(it, jt, kt, nsub, nmesh, N, seed)

        # 2. Interpolate coarse noise to tile grid, compute residual
        #    Coarse noise has variance 1/cell (normalized by 1/√block³).
        #    Interpolated to fine grid and scaled back: each fine cell gets
        #    the coarse cell's contribution = coarse_val * √block³ / block³ = coarse_val / √block³.
        #    But fine noise has variance 1/cell. The residual ξ₁ = ξ_fine - spread(ξ_coarse)
        #    has zero coarse-cell mean (Hoffman-Ribak property).
        # The coarse noise (from _downsample_noise) has unit variance per cell,
        # normalized by 1/√block³.  Each coarse cell represents the mean of
        # block³ fine cells times √block³.  To get the actual cell-average of
        # fine noise (which has variance 1/block³ per coarse cell), divide by
        # √block³.  This gives the block-mean of fine noise, which is what we
        # subtract from the fine noise to get the Hoffman-Ribak residual.
        coarse_noise_interp = _interpolate_to_tile(coarse_noise_cached,
            it, jt, kt, nsub, nmesh, N, M)
        block_mean_interp = coarse_noise_interp ./ Float32(sqrt(Float64(block)^3))
        residual = tile_noise .- block_mean_interp

        tile_noise = nothing; coarse_noise_interp = nothing; block_mean_interp = nothing

        # 3. δ_self via isolated FFT of residual noise
        delta_self = _isolated_convolve(residual, pk, boxsize_local, nmesh)

        # 4. δ_coarse interpolated from global coarse field
        delta_long = _interpolate_to_tile(delta_coarse, it, jt, kt, nsub, nmesh, N, M)

        # 5. Combined tile field
        delta_tile = delta_self .+ delta_long
        delta_self = nothing; delta_long = nothing

        # 1LPT displacements (same noise-level split)
        psi_tile = Vector{Array{Float32,3}}(undef, 3)
        for dim in 1:3
            psi_self = _isolated_convolve(residual, pk, boxsize_local, nmesh;
                                           kernel_fn=_kernel_1lpt(dim))
            psi_long = _interpolate_to_tile(psi_coarse[dim], it, jt, kt, nsub, nmesh, N, M)
            psi_tile[dim] = psi_self .+ psi_long
        end

        residual = nothing  # free after all isolated convolutions done

        # 2LPT: compute from combined delta on tile grid (periodic FFT OK here
        # since 2LPT is a second-order correction dominated by short modes)
        psi2_tile = nothing
        if ilpt >= 2
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

        # Laplacian (compute from combined delta on tile grid, simpler than split)
        lapd_tile = nothing
        if ioutshear >= 1
            delta_tile_k_tmp = rfft(delta_tile)
            _apply_kernel_inplace!(delta_tile_k_tmp, nmesh, boxsize_local, _kernel_laplacian())
            lapd_tile = irfft(delta_tile_k_tmp, nmesh)
            delta_tile_k_tmp = nothing
        end

        GC.gc()

        # ---- Phase 2: Peak finding ----
        tile_masks  = zeros(Int8, nmesh, nmesh, nmesh)
        tile_peaks  = PeakCandidate[]
        tile_Rf     = Float64[]
        tile_FcRf   = Float32[]
        tile_d2Rf   = Float32[]

        delta_tile_k = rfft(delta_tile)

        for ic in 1:length(filters)
            Rf = filters[ic][3]
            delta_s = smooth_field(delta_tile_k, nmesh, boxsize_local, Rf, wsmooth)

            xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)
            fcrit_tile = fcrit
            if ievol == 1
                z_tile = peak_redshift(obs[1], obs[2], obs[3], xbx, ybx, zbx, chi2z)
                fcrit_tile = Float32(fsc_of_z(z_tile, growth_tables))
            end

            new_peaks = find_peaks(delta_s, tile_masks, xbx, ybx, zbx,
                                    alatt, nbuff, fcrit_tile, Rf)
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

        # ---- Phase 3-4: Shell analysis ----
        isempty(tile_peaks) && continue

        fill!(tile_masks, 0)
        pg = PeakGrid(delta_tile,
                        psi_tile[1], psi_tile[2], psi_tile[3],
                        psi2_tile !== nothing ? psi2_tile[1] : nothing,
                        psi2_tile !== nothing ? psi2_tile[2] : nothing,
                        psi2_tile !== nothing ? psi2_tile[3] : nothing,
                        tile_masks, (nmesh, nmesh, nmesh),
                        lapd_tile)

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

            result = analyse_peak(pg, peak.ipp, alatt, ir2min, ZZon_pk, Rf, ct, shells;
                                   nbuff=nbuff, growth_tables=growth_tables,
                                   rmax2rs=cfg.rmax2rs)

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

        n_halos = ioutshear >= 1 ? length(halos_ext) : length(halos_basic)
        verbose && @info "  Tile $ti/$(length(tile_ids)) ($it,$jt,$kt): $(length(tile_peaks)) peaks, $n_halos total halos"
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

end # module MultiResolution
