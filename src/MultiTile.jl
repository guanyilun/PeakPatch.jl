module MultiTile

import ..Cosmology: CosmologyParams, Dlinear_tables, Dlinear_ab, chi,
    ChiToZTable, build_chi_to_z, chi_to_z, peak_redshift
import ..PowerSpectrum: load_pk, load_pk_nongaussian
import ..NonGaussian: apply_fnl_correlated!, apply_fnl_uncorrelated!
import ..RandomField: generate_grf, generate_grf_lcg
import ..LPT: displacements_1lpt, displacements_2lpt,
    displacement_1lpt_component, displacement_2lpt_component,
    compute_src2_k, compute_laplacian_field
import ..Filters: read_filterbank, smooth_field
import ..PeakFind: PeakCandidate, find_peaks
import ..RadialShell: ShellCell, PeakGrid, PeakResult, no_collapse,
    precompute_shells, analyse_peak, fsc_of_z
import ..Parameters: PipelineConfig
import ..Catalog: HaloRecord, ExtHaloRecord, write_pksc, read_pksc
import ..CollapseTable: CollapseTableInterp, read_homeltab

using FFTW
using StaticArrays

"""Convert 1-based flat index to (i,j,k) for column-major n×n×n array."""
function _ipp_to_ijk(ipp::Int, n::Int)
    k = (ipp - 1) ÷ (n * n) + 1
    rem = (ipp - 1) % (n * n)
    j = rem ÷ n + 1
    i = rem % n + 1
    return (i, j, k)
end

"""Compute k²-weighted field: ifft(k² × δ_k) / n³.

Matches the Fortran convention where the 'Laplacian' field uses positive k²
(no sign flip), used for the d2F diagnostic in the extended catalog.
"""
function _compute_laplacian(delta_k, n::Int, boxsize::Float64)
    Tf = real(eltype(delta_k))
    dk = Tf(2π / boxsize)
    nk = n ÷ 2 + 1
    nyq = n ÷ 2

    lapd_k = similar(delta_k)
    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = Tf(ix - 1) * dk
        ky = Tf(iy <= nyq + 1 ? iy - 1 : iy - 1 - n) * dk
        kz = Tf(iz <= nyq + 1 ? iz - 1 : iz - 1 - n) * dk
        k2 = kx^2 + ky^2 + kz^2
        lapd_k[ix, iy, iz] = k2 * delta_k[ix, iy, iz]
    end
    lapd = irfft(lapd_k, n) ./ Tf(n)^3
    return Tf.(lapd)
end

"""
    extract_tile(field, it, jt, kt, nsub, nmesh) -> Array

Extract an nmesh³ tile cube from the global N³ field.

Tile (it,jt,kt) occupies global cells (it-1)*nsub+1 : (it-1)*nsub+nmesh
in each dimension.  The global grid N = nsub*ntile + 2*nbuff is sized so
that all tiles (including their nbuff-wide buffer zones) fit within [1, N]
without wrap-around.  FFT periodicity ensures correct buffer physics at the
domain boundaries.
"""
function extract_tile(field::Array{T,3}, it::Int, jt::Int, kt::Int,
                      nsub::Int, nmesh::Int) where T
    i1 = (it - 1) * nsub + 1
    j1 = (jt - 1) * nsub + 1
    k1 = (kt - 1) * nsub + 1
    return field[i1:i1+nmesh-1, j1:j1+nmesh-1, k1:k1+nmesh-1]
end

"""
    tile_center(it, jt, kt, ntile, dcore_box) -> (xbx, ybx, zbx)

Physical center coordinates of tile (it,jt,kt), matching Fortran's
`xBh = cenx + (k1 - cen1) * dcore_box` with cenx=0, cen1=(ntile+1)/2.
"""
function tile_center(it::Int, jt::Int, kt::Int, ntile::Int, dcore_box::Float64)
    cen = (ntile + 1) / 2
    return (
        (it - cen) * dcore_box,
        (jt - cen) * dcore_box,
        (kt - cen) * dcore_box,
    )
end

"""
    run_multitile(cfg; ntile, seed=42, verbose=false, ...) -> Vector{HaloRecord/ExtHaloRecord}

Serial multi-tile driver.  Generates a Gaussian random field on the full
`(nsub*ntile + 2*nbuff)³` periodic grid, computes LPT displacements globally,
then processes each of the `ntile³` tiles sequentially using the existing
single-tile peak finder and shell analysis.

# Grid geometry (matches Fortran hpkvd)
- `nmesh = cfg.n`  (tile size including buffer)
- `nsub  = nmesh - 2*nbuff`  (core tile size, no overlap)
- `N     = nsub*ntile + 2*nbuff`  (full FFT grid per dimension)
- `boxsize_full = N * alatt`  where `alatt = cfg.boxsize / nmesh`
- For `ntile=1`:  `N = nmesh`, `boxsize_full = cfg.boxsize` (identical to `run_tile`)

Tile (it,jt,kt) occupies global cells `(it-1)*nsub+1` to `(it-1)*nsub+nmesh`
in each dimension.  Tiles overlap by `2*nbuff` cells with neighbours; peaks
are only found in the `nsub³` core (skipping `nbuff` cells at each edge), so
no deduplication is needed.
"""
function run_multitile(cfg::PipelineConfig; ntile::Int, seed::Integer=42,
                       verbose::Bool=false, use_lcg::Bool=false,
                       fortran_compat::Bool=false)
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

    fcrit_val = fsc_of_z(z_out, growth_tables)
    _, _, D_out = Dlinear_ab(a_out, growth_tables)

    Rfclmax = filters[1][3]
    nhunt = min(nbuff - 1, floor(Int, Rfclmax * 1.75 / alatt))
    shells = precompute_shells(nhunt)

    Omnr = cosmo.Om
    vTHvir0 = 100.0 * cosmo.h * sqrt(Omnr)
    ioutshear = cfg.ioutshear
    wsmooth = cfg.wsmooth

    # Lightcone mode
    ievol = cfg.ievol
    obs = (cfg.cenx, cfg.ceny, cfg.cenz)
    z_max = cfg.z_max
    chi2z = ievol == 1 ? build_chi_to_z(cosmo; z_max=z_max + 1.0) : nothing

    verbose && @info "Phase 0: N=$N, box=$(round(boxsize_full;digits=2)), ntile=$ntile, nsub=$nsub, nmesh=$nmesh, fcrit=$fcrit_val$(ievol == 1 ? ", lightcone mode" : "")"

    # ---- Phase 1: Field generation on full grid ----
    pk = load_pk(cfg.pkfile)
    delta_full = if use_lcg
        generate_grf_lcg(N, pk, boxsize_full, seed; fortran_compat=fortran_compat)
    else
        generate_grf(N, pk, boxsize_full, seed; fortran_compat=fortran_compat)
    end
    Tf = eltype(delta_full)  # field precision
    fcrit = Tf(fcrit_val)
    delta_k_full = rfft(delta_full)

    # Non-Gaussian corrections (modes 1-2)
    if cfg.NonGauss in (1, 2)
        _, tf_interp = load_pk_nongaussian(cfg.pkfile)
        if cfg.NonGauss == 1
            apply_fnl_correlated!(delta_full, delta_k_full, pk, tf_interp,
                                  N, boxsize_full, cfg.fNL)
            delta_k_full = rfft(delta_full)
        else
            apply_fnl_uncorrelated!(delta_full, delta_k_full, pk, tf_interp,
                                    N, boxsize_full, cfg.fNL, seed;
                                    use_lcg=use_lcg)
        end
        verbose && @info "  Applied fNL=$(cfg.fNL) (mode $(cfg.NonGauss))"
    end

    psi_x_full, psi_y_full, psi_z_full = displacements_1lpt(delta_k_full, N, boxsize_full)

    ilpt = cfg.ilpt
    psi2_x_full = psi2_y_full = psi2_z_full = nothing
    if ilpt >= 2
        psi2_x_full, psi2_y_full, psi2_z_full = displacements_2lpt(delta_k_full, N, boxsize_full)
    end

    lapd_full = nothing
    if ioutshear >= 1
        lapd_full = _compute_laplacian(delta_k_full, N, boxsize_full)
    end

    verbose && @info "Phase 1 done: full grid fields generated"

    # ---- Phase 2: Multi-scale peak finding per tile ----
    tile_ids = NTuple{3,Int}[]
    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        push!(tile_ids, (it, jt, kt))
    end

    # Lightcone: prune tiles beyond maximum_redshift
    if ievol == 1
        chi_max = chi(z_max, cosmo)
        half = dcore_box / 2.0
        ntiles_before = length(tile_ids)
        filter!(tile_ids) do tid
            it, jt, kt = tid
            xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)
            dx = max(abs(xbx - obs[1]) - half, 0.0)
            dy = max(abs(ybx - obs[2]) - half, 0.0)
            dz = max(abs(zbx - obs[3]) - half, 0.0)
            sqrt(dx^2 + dy^2 + dz^2) <= chi_max
        end
        verbose && @info "  Lightcone: $(length(tile_ids)) of $ntiles_before tiles within z_max=$z_max"
    end

    # Per-tile state
    tile_masks      = Dict(tid => zeros(Int8, nmesh, nmesh, nmesh) for tid in tile_ids)
    tile_peaks      = Dict(tid => PeakCandidate[] for tid in tile_ids)
    tile_peak_Rf    = Dict(tid => Float64[] for tid in tile_ids)
    tile_peak_FcRf  = Dict(tid => Tf[] for tid in tile_ids)
    tile_peak_d2Rf  = Dict(tid => Tf[] for tid in tile_ids)

    for ic in 1:length(filters)
        Rf = filters[ic][3]

        # Smooth globally, then extract per tile
        delta_s_full = smooth_field(delta_k_full, N, boxsize_full, Rf, wsmooth;
                                    fortran_compat=fortran_compat)

        lapd_s_full = nothing
        if ioutshear >= 1
            lapd_s_full = smooth_field(delta_k_full, N, boxsize_full, Rf, 3;
                                       fortran_compat=fortran_compat)
        end

        filter_peak_count = 0
        for tid in tile_ids
            it, jt, kt = tid
            delta_s_tile = extract_tile(Tf.(delta_s_full), it, jt, kt, nsub, nmesh)
            xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)

            # Per-tile fcrit in lightcone mode
            fcrit_tile = fcrit
            if ievol == 1
                z_tile = peak_redshift(obs[1], obs[2], obs[3], xbx, ybx, zbx, chi2z)
                fcrit_tile = Tf(fsc_of_z(z_tile, growth_tables))
            end

            new_peaks = find_peaks(delta_s_tile, tile_masks[tid],
                                   xbx, ybx, zbx, alatt, nbuff, fcrit_tile, Rf)

            append!(tile_peaks[tid], new_peaks)
            append!(tile_peak_Rf[tid], fill(Rf, length(new_peaks)))

            # Store per-peak filter-scale values for extended catalog
            if !isempty(new_peaks)
                lapd_s_tile = nothing
                if lapd_s_full !== nothing
                    lapd_s_tile = extract_tile(Tf.(lapd_s_full), it, jt, kt, nsub, nmesh)
                end

                for pk in new_peaks
                    i, j, k = _ipp_to_ijk(pk.ipp, nmesh)
                    push!(tile_peak_FcRf[tid], delta_s_tile[i, j, k])
                    push!(tile_peak_d2Rf[tid], lapd_s_tile !== nothing ? lapd_s_tile[i, j, k] : zero(Tf))
                end
            end

            filter_peak_count += length(new_peaks)
        end

        verbose && @info "  Filter $ic: Rf=$Rf, $filter_peak_count peaks"
    end

    total_peaks = sum(length(tile_peaks[tid]) for tid in tile_ids)
    verbose && @info "Phase 2 done: $total_peaks total peaks across $(ntile^3) tiles"

    # ---- Phase 3-4: Shell analysis per tile ----
    halos_basic = HaloRecord[]
    halos_ext = ExtHaloRecord[]

    for tid in tile_ids
        it, jt, kt = tid
        peaks = tile_peaks[tid]
        isempty(peaks) && continue

        # Extract all fields for this tile
        delta_tile = extract_tile(delta_full, it, jt, kt, nsub, nmesh)
        psi_x_tile = extract_tile(psi_x_full, it, jt, kt, nsub, nmesh)
        psi_y_tile = extract_tile(psi_y_full, it, jt, kt, nsub, nmesh)
        psi_z_tile = extract_tile(psi_z_full, it, jt, kt, nsub, nmesh)

        psi2_x_tile = psi2_y_tile = psi2_z_tile = nothing
        if psi2_x_full !== nothing
            psi2_x_tile = extract_tile(psi2_x_full, it, jt, kt, nsub, nmesh)
            psi2_y_tile = extract_tile(psi2_y_full, it, jt, kt, nsub, nmesh)
            psi2_z_tile = extract_tile(psi2_z_full, it, jt, kt, nsub, nmesh)
        end

        lapd_tile = nothing
        if lapd_full !== nothing
            lapd_tile = extract_tile(lapd_full, it, jt, kt, nsub, nmesh)
        end

        # Reset mask for shell analysis (was used for peak-finding exclusion)
        mask_tile = tile_masks[tid]
        fill!(mask_tile, 0)

        pg = PeakGrid(delta_tile, psi_x_tile, psi_y_tile, psi_z_tile,
                       psi2_x_tile, psi2_y_tile, psi2_z_tile,
                       mask_tile, (nmesh, nmesh, nmesh), lapd_tile)

        # Pre-compute per-peak ZZon for lightcone mode
        npeaks_tile = length(peaks)
        tile_results = Vector{PeakResult}(undef, npeaks_tile)
        tile_Rfs = tile_peak_Rf[tid]

        ZZon_tile = Vector{Float64}(undef, npeaks_tile)
        for idx in 1:npeaks_tile
            if ievol == 1
                pk = peaks[idx]
                z_pk = peak_redshift(obs[1], obs[2], obs[3], pk.x, pk.y, pk.z, chi2z)
                ZZon_tile[idx] = 1.0 + z_pk
            else
                ZZon_tile[idx] = ZZon
            end
        end

        Threads.@threads for idx in 1:npeaks_tile
            Rf = tile_Rfs[idx]
            ir2min = min(floor(Int, (1.75 * Rf / alatt)^2),
                         floor(Int, (40.0 / alatt - 1)^2))
            tile_results[idx] = analyse_peak(pg, peaks[idx].ipp, alatt, ir2min,
                                             ZZon_tile[idx], Rf, ct, shells;
                                             nbuff=nbuff,
                                             growth_tables=growth_tables,
                                             rmax2rs=cfg.rmax2rs,
                                             fortran_compat=fortran_compat)
        end

        # Collect halos sequentially
        for idx in 1:npeaks_tile
            result = tile_results[idx]
            result.RTHL <= 0 && continue

            # Lightcone: skip peaks beyond maximum redshift
            z_pk = ZZon_tile[idx] - 1.0
            if ievol == 1 && z_pk > z_max
                continue
            end

            # Per-peak growth factor and velocity scaling
            a_pk = 1.0 / ZZon_tile[idx]
            _, _, D_pk = Dlinear_ab(a_pk, growth_tables)

            peak = peaks[idx]
            Rf = tile_Rfs[idx]
            RTHL_phys = Float32(result.RTHL * alatt)
            Sbar_vel = result.Sbar .* D_pk
            Sbar2_vel = if psi2_x_full !== nothing
                Om_a = Omnr * a_pk^3 / (Omnr * a_pk^3 + cosmo.OL)
                result.Sbar2 .* (-(-3.0/7.0 * Om_a^(-1.0/143) * D_pk^2))
            else
                @SVector zeros(3)
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
                    Float32(result.d2F),
                    Float32(result.zvir_half),
                    Float32(result.gradpk[1]), Float32(result.gradpk[2]), Float32(result.gradpk[3]),
                    Float32(result.gradpkf[1]), Float32(result.gradpkf[2]), Float32(result.gradpkf[3]),
                    Float32(Rf),
                    tile_peak_FcRf[tid][idx],
                    tile_peak_d2Rf[tid][idx],
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

        tile_halos = ioutshear >= 1 ? length(halos_ext) : length(halos_basic)
        verbose && @info "  Tile $tid: $(length(peaks)) peaks → $tile_halos halos (cumulative)"
    end

    halos = ioutshear >= 1 ? halos_ext : halos_basic
    verbose && @info "Phase 4 done: $(length(halos)) halos from $(ntile^3) tiles"

    # ---- Phase 5: Write catalog ----
    if ioutshear >= 1
        RTHLmax = isempty(halos_ext) ? Float32(0) : maximum(h.RTHL for h in halos_ext)
        write_pksc(cfg.fileout, halos_ext, RTHLmax, Float32(z_out))
    else
        RTHLmax = isempty(halos_basic) ? Float32(0) : maximum(h.RTHL for h in halos_basic)
        write_pksc(cfg.fileout, halos_basic, RTHLmax, Float32(z_out))
    end

    verbose && @info "Phase 5 done: wrote $(length(halos)) halos to $(cfg.fileout)"
    return halos
end

"""
    _extract_one_tile(field, tid, nsub, nmesh, Tf) -> Array{Tf,3}

Extract a single tile sub-array from a full-grid field.
"""
function _extract_one_tile(field::AbstractArray{<:Any,3}, tid::NTuple{3,Int},
                           nsub::Int, nmesh::Int, Tf::Type)
    it, jt, kt = tid
    return Tf.(extract_tile(field, it, jt, kt, nsub, nmesh))
end

"""
    run_multitile_lowmem(cfg; ntile, seed=42, verbose=false, ...) -> Vector{HaloRecord/ExtHaloRecord}

Memory-efficient multi-tile driver. Same physics as `run_multitile`, but
computes displacement fields per-tile on the fly during shell analysis,
rather than storing all fields on the full grid simultaneously.

Only `delta_k_full` and `src2_k` (2LPT source) are kept in memory. For each
tile's shell analysis, displacement fields are recomputed from these via
IRFFT and tile extraction. This trades extra FFTs for dramatically less memory.

Memory budget (N=1760 Float64, ntile=4, ~40 active tiles, nmesh=512):
- Standard `run_multitile`: ~500 GB peak (9 full-grid arrays)
- This function: ~140 GB peak (delta_k + src2_k + 1 temp + 8 tile arrays)
- Extra cost: ntile_active × 8 full-grid IRFFTs (~312 for 39 tiles)
"""
function run_multitile_lowmem(cfg::PipelineConfig; ntile::Int, seed::Integer=42,
                               verbose::Bool=false, use_lcg::Bool=false,
                               fortran_compat::Bool=false)
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

    fcrit_val = fsc_of_z(z_out, growth_tables)
    _, _, D_out = Dlinear_ab(a_out, growth_tables)

    Rfclmax = filters[1][3]
    nhunt = min(nbuff - 1, floor(Int, Rfclmax * 1.75 / alatt))
    shells = precompute_shells(nhunt)

    Omnr = cosmo.Om
    vTHvir0 = 100.0 * cosmo.h * sqrt(Omnr)
    ioutshear = cfg.ioutshear
    wsmooth = cfg.wsmooth
    ilpt = cfg.ilpt

    ievol = cfg.ievol
    obs = (cfg.cenx, cfg.ceny, cfg.cenz)
    z_max = cfg.z_max
    chi2z = ievol == 1 ? build_chi_to_z(cosmo; z_max=z_max + 1.0) : nothing

    verbose && @info "Phase 0 (lowmem): N=$N, box=$(round(boxsize_full;digits=2)), ntile=$ntile, nsub=$nsub, nmesh=$nmesh$(ievol == 1 ? ", lightcone mode" : "")"

    # ---- Build tile list ----
    tile_ids = NTuple{3,Int}[]
    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        push!(tile_ids, (it, jt, kt))
    end

    if ievol == 1
        chi_max = chi(z_max, cosmo)
        half = dcore_box / 2.0
        ntiles_before = length(tile_ids)
        filter!(tile_ids) do tid
            it, jt, kt = tid
            xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)
            dx = max(abs(xbx - obs[1]) - half, 0.0)
            dy = max(abs(ybx - obs[2]) - half, 0.0)
            dz = max(abs(zbx - obs[3]) - half, 0.0)
            sqrt(dx^2 + dy^2 + dz^2) <= chi_max
        end
        verbose && @info "  Lightcone: $(length(tile_ids)) of $ntiles_before tiles within z_max=$z_max"
    end

    # ---- Phase 1a: Generate delta_k_full (kept throughout) ----
    pk = load_pk(cfg.pkfile)
    delta_full = if use_lcg
        generate_grf_lcg(N, pk, boxsize_full, seed; fortran_compat=fortran_compat)
    else
        generate_grf(N, pk, boxsize_full, seed; fortran_compat=fortran_compat)
    end
    Tf = eltype(delta_full)
    fcrit = Tf(fcrit_val)
    delta_k_full = rfft(delta_full)

    # Non-Gaussian corrections
    if cfg.NonGauss in (1, 2)
        _, tf_interp = load_pk_nongaussian(cfg.pkfile)
        if cfg.NonGauss == 1
            apply_fnl_correlated!(delta_full, delta_k_full, pk, tf_interp,
                                  N, boxsize_full, cfg.fNL)
            delta_k_full = rfft(delta_full)
        else
            apply_fnl_uncorrelated!(delta_full, delta_k_full, pk, tf_interp,
                                    N, boxsize_full, cfg.fNL, seed;
                                    use_lcg=use_lcg)
        end
        verbose && @info "  Applied fNL=$(cfg.fNL) (mode $(cfg.NonGauss))"
    end

    # Free delta_full — will recompute per-tile via IRFFT(delta_k) later
    delta_full = nothing
    GC.gc()

    verbose && @info "Phase 1a done: delta_k_full in memory"

    # ---- Phase 1b: Compute 2LPT source incrementally (kept throughout) ----
    src2_k = nothing
    if ilpt >= 2
        verbose && @info "  Computing 2LPT source incrementally..."
        src2_k = compute_src2_k(delta_k_full, N, boxsize_full)
        verbose && @info "  2LPT source computed"
    end

    verbose && @info "Phase 1 done (lowmem): delta_k + src2_k in memory"

    # ---- Phase 2: Multi-scale peak finding ----
    tile_masks      = Dict(tid => zeros(Int8, nmesh, nmesh, nmesh) for tid in tile_ids)
    tile_peaks      = Dict(tid => PeakCandidate[] for tid in tile_ids)
    tile_peak_Rf    = Dict(tid => Float64[] for tid in tile_ids)
    tile_peak_FcRf  = Dict(tid => Tf[] for tid in tile_ids)
    tile_peak_d2Rf  = Dict(tid => Tf[] for tid in tile_ids)

    for ic in 1:length(filters)
        Rf = filters[ic][3]

        delta_s_full = smooth_field(delta_k_full, N, boxsize_full, Rf, wsmooth;
                                    fortran_compat=fortran_compat)

        lapd_s_full = nothing
        if ioutshear >= 1
            lapd_s_full = smooth_field(delta_k_full, N, boxsize_full, Rf, 3;
                                       fortran_compat=fortran_compat)
        end

        filter_peak_count = 0
        for tid in tile_ids
            it, jt, kt = tid
            delta_s_tile = extract_tile(Tf.(delta_s_full), it, jt, kt, nsub, nmesh)
            xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)

            fcrit_tile = fcrit
            if ievol == 1
                z_tile = peak_redshift(obs[1], obs[2], obs[3], xbx, ybx, zbx, chi2z)
                fcrit_tile = Tf(fsc_of_z(z_tile, growth_tables))
            end

            new_peaks = find_peaks(delta_s_tile, tile_masks[tid],
                                   xbx, ybx, zbx, alatt, nbuff, fcrit_tile, Rf)

            append!(tile_peaks[tid], new_peaks)
            append!(tile_peak_Rf[tid], fill(Rf, length(new_peaks)))

            if !isempty(new_peaks)
                lapd_s_tile = nothing
                if lapd_s_full !== nothing
                    lapd_s_tile = extract_tile(Tf.(lapd_s_full), it, jt, kt, nsub, nmesh)
                end
                for pk in new_peaks
                    i, j, k = _ipp_to_ijk(pk.ipp, nmesh)
                    push!(tile_peak_FcRf[tid], delta_s_tile[i, j, k])
                    push!(tile_peak_d2Rf[tid], lapd_s_tile !== nothing ? lapd_s_tile[i, j, k] : zero(Tf))
                end
            end

            filter_peak_count += length(new_peaks)
        end

        verbose && @info "  Filter $ic: Rf=$Rf, $filter_peak_count peaks"
    end

    total_peaks = sum(length(tile_peaks[tid]) for tid in tile_ids)
    verbose && @info "Phase 2 done: $total_peaks total peaks across $(length(tile_ids)) tiles"

    # ---- Phase 3-4: Shell analysis — compute fields per tile on the fly ----
    # Only delta_k_full and src2_k are in memory. For each tile, we IRFFT
    # the full grid with the appropriate kernel, extract the tile, then free.
    halos_basic = HaloRecord[]
    halos_ext = ExtHaloRecord[]

    for (itile, tid) in enumerate(tile_ids)
        it, jt, kt = tid
        peaks = tile_peaks[tid]
        isempty(peaks) && continue

        verbose && @info "  Tile $tid ($itile/$(length(tile_ids))): computing fields..."

        # Compute delta tile via IRFFT of delta_k
        delta_full_tmp = irfft(copy(delta_k_full), N)
        delta_tile = _extract_one_tile(delta_full_tmp, tid, nsub, nmesh, Tf)
        delta_full_tmp = nothing

        # 1LPT displacement tiles
        psi_x_tile = _extract_one_tile(
            displacement_1lpt_component(delta_k_full, N, boxsize_full, 1), tid, nsub, nmesh, Tf)
        psi_y_tile = _extract_one_tile(
            displacement_1lpt_component(delta_k_full, N, boxsize_full, 2), tid, nsub, nmesh, Tf)
        psi_z_tile = _extract_one_tile(
            displacement_1lpt_component(delta_k_full, N, boxsize_full, 3), tid, nsub, nmesh, Tf)

        # 2LPT displacement tiles
        psi2_x_tile = psi2_y_tile = psi2_z_tile = nothing
        if src2_k !== nothing
            psi2_x_tile = _extract_one_tile(
                displacement_2lpt_component(src2_k, N, boxsize_full, 1), tid, nsub, nmesh, Tf)
            psi2_y_tile = _extract_one_tile(
                displacement_2lpt_component(src2_k, N, boxsize_full, 2), tid, nsub, nmesh, Tf)
            psi2_z_tile = _extract_one_tile(
                displacement_2lpt_component(src2_k, N, boxsize_full, 3), tid, nsub, nmesh, Tf)
        end

        # Laplacian tile
        lapd_tile = nothing
        if ioutshear >= 1
            lapd_tile = _extract_one_tile(
                compute_laplacian_field(delta_k_full, N, boxsize_full), tid, nsub, nmesh, Tf)
        end

        GC.gc()

        # Shell analysis (same as run_multitile)
        mask_tile = tile_masks[tid]
        fill!(mask_tile, 0)

        pg = PeakGrid(delta_tile, psi_x_tile, psi_y_tile, psi_z_tile,
                       psi2_x_tile, psi2_y_tile, psi2_z_tile,
                       mask_tile, (nmesh, nmesh, nmesh), lapd_tile)

        npeaks_tile = length(peaks)
        tile_results = Vector{PeakResult}(undef, npeaks_tile)
        tile_Rfs = tile_peak_Rf[tid]

        ZZon_tile = Vector{Float64}(undef, npeaks_tile)
        for idx in 1:npeaks_tile
            if ievol == 1
                pk = peaks[idx]
                z_pk = peak_redshift(obs[1], obs[2], obs[3], pk.x, pk.y, pk.z, chi2z)
                ZZon_tile[idx] = 1.0 + z_pk
            else
                ZZon_tile[idx] = ZZon
            end
        end

        Threads.@threads for idx in 1:npeaks_tile
            Rf = tile_Rfs[idx]
            ir2min = min(floor(Int, (1.75 * Rf / alatt)^2),
                         floor(Int, (40.0 / alatt - 1)^2))
            tile_results[idx] = analyse_peak(pg, peaks[idx].ipp, alatt, ir2min,
                                             ZZon_tile[idx], Rf, ct, shells;
                                             nbuff=nbuff,
                                             growth_tables=growth_tables,
                                             rmax2rs=cfg.rmax2rs,
                                             fortran_compat=fortran_compat)
        end

        # Collect halos
        for idx in 1:npeaks_tile
            result = tile_results[idx]
            result.RTHL <= 0 && continue

            z_pk = ZZon_tile[idx] - 1.0
            if ievol == 1 && z_pk > z_max
                continue
            end

            a_pk = 1.0 / ZZon_tile[idx]
            _, _, D_pk = Dlinear_ab(a_pk, growth_tables)

            peak = peaks[idx]
            Rf = tile_Rfs[idx]
            RTHL_phys = Float32(result.RTHL * alatt)
            Sbar_vel = result.Sbar .* D_pk
            Sbar2_vel = if src2_k !== nothing
                Om_a = Omnr * a_pk^3 / (Omnr * a_pk^3 + cosmo.OL)
                result.Sbar2 .* (-(-3.0/7.0 * Om_a^(-1.0/143) * D_pk^2))
            else
                @SVector zeros(3)
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
                    Float32(result.d2F),
                    Float32(result.zvir_half),
                    Float32(result.gradpk[1]), Float32(result.gradpk[2]), Float32(result.gradpk[3]),
                    Float32(result.gradpkf[1]), Float32(result.gradpkf[2]), Float32(result.gradpkf[3]),
                    Float32(Rf),
                    tile_peak_FcRf[tid][idx],
                    tile_peak_d2Rf[tid][idx],
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

        tile_halos = ioutshear >= 1 ? length(halos_ext) : length(halos_basic)
        verbose && @info "  Tile $tid: $(length(peaks)) peaks → $tile_halos halos (cumulative)"
    end

    halos = ioutshear >= 1 ? halos_ext : halos_basic
    verbose && @info "Phase 4 done: $(length(halos)) halos from $(length(tile_ids)) tiles"

    # ---- Phase 5: Write catalog ----
    if ioutshear >= 1
        RTHLmax = isempty(halos_ext) ? Float32(0) : maximum(h.RTHL for h in halos_ext)
        write_pksc(cfg.fileout, halos_ext, RTHLmax, Float32(z_out))
    else
        RTHLmax = isempty(halos_basic) ? Float32(0) : maximum(h.RTHL for h in halos_basic)
        write_pksc(cfg.fileout, halos_basic, RTHLmax, Float32(z_out))
    end

    verbose && @info "Phase 5 done: wrote $(length(halos)) halos to $(cfg.fileout)"
    return halos
end

end # module MultiTile
