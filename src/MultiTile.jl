module MultiTile

import ..Cosmology: CosmologyParams, Dlinear_tables, Dlinear_ab
import ..PowerSpectrum: load_pk, load_pk_nongaussian
import ..NonGaussian: apply_fnl_correlated!, apply_fnl_uncorrelated!
import ..RandomField: generate_grf, generate_grf_lcg
import ..LPT: displacements_1lpt, displacements_2lpt
import ..Filters: read_filterbank, smooth_field
import ..PeakFind: PeakCandidate, find_peaks
import ..RadialShell: ShellCell, PeakGrid, PeakResult, no_collapse,
    precompute_shells, analyse_peak, fsc_of_z
import ..Parameters: SimParams
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
    run_multitile(sp; ntile, seed=42, verbose=false, ...) -> Vector{HaloRecord/ExtHaloRecord}

Serial multi-tile driver.  Generates a Gaussian random field on the full
`(nsub*ntile + 2*nbuff)³` periodic grid, computes LPT displacements globally,
then processes each of the `ntile³` tiles sequentially using the existing
single-tile peak finder and shell analysis.

# Grid geometry (matches Fortran hpkvd)
- `nmesh = sp.nlx`  (tile size including buffer)
- `nsub  = nmesh - 2*nbuff`  (core tile size, no overlap)
- `N     = nsub*ntile + 2*nbuff`  (full FFT grid per dimension)
- `boxsize_full = N * alatt`  where `alatt = dL_box / nmesh`
- For `ntile=1`:  `N = nmesh`, `boxsize_full = dL_box` (identical to `run_tile`)

Tile (it,jt,kt) occupies global cells `(it-1)*nsub+1` to `(it-1)*nsub+nmesh`
in each dimension.  Tiles overlap by `2*nbuff` cells with neighbours; peaks
are only found in the `nsub³` core (skipping `nbuff` cells at each edge), so
no deduplication is needed.
"""
function run_multitile(sp::SimParams; ntile::Int, seed::Integer=42,
                       verbose::Bool=false, use_lcg::Bool=false,
                       fortran_compat::Bool=false)
    # ---- Geometry ----
    nmesh = Int(sp.nlx)
    nbuff = Int(sp.nbuff)
    nsub = nmesh - 2 * nbuff
    N = nsub * ntile + 2 * nbuff
    alatt = Float64(sp.dL_box) / nmesh
    boxsize_full = N * alatt
    dcore_box = nsub * alatt

    # ---- Phase 0: Initialization ----
    Om_total = Float64(sp.Omx) + Float64(sp.OmB)
    cosmo = CosmologyParams(Om_total, Float64(sp.OmB), Float64(sp.Omvac),
                            Float64(sp.h), 0.965, 0.808)
    growth_tables = Dlinear_tables(cosmo)

    ct_array, ct_params = read_homeltab(sp.TabInterpFile)
    ct = CollapseTableInterp(ct_array, ct_params)

    filters = read_filterbank(sp.filterfile)
    sort!(filters; by=f -> -f[3])

    z_out = Float64(sp.global_redshift)
    a_out = 1.0 / (1.0 + z_out)
    ZZon = 1.0 + z_out

    fcrit_val = fsc_of_z(z_out, growth_tables)
    _, _, D_out = Dlinear_ab(a_out, growth_tables)

    Rfclmax = filters[1][3]
    nhunt = min(nbuff - 1, floor(Int, Rfclmax * 1.75 / alatt))
    shells = precompute_shells(nhunt)

    Omnr = cosmo.Om
    vTHvir0 = 100.0 * cosmo.h * sqrt(Omnr)
    ioutshear = Int(sp.ioutshear)
    wsmooth = Int(sp.wsmooth)

    verbose && @info "Phase 0: N=$N, box=$(round(boxsize_full;digits=2)), ntile=$ntile, nsub=$nsub, nmesh=$nmesh, fcrit=$fcrit_val"

    # ---- Phase 1: Field generation on full grid ----
    pk = load_pk(sp.pkfile)
    delta_full = if use_lcg
        generate_grf_lcg(N, pk, boxsize_full, seed; fortran_compat=fortran_compat)
    else
        generate_grf(N, pk, boxsize_full, seed; fortran_compat=fortran_compat)
    end
    Tf = eltype(delta_full)  # field precision
    fcrit = Tf(fcrit_val)
    delta_k_full = rfft(delta_full)

    # Non-Gaussian corrections (modes 1-2)
    NonGauss = Int(sp.NonGauss)
    if NonGauss in (1, 2)
        _, tf_interp = load_pk_nongaussian(sp.pkfile)
        if NonGauss == 1
            apply_fnl_correlated!(delta_full, delta_k_full, pk, tf_interp,
                                  N, boxsize_full, Float64(sp.fNL))
            delta_k_full = rfft(delta_full)
        else
            apply_fnl_uncorrelated!(delta_full, delta_k_full, pk, tf_interp,
                                    N, boxsize_full, Float64(sp.fNL), seed;
                                    use_lcg=use_lcg)
        end
        verbose && @info "  Applied fNL=$(sp.fNL) (mode $NonGauss)"
    end

    psi_x_full, psi_y_full, psi_z_full = displacements_1lpt(delta_k_full, N, boxsize_full)

    ilpt = Int(sp.ilpt)
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

            new_peaks = find_peaks(delta_s_tile, tile_masks[tid],
                                   xbx, ybx, zbx, alatt, nbuff, fcrit, Rf)

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

        # Analyse peaks in parallel within tile
        npeaks_tile = length(peaks)
        tile_results = Vector{PeakResult}(undef, npeaks_tile)
        tile_Rfs = tile_peak_Rf[tid]
        rmax2rs_f = Float64(sp.rmax2rs)

        Threads.@threads for idx in 1:npeaks_tile
            Rf = tile_Rfs[idx]
            ir2min = min(floor(Int, (1.75 * Rf / alatt)^2),
                         floor(Int, (40.0 / alatt - 1)^2))
            tile_results[idx] = analyse_peak(pg, peaks[idx].ipp, alatt, ir2min,
                                             ZZon, Rf, ct, shells;
                                             nbuff=nbuff,
                                             growth_tables=growth_tables,
                                             rmax2rs=rmax2rs_f,
                                             fortran_compat=fortran_compat)
        end

        # Pre-compute 2LPT velocity factor
        Om_a_factor = if psi2_x_full !== nothing
            Om_a = Omnr * a_out^3 / (Omnr * a_out^3 + cosmo.OL)
            -(-3.0/7.0 * Om_a^(-1.0/143) * D_out^2)
        else
            0.0
        end

        # Collect halos sequentially
        for idx in 1:npeaks_tile
            result = tile_results[idx]
            result.RTHL <= 0 && continue

            peak = peaks[idx]
            Rf = tile_Rfs[idx]
            RTHL_phys = Float32(result.RTHL * alatt)
            Sbar_vel = result.Sbar .* D_out
            Sbar2_vel = psi2_x_full !== nothing ? result.Sbar2 .* Om_a_factor : @SVector zeros(3)

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
        write_pksc(sp.fileout, halos_ext, RTHLmax, Float32(z_out))
    else
        RTHLmax = isempty(halos_basic) ? Float32(0) : maximum(h.RTHL for h in halos_basic)
        write_pksc(sp.fileout, halos_basic, RTHLmax, Float32(z_out))
    end

    verbose && @info "Phase 5 done: wrote $(length(halos)) halos to $(sp.fileout)"
    return halos
end

end # module MultiTile
