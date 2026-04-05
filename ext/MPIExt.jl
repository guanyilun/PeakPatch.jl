module MPIExt

using PeakPatch
using MPI
using PencilFFTs
using PencilFFTs.Transforms: RFFT
using PencilArrays
using FFTW

# ============================================================
# Helpers
# ============================================================

"""Factor nranks into (p1, p2) for 2D pencil decomposition."""
function _factor_procs(nranks::Int)
    p1 = isqrt(nranks)
    while nranks % p1 != 0
        p1 -= 1
    end
    return (p1, nranks ÷ p1)
end

"""Block-assign ntile³ tiles to `nranks` MPI ranks.
Returns the list of (it,jt,kt) tuples assigned to `rank` (0-based)."""
function _assign_tiles(ntile::Int, nranks::Int, rank::Int)
    total = ntile^3
    tiles_per_rank = cld(total, nranks)
    start_idx = rank * tiles_per_rank + 1
    end_idx = min((rank + 1) * tiles_per_rank, total)

    tile_ids = NTuple{3,Int}[]
    idx = 0
    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        idx += 1
        if start_idx <= idx <= end_idx
            push!(tile_ids, (it, jt, kt))
        end
    end
    return tile_ids
end

"""Fill a PencilArray from a full N³ array (all ranks must have `full_arr`)."""
function _fill_pencil!(pencil_arr, full_arr)
    gv = PencilArrays.global_view(pencil_arr)
    pen = PencilArrays.pencil(pencil_arr)
    ranges = PencilArrays.range_local(pen)
    for k in ranges[3], j in ranges[2], i in ranges[1]
        gv[i, j, k] = full_arr[i, j, k]
    end
    return pencil_arr
end

"""Extract tile subcubes from a pencil-distributed real-space field.

Each rank receives only the nmesh³ subcubes for its own tiles via point-to-point
MPI (Isend/Irecv), instead of gathering the full N³ field.  Memory per rank is
O(ntile_local × nmesh³) instead of O(N³).

Returns `Dict{NTuple{3,Int}, Array{T,3}}` mapping tile IDs to nmesh³ arrays.
"""
function _extract_my_tiles(pencil_arr, my_tiles::Vector{NTuple{3,Int}},
                           ntile::Int, nsub::Int, nmesh::Int,
                           ::Type{T}, comm::MPI.Comm) where T
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    pen = PencilArrays.pencil(pencil_arr)
    local_ranges = PencilArrays.range_local(pen)
    gv = PencilArrays.global_view(pencil_arr)

    # --- Share pencil ranges across all ranks ---
    my_flat = Int32[first(local_ranges[1]), last(local_ranges[1]),
                    first(local_ranges[2]), last(local_ranges[2]),
                    first(local_ranges[3]), last(local_ranges[3])]
    all_flat = MPI.Allgather(my_flat, comm)
    all_ranges_mat = reshape(all_flat, 6, nranks)

    # --- Build tile owner map (deterministic, all ranks agree) ---
    tile_owner = Dict{NTuple{3,Int}, Int}()
    for r in 0:nranks-1
        for tid in _assign_tiles(ntile, nranks, r)
            tile_owner[tid] = r
        end
    end

    # --- Canonical tile list (consistent ordering → consistent MPI tags) ---
    all_tiles = NTuple{3,Int}[]
    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        push!(all_tiles, (it, jt, kt))
    end

    # --- Initialize output arrays ---
    result = Dict{NTuple{3,Int}, Array{T,3}}()
    for tid in my_tiles
        result[tid] = Array{T,3}(undef, nmesh, nmesh, nmesh)
    end

    # --- Compute overlaps, post sends and receives ---
    send_reqs = MPI.Request[]
    send_bufs = Vector{T}[]
    recv_reqs = MPI.Request[]
    recv_bufs = Vector{T}[]
    recv_meta = Tuple{NTuple{3,Int}, UnitRange{Int}, UnitRange{Int}, UnitRange{Int}}[]

    for (tidx, tid) in enumerate(all_tiles)
        it, jt, kt = tid
        gi1 = (it-1)*nsub + 1; gi2 = gi1 + nmesh - 1
        gj1 = (jt-1)*nsub + 1; gj2 = gj1 + nmesh - 1
        gk1 = (kt-1)*nsub + 1; gk2 = gk1 + nmesh - 1
        owner = tile_owner[tid]

        # Overlap between my pencil slice and this tile
        i_ovlp = max(gi1, first(local_ranges[1])):min(gi2, last(local_ranges[1]))
        j_ovlp = max(gj1, first(local_ranges[2])):min(gj2, last(local_ranges[2]))
        k_ovlp = max(gk1, first(local_ranges[3])):min(gk2, last(local_ranges[3]))
        has_data = !isempty(i_ovlp) && !isempty(j_ovlp) && !isempty(k_ovlp)

        if has_data && owner == rank
            # Local copy: I own both the data and the tile
            tile = result[tid]
            for gk in k_ovlp, gj in j_ovlp, gi in i_ovlp
                tile[gi - gi1 + 1, gj - gj1 + 1, gk - gk1 + 1] = T(gv[gi, gj, gk])
            end
        elseif has_data && owner != rank
            # Pack and send to tile owner
            buf = Vector{T}(undef, length(i_ovlp) * length(j_ovlp) * length(k_ovlp))
            idx = 0
            for gk in k_ovlp, gj in j_ovlp, gi in i_ovlp
                idx += 1
                buf[idx] = T(gv[gi, gj, gk])
            end
            push!(send_bufs, buf)
            push!(send_reqs, MPI.Isend(buf, owner, tidx, comm))
        end

        if owner == rank
            # Post receives from remote ranks that have overlapping data
            for src in 0:nranks-1
                src == rank && continue
                off = src * 6
                sr = (Int(all_ranges_mat[1, src+1]):Int(all_ranges_mat[2, src+1]),
                      Int(all_ranges_mat[3, src+1]):Int(all_ranges_mat[4, src+1]),
                      Int(all_ranges_mat[5, src+1]):Int(all_ranges_mat[6, src+1]))
                si = max(gi1, first(sr[1])):min(gi2, last(sr[1]))
                sj = max(gj1, first(sr[2])):min(gj2, last(sr[2]))
                sk = max(gk1, first(sr[3])):min(gk2, last(sr[3]))

                if !isempty(si) && !isempty(sj) && !isempty(sk)
                    buf = Vector{T}(undef, length(si) * length(sj) * length(sk))
                    push!(recv_bufs, buf)
                    push!(recv_reqs, MPI.Irecv!(buf, src, tidx, comm))
                    push!(recv_meta, (tid,
                                      si .- gi1 .+ 1,
                                      sj .- gj1 .+ 1,
                                      sk .- gk1 .+ 1))
                end
            end
        end
    end

    # --- Wait for receives, unpack into tile arrays ---
    MPI.Waitall(recv_reqs)
    for (idx, (tid, di_rng, dj_rng, dk_rng)) in enumerate(recv_meta)
        tile = result[tid]
        buf = recv_bufs[idx]
        buf_idx = 0
        for dk in dk_rng, dj in dj_rng, di in di_rng
            buf_idx += 1
            tile[di, dj, dk] = buf[buf_idx]
        end
    end

    # --- Wait for sends (buffers must stay alive until completion) ---
    MPI.Waitall(send_reqs)

    return result
end

"""Convert 1-based flat index to (i,j,k) for column-major n×n×n array."""
function _ipp_to_ijk(ipp::Int, n::Int)
    k = (ipp - 1) ÷ (n * n) + 1
    rem = (ipp - 1) % (n * n)
    j = rem ÷ n + 1
    i = rem % n + 1
    return (i, j, k)
end

# ============================================================
# Distributed k-space operations
# ============================================================

"""Multiply PencilArray delta_k by √(P(k) dk³ N³) in-place (distributed)."""
function _distributed_convolve_pk!(delta_k_pencil, pk, N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = FFTW.rfftfreq(N, N * dk)
    ky_arr = FFTW.fftfreq(N, N * dk)
    kz_arr = FFTW.fftfreq(N, N * dk)

    pen = PencilArrays.pencil(delta_k_pencil)
    ranges = PencilArrays.range_local(pen)
    gv = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0
            gv[gix, giy, giz] = 0.0
            continue
        end
        k = sqrt(k2)
        amp = sqrt(pk(k) * dk^3 * N^3)
        gv[gix, giy, giz] *= amp
    end
end

"""Apply a k-space window function to a copy of delta_k (distributed).
Returns a new PencilArray with window applied."""
function _distributed_apply_window(delta_k_pencil, N::Int, boxsize::Float64,
                                   Rf::Float64, wsmooth::Int)
    dk = 2π / boxsize
    kx_arr = FFTW.rfftfreq(N, N * dk)
    ky_arr = FFTW.fftfreq(N, N * dk)
    kz_arr = FFTW.fftfreq(N, N * dk)

    smoothed_k = similar(delta_k_pencil)
    smoothed_k .= delta_k_pencil

    pen = PencilArrays.pencil(smoothed_k)
    ranges = PencilArrays.range_local(pen)
    gv = PencilArrays.global_view(smoothed_k)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k = sqrt(kx^2 + ky^2 + kz^2)
        fkR = k * Rf
        fkR == 0.0 && continue
        if wsmooth == 0
            w = PeakPatch.gaussian_window_fortran(fkR)
        elseif wsmooth == 3
            w = PeakPatch.gaussian_window_fortran(fkR) * k^2
        else
            w = PeakPatch.tophat_window(fkR)
        end
        gv[gix, giy, giz] *= w
    end
    return smoothed_k
end

"""Compute 1LPT displacement component `dim` (1=x,2=y,3=z) from delta_k.
Returns a k-space PencilArray ready for inverse FFT."""
function _distributed_1lpt_component(delta_k_pencil, dim::Int,
                                     N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = Float64.(FFTW.rfftfreq(N, N * dk))
    ky_arr = Float64.(FFTW.fftfreq(N, N * dk))
    kz_arr = Float64.(FFTW.fftfreq(N, N * dk))
    nyq = N ÷ 2
    nk = N ÷ 2 + 1

    psi_k = similar(delta_k_pencil)
    pen = PencilArrays.pencil(psi_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(psi_k)
    gv_in = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        # Zero DC and Nyquist modes
        if k2 == 0.0 || gix == nk || giy == nyq + 1 || giz == nyq + 1
            gv_out[gix, giy, giz] = 0.0
            continue
        end
        ki = dim == 1 ? kx : dim == 2 ? ky : kz
        gv_out[gix, giy, giz] = im * ki / k2 * gv_in[gix, giy, giz]
    end
    return psi_k
end

"""Compute 2LPT second-derivative potential component phi_ij(k) = -ki*kj/k² δ(k).
`di`, `dj` are dimension indices (1=x, 2=y, 3=z)."""
function _distributed_phi_ij_component(delta_k_pencil, di::Int, dj::Int,
                                       N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = Float64.(FFTW.rfftfreq(N, N * dk))
    ky_arr = Float64.(FFTW.fftfreq(N, N * dk))
    kz_arr = Float64.(FFTW.fftfreq(N, N * dk))
    nyq = N ÷ 2
    nk = N ÷ 2 + 1

    phi_k = similar(delta_k_pencil)
    pen = PencilArrays.pencil(phi_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(phi_k)
    gv_in = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0 || gix == nk || giy == nyq + 1 || giz == nyq + 1
            gv_out[gix, giy, giz] = 0.0
            continue
        end
        ki = di == 1 ? kx : di == 2 ? ky : kz
        kj = dj == 1 ? kx : dj == 2 ? ky : kz
        gv_out[gix, giy, giz] = -ki * kj / k2 * gv_in[gix, giy, giz]
    end
    return phi_k
end

"""Compute 2LPT displacement component from src2_k.
psi2_i(k) = -im * ki / k² * src2(k)."""
function _distributed_2lpt_component(src2_k_pencil, dim::Int,
                                     N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = Float64.(FFTW.rfftfreq(N, N * dk))
    ky_arr = Float64.(FFTW.fftfreq(N, N * dk))
    kz_arr = Float64.(FFTW.fftfreq(N, N * dk))
    nyq = N ÷ 2
    nk = N ÷ 2 + 1

    psi2_k = similar(src2_k_pencil)
    pen = PencilArrays.pencil(psi2_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(psi2_k)
    gv_in = PencilArrays.global_view(src2_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0 || gix == nk || giy == nyq + 1 || giz == nyq + 1
            gv_out[gix, giy, giz] = 0.0
            continue
        end
        ki = dim == 1 ? kx : dim == 2 ? ky : kz
        gv_out[gix, giy, giz] = -im * ki / k2 * gv_in[gix, giy, giz]
    end
    return psi2_k
end

"""Compute k²-weighted field in k-space (for Laplacian diagnostic)."""
function _distributed_laplacian_k(delta_k_pencil, N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = FFTW.rfftfreq(N, N * dk)
    ky_arr = FFTW.fftfreq(N, N * dk)
    kz_arr = FFTW.fftfreq(N, N * dk)

    lapd_k = similar(delta_k_pencil)
    pen = PencilArrays.pencil(lapd_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(lapd_k)
    gv_in = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        gv_out[gix, giy, giz] = k2 * gv_in[gix, giy, giz]
    end
    return lapd_k
end

# ============================================================
# Main distributed driver
# ============================================================

function PeakPatch.run_multitile_mpi(sp::PeakPatch.SimParams;
        ntile::Int, seed::Integer=42, verbose::Bool=false,
        comm::MPI.Comm=MPI.COMM_WORLD)

    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    # ---- Geometry (all ranks, no comm) ----
    nmesh  = Int(sp.nlx)
    nbuff  = Int(sp.nbuff)
    nsub   = nmesh - 2 * nbuff
    N      = nsub * ntile + 2 * nbuff
    alatt  = Float64(sp.dL_box) / nmesh
    boxsize_full = N * alatt
    dcore_box = nsub * alatt

    # ---- Phase 0: Initialization (all ranks, no comm) ----
    Om_total = Float64(sp.Omx) + Float64(sp.OmB)
    cosmo = PeakPatch.CosmologyParams(Om_total, Float64(sp.OmB), Float64(sp.Omvac),
                                       Float64(sp.h), 0.965, 0.808)
    growth_tables = PeakPatch.Dlinear_tables(cosmo)

    ct_array, ct_params = PeakPatch.read_homeltab(sp.TabInterpFile)
    ct = PeakPatch.CollapseTableInterp(ct_array, ct_params)

    filters = PeakPatch.read_filterbank(sp.filterfile)
    sort!(filters; by=f -> -f[3])

    z_out  = Float64(sp.global_redshift)
    a_out  = 1.0 / (1.0 + z_out)
    ZZon   = 1.0 + z_out
    fcrit  = Float32(PeakPatch.fsc_of_z(z_out, growth_tables))
    _, _, D_out = PeakPatch.Dlinear_ab(a_out, growth_tables)

    Rfclmax = filters[1][3]
    nhunt  = min(nbuff - 1, floor(Int, Rfclmax * 1.75 / alatt))
    shells = PeakPatch.precompute_shells(nhunt)

    Omnr    = cosmo.Om
    vTHvir0 = 100.0 * cosmo.h * sqrt(Omnr)
    ioutshear = Int(sp.ioutshear)
    wsmooth   = Int(sp.wsmooth)
    ilpt      = Int(sp.ilpt)

    NonGauss = Int(sp.NonGauss)
    NonGauss != 0 && error("run_multitile_mpi: NonGauss=$NonGauss not yet supported (only 0)")

    # Lightcone mode
    ievol = Int(sp.ievol)
    obs = (Float64(sp.cenx), Float64(sp.ceny), Float64(sp.cenz))
    z_max = Float64(sp.maximum_redshift)
    chi2z = ievol == 1 ? PeakPatch.build_chi_to_z(cosmo; z_max=z_max + 1.0) : nothing

    my_tiles = _assign_tiles(ntile, nranks, rank)

    # Lightcone: prune tiles beyond maximum_redshift
    if ievol == 1
        chi_max = PeakPatch.chi(z_max, cosmo)
        half = dcore_box / 2.0
        filter!(my_tiles) do tid
            it, jt, kt = tid
            xbx, ybx, zbx = PeakPatch.tile_center(it, jt, kt, ntile, dcore_box)
            dx = max(abs(xbx - obs[1]) - half, 0.0)
            dy = max(abs(ybx - obs[2]) - half, 0.0)
            dz = max(abs(zbx - obs[3]) - half, 0.0)
            sqrt(dx^2 + dy^2 + dz^2) <= chi_max
        end
    end

    rank == 0 && verbose && @info "Phase 0: N=$N, ntile=$ntile, nranks=$nranks, tiles/rank=$(length(my_tiles))$(ievol == 1 ? ", lightcone mode" : "")"

    # ---- PencilFFT setup ----
    proc_dims = _factor_procs(nranks)
    plan = PencilFFTPlan((N, N, N), RFFT(), proc_dims, comm)

    # ---- Phase 1a: Distributed noise generation (Threefry counter-based RNG) ----
    pk = PeakPatch.load_pk(sp.pkfile)

    noise_pencil = PencilFFTs.allocate_input(plan)
    pen = PencilArrays.pencil(noise_pencil)
    ranges = PencilArrays.range_local(pen)
    PeakPatch.fill_noise_threefry_region!(parent(noise_pencil), ranges, N, seed)

    rank == 0 && verbose && @info "Phase 1a: distributed noise generated (Threefry)"

    # ---- Phase 1b: Distributed forward FFT + P(k) convolution ----
    delta_k_pencil = plan * noise_pencil
    noise_pencil = nothing; GC.gc()

    _distributed_convolve_pk!(delta_k_pencil, pk, N, boxsize_full)

    # Save delta_k for reuse in smoothing loop and LPT
    delta_k_saved = similar(delta_k_pencil)
    delta_k_saved .= delta_k_pencil

    rank == 0 && verbose && @info "Phase 1b: forward FFT + P(k) done"

    # ---- Phase 1c: delta tiles (distributed inverse FFT → tile extraction) ----
    delta_real_pencil = plan \ delta_k_pencil
    delta_tiles = _extract_my_tiles(delta_real_pencil, my_tiles, ntile, nsub, nmesh,
                                    Float32, comm)
    delta_real_pencil = nothing; GC.gc()

    # ---- Phase 1d: 1LPT displacement tiles (distributed) ----
    psi_tiles = Dict{Int, Dict{NTuple{3,Int}, Array{Float32,3}}}()
    for dim in 1:3
        psi_k = _distributed_1lpt_component(delta_k_saved, dim, N, boxsize_full)
        psi_real = plan \ psi_k
        psi_tiles[dim] = _extract_my_tiles(psi_real, my_tiles, ntile, nsub, nmesh,
                                           Float32, comm)
    end

    rank == 0 && verbose && @info "Phase 1d: 1LPT done"

    # ---- Phase 1e: 2LPT displacement tiles (distributed) ----
    psi2_tiles = nothing
    if ilpt >= 2
        # Compute phi_ij in k-space → inverse FFT → real-space PencilArrays
        phi_pairs = [(1,1), (2,2), (3,3), (1,2), (1,3), (2,3)]
        phi_real = Dict{Tuple{Int,Int}, Any}()
        for (di, dj) in phi_pairs
            phi_k = _distributed_phi_ij_component(delta_k_saved, di, dj, N, boxsize_full)
            phi_real[(di,dj)] = plan \ phi_k
        end

        # Compute src2 in distributed real space
        src2_pencil = PencilFFTs.allocate_input(plan)
        p11 = parent(phi_real[(1,1)]); p22 = parent(phi_real[(2,2)]); p33 = parent(phi_real[(3,3)])
        p12 = parent(phi_real[(1,2)]); p13 = parent(phi_real[(1,3)]); p23 = parent(phi_real[(2,3)])
        parent(src2_pencil) .= p11 .* p22 .- p12 .^ 2 .+
                                p11 .* p33 .- p13 .^ 2 .+
                                p22 .* p33 .- p23 .^ 2
        phi_real = nothing; GC.gc()

        # Forward FFT src2
        src2_k_pencil = plan * src2_pencil
        src2_pencil = nothing; GC.gc()

        # Compute psi2_i tiles from src2_k
        psi2_tiles = Dict{Int, Dict{NTuple{3,Int}, Array{Float32,3}}}()
        for dim in 1:3
            psi2_k = _distributed_2lpt_component(src2_k_pencil, dim, N, boxsize_full)
            psi2_real = plan \ psi2_k
            psi2_tiles[dim] = _extract_my_tiles(psi2_real, my_tiles, ntile, nsub, nmesh,
                                                Float32, comm)
        end
        src2_k_pencil = nothing; GC.gc()

        rank == 0 && verbose && @info "Phase 1e: 2LPT done"
    end

    # ---- Phase 1f: Laplacian tiles (distributed) ----
    lapd_tiles = nothing
    if ioutshear >= 1
        lapd_k = _distributed_laplacian_k(delta_k_saved, N, boxsize_full)
        lapd_real = plan \ lapd_k
        lapd_tiles_f64 = _extract_my_tiles(lapd_real, my_tiles, ntile, nsub, nmesh,
                                           Float64, comm)
        # Match serial: divide by N³ (serial _compute_laplacian does irfft / n³)
        lapd_tiles = Dict{NTuple{3,Int}, Array{Float32,3}}()
        for (tid, arr) in lapd_tiles_f64
            lapd_tiles[tid] = Float32.(arr ./ Float64(N)^3)
        end
        lapd_tiles_f64 = nothing; GC.gc()

        rank == 0 && verbose && @info "Phase 1f: Laplacian done"
    end

    # ---- Phase 2: Multi-scale peak finding ----
    tile_masks    = Dict(tid => zeros(Int8, nmesh, nmesh, nmesh) for tid in my_tiles)
    tile_peaks    = Dict(tid => PeakPatch.PeakCandidate[] for tid in my_tiles)
    tile_peak_Rf  = Dict(tid => Float64[] for tid in my_tiles)
    tile_peak_FcRf = Dict(tid => Float32[] for tid in my_tiles)
    tile_peak_d2Rf = Dict(tid => Float32[] for tid in my_tiles)

    for ic in 1:length(filters)
        Rf = filters[ic][3]

        # Distributed smoothing: k-space window + inverse FFT → tile extraction
        smoothed_k = _distributed_apply_window(delta_k_saved, N, boxsize_full, Rf, wsmooth)
        smoothed_real = plan \ smoothed_k
        delta_s_tiles = _extract_my_tiles(smoothed_real, my_tiles, ntile, nsub, nmesh,
                                          Float32, comm)

        lapd_s_tiles = nothing
        if ioutshear >= 1
            lapd_s_k = _distributed_apply_window(delta_k_saved, N, boxsize_full, Rf, 3)
            lapd_s_real = plan \ lapd_s_k
            lapd_s_tiles = _extract_my_tiles(lapd_s_real, my_tiles, ntile, nsub, nmesh,
                                             Float32, comm)
        end

        filter_peak_count = 0
        for tid in my_tiles
            it, jt, kt = tid
            delta_s_tile = delta_s_tiles[tid]
            xbx, ybx, zbx = PeakPatch.tile_center(it, jt, kt, ntile, dcore_box)

            # Per-tile fcrit in lightcone mode
            fcrit_tile = fcrit
            if ievol == 1
                z_tile = PeakPatch.peak_redshift(obs[1], obs[2], obs[3], xbx, ybx, zbx, chi2z)
                fcrit_tile = Float32(PeakPatch.fsc_of_z(z_tile, growth_tables))
            end

            new_peaks = PeakPatch.find_peaks(delta_s_tile, tile_masks[tid],
                                              xbx, ybx, zbx, alatt, nbuff, fcrit_tile, Rf)

            append!(tile_peaks[tid], new_peaks)
            append!(tile_peak_Rf[tid], fill(Rf, length(new_peaks)))

            if !isempty(new_peaks)
                lapd_s_tile = lapd_s_tiles !== nothing ? lapd_s_tiles[tid] : nothing
                for pk in new_peaks
                    i, j, k = _ipp_to_ijk(pk.ipp, nmesh)
                    push!(tile_peak_FcRf[tid], delta_s_tile[i, j, k])
                    push!(tile_peak_d2Rf[tid],
                          lapd_s_tile !== nothing ? lapd_s_tile[i, j, k] : 0.0f0)
                end
            end
            filter_peak_count += length(new_peaks)
        end

        rank == 0 && verbose && @info "  Filter $ic: Rf=$(round(Rf;digits=3)), peaks=$filter_peak_count"
    end

    delta_k_saved = nothing; GC.gc()

    total_peaks = sum(length(tile_peaks[tid]) for tid in my_tiles; init=0)
    rank == 0 && verbose && @info "Phase 2 done: $total_peaks peaks on rank $rank"

    # ---- Phase 3-4: Shell analysis per tile (local, no comm) ----
    halos_basic = PeakPatch.HaloRecord[]
    halos_ext   = PeakPatch.ExtHaloRecord[]

    for tid in my_tiles
        it, jt, kt = tid
        peaks = tile_peaks[tid]
        isempty(peaks) && continue

        delta_tile = delta_tiles[tid]
        psi_x_tile = psi_tiles[1][tid]
        psi_y_tile = psi_tiles[2][tid]
        psi_z_tile = psi_tiles[3][tid]

        psi2_x_tile = psi2_y_tile = psi2_z_tile = nothing
        if psi2_tiles !== nothing
            psi2_x_tile = psi2_tiles[1][tid]
            psi2_y_tile = psi2_tiles[2][tid]
            psi2_z_tile = psi2_tiles[3][tid]
        end

        lapd_tile = lapd_tiles !== nothing ? lapd_tiles[tid] : nothing

        mask_tile = tile_masks[tid]
        fill!(mask_tile, 0)

        pg = PeakPatch.PeakGrid(delta_tile, psi_x_tile, psi_y_tile, psi_z_tile,
                                  psi2_x_tile, psi2_y_tile, psi2_z_tile,
                                  mask_tile, (nmesh, nmesh, nmesh), lapd_tile)

        for idx in 1:length(peaks)
            peak = peaks[idx]
            Rf = tile_peak_Rf[tid][idx]
            ir2min = min(floor(Int, (1.75 * Rf / alatt)^2),
                         floor(Int, (40.0 / alatt - 1)^2))

            # Per-peak ZZon in lightcone mode
            ZZon_pk = ZZon
            if ievol == 1
                z_pk = PeakPatch.peak_redshift(obs[1], obs[2], obs[3], peak.x, peak.y, peak.z, chi2z)
                ZZon_pk = 1.0 + z_pk
                if z_pk > z_max
                    continue
                end
            end

            result = PeakPatch.analyse_peak(pg, peak.ipp, alatt, ir2min, ZZon_pk, Rf, ct, shells;
                                             nbuff=nbuff, growth_tables=growth_tables,
                                             rmax2rs=Float64(sp.rmax2rs))

            if result.RTHL > 0
                # Per-peak growth factor and velocity scaling
                a_pk = 1.0 / ZZon_pk
                _, _, D_pk = PeakPatch.Dlinear_ab(a_pk, growth_tables)

                RTHL_phys = Float32(result.RTHL * alatt)
                Sbar_vel = result.Sbar .* D_pk

                Sbar2_vel = zeros(3)
                if psi2_tiles !== nothing
                    Om_a = Omnr * a_pk^3 / (Omnr * a_pk^3 + cosmo.OL)
                    Sbar2_vel = -result.Sbar2 .* (-3.0/7.0 * Om_a^(-1.0/143) * D_pk^2)
                end

                if ioutshear >= 1
                    sm = result.strain_mat
                    push!(halos_ext, PeakPatch.ExtHaloRecord(
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
                    push!(halos_basic, PeakPatch.HaloRecord(
                        Float32(peak.x), Float32(peak.y), Float32(peak.z),
                        Float32(Sbar_vel[1]), Float32(Sbar_vel[2]), Float32(Sbar_vel[3]),
                        RTHL_phys,
                        Float32(Sbar2_vel[1]), Float32(Sbar2_vel[2]), Float32(Sbar2_vel[3]),
                        Float32(result.Fbarx)
                    ))
                end
            end
        end
    end

    local_halos = ioutshear >= 1 ? halos_ext : halos_basic

    rank == 0 && verbose && @info "Phase 3-4 done: $(length(local_halos)) halos on rank $rank"

    # ---- Phase 5: Gather halos to rank 0 ----
    nfields = ioutshear >= 1 ? fieldcount(PeakPatch.ExtHaloRecord) : fieldcount(PeakPatch.HaloRecord)
    T_halo = ioutshear >= 1 ? PeakPatch.ExtHaloRecord : PeakPatch.HaloRecord

    # Serialize local halos as flat Float32 vector
    local_data = Float32[]
    for h in local_halos
        for fn in fieldnames(typeof(h))
            push!(local_data, getfield(h, fn))
        end
    end

    # Gather counts from all ranks
    local_count = Int32[length(local_data)]
    all_counts = MPI.Allgather(local_count, comm)
    all_counts_vec = vec(all_counts)

    # Gather data to rank 0
    if rank == 0
        total_floats = sum(all_counts_vec)
        recvbuf = Vector{Float32}(undef, total_floats)
        MPI.Gatherv!(local_data, MPI.VBuffer(recvbuf, all_counts_vec), 0, comm)

        # Deserialize halos
        all_halos = T_halo[]
        pos = 1
        while pos + nfields - 1 <= length(recvbuf)
            vals = ntuple(i -> recvbuf[pos + i - 1], nfields)
            push!(all_halos, T_halo(vals...))
            pos += nfields
        end

        verbose && @info "Phase 5: gathered $(length(all_halos)) halos on rank 0"
        return all_halos
    else
        MPI.Gatherv!(local_data, nothing, 0, comm)
        return T_halo[]
    end
end

end # module MPIExt
