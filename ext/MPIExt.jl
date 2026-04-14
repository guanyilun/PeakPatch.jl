module MPIExt

using PeakPatch
using MPI
using PencilFFTs
using PencilFFTs.Transforms: RFFT, RFFT!
using PencilArrays
using PencilArrays: Transpositions
using FFTW

# ============================================================
# MPI large-count workaround
# ============================================================
# MPI.Buffer.count is Cint (Int32), limiting messages to 2^31-1 elements.
# At 6144³ with few ranks, pencil transposes exceed this.  Use MPI derived
# contiguous datatypes so that count = nelements ÷ chunk fits in Int32.

"""Find smallest divisor of `n` such that `n ÷ divisor ≤ typemax(Cint)`."""
function _large_count_chunk(n::Integer)
    max_count = Int(typemax(Cint))
    n ≤ max_count && return 1
    m = cld(n, max_count)
    while n % m != 0
        m += 1
    end
    return m
end

# Override PencilArrays' internal mpi_buffer to handle large messages.
# The original (Transpositions.jl:272) does MPI.Buffer(view) which calls
# Cint(length(view)) and crashes for >2^31 elements.
# Deferred to __init__ to avoid "method overwriting during precompilation" error.
function _mpi_buffer_large(buf::AbstractArray, off, len)
    inds = (off + 1):(off + len)
    v = view(buf, inds)
    if len ≤ typemax(Cint)
        return MPI.Buffer(v)
    end
    chunk = _large_count_chunk(len)
    bigtype = MPI.Types.create_contiguous(Cint(chunk), MPI.Datatype(eltype(buf)))
    MPI.Types.commit!(bigtype)
    return MPI.Buffer(v, Cint(len ÷ chunk), bigtype)
end

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

When `global_tile_list` is provided, only those tiles participate in communication
(for batched extraction). All ranks must pass the same `global_tile_list`.

Returns `Dict{NTuple{3,Int}, Array{T,3}}` mapping tile IDs to nmesh³ arrays.
"""
function _extract_my_tiles(pencil_arr, my_tiles::Vector{NTuple{3,Int}},
                           ntile::Int, nsub::Int, nmesh::Int,
                           ::Type{T}, comm::MPI.Comm;
                           global_tile_list::Union{Nothing, Vector{NTuple{3,Int}}}=nothing) where T
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
    if global_tile_list === nothing
        all_tiles = NTuple{3,Int}[]
        for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
            push!(all_tiles, (it, jt, kt))
        end
    else
        all_tiles = global_tile_list
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
# In-place k-space kernel helpers (write into pre-allocated output)
# ============================================================
# These avoid allocating a ~155 GB temporary PencilArray per call.

"""Apply 1LPT kernel im*ki/k² to delta_k, writing result into output_k."""
function _apply_1lpt_kernel!(output_k, delta_k_pencil, dim::Int, N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = Float64.(FFTW.rfftfreq(N, N * dk))
    ky_arr = Float64.(FFTW.fftfreq(N, N * dk))
    kz_arr = Float64.(FFTW.fftfreq(N, N * dk))
    nyq = N ÷ 2; nk = N ÷ 2 + 1

    pen = PencilArrays.pencil(output_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(output_k)
    gv_in = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0 || gix == nk || giy == nyq + 1 || giz == nyq + 1
            gv_out[gix, giy, giz] = 0
            continue
        end
        ki = dim == 1 ? kx : dim == 2 ? ky : kz
        gv_out[gix, giy, giz] = im * ki / k2 * gv_in[gix, giy, giz]
    end
end

"""Apply 2LPT kernel -im*ki/k² to src2_k, writing result into output_k."""
function _apply_2lpt_kernel!(output_k, src2_k_pencil, dim::Int, N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = Float64.(FFTW.rfftfreq(N, N * dk))
    ky_arr = Float64.(FFTW.fftfreq(N, N * dk))
    kz_arr = Float64.(FFTW.fftfreq(N, N * dk))
    nyq = N ÷ 2; nk = N ÷ 2 + 1

    pen = PencilArrays.pencil(output_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(output_k)
    gv_in = PencilArrays.global_view(src2_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0 || gix == nk || giy == nyq + 1 || giz == nyq + 1
            gv_out[gix, giy, giz] = 0
            continue
        end
        ki = dim == 1 ? kx : dim == 2 ? ky : kz
        gv_out[gix, giy, giz] = -im * ki / k2 * gv_in[gix, giy, giz]
    end
end

"""Apply phi_ij kernel -ki*kj/k² to delta_k, writing result into output_k."""
function _apply_phi_ij_kernel!(output_k, delta_k_pencil, di::Int, dj::Int,
                               N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = Float64.(FFTW.rfftfreq(N, N * dk))
    ky_arr = Float64.(FFTW.fftfreq(N, N * dk))
    kz_arr = Float64.(FFTW.fftfreq(N, N * dk))
    nyq = N ÷ 2; nk = N ÷ 2 + 1

    pen = PencilArrays.pencil(output_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(output_k)
    gv_in = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        if k2 == 0.0 || gix == nk || giy == nyq + 1 || giz == nyq + 1
            gv_out[gix, giy, giz] = 0
            continue
        end
        ki = di == 1 ? kx : di == 2 ? ky : kz
        kj = dj == 1 ? kx : dj == 2 ? ky : kz
        gv_out[gix, giy, giz] = -ki * kj / k2 * gv_in[gix, giy, giz]
    end
end

"""Apply smoothing window to delta_k, writing result into output_k."""
function _apply_window_kernel!(output_k, delta_k_pencil, N::Int, boxsize::Float64,
                               Rf::Float64, wsmooth::Int)
    dk = 2π / boxsize
    kx_arr = FFTW.rfftfreq(N, N * dk)
    ky_arr = FFTW.fftfreq(N, N * dk)
    kz_arr = FFTW.fftfreq(N, N * dk)

    pen = PencilArrays.pencil(output_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(output_k)
    gv_in = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k = sqrt(kx^2 + ky^2 + kz^2)
        fkR = k * Rf
        if fkR == 0.0
            gv_out[gix, giy, giz] = gv_in[gix, giy, giz]
            continue
        end
        if wsmooth == 0
            w = PeakPatch.gaussian_window_fortran(fkR)
        elseif wsmooth == 3
            w = PeakPatch.gaussian_window_fortran(fkR) * k^2
        else
            w = PeakPatch.tophat_window(fkR)
        end
        gv_out[gix, giy, giz] = w * gv_in[gix, giy, giz]
    end
end

"""Apply Laplacian kernel k² to delta_k, writing result into output_k."""
function _apply_laplacian_kernel!(output_k, delta_k_pencil, N::Int, boxsize::Float64)
    dk = 2π / boxsize
    kx_arr = FFTW.rfftfreq(N, N * dk)
    ky_arr = FFTW.fftfreq(N, N * dk)
    kz_arr = FFTW.fftfreq(N, N * dk)

    pen = PencilArrays.pencil(output_k)
    ranges = PencilArrays.range_local(pen)
    gv_out = PencilArrays.global_view(output_k)
    gv_in = PencilArrays.global_view(delta_k_pencil)

    for giz in ranges[3], giy in ranges[2], gix in ranges[1]
        kx = kx_arr[gix]; ky = ky_arr[giy]; kz = kz_arr[giz]
        k2 = kx^2 + ky^2 + kz^2
        gv_out[gix, giy, giz] = k2 * gv_in[gix, giy, giz]
    end
end

# ============================================================
# Shell analysis helper (shared by standard and lowmem paths)
# ============================================================

"""Run shell analysis for peaks in one tile. Appends halos to output vectors."""
function _analyse_tile_shells!(halos_basic, halos_ext,
        tid, peaks, peak_Rfs, peak_FcRfs, peak_d2Rfs, mask_tile,
        delta_tile, psi_x, psi_y, psi_z,
        psi2_x, psi2_y, psi2_z, lapd_tile,
        nmesh, nbuff, alatt, ntile, dcore_box,
        obs, ievol, z_max, chi2z, growth_tables, ZZon,
        Omnr, cosmo, ct, shells, cfg, ioutshear, ilpt)

    isempty(peaks) && return

    it, jt, kt = tid
    fill!(mask_tile, 0)

    pg = PeakPatch.PeakGrid(delta_tile, psi_x, psi_y, psi_z,
                              psi2_x, psi2_y, psi2_z,
                              mask_tile, (nmesh, nmesh, nmesh), lapd_tile)

    for idx in 1:length(peaks)
        peak = peaks[idx]
        Rf = peak_Rfs[idx]
        ir2min = min(floor(Int, (1.75 * Rf / alatt)^2),
                     floor(Int, (40.0 / alatt - 1)^2))

        ZZon_pk = ZZon
        if ievol == 1
            z_pk = PeakPatch.peak_redshift(obs[1], obs[2], obs[3], peak.x, peak.y, peak.z, chi2z)
            ZZon_pk = 1.0 + z_pk
            z_pk > z_max && continue
        end

        result = PeakPatch.analyse_peak(pg, peak.ipp, alatt, ir2min, ZZon_pk, Rf, ct, shells;
                                         nbuff=nbuff, growth_tables=growth_tables,
                                         rmax2rs=cfg.rmax2rs)

        if result.RTHL > 0
            a_pk = 1.0 / ZZon_pk
            _, _, D_pk = PeakPatch.Dlinear_ab(a_pk, growth_tables)

            RTHL_phys = Float32(result.RTHL * alatt)
            Sbar_vel = result.Sbar .* D_pk

            Sbar2_vel = zeros(3)
            if ilpt >= 2
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
                    peak_FcRfs[idx],
                    peak_d2Rfs[idx],
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

# ============================================================
# Main distributed driver
# ============================================================

function PeakPatch.run_multitile_mpi(cfg::PeakPatch.PipelineConfig;
        ntile::Int, seed::Integer=42, verbose::Bool=false,
        comm::MPI.Comm=MPI.COMM_WORLD, lowmem::Bool=false,
        available_gb::Float64=900.0)

    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    # ---- Geometry (all ranks, no comm) ----
    nmesh  = cfg.n
    nbuff  = cfg.nbuff
    nsub   = nmesh - 2 * nbuff
    N      = nsub * ntile + 2 * nbuff
    alatt  = cfg.boxsize / nmesh
    boxsize_full = N * alatt
    dcore_box = nsub * alatt

    # ---- Phase 0: Initialization (all ranks, no comm) ----
    Om_total = cfg.Omx + cfg.OmB
    cosmo = PeakPatch.CosmologyParams(Om_total, cfg.OmB, cfg.Omvac,
                                       cfg.h, 0.965, 0.808)
    growth_tables = PeakPatch.Dlinear_tables(cosmo)

    ct_array, ct_params = PeakPatch.read_homeltab(cfg.tabfile)
    ct = PeakPatch.CollapseTableInterp(ct_array, ct_params)

    filters = PeakPatch.read_filterbank(cfg.filterfile)
    sort!(filters; by=f -> -f[3])

    z_out  = cfg.z_out
    a_out  = 1.0 / (1.0 + z_out)
    ZZon   = 1.0 + z_out
    fcrit  = Float32(PeakPatch.fsc_of_z(z_out, growth_tables))
    _, _, D_out = PeakPatch.Dlinear_ab(a_out, growth_tables)

    Rfclmax = filters[1][3]
    nhunt  = min(nbuff - 1, floor(Int, Rfclmax * 1.75 / alatt))
    shells = PeakPatch.precompute_shells(nhunt)

    Omnr    = cosmo.Om
    vTHvir0 = 100.0 * cosmo.h * sqrt(Omnr)
    ioutshear = cfg.ioutshear
    wsmooth   = cfg.wsmooth
    ilpt      = cfg.ilpt

    cfg.NonGauss != 0 && error("run_multitile_mpi: NonGauss=$(cfg.NonGauss) not yet supported (only 0)")

    # Lightcone mode
    ievol = cfg.ievol
    obs = (cfg.cenx, cfg.ceny, cfg.cenz)
    z_max = cfg.z_max
    chi2z = ievol == 1 ? PeakPatch.build_chi_to_z(cosmo; z_max=z_max + 1.0) : nothing

    my_tiles = _assign_tiles(ntile, nranks, rank)

    # Lightcone: prune tiles beyond maximum_redshift
    chi_max = 0.0
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

    rank == 0 && verbose && @info "Phase 0: N=$N, ntile=$ntile, nranks=$nranks, tiles/rank=$(length(my_tiles))$(ievol == 1 ? ", lightcone mode" : "")$(lowmem ? ", lowmem" : "")"

    # ---- PencilFFT setup (in-place RFFT!, Float32) ----
    # In-place RFFT! saves ~310 GB/node by sharing a single data buffer for
    # real and complex views, eliminating ibuf/obuf intermediate arrays.
    nthreads = Threads.nthreads()
    FFTW.set_num_threads(nthreads)
    rank == 0 && verbose && @info "FFTW threads: $nthreads"
    proc_dims = _factor_procs(nranks)
    _plan_tmp = PencilFFTPlan((N, N, N), RFFT(), proc_dims, comm)
    _noise_tmp = PencilFFTs.allocate_input(_plan_tmp)
    _pen = PencilArrays.pencil(_noise_tmp)
    noise_pencil = PencilArray{Float32}(undef, _pen)
    _plan_tmp = nothing; _noise_tmp = nothing
    plan = PencilFFTPlan(noise_pencil, RFFT!())
    noise_pencil = nothing
    # A is a ManyPencilArrayRFFT! — first(A) = real view, last(A) = complex k-space view.
    # All views share the same memory buffer (~155 GB at 6144³ Float32, 6 ranks).
    A = PencilFFTs.allocate_input(plan)

    # ---- Phase 1a: Distributed noise generation (Threefry counter-based RNG) ----
    pk = PeakPatch.load_pk(cfg.pkfile)

    pen = PencilArrays.pencil(first(A))
    ranges = PencilArrays.range_local(pen)
    PeakPatch.fill_noise_threefry_region!(parent(first(A)), ranges, N, seed)

    rank == 0 && verbose && @info "Phase 1a: distributed noise generated (Threefry, Float32)"

    # ---- Phase 1b: Distributed forward FFT + P(k) convolution ----
    plan * A
    # last(A) now contains delta_k in the final pencil configuration
    _distributed_convolve_pk!(last(A), pk, N, boxsize_full)

    # Save delta_k for reuse (separate allocation; last(A) will be overwritten by IFFTs)
    delta_k_saved = similar(last(A))
    delta_k_saved .= last(A)

    rank == 0 && verbose && @info "Phase 1b: forward FFT + P(k) done"

    # ================================================================
    # Phase 1 field extraction: standard vs lowmem
    # ================================================================
    # All backward FFTs use in-place: write kernel into last(A), plan \ A,
    # read result from first(A).  first(A) is overwritten each time.

    delta_tiles = nothing
    psi_tiles = nothing
    psi2_tiles = nothing
    lapd_tiles = nothing
    src2_k_pencil = nothing

    if !lowmem
        # ---- Standard: extract all tiles during Phase 1 ----

        # Phase 1c: delta tiles (IFFT of delta_k)
        last(A) .= delta_k_saved
        plan \ A
        delta_tiles = _extract_my_tiles(first(A), my_tiles, ntile, nsub, nmesh,
                                        Float32, comm)

        # Phase 1d: 1LPT displacement tiles
        psi_tiles = Dict{Int, Dict{NTuple{3,Int}, Array{Float32,3}}}()
        for dim in 1:3
            _apply_1lpt_kernel!(last(A), delta_k_saved, dim, N, boxsize_full)
            plan \ A
            psi_tiles[dim] = _extract_my_tiles(first(A), my_tiles, ntile, nsub, nmesh,
                                               Float32, comm)
        end

        rank == 0 && verbose && @info "Phase 1d: 1LPT done"

        # Phase 1e: 2LPT displacement tiles (trace identity for src2)
        if ilpt >= 2
            # src2 = (delta² - Σ phi_ii²) / 2 - Σ_{i<j} phi_ij²
            # Uses trace identity: Σ_{i<j} phi_ii*phi_jj = ((Σ phi_ii)² - Σ phi_ii²)/2
            # where Σ phi_ii = -delta (Poisson equation).
            # This avoids holding 2 phi arrays simultaneously (saves ~155 GB).
            src2_pencil = PencilArray{Float32}(undef, PencilArrays.pencil(first(A)))
            parent(src2_pencil) .= 0.0f0

            # delta² / 2  (trace squared = delta²)
            last(A) .= delta_k_saved
            plan \ A
            parent(src2_pencil) .+= parent(first(A)) .^ 2 .* 0.5f0

            # - phi_ii² / 2 for each diagonal
            for d in 1:3
                _apply_phi_ij_kernel!(last(A), delta_k_saved, d, d, N, boxsize_full)
                plan \ A
                parent(src2_pencil) .-= parent(first(A)) .^ 2 .* 0.5f0
            end

            # - phi_ij² for each off-diagonal
            for (di, dj) in ((1,2), (1,3), (2,3))
                _apply_phi_ij_kernel!(last(A), delta_k_saved, di, dj, N, boxsize_full)
                plan \ A
                parent(src2_pencil) .-= parent(first(A)) .^ 2
            end

            # Forward FFT src2 to get src2_k
            parent(first(A)) .= parent(src2_pencil)
            src2_pencil = nothing; GC.gc()
            plan * A
            src2_k_local = similar(last(A))
            src2_k_local .= last(A)

            # Extract psi2 tiles from src2_k
            psi2_tiles = Dict{Int, Dict{NTuple{3,Int}, Array{Float32,3}}}()
            for dim in 1:3
                _apply_2lpt_kernel!(last(A), src2_k_local, dim, N, boxsize_full)
                plan \ A
                psi2_tiles[dim] = _extract_my_tiles(first(A), my_tiles, ntile, nsub, nmesh,
                                                    Float32, comm)
            end
            src2_k_local = nothing; GC.gc()

            rank == 0 && verbose && @info "Phase 1e: 2LPT done (trace identity src2)"
        end

        # Phase 1f: Laplacian tiles
        if ioutshear >= 1
            _apply_laplacian_kernel!(last(A), delta_k_saved, N, boxsize_full)
            plan \ A
            lapd_tiles = _extract_my_tiles(first(A), my_tiles, ntile, nsub, nmesh,
                                           Float32, comm)
            scale = Float32(1.0 / Float64(N)^3)
            for (_, arr) in lapd_tiles
                arr .*= scale
            end

            rank == 0 && verbose && @info "Phase 1f: Laplacian done"
        end

    else
        # ---- Lowmem: compute src2_k only, defer tile extraction to Phase 3-4 ----

        if ilpt >= 2
            src2_pencil = PencilArray{Float32}(undef, PencilArrays.pencil(first(A)))
            parent(src2_pencil) .= 0.0f0

            # Trace identity (same as standard path)
            last(A) .= delta_k_saved
            plan \ A
            parent(src2_pencil) .+= parent(first(A)) .^ 2 .* 0.5f0

            for d in 1:3
                _apply_phi_ij_kernel!(last(A), delta_k_saved, d, d, N, boxsize_full)
                plan \ A
                parent(src2_pencil) .-= parent(first(A)) .^ 2 .* 0.5f0
            end

            for (di, dj) in ((1,2), (1,3), (2,3))
                _apply_phi_ij_kernel!(last(A), delta_k_saved, di, dj, N, boxsize_full)
                plan \ A
                parent(src2_pencil) .-= parent(first(A)) .^ 2
            end

            parent(first(A)) .= parent(src2_pencil)
            src2_pencil = nothing; GC.gc()
            plan * A
            src2_k_pencil = similar(last(A))
            src2_k_pencil .= last(A)

            rank == 0 && verbose && @info "Phase 1e: src2_k computed (lowmem, trace identity)"
        end
    end

    # ================================================================
    # Phase 2: Multi-scale peak finding (shared by both modes)
    # ================================================================

    tile_masks    = Dict(tid => zeros(Int8, nmesh, nmesh, nmesh) for tid in my_tiles)
    tile_peaks    = Dict(tid => PeakPatch.PeakCandidate[] for tid in my_tiles)
    tile_peak_Rf  = Dict(tid => Float64[] for tid in my_tiles)
    tile_peak_FcRf = Dict(tid => Float32[] for tid in my_tiles)
    tile_peak_d2Rf = Dict(tid => Float32[] for tid in my_tiles)

    for ic in 1:length(filters)
        Rf = filters[ic][3]

        _apply_window_kernel!(last(A), delta_k_saved, N, boxsize_full, Rf, wsmooth)
        plan \ A
        delta_s_tiles = _extract_my_tiles(first(A), my_tiles, ntile, nsub, nmesh,
                                          Float32, comm)

        lapd_s_tiles = nothing
        if ioutshear >= 1
            _apply_window_kernel!(last(A), delta_k_saved, N, boxsize_full, Rf, 3)
            plan \ A
            lapd_s_tiles = _extract_my_tiles(first(A), my_tiles, ntile, nsub, nmesh,
                                             Float32, comm)
        end

        filter_peak_count = 0
        for tid in my_tiles
            it, jt, kt = tid
            delta_s_tile = delta_s_tiles[tid]
            xbx, ybx, zbx = PeakPatch.tile_center(it, jt, kt, ntile, dcore_box)

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

    if !lowmem
        delta_k_saved = nothing; GC.gc()
    end

    total_peaks = sum(length(tile_peaks[tid]) for tid in my_tiles; init=0)
    rank == 0 && verbose && @info "Phase 2 done: $total_peaks peaks on rank $rank"

    # ================================================================
    # Phase 3-4: Shell analysis
    # ================================================================

    halos_basic = PeakPatch.HaloRecord[]
    halos_ext   = PeakPatch.ExtHaloRecord[]

    if !lowmem
        # ---- Standard: use pre-extracted tiles ----
        for tid in my_tiles
            _analyse_tile_shells!(halos_basic, halos_ext,
                tid, tile_peaks[tid], tile_peak_Rf[tid],
                tile_peak_FcRf[tid], tile_peak_d2Rf[tid], tile_masks[tid],
                delta_tiles[tid],
                psi_tiles[1][tid], psi_tiles[2][tid], psi_tiles[3][tid],
                psi2_tiles !== nothing ? psi2_tiles[1][tid] : nothing,
                psi2_tiles !== nothing ? psi2_tiles[2][tid] : nothing,
                psi2_tiles !== nothing ? psi2_tiles[3][tid] : nothing,
                lapd_tiles !== nothing ? lapd_tiles[tid] : nothing,
                nmesh, nbuff, alatt, ntile, dcore_box,
                obs, ievol, z_max, chi2z, growth_tables, ZZon,
                Omnr, cosmo, ct, shells, cfg, ioutshear, ilpt)
        end
    else
        # ---- Lowmem: batched tile extraction + shell analysis ----

        # Build global active tile list (deterministic, all ranks agree)
        all_active_tiles = NTuple{3,Int}[]
        half_lc = dcore_box / 2.0
        for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
            tid = (it, jt, kt)
            if ievol == 1
                xbx, ybx, zbx = PeakPatch.tile_center(it, jt, kt, ntile, dcore_box)
                dx = max(abs(xbx - obs[1]) - half_lc, 0.0)
                dy = max(abs(ybx - obs[2]) - half_lc, 0.0)
                dz = max(abs(zbx - obs[3]) - half_lc, 0.0)
                sqrt(dx^2 + dy^2 + dz^2) <= chi_max || continue
            end
            push!(all_active_tiles, tid)
        end

        # Compute batch size from memory budget
        # Peak during IFFT: delta_k + src2_k + field_k + field_real + (n_fields-1) batch tiles
        # = (n_stored + 2) × dk_bytes + (n_fields-1) × B × tile_bytes
        local_dk_bytes = sizeof(parent(delta_k_saved))
        tile_field_bytes = nmesh^3 * sizeof(Float32)
        n_stored = src2_k_pencil !== nothing ? 2 : 1  # delta_k + maybe src2_k
        n_tile_fields = 1 + 3 + (ilpt >= 2 ? 3 : 0) + (ioutshear >= 1 ? 1 : 0)
        overhead_bytes = (n_stored + 2) * local_dk_bytes
        avail_for_tiles = available_gb * 1e9 - overhead_bytes
        tiles_per_batch = max(1, floor(Int, avail_for_tiles / ((n_tile_fields - 1) * tile_field_bytes)))

        my_tile_set = Set(my_tiles)
        n_batches = cld(length(all_active_tiles), tiles_per_batch)

        rank == 0 && verbose && @info "Phase 3-4: lowmem batched, $(length(all_active_tiles)) active tiles, batch_size=$tiles_per_batch, $n_batches batches"

        for bi in 1:n_batches
            batch_start = (bi - 1) * tiles_per_batch + 1
            batch_end = min(bi * tiles_per_batch, length(all_active_tiles))
            batch_global = all_active_tiles[batch_start:batch_end]
            my_batch = [t for t in batch_global if t in my_tile_set]

            rank == 0 && verbose && @info "  Batch $bi/$n_batches: $(length(batch_global)) tiles ($(length(my_batch)) local)"

            # Extract delta tiles
            last(A) .= delta_k_saved
            plan \ A
            batch_delta = _extract_my_tiles(first(A), my_batch, ntile, nsub, nmesh,
                                            Float32, comm; global_tile_list=batch_global)

            # Extract 1LPT displacement tiles
            batch_psi = Dict{Int, Dict{NTuple{3,Int}, Array{Float32,3}}}()
            for dim in 1:3
                _apply_1lpt_kernel!(last(A), delta_k_saved, dim, N, boxsize_full)
                plan \ A
                batch_psi[dim] = _extract_my_tiles(first(A), my_batch, ntile, nsub, nmesh,
                                                   Float32, comm; global_tile_list=batch_global)
            end

            # Extract 2LPT displacement tiles
            batch_psi2 = Dict{Int, Dict{NTuple{3,Int}, Array{Float32,3}}}()
            if ilpt >= 2 && src2_k_pencil !== nothing
                for dim in 1:3
                    _apply_2lpt_kernel!(last(A), src2_k_pencil, dim, N, boxsize_full)
                    plan \ A
                    batch_psi2[dim] = _extract_my_tiles(first(A), my_batch, ntile, nsub, nmesh,
                                                        Float32, comm; global_tile_list=batch_global)
                end
            end

            # Extract Laplacian tiles
            batch_lapd = nothing
            if ioutshear >= 1
                _apply_laplacian_kernel!(last(A), delta_k_saved, N, boxsize_full)
                plan \ A
                batch_lapd = _extract_my_tiles(first(A), my_batch, ntile, nsub, nmesh,
                                               Float32, comm; global_tile_list=batch_global)
                scale = Float32(1.0 / Float64(N)^3)
                for (_, arr) in batch_lapd
                    arr .*= scale
                end
            end

            # Shell analysis for this batch
            for tid in my_batch
                _analyse_tile_shells!(halos_basic, halos_ext,
                    tid, tile_peaks[tid], tile_peak_Rf[tid],
                    tile_peak_FcRf[tid], tile_peak_d2Rf[tid], tile_masks[tid],
                    batch_delta[tid],
                    batch_psi[1][tid], batch_psi[2][tid], batch_psi[3][tid],
                    haskey(batch_psi2, 1) ? batch_psi2[1][tid] : nothing,
                    haskey(batch_psi2, 2) ? batch_psi2[2][tid] : nothing,
                    haskey(batch_psi2, 3) ? batch_psi2[3][tid] : nothing,
                    batch_lapd !== nothing ? batch_lapd[tid] : nothing,
                    nmesh, nbuff, alatt, ntile, dcore_box,
                    obs, ievol, z_max, chi2z, growth_tables, ZZon,
                    Omnr, cosmo, ct, shells, cfg, ioutshear, ilpt)
            end

            # Free batch data
            batch_delta = nothing; batch_psi = nothing
            batch_psi2 = nothing; batch_lapd = nothing; GC.gc()
        end

        delta_k_saved = nothing; src2_k_pencil = nothing; GC.gc()
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

# Install mpi_buffer override at module load time (not precompile time)
function __init__()
    @eval Transpositions.mpi_buffer(buf::AbstractArray, off, len) =
        _mpi_buffer_large(buf, off, len)
end

end # module MPIExt
