module CUDAExt

# GPU port of the PeakPatch shell analysis, starting with the block-per-peak
# delta-only gather primitive (Phase 2, step 5 of docs/gpu_implementation_plan.org).
#
# This extension is loaded automatically when the user does `using CUDA` in a
# session that also has PeakPatch available. The CPU fallback (analyse_peak_gpu
# in src/ShellAnalysisGPU.jl) remains the reference implementation.

using PeakPatch
using PeakPatch: ShellTables
using CUDA

# Int32 constants used in kernels (avoids Int64 promotion surprises).
const _I0 = Int32(0)
const _I1 = Int32(1)
const _I5 = Int32(5)
const _I31 = Int32(31)

# ============================================================
# GPU-resident shell offset tables
# ============================================================

struct ShellTablesGPU
    offsets_di::CuVector{Int32}
    offsets_dj::CuVector{Int32}
    offsets_dk::CuVector{Int32}
    shell_start::CuVector{Int32}
    shell_count::CuVector{Int32}
    shell_r2::CuVector{Int32}
    nshells::Int
end

"""
    ShellTablesGPU(stab::ShellTables) -> ShellTablesGPU

Upload a CPU `ShellTables` to GPU memory. Allocates ~7 Int32 arrays totalling
O(rmax^3) elements (~6 MB for rmax=50).
"""
function ShellTablesGPU(stab::ShellTables)
    return ShellTablesGPU(
        CuArray{Int32}(stab.offsets_di),
        CuArray{Int32}(stab.offsets_dj),
        CuArray{Int32}(stab.offsets_dk),
        CuArray{Int32}(stab.shell_start),
        CuArray{Int32}(stab.shell_count),
        CuArray{Int32}(stab.shell_r2),
        stab.nshells,
    )
end

# ============================================================
# SPH kernel lookup table (GPU port of atab4)
# ============================================================
#
# The CPU code caches a 1101-element Float64 lookup table in `_AKK_TAB`.
# For GPU we upload a Float32 version once per extension load (module
# __init__) and expose as `_AKK_GPU`. Callers reference the global.

const _AKK_CPU_F32 = Float32[x for x in PeakPatch.atab4()]
const _AKK_GPU = Ref{CuArray{Float32,1}}()

function __init__()
    _AKK_GPU[] = CuArray(_AKK_CPU_F32)
    return
end

# ============================================================
# Cached cuFFT plans (one per (n, device) pair)
# ============================================================
#
# cuFFT plan creation is non-trivial (workspace alloc + twiddle factor build).
# Plans are size-only and can be reused across any contiguous CuArray of the
# same shape, so we memoize them keyed on (cube size, device id). Each
# multi-GPU worker task binds itself to a CUDA device before invoking these
# helpers (see `set_cuda_device!`), so the device id captured here is the
# task-local one.
#
# Used by isolated_convolve_gpu (n2 = 2*nmesh), compute_2lpt_gpu /
# compute_laplacian_gpu / peak_find_tile_gpu (all n = nmesh).

const _CUFFT_PLAN_LOCK = ReentrantLock()
const _CUFFT_FWD_CACHE = Dict{Tuple{Int,Int},Any}()  # (n, devid) -> plan_rfft
const _CUFFT_INV_CACHE = Dict{Tuple{Int,Int},Any}()  # (n, devid) -> plan_irfft

@inline function _cufft_fwd_plan(arr::CuArray{Float32,3})
    n = size(arr, 1)
    devid = CUDA.deviceid(CUDA.device())
    key = (n, devid)
    lock(_CUFFT_PLAN_LOCK) do
        get!(_CUFFT_FWD_CACHE, key) do
            CUDA.CUFFT.plan_rfft(arr)
        end
    end
end

@inline function _cufft_inv_plan(arr::CuArray{ComplexF32,3}, n::Integer)
    n_int = Int(n)
    devid = CUDA.deviceid(CUDA.device())
    key = (n_int, devid)
    lock(_CUFFT_PLAN_LOCK) do
        get!(_CUFFT_INV_CACHE, key) do
            CUDA.CUFFT.plan_irfft(arr, n_int)
        end
    end
end

@inline _cufft_rfft(arr::CuArray{Float32,3})            = _cufft_fwd_plan(arr) * arr
@inline _cufft_irfft(arr::CuArray{ComplexF32,3}, n::Integer) = _cufft_inv_plan(arr, n) * arr

# ============================================================
# hRinteg (GPU port) — piecewise polynomial used by kernel_strain
# ============================================================

@inline function _hRinteg_f32(u0::Float32, u12::Float32)
    x0 = abs(u0)
    if x0 >= 2.0f0
        return 0.0f0
    end

    one7 = 1.0f0 / 7.0f0
    result = 1.6f0 * u12 - 6.4f0 * one7

    x02 = x0 * x0
    x04 = x02 * x02

    if x0 >= 1.0f0
        h1 = x04 * (2.0f0 - x0 * (2.4f0 - x0 * (1.0f0 - x0 * one7)))
        h2 = x02 * (4.0f0 - x0 * (4.0f0 - x0 * (1.5f0 - x0 * 0.2f0)))
        result = result - h2 * u12 + h1
    else
        result = result - 0.2f0 * (u12 - one7)
        if x0 > 0.0f0
            h1 = x04 * (1.0f0 - x02 * (1.0f0 - 3.0f0 * x0 * one7))
            h2 = x02 * (2.0f0 - x02 * (1.5f0 - 0.6f0 * x0))
            result = result - h2 * u12 + h1
        end
    end
    return result
end

# ============================================================
# kernel_strain (GPU port)
# ============================================================
#
# Device-side radial integration. Operates on per-shell profile arrays
# (rad, Gshell, Gshellf, SRshell) for a single peak — same layout as the
# gather kernel outputs: (3, nshells), (3, 3, nshells) etc., indexed in
# the *shell* dimension from `mlow` to `mupp`, centered on `mp`.
#
# Returns (Ebar_11, Ebar_12, Ebar_13, Ebar_22, Ebar_23, Ebar_33,
#          grad_x, grad_y, grad_z, gradf_x, gradf_y, gradf_z)
# as a tuple of Float32 — 12 scalars. Ebar is symmetric so we return
# only the 6 upper-triangle entries.
#
# Two branches mirror CPU:
#   mp > 1  → hRinteg weighting, scale = wRnor / rad[mp]
#   mp == 1 → akk SPH-kernel table lookup, scale = aRnor

@inline function _kernel_strain_f32(rad, Gshell, Gshellf, SRshell, akk_tab,
                                    mp::Int32, mlow::Int32, mupp::Int32,
                                    wRnor::Float32, aRnor::Float32,
                                    hlatt_1::Float32, hlatt_2::Float32)
    E11 = 0.0f0; E12 = 0.0f0; E13 = 0.0f0
    E22 = 0.0f0; E23 = 0.0f0; E33 = 0.0f0
    gx = 0.0f0; gy = 0.0f0; gz = 0.0f0
    gfx = 0.0f0; gfy = 0.0f0; gfz = 0.0f0

    if mp > _I1
        @inbounds rmp = rad[mp]
        rmp2 = rmp * rmp
        m1 = mlow
        while m1 <= mupp
            @inbounds rm1 = rad[m1]
            u0 = hlatt_1 * (rmp - rm1)
            u02 = u0 * u0
            if u02 < 4.0f0
                u12 = hlatt_2 * (rmp2 + rm1 * rm1)
                denom = (u12 - u02)
                wt = _hRinteg_f32(u0, u12) / (denom * denom)
                @inbounds Gx = Gshell[_I1, m1]
                @inbounds Gy = Gshell[Int32(2), m1]
                @inbounds Gz = Gshell[Int32(3), m1]
                @inbounds Gfx_ = Gshellf[_I1, m1]
                @inbounds Gfy_ = Gshellf[Int32(2), m1]
                @inbounds Gfz_ = Gshellf[Int32(3), m1]
                gx  += wt * Gx;   gy  += wt * Gy;   gz  += wt * Gz
                gfx += wt * Gfx_; gfy += wt * Gfy_; gfz += wt * Gfz_
                # Symmetrised SR: 0.5 * wt * (SR[L,K] + SR[K,L])
                @inbounds SR11 = SRshell[_I1, _I1, m1]
                @inbounds SR22 = SRshell[Int32(2), Int32(2), m1]
                @inbounds SR33 = SRshell[Int32(3), Int32(3), m1]
                @inbounds SR12 = SRshell[_I1, Int32(2), m1]
                @inbounds SR21 = SRshell[Int32(2), _I1, m1]
                @inbounds SR13 = SRshell[_I1, Int32(3), m1]
                @inbounds SR31 = SRshell[Int32(3), _I1, m1]
                @inbounds SR23 = SRshell[Int32(2), Int32(3), m1]
                @inbounds SR32 = SRshell[Int32(3), Int32(2), m1]
                E11 -= wt * SR11   # 0.5*wt*(SR11+SR11) = wt*SR11
                E22 -= wt * SR22
                E33 -= wt * SR33
                E12 -= 0.5f0 * wt * (SR12 + SR21)
                E13 -= 0.5f0 * wt * (SR13 + SR31)
                E23 -= 0.5f0 * wt * (SR23 + SR32)
            end
            m1 += _I1
        end
        scale = wRnor / rmp
        E11 *= scale; E12 *= scale; E13 *= scale
        E22 *= scale; E23 *= scale; E33 *= scale
        gx *= scale; gy *= scale; gz *= scale
        gfx *= scale; gfy *= scale; gfz *= scale
    else
        # mp == 1: akk SPH-kernel table
        m1 = mlow
        while m1 <= mupp
            @inbounds rm1 = rad[m1]
            u0 = hlatt_1 * rm1
            con = 100.0f0 * u0 * u0
            icon = Int32(floor(con))
            diff = con - Float32(icon)
            i2c = icon + _I1  # 1-based index
            # Guard against overflow (akk has length 1101)
            if i2c >= Int32(1) && i2c < Int32(1101)
                @inbounds a0 = akk_tab[i2c]
                @inbounds a1 = akk_tab[i2c + _I1]
                aww = a0 + diff * (a1 - a0)
                wt = 0.5f0 * aww / rm1
                @inbounds Gx = Gshell[_I1, m1]
                @inbounds Gy = Gshell[Int32(2), m1]
                @inbounds Gz = Gshell[Int32(3), m1]
                @inbounds Gfx_ = Gshellf[_I1, m1]
                @inbounds Gfy_ = Gshellf[Int32(2), m1]
                @inbounds Gfz_ = Gshellf[Int32(3), m1]
                gx  += wt * Gx;   gy  += wt * Gy;   gz  += wt * Gz
                gfx += wt * Gfx_; gfy += wt * Gfy_; gfz += wt * Gfz_
                @inbounds SR11 = SRshell[_I1, _I1, m1]
                @inbounds SR22 = SRshell[Int32(2), Int32(2), m1]
                @inbounds SR33 = SRshell[Int32(3), Int32(3), m1]
                @inbounds SR12 = SRshell[_I1, Int32(2), m1]
                @inbounds SR21 = SRshell[Int32(2), _I1, m1]
                @inbounds SR13 = SRshell[_I1, Int32(3), m1]
                @inbounds SR31 = SRshell[Int32(3), _I1, m1]
                @inbounds SR23 = SRshell[Int32(2), Int32(3), m1]
                @inbounds SR32 = SRshell[Int32(3), Int32(2), m1]
                E11 -= wt * SR11
                E22 -= wt * SR22
                E33 -= wt * SR33
                E12 -= 0.5f0 * wt * (SR12 + SR21)
                E13 -= 0.5f0 * wt * (SR13 + SR31)
                E23 -= 0.5f0 * wt * (SR23 + SR32)
            end
            m1 += _I1
        end
        E11 *= aRnor; E12 *= aRnor; E13 *= aRnor
        E22 *= aRnor; E23 *= aRnor; E33 *= aRnor
        gx *= aRnor; gy *= aRnor; gz *= aRnor
        gfx *= aRnor; gfy *= aRnor; gfz *= aRnor
    end
    return (E11, E12, E13, E22, E23, E33,
            gx, gy, gz, gfx, gfy, gfz)
end

# ============================================================
# Trilinear interpolation of a 3D Float32 table (GPU port of CollapseTable)
# ============================================================
#
# CPU path uses `Interpolations.jl` with `BSpline(Linear())` and
# `extrapolate(sitp, out_val)`. On GPU we do manual trilinear lookup
# with explicit in-range clamp; out-of-range queries return `out_val`.
# Grid axes are evenly-spaced (`range(X1, X2, Nx)`), so we use
# multiplicative step inverses rather than per-axis searching.

@inline function _interp3_trilinear_f32(table, x::Float32, y::Float32, z::Float32,
                                         x1::Float32, y1::Float32, z1::Float32,
                                         dx_inv::Float32, dy_inv::Float32, dz_inv::Float32,
                                         nx::Int32, ny::Int32, nz::Int32,
                                         out_val::Float32)
    fx = (x - x1) * dx_inv
    fy = (y - y1) * dy_inv
    fz = (z - z1) * dz_inv

    # Out-of-range → return fallback (matches extrapolate(sitp, out_val))
    if (fx < 0.0f0) | (fy < 0.0f0) | (fz < 0.0f0) |
       (fx > Float32(nx - _I1)) | (fy > Float32(ny - _I1)) | (fz > Float32(nz - _I1))
        return out_val
    end

    # 1-based indices of lower corner; clamp so ix+1 ≤ nx
    ix = Int32(floor(fx)) + _I1
    iy = Int32(floor(fy)) + _I1
    iz = Int32(floor(fz)) + _I1
    ix = ix > (nx - _I1) ? (nx - _I1) : ix
    iy = iy > (ny - _I1) ? (ny - _I1) : iy
    iz = iz > (nz - _I1) ? (nz - _I1) : iz

    tx = fx - Float32(ix - _I1)
    ty = fy - Float32(iy - _I1)
    tz = fz - Float32(iz - _I1)
    ox = 1.0f0 - tx; oy = 1.0f0 - ty; oz = 1.0f0 - tz

    @inbounds v000 = table[ix,       iy,       iz]
    @inbounds v100 = table[ix + _I1, iy,       iz]
    @inbounds v010 = table[ix,       iy + _I1, iz]
    @inbounds v110 = table[ix + _I1, iy + _I1, iz]
    @inbounds v001 = table[ix,       iy,       iz + _I1]
    @inbounds v101 = table[ix + _I1, iy,       iz + _I1]
    @inbounds v011 = table[ix,       iy + _I1, iz + _I1]
    @inbounds v111 = table[ix + _I1, iy + _I1, iz + _I1]

    c00 = v000 * ox + v100 * tx
    c10 = v010 * ox + v110 * tx
    c01 = v001 * ox + v101 * tx
    c11 = v011 * ox + v111 * tx
    c0 = c00 * oy + c10 * ty
    c1 = c01 * oy + c11 * ty
    return c0 * oz + c1 * tz
end

# ============================================================
# 3×3 symmetric eigenvalues via Cardano (trigonometric form)
# ============================================================
#
# GPU port of `RadialShell.get_evals`. The CPU path uses LAPACK via
# `eigvals(Symmetric(mat))`; on-device we use the classical analytic
# Cardano formula for the depressed cubic, which is stable for 3×3
# real symmetric matrices and runs in constant ops (no iteration).
#
# Reference: Smith (1961) "Eigenvalues of a Symmetric 3×3 Matrix",
# Comm. ACM 4(4). Also https://en.wikipedia.org/wiki/Eigenvalue_algorithm
# #3x3_matrices
#
# Returns a NTuple{3,Float32} of eigenvalues sorted ascending. The
# caller gets `Lam[1] <= Lam[2] <= Lam[3]` just like the CPU path.
# Input is treated as symmetric via (mat + mat') / 2 style averaging on
# the fly; caller must pass only the 6 unique entries.

@inline function _eig3_symmetric_f32(a11::Float32, a22::Float32, a33::Float32,
                                     a12::Float32, a13::Float32, a23::Float32)
    # Trace/3 and the trace-free invariants
    q = (a11 + a22 + a33) * (1.0f0 / 3.0f0)
    # Subtract isotropic component
    b11 = a11 - q; b22 = a22 - q; b33 = a33 - q
    # p² = (1/6) * sum of squares of a unit-shifted matrix
    p2 = (b11*b11 + b22*b22 + b33*b33 + 2.0f0*(a12*a12 + a13*a13 + a23*a23)) *
         (1.0f0 / 6.0f0)

    if p2 <= 1.0f-30  # matrix is (nearly) a multiple of I
        return (q, q, q)
    end

    p = sqrt(p2)
    inv_p = 1.0f0 / p
    # B = (A - qI) / p has unit-norm structure; det(B)/2 = r ∈ [-1, 1]
    c11 = b11 * inv_p; c22 = b22 * inv_p; c33 = b33 * inv_p
    c12 = a12 * inv_p; c13 = a13 * inv_p; c23 = a23 * inv_p

    # det(B) / 2 using cofactor expansion
    det_b = c11*(c22*c33 - c23*c23) -
            c12*(c12*c33 - c23*c13) +
            c13*(c12*c23 - c22*c13)
    r = det_b * 0.5f0
    # Clamp for safety against Float32 round-off crossing ±1
    r = r < -1.0f0 ? -1.0f0 : (r > 1.0f0 ? 1.0f0 : r)

    phi = acos(r) * (1.0f0 / 3.0f0)

    # Eigenvalues, ordered largest / smallest via cos offsets
    eig_max = q + 2.0f0 * p * cos(phi)
    eig_min = q + 2.0f0 * p * cos(phi + (2.0f0 * Float32(pi) / 3.0f0))
    eig_mid = 3.0f0 * q - eig_min - eig_max  # from trace: sum = 3q

    # Ascending sort
    return (eig_min, eig_mid, eig_max)
end

# ============================================================
# Block reduction primitives
# ============================================================

# Warp-level sum reduction via shuffle-down. Assumes warp size 32 (NVIDIA).
@inline function warp_reduce_sum_f32(val::Float32)
    val += CUDA.shfl_down_sync(0xffffffff, val, 16)
    val += CUDA.shfl_down_sync(0xffffffff, val, 8)
    val += CUDA.shfl_down_sync(0xffffffff, val, 4)
    val += CUDA.shfl_down_sync(0xffffffff, val, 2)
    val += CUDA.shfl_down_sync(0xffffffff, val, 1)
    return val
end

@inline function warp_reduce_sum_i32(val::Int32)
    val += CUDA.shfl_down_sync(0xffffffff, val, 16)
    val += CUDA.shfl_down_sync(0xffffffff, val, 8)
    val += CUDA.shfl_down_sync(0xffffffff, val, 4)
    val += CUDA.shfl_down_sync(0xffffffff, val, 2)
    val += CUDA.shfl_down_sync(0xffffffff, val, 1)
    return val
end

# Block-level sum reduction: warp-shuffle within each warp, then gather
# partials through shared memory and let warp 0 reduce again.
# Returns the reduced value on thread 1 (valid only there).
@inline function block_reduce_sum_f32(val::Float32, sdata)
    tid = threadIdx().x
    lane = (tid - _I1) & _I31
    warp = (tid - _I1) >> _I5
    v = warp_reduce_sum_f32(val)
    if lane == _I0
        @inbounds sdata[warp + _I1] = v
    end
    sync_threads()
    nwarps = (blockDim().x + _I31) >> _I5
    v2 = (tid - _I1) < nwarps ? (@inbounds sdata[tid]) : 0.0f0
    if warp == _I0
        v2 = warp_reduce_sum_f32(v2)
    end
    return v2
end

@inline function block_reduce_sum_i32(val::Int32, sdata)
    tid = threadIdx().x
    lane = (tid - _I1) & _I31
    warp = (tid - _I1) >> _I5
    v = warp_reduce_sum_i32(val)
    if lane == _I0
        @inbounds sdata[warp + _I1] = v
    end
    sync_threads()
    nwarps = (blockDim().x + _I31) >> _I5
    v2 = (tid - _I1) < nwarps ? (@inbounds sdata[tid]) : _I0
    if warp == _I0
        v2 = warp_reduce_sum_i32(v2)
    end
    return v2
end

# ============================================================
# Block-per-peak shell gather kernel (delta-only)
# ============================================================

function _shell_fbar_kernel!(Fshell_out, nshell_out,
                             delta, peaks_i, peaks_j, peaks_k,
                             off_di, off_dj, off_dk,
                             shell_start, shell_count,
                             nshells::Int32, n1::Int32, n2::Int32, n3::Int32)
    peak_id = blockIdx().x
    tid = threadIdx().x
    bdim = blockDim().x

    # Shared-memory scratch pad: 32 Float32 slots followed by 32 Int32 slots.
    # Accommodates blocks up to 1024 threads (32 warps).
    sdata_f = CuDynamicSharedArray(Float32, 32)
    sdata_i = CuDynamicSharedArray(Int32, 32, 32 * sizeof(Float32))

    @inbounds ci = peaks_i[peak_id]
    @inbounds cj = peaks_j[peak_id]
    @inbounds ck = peaks_k[peak_id]

    s = _I1
    while s <= nshells
        @inbounds s0 = shell_start[s]
        @inbounds nc = shell_count[s]

        local_sum = 0.0f0
        local_nsh = _I0

        # Cooperative strided gather over cells in this shell
        c = tid
        while c <= nc
            idx = s0 + c - _I1
            @inbounds di = off_di[idx]
            @inbounds dj = off_dj[idx]
            @inbounds dk = off_dk[idx]

            iv1 = ci + di
            iv2 = cj + dj
            iv3 = ck + dk

            if (_I1 <= iv1) & (iv1 <= n1) & (_I1 <= iv2) & (iv2 <= n2) & (_I1 <= iv3) & (iv3 <= n3)
                @inbounds local_sum += delta[iv1, iv2, iv3]
                local_nsh += _I1
            end
            c += bdim
        end

        total_sum = block_reduce_sum_f32(local_sum, sdata_f)
        sync_threads()
        total_nsh = block_reduce_sum_i32(local_nsh, sdata_i)

        if tid == _I1
            @inbounds Fshell_out[s, peak_id] = total_sum
            @inbounds nshell_out[s, peak_id] = total_nsh
        end
        sync_threads()  # ensure sdata scratch is free for next shell
        s += _I1
    end
    return
end

# ============================================================
# Host entry point
# ============================================================

function PeakPatch.shell_fbar_gather_gpu(delta_h::AbstractArray{<:Real,3},
                                         peaks_i::AbstractVector{<:Integer},
                                         peaks_j::AbstractVector{<:Integer},
                                         peaks_k::AbstractVector{<:Integer},
                                         stab::ShellTables;
                                         threads::Int=128)
    @assert length(peaks_i) == length(peaks_j) == length(peaks_k)
    npeaks = length(peaks_i)
    nshells = stab.nshells

    # Upload inputs
    delta_d = CuArray{Float32}(delta_h)
    pi_d = CuArray{Int32}(peaks_i)
    pj_d = CuArray{Int32}(peaks_j)
    pk_d = CuArray{Int32}(peaks_k)
    stab_d = ShellTablesGPU(stab)

    Fshell_d = CUDA.zeros(Float32, nshells, npeaks)
    nshell_d = CUDA.zeros(Int32, nshells, npeaks)

    n1, n2, n3 = size(delta_d)
    shmem = 32 * sizeof(Float32) + 32 * sizeof(Int32)

    @cuda threads=threads blocks=npeaks shmem=shmem _shell_fbar_kernel!(
        Fshell_d, nshell_d,
        delta_d, pi_d, pj_d, pk_d,
        stab_d.offsets_di, stab_d.offsets_dj, stab_d.offsets_dk,
        stab_d.shell_start, stab_d.shell_count,
        Int32(nshells), Int32(n1), Int32(n2), Int32(n3),
    )

    return Array(Fshell_d), Array(nshell_d)
end

# ============================================================
# Block-per-peak shell gather kernel (delta + ψ×3 displacement fields)
# ============================================================
#
# Extends the delta-only kernel to also accumulate per-shell:
#   Sshell[L]  = Σ eta_L[cell]              (L ∈ {x,y,z})
#   Gshell[L]  = Σ delta[cell] * d_L        (un-normalised; divide by rad[s] at caller)
#
# Per-thread local state: 8 Float32 accumulators + 1 Int32 counter. Each
# cell gather does 4 Float32 reads (delta, etax, etay, etaz) and 3 Int32
# offset reads. Eight block reductions per shell (7 Float32 + 1 Int32).

function _shell_gather_psi_kernel!(Fshell_out, nshell_out, Sshell_out, Gshell_out,
                                   delta, etax, etay, etaz,
                                   peaks_i, peaks_j, peaks_k,
                                   off_di, off_dj, off_dk,
                                   shell_start, shell_count,
                                   nshells::Int32, n1::Int32, n2::Int32, n3::Int32)
    peak_id = blockIdx().x
    tid = threadIdx().x
    bdim = blockDim().x

    # Shared scratchpad: 32 Float32 for f32 reduce, 32 Int32 for i32 reduce.
    sdata_f = CuDynamicSharedArray(Float32, 32)
    sdata_i = CuDynamicSharedArray(Int32, 32, 32 * sizeof(Float32))

    @inbounds ci = peaks_i[peak_id]
    @inbounds cj = peaks_j[peak_id]
    @inbounds ck = peaks_k[peak_id]

    s = _I1
    while s <= nshells
        @inbounds s0 = shell_start[s]
        @inbounds nc = shell_count[s]

        local_F = 0.0f0
        local_Sx = 0.0f0; local_Sy = 0.0f0; local_Sz = 0.0f0
        local_Gx = 0.0f0; local_Gy = 0.0f0; local_Gz = 0.0f0
        local_n = _I0

        c = tid
        while c <= nc
            idx = s0 + c - _I1
            @inbounds di = off_di[idx]
            @inbounds dj = off_dj[idx]
            @inbounds dk = off_dk[idx]

            iv1 = ci + di
            iv2 = cj + dj
            iv3 = ck + dk

            if (_I1 <= iv1) & (iv1 <= n1) & (_I1 <= iv2) & (iv2 <= n2) & (_I1 <= iv3) & (iv3 <= n3)
                @inbounds d_val = delta[iv1, iv2, iv3]
                @inbounds ex = etax[iv1, iv2, iv3]
                @inbounds ey = etay[iv1, iv2, iv3]
                @inbounds ez = etaz[iv1, iv2, iv3]

                local_F += d_val
                local_Sx += ex
                local_Sy += ey
                local_Sz += ez
                # Cast Int32 offsets to Float32 for multiplication
                fdi = Float32(di); fdj = Float32(dj); fdk = Float32(dk)
                local_Gx += d_val * fdi
                local_Gy += d_val * fdj
                local_Gz += d_val * fdk
                local_n += _I1
            end
            c += bdim
        end

        # 7 Float32 block reductions + 1 Int32. Sync between each reuse of sdata.
        F_tot  = block_reduce_sum_f32(local_F,  sdata_f); sync_threads()
        Sx_tot = block_reduce_sum_f32(local_Sx, sdata_f); sync_threads()
        Sy_tot = block_reduce_sum_f32(local_Sy, sdata_f); sync_threads()
        Sz_tot = block_reduce_sum_f32(local_Sz, sdata_f); sync_threads()
        Gx_tot = block_reduce_sum_f32(local_Gx, sdata_f); sync_threads()
        Gy_tot = block_reduce_sum_f32(local_Gy, sdata_f); sync_threads()
        Gz_tot = block_reduce_sum_f32(local_Gz, sdata_f); sync_threads()
        n_tot  = block_reduce_sum_i32(local_n,  sdata_i)

        if tid == _I1
            @inbounds Fshell_out[s, peak_id] = F_tot
            @inbounds nshell_out[s, peak_id] = n_tot
            @inbounds Sshell_out[_I1, s, peak_id] = Sx_tot
            @inbounds Sshell_out[Int32(2), s, peak_id] = Sy_tot
            @inbounds Sshell_out[Int32(3), s, peak_id] = Sz_tot
            @inbounds Gshell_out[_I1, s, peak_id] = Gx_tot
            @inbounds Gshell_out[Int32(2), s, peak_id] = Gy_tot
            @inbounds Gshell_out[Int32(3), s, peak_id] = Gz_tot
        end
        sync_threads()
        s += _I1
    end
    return
end

function PeakPatch.shell_gather_psi_gpu(delta_h::AbstractArray{<:Real,3},
                                        etax_h::AbstractArray{<:Real,3},
                                        etay_h::AbstractArray{<:Real,3},
                                        etaz_h::AbstractArray{<:Real,3},
                                        peaks_i::AbstractVector{<:Integer},
                                        peaks_j::AbstractVector{<:Integer},
                                        peaks_k::AbstractVector{<:Integer},
                                        stab::ShellTables;
                                        threads::Int=128)
    @assert size(delta_h) == size(etax_h) == size(etay_h) == size(etaz_h)
    @assert length(peaks_i) == length(peaks_j) == length(peaks_k)
    npeaks = length(peaks_i)
    nshells = stab.nshells

    delta_d = CuArray{Float32}(delta_h)
    etax_d  = CuArray{Float32}(etax_h)
    etay_d  = CuArray{Float32}(etay_h)
    etaz_d  = CuArray{Float32}(etaz_h)
    pi_d = CuArray{Int32}(peaks_i)
    pj_d = CuArray{Int32}(peaks_j)
    pk_d = CuArray{Int32}(peaks_k)
    stab_d = ShellTablesGPU(stab)

    Fshell_d = CUDA.zeros(Float32, nshells, npeaks)
    nshell_d = CUDA.zeros(Int32, nshells, npeaks)
    Sshell_d = CUDA.zeros(Float32, 3, nshells, npeaks)
    Gshell_d = CUDA.zeros(Float32, 3, nshells, npeaks)

    n1, n2, n3 = size(delta_d)
    shmem = 32 * sizeof(Float32) + 32 * sizeof(Int32)

    @cuda threads=threads blocks=npeaks shmem=shmem _shell_gather_psi_kernel!(
        Fshell_d, nshell_d, Sshell_d, Gshell_d,
        delta_d, etax_d, etay_d, etaz_d,
        pi_d, pj_d, pk_d,
        stab_d.offsets_di, stab_d.offsets_dj, stab_d.offsets_dk,
        stab_d.shell_start, stab_d.shell_count,
        Int32(nshells), Int32(n1), Int32(n2), Int32(n3),
    )

    return Array(Fshell_d), Array(nshell_d), Array(Sshell_d), Array(Gshell_d)
end

# ============================================================
# Block-per-peak shell gather kernel (delta + ψ×3 + strain tensor SR)
# ============================================================
#
# Extends `_shell_gather_psi_kernel!` with the 9-element raw strain tensor
#   SRshell[L,K] = Σ eta_L[cell] * d_K
# Not symmetric per-shell; symmetry is applied later via `Symmetric(Ebar)`.
#
# Per-thread accumulator count: 8 Float32 (F + 3 S + 3 G) + 9 Float32 (SR)
#   + 1 Int32 = 17 Float32 + 1 Int32. Still well within register budget on
#   Ampere (255 regs/thread max).
#
# Reduction count per shell: 16 Float32 + 1 Int32 = 17 block reductions.

function _shell_gather_strain_kernel!(Fshell_out, nshell_out,
                                      Sshell_out, Gshell_out, SRshell_out,
                                      delta, etax, etay, etaz,
                                      peaks_i, peaks_j, peaks_k,
                                      off_di, off_dj, off_dk,
                                      shell_start, shell_count,
                                      nshells::Int32, n1::Int32, n2::Int32, n3::Int32)
    peak_id = blockIdx().x
    tid = threadIdx().x
    bdim = blockDim().x

    sdata_f = CuDynamicSharedArray(Float32, 32)
    sdata_i = CuDynamicSharedArray(Int32, 32, 32 * sizeof(Float32))

    @inbounds ci = peaks_i[peak_id]
    @inbounds cj = peaks_j[peak_id]
    @inbounds ck = peaks_k[peak_id]

    s = _I1
    while s <= nshells
        @inbounds s0 = shell_start[s]
        @inbounds nc = shell_count[s]

        local_F = 0.0f0
        local_Sx = 0.0f0; local_Sy = 0.0f0; local_Sz = 0.0f0
        local_Gx = 0.0f0; local_Gy = 0.0f0; local_Gz = 0.0f0
        # Strain: 9 accumulators for SR[L,K], L ∈ {x,y,z}, K ∈ {x,y,z}
        local_SRxx = 0.0f0; local_SRxy = 0.0f0; local_SRxz = 0.0f0
        local_SRyx = 0.0f0; local_SRyy = 0.0f0; local_SRyz = 0.0f0
        local_SRzx = 0.0f0; local_SRzy = 0.0f0; local_SRzz = 0.0f0
        local_n = _I0

        c = tid
        while c <= nc
            idx = s0 + c - _I1
            @inbounds di = off_di[idx]
            @inbounds dj = off_dj[idx]
            @inbounds dk = off_dk[idx]

            iv1 = ci + di
            iv2 = cj + dj
            iv3 = ck + dk

            if (_I1 <= iv1) & (iv1 <= n1) & (_I1 <= iv2) & (iv2 <= n2) & (_I1 <= iv3) & (iv3 <= n3)
                @inbounds d_val = delta[iv1, iv2, iv3]
                @inbounds ex = etax[iv1, iv2, iv3]
                @inbounds ey = etay[iv1, iv2, iv3]
                @inbounds ez = etaz[iv1, iv2, iv3]

                local_F += d_val
                local_Sx += ex; local_Sy += ey; local_Sz += ez

                fdi = Float32(di); fdj = Float32(dj); fdk = Float32(dk)
                local_Gx += d_val * fdi
                local_Gy += d_val * fdj
                local_Gz += d_val * fdk

                local_SRxx += ex * fdi; local_SRxy += ex * fdj; local_SRxz += ex * fdk
                local_SRyx += ey * fdi; local_SRyy += ey * fdj; local_SRyz += ey * fdk
                local_SRzx += ez * fdi; local_SRzy += ez * fdj; local_SRzz += ez * fdk

                local_n += _I1
            end
            c += bdim
        end

        # 16 Float32 reductions + 1 Int32 reduction per shell.
        F_tot  = block_reduce_sum_f32(local_F,  sdata_f); sync_threads()
        Sx_tot = block_reduce_sum_f32(local_Sx, sdata_f); sync_threads()
        Sy_tot = block_reduce_sum_f32(local_Sy, sdata_f); sync_threads()
        Sz_tot = block_reduce_sum_f32(local_Sz, sdata_f); sync_threads()
        Gx_tot = block_reduce_sum_f32(local_Gx, sdata_f); sync_threads()
        Gy_tot = block_reduce_sum_f32(local_Gy, sdata_f); sync_threads()
        Gz_tot = block_reduce_sum_f32(local_Gz, sdata_f); sync_threads()
        SRxx_tot = block_reduce_sum_f32(local_SRxx, sdata_f); sync_threads()
        SRxy_tot = block_reduce_sum_f32(local_SRxy, sdata_f); sync_threads()
        SRxz_tot = block_reduce_sum_f32(local_SRxz, sdata_f); sync_threads()
        SRyx_tot = block_reduce_sum_f32(local_SRyx, sdata_f); sync_threads()
        SRyy_tot = block_reduce_sum_f32(local_SRyy, sdata_f); sync_threads()
        SRyz_tot = block_reduce_sum_f32(local_SRyz, sdata_f); sync_threads()
        SRzx_tot = block_reduce_sum_f32(local_SRzx, sdata_f); sync_threads()
        SRzy_tot = block_reduce_sum_f32(local_SRzy, sdata_f); sync_threads()
        SRzz_tot = block_reduce_sum_f32(local_SRzz, sdata_f); sync_threads()
        n_tot    = block_reduce_sum_i32(local_n,  sdata_i)

        if tid == _I1
            @inbounds Fshell_out[s, peak_id] = F_tot
            @inbounds nshell_out[s, peak_id] = n_tot
            @inbounds Sshell_out[_I1, s, peak_id] = Sx_tot
            @inbounds Sshell_out[Int32(2), s, peak_id] = Sy_tot
            @inbounds Sshell_out[Int32(3), s, peak_id] = Sz_tot
            @inbounds Gshell_out[_I1, s, peak_id] = Gx_tot
            @inbounds Gshell_out[Int32(2), s, peak_id] = Gy_tot
            @inbounds Gshell_out[Int32(3), s, peak_id] = Gz_tot
            @inbounds SRshell_out[_I1, _I1, s, peak_id] = SRxx_tot
            @inbounds SRshell_out[_I1, Int32(2), s, peak_id] = SRxy_tot
            @inbounds SRshell_out[_I1, Int32(3), s, peak_id] = SRxz_tot
            @inbounds SRshell_out[Int32(2), _I1, s, peak_id] = SRyx_tot
            @inbounds SRshell_out[Int32(2), Int32(2), s, peak_id] = SRyy_tot
            @inbounds SRshell_out[Int32(2), Int32(3), s, peak_id] = SRyz_tot
            @inbounds SRshell_out[Int32(3), _I1, s, peak_id] = SRzx_tot
            @inbounds SRshell_out[Int32(3), Int32(2), s, peak_id] = SRzy_tot
            @inbounds SRshell_out[Int32(3), Int32(3), s, peak_id] = SRzz_tot
        end
        sync_threads()
        s += _I1
    end
    return
end

function PeakPatch.shell_gather_strain_gpu(delta_h::AbstractArray{<:Real,3},
                                           etax_h::AbstractArray{<:Real,3},
                                           etay_h::AbstractArray{<:Real,3},
                                           etaz_h::AbstractArray{<:Real,3},
                                           peaks_i::AbstractVector{<:Integer},
                                           peaks_j::AbstractVector{<:Integer},
                                           peaks_k::AbstractVector{<:Integer},
                                           stab::ShellTables;
                                           threads::Int=128)
    @assert size(delta_h) == size(etax_h) == size(etay_h) == size(etaz_h)
    @assert length(peaks_i) == length(peaks_j) == length(peaks_k)
    npeaks = length(peaks_i)
    nshells = stab.nshells

    delta_d = CuArray{Float32}(delta_h)
    etax_d  = CuArray{Float32}(etax_h)
    etay_d  = CuArray{Float32}(etay_h)
    etaz_d  = CuArray{Float32}(etaz_h)
    pi_d = CuArray{Int32}(peaks_i)
    pj_d = CuArray{Int32}(peaks_j)
    pk_d = CuArray{Int32}(peaks_k)
    stab_d = ShellTablesGPU(stab)

    Fshell_d  = CUDA.zeros(Float32, nshells, npeaks)
    nshell_d  = CUDA.zeros(Int32,   nshells, npeaks)
    Sshell_d  = CUDA.zeros(Float32, 3, nshells, npeaks)
    Gshell_d  = CUDA.zeros(Float32, 3, nshells, npeaks)
    SRshell_d = CUDA.zeros(Float32, 3, 3, nshells, npeaks)

    n1, n2, n3 = size(delta_d)
    shmem = 32 * sizeof(Float32) + 32 * sizeof(Int32)

    @cuda threads=threads blocks=npeaks shmem=shmem _shell_gather_strain_kernel!(
        Fshell_d, nshell_d, Sshell_d, Gshell_d, SRshell_d,
        delta_d, etax_d, etay_d, etaz_d,
        pi_d, pj_d, pk_d,
        stab_d.offsets_di, stab_d.offsets_dj, stab_d.offsets_dk,
        stab_d.shell_start, stab_d.shell_count,
        Int32(nshells), Int32(n1), Int32(n2), Int32(n3),
    )

    return Array(Fshell_d), Array(nshell_d), Array(Sshell_d),
           Array(Gshell_d), Array(SRshell_d)
end

# ============================================================
# Block-per-peak shell gather kernel (all fields: + 2LPT + Gf)
# ============================================================
#
# Full-fidelity port of the CPU gather in `analyse_peak_gpu`
# (`src/ShellAnalysisGPU.jl:196-232`).  Adds to `_shell_gather_strain_kernel!`:
#   S2shell[L]  = Σ eta2_L[cell]                           (2LPT; pass zeros to skip)
#   Gfshell[L]  = Σ delta[iv3, iv1, iv2] * d_L             (transposed gradient)
#
# Per-cell memory reads: 8 Float32 (delta × 2 reads, 3 eta, 3 eta2) + 3 Int32 offsets.
# Per-thread accumulators: 22 Float32 + 1 Int32.

function _shell_gather_full_kernel!(Fshell_out, nshell_out,
                                    Sshell_out, S2shell_out, Gshell_out,
                                    Gfshell_out, SRshell_out, lapdshell_out,
                                    delta, etax, etay, etaz,
                                    eta2x, eta2y, eta2z, lapd,
                                    peaks_i, peaks_j, peaks_k,
                                    off_di, off_dj, off_dk,
                                    shell_start, shell_count,
                                    nshells::Int32, n1::Int32, n2::Int32, n3::Int32)
    peak_id = blockIdx().x
    tid = threadIdx().x
    bdim = blockDim().x

    sdata_f = CuDynamicSharedArray(Float32, 32)
    sdata_i = CuDynamicSharedArray(Int32, 32, 32 * sizeof(Float32))

    @inbounds ci = peaks_i[peak_id]
    @inbounds cj = peaks_j[peak_id]
    @inbounds ck = peaks_k[peak_id]

    s = _I1
    while s <= nshells
        @inbounds s0 = shell_start[s]
        @inbounds nc = shell_count[s]

        local_F = 0.0f0
        local_Sx  = 0.0f0; local_Sy  = 0.0f0; local_Sz  = 0.0f0
        local_S2x = 0.0f0; local_S2y = 0.0f0; local_S2z = 0.0f0
        local_Gx  = 0.0f0; local_Gy  = 0.0f0; local_Gz  = 0.0f0
        local_Gfx = 0.0f0; local_Gfy = 0.0f0; local_Gfz = 0.0f0
        local_SRxx = 0.0f0; local_SRxy = 0.0f0; local_SRxz = 0.0f0
        local_SRyx = 0.0f0; local_SRyy = 0.0f0; local_SRyz = 0.0f0
        local_SRzx = 0.0f0; local_SRzy = 0.0f0; local_SRzz = 0.0f0
        local_lapd = 0.0f0
        local_n = _I0

        c = tid
        while c <= nc
            idx = s0 + c - _I1
            @inbounds di = off_di[idx]
            @inbounds dj = off_dj[idx]
            @inbounds dk = off_dk[idx]

            iv1 = ci + di
            iv2 = cj + dj
            iv3 = ck + dk

            if (_I1 <= iv1) & (iv1 <= n1) & (_I1 <= iv2) & (iv2 <= n2) & (_I1 <= iv3) & (iv3 <= n3)
                @inbounds d_val = delta[iv1, iv2, iv3]
                # Transposed read for Gfshell: delta[iv3, iv1, iv2]
                # Note: transpose indices are still in-range because n1=n2=n3 in
                # practice (cubic grid); we already bounds-checked all three.
                @inbounds df_val = delta[iv3, iv1, iv2]
                @inbounds ex = etax[iv1, iv2, iv3]
                @inbounds ey = etay[iv1, iv2, iv3]
                @inbounds ez = etaz[iv1, iv2, iv3]
                @inbounds e2x = eta2x[iv1, iv2, iv3]
                @inbounds e2y = eta2y[iv1, iv2, iv3]
                @inbounds e2z = eta2z[iv1, iv2, iv3]

                local_F += d_val
                local_Sx  += ex;  local_Sy  += ey;  local_Sz  += ez
                local_S2x += e2x; local_S2y += e2y; local_S2z += e2z

                fdi = Float32(di); fdj = Float32(dj); fdk = Float32(dk)
                local_Gx  += d_val  * fdi; local_Gy  += d_val  * fdj; local_Gz  += d_val  * fdk
                local_Gfx += df_val * fdi; local_Gfy += df_val * fdj; local_Gfz += df_val * fdk

                local_SRxx += ex * fdi; local_SRxy += ex * fdj; local_SRxz += ex * fdk
                local_SRyx += ey * fdi; local_SRyy += ey * fdj; local_SRyz += ey * fdk
                local_SRzx += ez * fdi; local_SRzy += ez * fdj; local_SRzz += ez * fdk

                @inbounds local_lapd += lapd[iv1, iv2, iv3]
                local_n += _I1
            end
            c += bdim
        end

        # 22 Float32 block reductions + 1 Int32 reduction per shell.
        F_tot   = block_reduce_sum_f32(local_F,   sdata_f); sync_threads()
        Sx_tot  = block_reduce_sum_f32(local_Sx,  sdata_f); sync_threads()
        Sy_tot  = block_reduce_sum_f32(local_Sy,  sdata_f); sync_threads()
        Sz_tot  = block_reduce_sum_f32(local_Sz,  sdata_f); sync_threads()
        S2x_tot = block_reduce_sum_f32(local_S2x, sdata_f); sync_threads()
        S2y_tot = block_reduce_sum_f32(local_S2y, sdata_f); sync_threads()
        S2z_tot = block_reduce_sum_f32(local_S2z, sdata_f); sync_threads()
        Gx_tot  = block_reduce_sum_f32(local_Gx,  sdata_f); sync_threads()
        Gy_tot  = block_reduce_sum_f32(local_Gy,  sdata_f); sync_threads()
        Gz_tot  = block_reduce_sum_f32(local_Gz,  sdata_f); sync_threads()
        Gfx_tot = block_reduce_sum_f32(local_Gfx, sdata_f); sync_threads()
        Gfy_tot = block_reduce_sum_f32(local_Gfy, sdata_f); sync_threads()
        Gfz_tot = block_reduce_sum_f32(local_Gfz, sdata_f); sync_threads()
        SRxx_tot = block_reduce_sum_f32(local_SRxx, sdata_f); sync_threads()
        SRxy_tot = block_reduce_sum_f32(local_SRxy, sdata_f); sync_threads()
        SRxz_tot = block_reduce_sum_f32(local_SRxz, sdata_f); sync_threads()
        SRyx_tot = block_reduce_sum_f32(local_SRyx, sdata_f); sync_threads()
        SRyy_tot = block_reduce_sum_f32(local_SRyy, sdata_f); sync_threads()
        SRyz_tot = block_reduce_sum_f32(local_SRyz, sdata_f); sync_threads()
        SRzx_tot = block_reduce_sum_f32(local_SRzx, sdata_f); sync_threads()
        SRzy_tot = block_reduce_sum_f32(local_SRzy, sdata_f); sync_threads()
        SRzz_tot = block_reduce_sum_f32(local_SRzz, sdata_f); sync_threads()
        lapd_tot = block_reduce_sum_f32(local_lapd, sdata_f); sync_threads()
        n_tot    = block_reduce_sum_i32(local_n,   sdata_i)

        if tid == _I1
            @inbounds Fshell_out[s, peak_id] = F_tot
            @inbounds nshell_out[s, peak_id] = n_tot
            @inbounds Sshell_out[_I1, s, peak_id] = Sx_tot
            @inbounds Sshell_out[Int32(2), s, peak_id] = Sy_tot
            @inbounds Sshell_out[Int32(3), s, peak_id] = Sz_tot
            @inbounds S2shell_out[_I1, s, peak_id] = S2x_tot
            @inbounds S2shell_out[Int32(2), s, peak_id] = S2y_tot
            @inbounds S2shell_out[Int32(3), s, peak_id] = S2z_tot
            @inbounds Gshell_out[_I1, s, peak_id] = Gx_tot
            @inbounds Gshell_out[Int32(2), s, peak_id] = Gy_tot
            @inbounds Gshell_out[Int32(3), s, peak_id] = Gz_tot
            @inbounds Gfshell_out[_I1, s, peak_id] = Gfx_tot
            @inbounds Gfshell_out[Int32(2), s, peak_id] = Gfy_tot
            @inbounds Gfshell_out[Int32(3), s, peak_id] = Gfz_tot
            @inbounds SRshell_out[_I1, _I1, s, peak_id] = SRxx_tot
            @inbounds SRshell_out[_I1, Int32(2), s, peak_id] = SRxy_tot
            @inbounds SRshell_out[_I1, Int32(3), s, peak_id] = SRxz_tot
            @inbounds SRshell_out[Int32(2), _I1, s, peak_id] = SRyx_tot
            @inbounds SRshell_out[Int32(2), Int32(2), s, peak_id] = SRyy_tot
            @inbounds SRshell_out[Int32(2), Int32(3), s, peak_id] = SRyz_tot
            @inbounds SRshell_out[Int32(3), _I1, s, peak_id] = SRzx_tot
            @inbounds SRshell_out[Int32(3), Int32(2), s, peak_id] = SRzy_tot
            @inbounds SRshell_out[Int32(3), Int32(3), s, peak_id] = SRzz_tot
            @inbounds lapdshell_out[s, peak_id] = lapd_tot
        end
        sync_threads()
        s += _I1
    end
    return
end

# Internal launcher: operates entirely on pre-allocated GPU arrays.
# Exposed so `analyse_peak_gpu_cuda` can drive the gather without a
# host→GPU→host round-trip on the per-shell profiles.
function _launch_shell_gather_full!(
        # GPU outputs (pre-allocated)
        Fshell_d, nshell_d, Sshell_d, S2shell_d,
        Gshell_d, Gfshell_d, SRshell_d, lapdshell_d,
        # GPU inputs (already uploaded)
        delta_d, etax_d, etay_d, etaz_d,
        eta2x_d, eta2y_d, eta2z_d, lapd_d,
        pi_d, pj_d, pk_d,
        stab_d::ShellTablesGPU;
        threads::Int=128)
    npeaks  = length(pi_d)
    nshells = stab_d.nshells
    n1, n2, n3 = size(delta_d)
    shmem = 32 * sizeof(Float32) + 32 * sizeof(Int32)
    @cuda threads=threads blocks=npeaks shmem=shmem _shell_gather_full_kernel!(
        Fshell_d, nshell_d, Sshell_d, S2shell_d, Gshell_d, Gfshell_d, SRshell_d,
        lapdshell_d,
        delta_d, etax_d, etay_d, etaz_d, eta2x_d, eta2y_d, eta2z_d, lapd_d,
        pi_d, pj_d, pk_d,
        stab_d.offsets_di, stab_d.offsets_dj, stab_d.offsets_dk,
        stab_d.shell_start, stab_d.shell_count,
        Int32(nshells), Int32(n1), Int32(n2), Int32(n3),
    )
    return nothing
end

function PeakPatch.shell_gather_full_gpu(delta_h::AbstractArray{<:Real,3},
                                         etax_h::AbstractArray{<:Real,3},
                                         etay_h::AbstractArray{<:Real,3},
                                         etaz_h::AbstractArray{<:Real,3},
                                         eta2x_h::AbstractArray{<:Real,3},
                                         eta2y_h::AbstractArray{<:Real,3},
                                         eta2z_h::AbstractArray{<:Real,3},
                                         peaks_i::AbstractVector{<:Integer},
                                         peaks_j::AbstractVector{<:Integer},
                                         peaks_k::AbstractVector{<:Integer},
                                         stab::ShellTables;
                                         threads::Int=128,
                                         lapd::Union{Nothing, AbstractArray{<:Real,3}}=nothing,
                                         return_lapd::Bool=false)
    @assert size(delta_h) == size(etax_h) == size(etay_h) == size(etaz_h) ==
            size(eta2x_h) == size(eta2y_h) == size(eta2z_h)
    @assert length(peaks_i) == length(peaks_j) == length(peaks_k)
    npeaks = length(peaks_i)
    nshells = stab.nshells

    delta_d = CuArray{Float32}(delta_h)
    etax_d  = CuArray{Float32}(etax_h)
    etay_d  = CuArray{Float32}(etay_h)
    etaz_d  = CuArray{Float32}(etaz_h)
    eta2x_d = CuArray{Float32}(eta2x_h)
    eta2y_d = CuArray{Float32}(eta2y_h)
    eta2z_d = CuArray{Float32}(eta2z_h)
    # lapd: upload if provided, else allocate zeros (no-op contribution)
    lapd_d = if lapd === nothing
        CUDA.zeros(Float32, size(delta_h)...)
    else
        @assert size(lapd) == size(delta_h)
        CuArray{Float32}(lapd)
    end
    pi_d = CuArray{Int32}(peaks_i)
    pj_d = CuArray{Int32}(peaks_j)
    pk_d = CuArray{Int32}(peaks_k)
    stab_d = ShellTablesGPU(stab)

    Fshell_d    = CUDA.zeros(Float32, nshells, npeaks)
    nshell_d    = CUDA.zeros(Int32,   nshells, npeaks)
    Sshell_d    = CUDA.zeros(Float32, 3, nshells, npeaks)
    S2shell_d   = CUDA.zeros(Float32, 3, nshells, npeaks)
    Gshell_d    = CUDA.zeros(Float32, 3, nshells, npeaks)
    Gfshell_d   = CUDA.zeros(Float32, 3, nshells, npeaks)
    SRshell_d   = CUDA.zeros(Float32, 3, 3, nshells, npeaks)
    lapdshell_d = CUDA.zeros(Float32, nshells, npeaks)

    _launch_shell_gather_full!(
        Fshell_d, nshell_d, Sshell_d, S2shell_d, Gshell_d, Gfshell_d, SRshell_d,
        lapdshell_d,
        delta_d, etax_d, etay_d, etaz_d, eta2x_d, eta2y_d, eta2z_d, lapd_d,
        pi_d, pj_d, pk_d, stab_d; threads=threads)

    # Backward-compatible return: existing tests unpack 7 values. Opt-in to
    # the 8th (lapdshell) via return_lapd=true for callers that need it.
    if return_lapd
        return (Array(Fshell_d), Array(nshell_d),
                Array(Sshell_d), Array(S2shell_d),
                Array(Gshell_d), Array(Gfshell_d),
                Array(SRshell_d), Array(lapdshell_d))
    else
        return (Array(Fshell_d), Array(nshell_d),
                Array(Sshell_d), Array(S2shell_d),
                Array(Gshell_d), Array(Gfshell_d),
                Array(SRshell_d))
    end
end

# ============================================================
# Batch 3×3 eigendecomposition (test harness for `_eig3_symmetric_f32`)
# ============================================================

function _eig3_batch_kernel!(eigvals_out, mats, n::Int32)
    idx = (blockIdx().x - _I1) * blockDim().x + threadIdx().x
    if idx <= n
        @inbounds a11 = mats[1, 1, idx]
        @inbounds a22 = mats[2, 2, idx]
        @inbounds a33 = mats[3, 3, idx]
        # Symmetrise off-diagonals: average (i,j) and (j,i) entries
        @inbounds a12 = 0.5f0 * (mats[1, 2, idx] + mats[2, 1, idx])
        @inbounds a13 = 0.5f0 * (mats[1, 3, idx] + mats[3, 1, idx])
        @inbounds a23 = 0.5f0 * (mats[2, 3, idx] + mats[3, 2, idx])
        e1, e2, e3 = _eig3_symmetric_f32(a11, a22, a33, a12, a13, a23)
        @inbounds eigvals_out[1, idx] = e1
        @inbounds eigvals_out[2, idx] = e2
        @inbounds eigvals_out[3, idx] = e3
    end
    return
end

function PeakPatch.eig3_symmetric_batch_gpu(mats_h::AbstractArray{<:Real,3};
                                            threads::Int=128)
    @assert size(mats_h, 1) == 3 && size(mats_h, 2) == 3
    n = size(mats_h, 3)
    mats_d = CuArray{Float32}(mats_h)
    out_d = CUDA.zeros(Float32, 3, n)
    blocks = cld(n, threads)
    @cuda threads=threads blocks=blocks _eig3_batch_kernel!(out_d, mats_d, Int32(n))
    return Array(out_d)
end

# ============================================================
# Batch trilinear interpolation (test harness for _interp3_trilinear_f32)
# ============================================================

function _interp3_batch_kernel!(out, table, xs, ys, zs, n::Int32,
                                 x1::Float32, y1::Float32, z1::Float32,
                                 dx_inv::Float32, dy_inv::Float32, dz_inv::Float32,
                                 nx::Int32, ny::Int32, nz::Int32,
                                 out_val::Float32)
    idx = (blockIdx().x - _I1) * blockDim().x + threadIdx().x
    if idx <= n
        @inbounds x = xs[idx]
        @inbounds y = ys[idx]
        @inbounds z = zs[idx]
        v = _interp3_trilinear_f32(table, x, y, z, x1, y1, z1,
                                    dx_inv, dy_inv, dz_inv,
                                    nx, ny, nz, out_val)
        @inbounds out[idx] = v
    end
    return
end

function PeakPatch.interp3_trilinear_gpu(table_h::AbstractArray{<:Real,3},
                                         X1::Real, X2::Real,
                                         Y1::Real, Y2::Real,
                                         Z1::Real, Z2::Real,
                                         xs_h::AbstractVector{<:Real},
                                         ys_h::AbstractVector{<:Real},
                                         zs_h::AbstractVector{<:Real},
                                         out_val::Real=-1.0;
                                         threads::Int=128)
    @assert length(xs_h) == length(ys_h) == length(zs_h)
    n = length(xs_h)
    nx, ny, nz = size(table_h)

    table_d = CuArray{Float32}(table_h)
    xs_d = CuArray{Float32}(xs_h); ys_d = CuArray{Float32}(ys_h); zs_d = CuArray{Float32}(zs_h)
    out_d = CUDA.zeros(Float32, n)

    # Step inverse = (nx - 1) / (X2 - X1) for range(X1, X2, nx)
    dx_inv = Float32((nx - 1) / (X2 - X1))
    dy_inv = Float32((ny - 1) / (Y2 - Y1))
    dz_inv = Float32((nz - 1) / (Z2 - Z1))

    blocks = cld(n, threads)
    @cuda threads=threads blocks=blocks _interp3_batch_kernel!(
        out_d, table_d, xs_d, ys_d, zs_d, Int32(n),
        Float32(X1), Float32(Y1), Float32(Z1),
        dx_inv, dy_inv, dz_inv,
        Int32(nx), Int32(ny), Int32(nz),
        Float32(out_val),
    )
    return Array(out_d)
end

# ============================================================
# kernel_strain test harness
# ============================================================

function _kernel_strain_test_kernel!(out, rad, Gshell, Gshellf, SRshell, akk,
                                     mp::Int32, mlow::Int32, mupp::Int32,
                                     wRnor::Float32, aRnor::Float32,
                                     hlatt_1::Float32, hlatt_2::Float32)
    # Run once on thread (1,1,1)
    if threadIdx().x == _I1 && blockIdx().x == _I1
        (E11, E12, E13, E22, E23, E33, gx, gy, gz, gfx, gfy, gfz) =
            _kernel_strain_f32(rad, Gshell, Gshellf, SRshell, akk,
                                mp, mlow, mupp, wRnor, aRnor, hlatt_1, hlatt_2)
        @inbounds out[1]  = E11
        @inbounds out[2]  = E12
        @inbounds out[3]  = E13
        @inbounds out[4]  = E22
        @inbounds out[5]  = E23
        @inbounds out[6]  = E33
        @inbounds out[7]  = gx
        @inbounds out[8]  = gy
        @inbounds out[9]  = gz
        @inbounds out[10] = gfx
        @inbounds out[11] = gfy
        @inbounds out[12] = gfz
    end
    return
end

function PeakPatch.kernel_strain_gpu(rad_h::AbstractVector{<:Real},
                                     Gshell_h::AbstractMatrix{<:Real},
                                     Gshellf_h::AbstractMatrix{<:Real},
                                     SRshell_h::AbstractArray{<:Real,3},
                                     mp::Integer, mlow::Integer, mupp::Integer,
                                     wRnor::Real, aRnor::Real,
                                     hlatt_1::Real, hlatt_2::Real)
    @assert size(Gshell_h, 1) == 3 && size(Gshellf_h, 1) == 3
    @assert size(SRshell_h, 1) == 3 && size(SRshell_h, 2) == 3

    rad_d     = CuArray{Float32}(rad_h)
    Gshell_d  = CuArray{Float32}(Gshell_h)
    Gshellf_d = CuArray{Float32}(Gshellf_h)
    SRshell_d = CuArray{Float32}(SRshell_h)
    out_d     = CUDA.zeros(Float32, 12)

    @cuda threads=1 blocks=1 _kernel_strain_test_kernel!(
        out_d, rad_d, Gshell_d, Gshellf_d, SRshell_d, _AKK_GPU[],
        Int32(mp), Int32(mlow), Int32(mupp),
        Float32(wRnor), Float32(aRnor),
        Float32(hlatt_1), Float32(hlatt_2),
    )
    out = Array(out_d)
    Ebar = Float32[out[1] out[2] out[3];
                   out[2] out[4] out[5];
                   out[3] out[5] out[6]]
    grad  = Float32[out[7],  out[8],  out[9]]
    gradf = Float32[out[10], out[11], out[12]]
    return (Ebar, grad, gradf)
end

# ============================================================
# Post-process kernel: phases 4 (fcrit) + 5 (inward walk) + 6 (RTHL interp)
# ============================================================
#
# Two-kernel pipeline for `analyse_peak_gpu_cuda`:
#   1. `_shell_gather_full_kernel!`  - parallel gather of all 7 field profiles
#   2. `_post_process_kernel!`       - sequential phases 4-6 on thread 0/block
#
# One block per peak; thread 0 runs the post-process, other threads are idle.
# Per-peak profile data lives in shared memory (loaded cooperatively from
# global), so thread 0's sequential reads hit SRAM not HBM.
#
# Phase 3 (mrf + gradpkrf) and phase 7 (Sbar, Srb, d2F, mask) are not yet
# implemented — they're trivially sequential once phases 4-6 work.
#
# Output: compact Float32 fields per peak — RTHL, Fbarx, e_v, p_v, strain_final
# (3×3 symmetric), eigenvalues of strain_final (ascending).
#
# RTHL = -1.0f0 signals no_collapse (same convention as CPU).

# Compile-time upper bound on shells per peak. rmax=30 gives ~30 shells;
# this bound of 200 handles rmax up to ~70 comfortably.
const _MAX_SHELLS_GPU = Int32(200)

# Device helper: write `no_collapse` zeros across all PeakResult outputs.
# Mirrors `RadialShell.no_collapse()` which zeros every field except
# RTHL = -1.0 and zvir_half = -1.0.
@inline function _write_no_collapse_outputs!(
        RTHL_out, Fbarx_out, e_v_out, p_v_out,
        strain_final_out, eigs_out,
        Srb_out, Sbar_out, Sbar2_out,
        gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
        zvir_half_out,
        peak_id)
    @inbounds RTHL_out[peak_id]  = -1.0f0
    @inbounds Fbarx_out[peak_id] = 0.0f0
    @inbounds e_v_out[peak_id]   = 0.0f0
    @inbounds p_v_out[peak_id]   = 0.0f0
    @inbounds Srb_out[peak_id]   = 0.0f0
    @inbounds d2F_out[peak_id]   = 0.0f0
    @inbounds zvir_half_out[peak_id] = -1.0f0
    for L in _I1:Int32(3)
        @inbounds Sbar_out[L, peak_id]     = 0.0f0
        @inbounds Sbar2_out[L, peak_id]    = 0.0f0
        @inbounds gradpk_out[L, peak_id]   = 0.0f0
        @inbounds gradpkf_out[L, peak_id]  = 0.0f0
        @inbounds gradpkrf_out[L, peak_id] = 0.0f0
        @inbounds eigs_out[L, peak_id]     = 0.0f0
        for K in _I1:Int32(3)
            @inbounds strain_final_out[L, K, peak_id] = 0.0f0
        end
    end
    return
end

function _post_process_kernel!(
    # Outputs (indexed by peak)
    RTHL_out,           # Vector{Float32} (npeaks)
    Fbarx_out,          # Vector{Float32}
    e_v_out,            # Vector{Float32}
    p_v_out,            # Vector{Float32}
    strain_final_out,   # Array{Float32,3} (3, 3, npeaks)
    eigs_out,           # Array{Float32,2} (3, npeaks)
    Srb_out,            # Vector{Float32}
    Sbar_out,           # Array{Float32,2} (3, npeaks)
    Sbar2_out,          # Array{Float32,2} (3, npeaks)
    gradpk_out,         # Array{Float32,2} (3, npeaks)
    gradpkf_out,        # Array{Float32,2} (3, npeaks)
    gradpkrf_out,       # Array{Float32,2} (3, npeaks)
    d2F_out,            # Vector{Float32} — Laplacian average within RTHL
    zvir_half_out,      # Vector{Float32} — formation redshift (collapse z at RTHL/2)
    mask_out,           # Array{Int8,3} (n1, n2, n3) — side-effect, cells within RTHL → 1
    # Inputs from gather (un-normalised)
    Fshell, nshell_in, Sshell_in, S2shell_in, Gshell_raw, Gfshell_raw, SRshell_raw,
    lapdshell_in,
    # Peak centers for mask update
    peaks_i, peaks_j, peaks_k,
    # Shell offset tables (also for mask update)
    off_di, off_dj, off_dk, shell_start, shell_count,
    # Shell geometry
    shell_r2,
    # Collapse table
    ct_table, ct_X1::Float32, ct_Y1::Float32, ct_Z1::Float32,
    ct_dxi::Float32, ct_dyi::Float32, ct_dzi::Float32,
    ct_nx::Int32, ct_ny::Int32, ct_nz::Int32, ct_out_val::Float32,
    # Kernel tables & scaling
    akk_tab,
    # Per-peak collapse state: ZZon_pp[i] = 1+z_collapse, fcrit_pp[i] the
    # corresponding Fbar threshold. For ievol==0 both are broadcast scalars;
    # for ievol==1 each peak carries its own lightcone redshift.
    ZZon_pp, fcrit_pp, Rfclvi_r2::Float32,
    wRnor::Float32, aRnor::Float32, hlatt_1::Float32, hlatt_2::Float32,
    ir2min::Int32,
    nshells::Int32, nshells_max::Int32,
    n1::Int32, n2::Int32, n3::Int32, nbuff::Int32, update_mask::Bool)

    peak_id = blockIdx().x
    @inbounds ZZon  = ZZon_pp[peak_id]
    @inbounds fcrit = fcrit_pp[peak_id]

    # Shared-memory per-peak profile arrays. Sized to compile-time upper
    # bound, though we only use the first `nshells_max` slots.
    shmem_rad   = CuDynamicSharedArray(Float32, _MAX_SHELLS_GPU, 0)
    shmem_Fbar  = CuDynamicSharedArray(Float32, _MAX_SHELLS_GPU,
                                        Int(_MAX_SHELLS_GPU) * sizeof(Float32))
    shmem_Gn    = CuDynamicSharedArray(Float32, (Int32(3), _MAX_SHELLS_GPU),
                                        Int(_MAX_SHELLS_GPU) * sizeof(Float32) * 2)
    shmem_Gfn   = CuDynamicSharedArray(Float32, (Int32(3), _MAX_SHELLS_GPU),
                                        Int(_MAX_SHELLS_GPU) * sizeof(Float32) * (2 + 3))
    shmem_SRn   = CuDynamicSharedArray(Float32, (Int32(3), Int32(3), _MAX_SHELLS_GPU),
                                        Int(_MAX_SHELLS_GPU) * sizeof(Float32) * (2 + 3 + 3))

    tid = threadIdx().x
    if tid != _I1
        return
    end

    # ============================================================
    # Build normalised per-shell profile for this peak
    # ============================================================
    # Shell 1 (center, r²=0): rad=0, special-case Fbar and zero G/Gf/SR
    @inbounds n1_cnt = nshell_in[1, peak_id]
    @inbounds F1_sum = Fshell[1, peak_id]
    dFbar_1 = n1_cnt > _I0 ? F1_sum / Float32(n1_cnt) : 0.0f0
    @inbounds shmem_rad[1]  = 0.0f0
    @inbounds shmem_Fbar[1] = dFbar_1
    @inbounds shmem_Gn[_I1, 1]  = 0.0f0
    @inbounds shmem_Gn[Int32(2), 1] = 0.0f0
    @inbounds shmem_Gn[Int32(3), 1] = 0.0f0
    @inbounds shmem_Gfn[_I1, 1] = 0.0f0
    @inbounds shmem_Gfn[Int32(2), 1] = 0.0f0
    @inbounds shmem_Gfn[Int32(3), 1] = 0.0f0
    @inbounds shmem_SRn[_I1, _I1, 1] = 0.0f0
    @inbounds shmem_SRn[_I1, Int32(2), 1] = 0.0f0
    @inbounds shmem_SRn[_I1, Int32(3), 1] = 0.0f0
    @inbounds shmem_SRn[Int32(2), _I1, 1] = 0.0f0
    @inbounds shmem_SRn[Int32(2), Int32(2), 1] = 0.0f0
    @inbounds shmem_SRn[Int32(2), Int32(3), 1] = 0.0f0
    @inbounds shmem_SRn[Int32(3), _I1, 1] = 0.0f0
    @inbounds shmem_SRn[Int32(3), Int32(2), 1] = 0.0f0
    @inbounds shmem_SRn[Int32(3), Int32(3), 1] = 0.0f0

    dFbarp = dFbar_1
    s = Int32(2)
    while s <= nshells_max
        @inbounds r2 = shell_r2[s]
        rad_s = sqrt(Float32(r2))
        @inbounds shmem_rad[s] = rad_s
        rad_s_inv = 1.0f0 / rad_s
        @inbounds n_s = nshell_in[s, peak_id]
        @inbounds F_s = Fshell[s, peak_id]
        dFbar_s = n_s > _I0 ? F_s / Float32(n_s) : 0.0f0
        @inbounds rad_sp = shmem_rad[s - _I1]
        rad3p = rad_sp * rad_sp * rad_sp
        rad3  = rad_s * rad_s * rad_s
        @inbounds Fbarp = shmem_Fbar[s - _I1]
        Fbar_s = (rad3p * Fbarp + 0.5f0 * (dFbarp + dFbar_s) * (rad3 - rad3p)) / rad3
        @inbounds shmem_Fbar[s] = Fbar_s
        dFbarp = dFbar_s
        # Normalise G, Gf, SR by rad_s (matches CPU Gshell[L][m] = local_G[L] / rad[m])
        for L in _I1:Int32(3)
            @inbounds shmem_Gn[L, s]  = Gshell_raw[L, s, peak_id]  * rad_s_inv
            @inbounds shmem_Gfn[L, s] = Gfshell_raw[L, s, peak_id] * rad_s_inv
            for K in _I1:Int32(3)
                @inbounds shmem_SRn[L, K, s] = SRshell_raw[L, K, s, peak_id] * rad_s_inv
            end
        end
        s += _I1
    end
    m = nshells_max

    # ============================================================
    # Phase 3: gradient at filter scale (gradpkrf)
    # ============================================================
    # Find mrf: shell index closest to Rfclvi/alatt in r² distance.
    # CPU analyse_peak_gpu (line 280-301) iterates over shells s=2..nshells_max
    # and tracks mrfi (increments once per distinct r² value). Since our
    # build_shell_tables assigns unique r² per shell, mrfi == s-1 after s=2.
    # We replicate this by: mrf defaults to 1; walk s=2..nshells_max; when
    # |r²(s) - rfclvi_r2| is a new minimum, mrf ← mrfi; increment mrfi.
    mrf = _I1
    rmrf = 10.0f0
    mrfi = _I1
    s_mrf = Int32(2)
    while s_mrf <= nshells_max
        @inbounds r2_s = shell_r2[s_mrf]
        dist_mrf = abs(Float32(r2_s) - Rfclvi_r2)
        if dist_mrf < rmrf
            mrf = mrfi
            rmrf = dist_mrf
        end
        if s_mrf < nshells_max
            mrfi += _I1
        end
        s_mrf += _I1
    end
    # Compute mlow_rf, mupp_rf
    @inbounds rlow_rf = max(shmem_rad[mrf] - 2.0f0, 0.0f0)
    mlow_rf = mrf
    if mrf > _I1
        m1_rf = mrf - _I1
        while m1_rf >= _I1
            @inbounds if shmem_rad[m1_rf] > rlow_rf
                mlow_rf = m1_rf
            end
            m1_rf -= _I1
        end
    end
    @inbounds rupp_rf = shmem_rad[mrf] + 2.0f0
    mupp_rf = mrf
    m1_rf = mrf + _I1
    while m1_rf <= m
        @inbounds if shmem_rad[m1_rf] < rupp_rf
            mupp_rf = m1_rf
        end
        m1_rf += _I1
    end
    # kernel_strain at mrf — CPU destructures `_, gradpkrf, _` meaning
    # gradpkrf = grad (the MIDDLE return). Our tuple is (E×6, gx/gy/gz, gfx/gfy/gfz)
    # so gradpkrf = (gx, gy, gz).
    (_, _, _, _, _, _,
     gradpkrf_x, gradpkrf_y, gradpkrf_z, _, _, _) = _kernel_strain_f32(
        shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
        mrf, mlow_rf, mupp_rf, wRnor, aRnor, hlatt_1, hlatt_2)

    # ============================================================
    # Phase 4a: reconstruct the initial m0/mupp tracking from gather loop
    # ============================================================
    m0 = _I1
    ifcrit = _I1
    rupp_init = sqrt(Float32(ir2min)) + 2.0f0
    mupp = _I1
    s = Int32(2)
    while s <= nshells_max
        @inbounds r2 = shell_r2[s]
        if ifcrit == _I1
            m0 = s
            if Float32(r2) > Float32(ir2min)
                @inbounds rupp_init = shmem_rad[m0] + 2.0f0
                ifcrit = _I0
                mupp = s
            end
        else
            @inbounds if shmem_rad[s] < rupp_init
                mupp = s
            end
        end
        s += _I1
    end

    # ============================================================
    # Phase 4b: find actual fcrit crossing near m0
    # ============================================================
    @inbounds Fbar_m0 = shmem_Fbar[m0]
    found_cross = false
    if Fbar_m0 >= fcrit
        # Search outward
        mm = m0
        while mm <= m
            @inbounds if shmem_Fbar[mm] < fcrit
                m0 = mm
                found_cross = true
                break
            end
            mm += _I1
        end
        if !found_cross
            mupp = m
        end
        if found_cross
            @inbounds rupp_m0 = shmem_rad[m0] + 2.0f0
            mupp = m0
            mm = m0 + _I1
            while mm <= m
                @inbounds if shmem_rad[mm] < rupp_m0
                    mupp = mm
                end
                mm += _I1
            end
        end
    else
        # Search inward from m0-1
        mstart = m0 > _I1 ? (m0 - _I1) : _I1
        mp = mstart
        while mp >= _I1
            @inbounds if shmem_Fbar[mp] >= fcrit
                found_cross = true
                break
            end
            m0 = mp
            mp -= _I1
        end
        if !found_cross
            _write_no_collapse_outputs!(
                RTHL_out, Fbarx_out, e_v_out, p_v_out,
                strain_final_out, eigs_out,
                Srb_out, Sbar_out, Sbar2_out,
                gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
                zvir_half_out, peak_id)
            return
        end
        @inbounds rupp_m0 = shmem_rad[m0] + 2.0f0
        mupp_new = m0
        mm = m0 + _I1
        while mm <= m
            @inbounds if shmem_rad[mm] < rupp_m0
                mupp_new = mm
            end
            mm += _I1
        end
        mupp = mupp_new
    end

    # ============================================================
    # Phase 5: strain at m0 + inward walk with zvir check
    # ============================================================
    mlow = m0
    @inbounds rlow = max(shmem_rad[m0] - 2.0f0, 0.0f0)
    if m0 > _I1
        m1 = m0 - _I1
        while m1 >= _I1
            @inbounds if shmem_rad[m1] > rlow
                mlow = m1
            end
            m1 -= _I1
        end
    end

    (E11_m0, E12_m0, E13_m0, E22_m0, E23_m0, E33_m0,
     gx_m0, gy_m0, gz_m0, gfx_m0, gfy_m0, gfz_m0) = _kernel_strain_f32(
        shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
        m0, mlow, mupp, wRnor, aRnor, hlatt_1, hlatt_2)
    (lam1_m0, lam2_m0, lam3_m0) = _eig3_symmetric_f32(
        E11_m0, E22_m0, E33_m0, E12_m0, E13_m0, E23_m0)
    Frho_m0_raw = lam1_m0 + lam2_m0 + lam3_m0
    e_v_m0 = Frho_m0_raw > 0.0f0 ? 0.5f0 * (lam3_m0 - lam1_m0) / Frho_m0_raw : 0.0f0
    p_v_m0 = Frho_m0_raw > 0.0f0 ? 0.5f0 * (lam3_m0 + lam1_m0 - 2.0f0 * lam2_m0) / Frho_m0_raw : 0.0f0
    Frhoh_m0 = Frho_m0_raw
    @inbounds Frho_m0_fbar = shmem_Fbar[m0]
    # CPU prototype uses -1.0f0 initial zvir1p (Fortran-compat uninitialized)
    zvir1p_m0 = -1.0f0

    # Running state (initialised at m0, updated during walk)
    Frhoh = Frhoh_m0
    Frho_val = Frho_m0_fbar
    e_v_curr = e_v_m0
    p_v_curr = p_v_m0
    E11_last = E11_m0; E12_last = E12_m0; E13_last = E13_m0
    E22_last = E22_m0; E23_last = E23_m0; E33_last = E33_m0
    # "Previous-frame" strain, updated when we finalise a non-collapsing shell
    E11_prev = E11_m0; E12_prev = E12_m0; E13_prev = E13_m0
    E22_prev = E22_m0; E23_prev = E23_m0; E33_prev = E33_m0
    # gradpk / gradpkf tracking (CPU analyse_peak_gpu:391-394)
    gx_last = gx_m0;  gy_last = gy_m0;  gz_last = gz_m0
    gfx_last = gfx_m0; gfy_last = gfy_m0; gfz_last = gfz_m0
    gx_prev = gx_m0;  gy_prev = gy_m0;  gz_prev = gz_m0
    gfx_prev = gfx_m0; gfy_prev = gfy_m0; gfz_prev = gfz_m0
    Frhpk = Frhoh_m0; Fnupk = Frho_m0_fbar; Fevpk = e_v_m0; Fpvpk = p_v_m0

    collapsed = false
    zvir1 = -1.0f0
    zvir1p = zvir1p_m0

    if zvir1p_m0 >= ZZon
        collapsed = true
    else
        zvir1 = zvir1p_m0
        if m0 == _I1
            _write_no_collapse_outputs!(
                RTHL_out, Fbarx_out, e_v_out, p_v_out,
                strain_final_out, eigs_out,
                Srb_out, Sbar_out, Sbar2_out,
                gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
                zvir_half_out, peak_id)
            return
        end

        mp = m0 - _I1
        while mp >= _I1
            @inbounds rupp_mp = shmem_rad[mp] + 2.0f0
            @inbounds rlow_mp = max(shmem_rad[mp] - 2.0f0, 0.0f0)

            @inbounds if shmem_rad[mupp] > rupp_mp
                mupp_new = mp
                m1 = mp + _I1
                while m1 <= mupp
                    @inbounds if shmem_rad[m1] < rupp_mp
                        mupp_new = m1
                    end
                    m1 += _I1
                end
                mupp = mupp_new
            end
            @inbounds if mlow > _I1 && shmem_rad[mlow - _I1] > rlow_mp
                mlownew = mlow
                m1 = mlow - _I1
                while m1 >= _I1
                    @inbounds if shmem_rad[m1] > rlow_mp
                        mlownew = m1
                    end
                    m1 -= _I1
                end
                mlow = mlownew
            end

            (E11_mp, E12_mp, E13_mp, E22_mp, E23_mp, E33_mp,
             gx_mp, gy_mp, gz_mp, gfx_mp, gfy_mp, gfz_mp) = _kernel_strain_f32(
                shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
                mp, mlow, mupp, wRnor, aRnor, hlatt_1, hlatt_2)
            (lam1_mp, lam2_mp, lam3_mp) = _eig3_symmetric_f32(
                E11_mp, E22_mp, E33_mp, E12_mp, E13_mp, E23_mp)
            Frho_mp_raw = lam1_mp + lam2_mp + lam3_mp
            e_v_mp = Frho_mp_raw > 0.0f0 ? 0.5f0 * (lam3_mp - lam1_mp) / Frho_mp_raw : 0.0f0
            p_v_mp = Frho_mp_raw > 0.0f0 ? 0.5f0 * (lam3_mp + lam1_mp - 2.0f0 * lam2_mp) / Frho_mp_raw : 0.0f0

            Frhoh_mp = Frho_mp_raw
            @inbounds Frho_mp_fbar = shmem_Fbar[mp]

            zvir1p_mp = -1.0f0
            if Frho_mp_fbar > 0.0f0
                poe_mp = e_v_mp < 1.0f-5 ? 0.0f0 : p_v_mp / e_v_mp
                zvir1p_mp = _interp3_trilinear_f32(
                    ct_table, log10(Frho_mp_fbar), e_v_mp, poe_mp,
                    ct_X1, ct_Y1, ct_Z1, ct_dxi, ct_dyi, ct_dzi,
                    ct_nx, ct_ny, ct_nz, ct_out_val)
            end

            if zvir1p_mp >= ZZon
                m0 = mp + _I1
                zvir1p = zvir1p_mp
                E11_last = E11_mp; E12_last = E12_mp; E13_last = E13_mp
                E22_last = E22_mp; E23_last = E23_mp; E33_last = E33_mp
                gx_last = gx_mp;   gy_last = gy_mp;   gz_last = gz_mp
                gfx_last = gfx_mp; gfy_last = gfy_mp; gfz_last = gfz_mp
                Frhoh = Frhoh_mp
                Frho_val = Frho_mp_fbar
                e_v_curr = e_v_mp
                p_v_curr = p_v_mp
                collapsed = true
                break
            end

            zvir1 = zvir1p_mp
            Frhpk = Frhoh_mp; Fnupk = Frho_mp_fbar; Fevpk = e_v_mp; Fpvpk = p_v_mp
            E11_prev = E11_mp; E12_prev = E12_mp; E13_prev = E13_mp
            E22_prev = E22_mp; E23_prev = E23_mp; E33_prev = E33_mp
            gx_prev = gx_mp;   gy_prev = gy_mp;   gz_prev = gz_mp
            gfx_prev = gfx_mp; gfy_prev = gfy_mp; gfz_prev = gfz_mp
            Frhoh = Frhoh_mp
            Frho_val = Frho_mp_fbar
            e_v_curr = e_v_mp
            p_v_curr = p_v_mp
            E11_last = E11_mp; E12_last = E12_mp; E13_last = E13_mp
            E22_last = E22_mp; E23_last = E23_mp; E33_last = E33_mp
            gx_last = gx_mp;   gy_last = gy_mp;   gz_last = gz_mp
            gfx_last = gfx_mp; gfy_last = gfy_mp; gfz_last = gfz_mp
            mp -= _I1
        end
    end

    if !collapsed
        _write_no_collapse_outputs!(
            RTHL_out, Fbarx_out, e_v_out, p_v_out,
            strain_final_out, eigs_out,
            Srb_out, Sbar_out, Sbar2_out,
            gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
            zvir_half_out, peak_id)
        return
    end

    # ============================================================
    # Phase 6: RTHL interpolation
    # ============================================================
    dZvir = zvir1p - zvir1
    RTHL = 0.0f0
    Fbarx = 0.0f0
    frac = 0.0f0
    if zvir1 > 0.0f0 && dZvir != 0.0f0
        @inbounds radmp = shmem_rad[m0 - _I1]
        @inbounds radm  = shmem_rad[m0]
        RTHL3 = radmp * radmp * radmp +
                (radm * radm * radm - radmp * radmp * radmp) * (zvir1p - ZZon) / dZvir
        RTHL = cbrt(RTHL3)
    else
        @inbounds RTHL = shmem_rad[m0 - _I1]
    end

    if RTHL <= 0.0f0
        _write_no_collapse_outputs!(
            RTHL_out, Fbarx_out, e_v_out, p_v_out,
            strain_final_out, eigs_out,
            Srb_out, Sbar_out, Sbar2_out,
            gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
            zvir_half_out, peak_id)
        return
    end

    @inbounds rad3p = shmem_rad[m0 - _I1] * shmem_rad[m0 - _I1] * shmem_rad[m0 - _I1]
    @inbounds rad3  = shmem_rad[m0] * shmem_rad[m0] * shmem_rad[m0]
    drad3 = rad3 - rad3p
    RTHL3 = RTHL * RTHL * RTHL
    @inbounds dFbar = shmem_Fbar[m0] - shmem_Fbar[m0 - _I1]

    if zvir1 > 0.0f0 && drad3 != 0.0f0
        frac = (RTHL3 - rad3p) / drad3
        @inbounds Fbarx = shmem_Fbar[m0 - _I1] + frac * dFbar
        Frhpk = Frhoh + frac * (Frhpk - Frhoh)
        Fnupk = Frho_val + frac * (Fnupk - Frho_val)
        Fevpk = e_v_curr + frac * (Fevpk - e_v_curr)
        Fpvpk = p_v_curr + frac * (Fpvpk - p_v_curr)
        # Strain interp: E_last + frac * (E_prev - E_last)
        E11_mat = E11_last + frac * (E11_prev - E11_last)
        E12_mat = E12_last + frac * (E12_prev - E12_last)
        E13_mat = E13_last + frac * (E13_prev - E13_last)
        E22_mat = E22_last + frac * (E22_prev - E22_last)
        E23_mat = E23_last + frac * (E23_prev - E23_last)
        E33_mat = E33_last + frac * (E33_prev - E33_last)
        # Gradpk / gradpkf interp (same frac as strain)
        gx_mat  = gx_last  + frac * (gx_prev  - gx_last)
        gy_mat  = gy_last  + frac * (gy_prev  - gy_last)
        gz_mat  = gz_last  + frac * (gz_prev  - gz_last)
        gfx_mat = gfx_last + frac * (gfx_prev - gfx_last)
        gfy_mat = gfy_last + frac * (gfy_prev - gfy_last)
        gfz_mat = gfz_last + frac * (gfz_prev - gfz_last)
    else
        @inbounds Fbarx = shmem_Fbar[m0 - _I1]
        Fnupk = Frho_val; Fevpk = e_v_curr; Fpvpk = p_v_curr
        E11_mat = E11_last; E12_mat = E12_last; E13_mat = E13_last
        E22_mat = E22_last; E23_mat = E23_last; E33_mat = E33_last
        gx_mat  = gx_prev;  gy_mat  = gy_prev;  gz_mat  = gz_prev
        gfx_mat = gfx_prev; gfy_mat = gfy_prev; gfz_mat = gfz_prev
    end

    # Normalise interpolated strain by its trace (tr = 3*Frho), rescale to Fbarx
    tr_E = E11_mat + E22_mat + E33_mat
    if tr_E != 0.0f0
        sc = Fbarx / tr_E
        E11n = E11_mat * sc; E12n = E12_mat * sc; E13n = E13_mat * sc
        E22n = E22_mat * sc; E23n = E23_mat * sc; E33n = E33_mat * sc
    else
        E11n = E11_mat; E12n = E12_mat; E13n = E13_mat
        E22n = E22_mat; E23n = E23_mat; E33n = E33_mat
    end
    (lamf1, lamf2, lamf3) = _eig3_symmetric_f32(E11n, E22n, E33n, E12n, E13n, E23n)
    Frhoc = lamf1 + lamf2 + lamf3
    e_v_f = Frhoc > 0.0f0 ? 0.5f0 * (lamf3 - lamf1) / Frhoc : 0.0f0
    p_v_f = Frhoc > 0.0f0 ? 0.5f0 * (lamf3 + lamf1 - 2.0f0 * lamf2) / Frhoc : 0.0f0

    # Final strain matrix = Ebar_last * (Fbarx / Frhoh), per CPU analyse_peak_gpu:509
    if Frhoh != 0.0f0
        sc2 = Fbarx / Frhoh
        E11f = E11_last * sc2; E12f = E12_last * sc2; E13f = E13_last * sc2
        E22f = E22_last * sc2; E23f = E23_last * sc2; E33f = E33_last * sc2
    else
        E11f = E11_last; E12f = E12_last; E13f = E13_last
        E22f = E22_last; E23f = E23_last; E33f = E33_last
    end

    # ============================================================
    # Phase 7: Sbar, Sbar2 (mean displacement over collapsed shells)
    #          and Srb (energy integral)
    # ============================================================
    # Sbar[L] = (1/nSbar) * Σ_{m1=1..m0-1} Sshell[L][m1]  (un-normalised sum
    # divided by total in-bounds cell count). Reads raw Sshell/S2shell from
    # global memory (we skipped staging these in shared mem to save bytes).
    Sbar_x  = 0.0f0; Sbar_y  = 0.0f0; Sbar_z  = 0.0f0
    Sbar2_x = 0.0f0; Sbar2_y = 0.0f0; Sbar2_z = 0.0f0
    nSbar = _I0
    m1 = _I1
    while m1 <= (m0 - _I1)
        @inbounds n_m1 = nshell_in[m1, peak_id]
        if n_m1 > _I0
            @inbounds Sbar_x += Sshell_in[_I1, m1, peak_id]
            @inbounds Sbar_y += Sshell_in[Int32(2), m1, peak_id]
            @inbounds Sbar_z += Sshell_in[Int32(3), m1, peak_id]
            @inbounds Sbar2_x += S2shell_in[_I1, m1, peak_id]
            @inbounds Sbar2_y += S2shell_in[Int32(2), m1, peak_id]
            @inbounds Sbar2_z += S2shell_in[Int32(3), m1, peak_id]
            nSbar += n_m1
        end
        m1 += _I1
    end
    if nSbar > _I0
        invN = 1.0f0 / Float32(nSbar)
        Sbar_x *= invN; Sbar_y *= invN; Sbar_z *= invN
        Sbar2_x *= invN; Sbar2_y *= invN; Sbar2_z *= invN
    end

    # Srb: energy integral ∫ Fbar(r) r^4 dr from 0 to RTHL, normalised by
    # Fbarx * RTHL^5. Trapezoidal in r^5 between successive shells.
    Srb = 0.0f0
    rad5p = 0.0f0
    if m0 > Int32(2)
        mp = Int32(2)
        while mp <= (m0 - _I1)
            @inbounds rad_mp = shmem_rad[mp]
            rad5 = rad_mp * rad_mp * rad_mp * rad_mp * rad_mp
            @inbounds Srb += 0.5f0 * (shmem_Fbar[mp - _I1] + shmem_Fbar[mp]) * (rad5 - rad5p)
            rad5p = rad5
            mp += _I1
        end
    end
    RTHL5 = RTHL3 * RTHL * RTHL  # RTHL^5
    if zvir1 > 0.0f0 && dZvir != 0.0f0
        @inbounds Srb += 0.5f0 * (shmem_Fbar[m0 - _I1] + Fbarx) * (RTHL5 - rad5p)
    end
    denom_Srb = Fbarx * RTHL5
    if denom_Srb > 0.0f0
        Srb /= denom_Srb
    end

    # ============================================================
    # Phase 7b: d2F Laplacian average over cells within RTHL
    # ============================================================
    # CPU (analyse_peak_gpu:549-564): iterate shells whose sqrt(r²) ≤ RTHL,
    # sum lapd over in-bounds cells, divide by total cell count.
    # We use the per-shell lapd sum from the gather kernel and the per-shell
    # in-bounds counts already materialised in nshell_in.
    d2F = 0.0f0
    nd2 = _I0
    s_lap = _I1
    while s_lap <= nshells
        @inbounds r2_s = shell_r2[s_lap]
        sqrt_r2 = sqrt(Float32(r2_s))
        if sqrt_r2 > RTHL
            break
        end
        @inbounds d2F += lapdshell_in[s_lap, peak_id]
        @inbounds nd2 += nshell_in[s_lap, peak_id]
        s_lap += _I1
    end
    if nd2 > _I0
        d2F /= Float32(nd2)
    end

    # ============================================================
    # Phase 8: zvir_half — formation redshift at RTHL/2
    # ============================================================
    # Matches RadialShell.analyse_peak:821-885. Evaluates the strain tensor
    # at 10 radii between RTHL/2 and RTHL (biased by (x)^(1/3) so we sample
    # denser near the smaller radii), runs eig3 + ct-table lookup at each,
    # and tracks max(zvir_rc - 1.0). Gated on RTHL ≥ 3.0 (hlatt = 1.0).
    # All inputs are already in shmem; no global-memory reads besides ct_table.
    zvir_half = -1.0f0
    if RTHL >= 3.0f0
        jj_int = Int32(1)
        while jj_int <= Int32(10)
            # rcur = RTHL * ((jj-1)/9 * 0.5 + 0.5)^(1/3)
            tfrac = (Float32(jj_int - _I1)) * (1.0f0 / 9.0f0) * 0.5f0 + 0.5f0
            rcur = RTHL * cbrt(tfrac)

            # ---- Find mupp_rc: smallest m1 in 2..nshells_max-1 with rad[m1] > rcur+2 ----
            rupp_rc = rcur + 2.0f0
            mupp_rc = _I1
            m1_u = Int32(2)
            while m1_u <= nshells_max - _I1
                @inbounds r_u = shmem_rad[m1_u]
                if r_u > rupp_rc
                    mupp_rc = m1_u
                    break
                elseif r_u == 0.0f0
                    mupp_rc = m1_u - _I1
                    break
                end
                mupp_rc = m1_u
                m1_u += _I1
            end

            # ---- Find mlow_rc: first m1 in 1..nshells_max-1 with rad[m1] > rlow_rc ----
            rlow_rc = rcur - 2.0f0
            if rlow_rc < 0.0f0
                rlow_rc = 0.0f0
            end
            mlow_rc = _I1
            m1_l = _I1
            while m1_l <= nshells_max - _I1
                @inbounds if shmem_rad[m1_l] > rlow_rc
                    mlow_rc = m1_l
                    break
                end
                m1_l += _I1
            end

            # ---- Find m0_rc: last m1 in 1..nshells_max-1 with rad[m1] ≤ rcur ----
            # (CPU: first m1 with rad[m1] > rcur, or m-1 if none)
            m0_rc = _I1
            m1_0 = _I1
            while m1_0 <= nshells_max - _I1
                @inbounds if shmem_rad[m1_0] > rcur
                    m0_rc = m1_0
                    break
                end
                m0_rc = m1_0
                m1_0 += _I1
            end

            # ---- kernel_strain at mp = m0_rc (hlatt_1 = hlatt_2 = 1.0) ----
            (E11_rc, E12_rc, E13_rc, E22_rc, E23_rc, E33_rc,
             _, _, _, _, _, _) = _kernel_strain_f32(
                shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
                m0_rc, mlow_rc, mupp_rc, wRnor, aRnor, 1.0f0, 1.0f0)

            # ---- Eigenvalues and normalised e_v, p_v ----
            (lam1_rc, lam2_rc, lam3_rc) = _eig3_symmetric_f32(
                E11_rc, E22_rc, E33_rc, E12_rc, E13_rc, E23_rc)
            Frho_rc = lam1_rc + lam2_rc + lam3_rc
            e_v_rc = 0.0f0; p_v_rc = 0.0f0
            if Frho_rc > 0.0f0
                e_v_rc = 0.5f0 * (lam3_rc - lam1_rc) / Frho_rc
                p_v_rc = 0.5f0 * (lam3_rc + lam1_rc - 2.0f0 * lam2_rc) / Frho_rc
            end

            @inbounds Frhoc_rc = shmem_Fbar[m0_rc]
            poe_rc = e_v_rc < 1.0f-5 ? 0.0f0 : p_v_rc / e_v_rc
            zvir_rc = -1.0f0
            # iflag_rc == 0 always on CPU (LAPACK) — skip the check on GPU too
            if Frhoc_rc > 0.0f0
                zvir_rc = _interp3_trilinear_f32(
                    ct_table, log10(Frhoc_rc), e_v_rc, poe_rc,
                    ct_X1, ct_Y1, ct_Z1, ct_dxi, ct_dyi, ct_dzi,
                    ct_nx, ct_ny, ct_nz, ct_out_val)
            end
            zrc_shifted = zvir_rc - 1.0f0
            if zrc_shifted > zvir_half
                zvir_half = zrc_shifted
            end

            jj_int += _I1
        end
    end

    # ============================================================
    # Write all outputs
    # ============================================================
    @inbounds RTHL_out[peak_id]  = RTHL
    @inbounds Fbarx_out[peak_id] = Fbarx
    @inbounds e_v_out[peak_id]   = e_v_f
    @inbounds p_v_out[peak_id]   = p_v_f
    @inbounds Srb_out[peak_id]   = Srb
    @inbounds d2F_out[peak_id]   = d2F
    @inbounds zvir_half_out[peak_id] = zvir_half
    @inbounds Sbar_out[_I1, peak_id]  = Sbar_x
    @inbounds Sbar_out[Int32(2), peak_id] = Sbar_y
    @inbounds Sbar_out[Int32(3), peak_id] = Sbar_z
    @inbounds Sbar2_out[_I1, peak_id] = Sbar2_x
    @inbounds Sbar2_out[Int32(2), peak_id] = Sbar2_y
    @inbounds Sbar2_out[Int32(3), peak_id] = Sbar2_z
    @inbounds gradpk_out[_I1, peak_id]  = gx_mat
    @inbounds gradpk_out[Int32(2), peak_id] = gy_mat
    @inbounds gradpk_out[Int32(3), peak_id] = gz_mat
    @inbounds gradpkf_out[_I1, peak_id]  = gfx_mat
    @inbounds gradpkf_out[Int32(2), peak_id] = gfy_mat
    @inbounds gradpkf_out[Int32(3), peak_id] = gfz_mat
    @inbounds gradpkrf_out[_I1, peak_id]  = gradpkrf_x
    @inbounds gradpkrf_out[Int32(2), peak_id] = gradpkrf_y
    @inbounds gradpkrf_out[Int32(3), peak_id] = gradpkrf_z
    @inbounds strain_final_out[1, 1, peak_id] = E11f
    @inbounds strain_final_out[1, 2, peak_id] = E12f
    @inbounds strain_final_out[1, 3, peak_id] = E13f
    @inbounds strain_final_out[2, 1, peak_id] = E12f
    @inbounds strain_final_out[2, 2, peak_id] = E22f
    @inbounds strain_final_out[2, 3, peak_id] = E23f
    @inbounds strain_final_out[3, 1, peak_id] = E13f
    @inbounds strain_final_out[3, 2, peak_id] = E23f
    @inbounds strain_final_out[3, 3, peak_id] = E33f
    # Eigenvalues of normalised (trace-matched) strain
    @inbounds eigs_out[1, peak_id] = lamf1
    @inbounds eigs_out[2, peak_id] = lamf2
    @inbounds eigs_out[3, peak_id] = lamf3

    # ============================================================
    # Phase 7c: mask update (side effect on shared mask grid)
    # ============================================================
    # Matches CPU analyse_peak_gpu:567-582. Iterate shells within RTHL,
    # write mask[cell] = 1 for in-bounds, non-nbuff cells.
    # Writes are idempotent (same value 1 from any peak). Races between
    # concurrent blocks writing same cell are harmless because the value
    # is always 1.
    if update_mask
        @inbounds ci = peaks_i[peak_id]
        @inbounds cj = peaks_j[peak_id]
        @inbounds ck = peaks_k[peak_id]
        s_mask = _I1
        while s_mask <= nshells
            @inbounds r2_s = shell_r2[s_mask]
            if sqrt(Float32(r2_s)) > RTHL
                break
            end
            @inbounds s0 = shell_start[s_mask]
            @inbounds nc = shell_count[s_mask]
            c_mask = _I1
            while c_mask <= nc
                idx = s0 + c_mask - _I1
                @inbounds di = off_di[idx]
                @inbounds dj = off_dj[idx]
                @inbounds dk = off_dk[idx]
                iv1 = ci + di; iv2 = cj + dj; iv3 = ck + dk
                if (_I1 <= iv1) & (iv1 <= n1) & (_I1 <= iv2) & (iv2 <= n2) & (_I1 <= iv3) & (iv3 <= n3)
                    if (iv1 > nbuff) & (iv1 <= n1 - nbuff) &
                       (iv2 > nbuff) & (iv2 <= n2 - nbuff) &
                       (iv3 > nbuff) & (iv3 <= n3 - nbuff)
                        @inbounds mask_out[iv1, iv2, iv3] = Int8(1)
                    end
                end
                c_mask += _I1
            end
            s_mask += _I1
        end
    end
    return
end

# ============================================================
# Fused gather + post-process kernel
# ============================================================
#
# Combines `_shell_gather_full_kernel!` and `_post_process_kernel!` into a
# single launch. Per-shell sums never touch global memory — the gather
# reductions stream directly into shared memory, and the sequential
# post-process on thread 0 reads from the same shared arrays.
#
# Savings vs the two-kernel path:
#   - Eliminates ~25 × nshells × npeaks × 4 B of per-tile global writes
#     (Fshell, nshell, Sshell, S2shell, Gshell, Gfshell, SRshell, lapdshell).
#   - At nshells=200, npeaks=5k, that's ~100 MB of per-batch global traffic
#     saved on both the write (gather) and the read (post-process) side.
#
# Shared-memory layout per block (one block = one peak), offsets in bytes:
#
#    0     shmem_rad     Float32 × 200    (800 B)
#  800     shmem_Fbar    Float32 × 200    (800 B)
# 1600     shmem_Gn      Float32 × 3 × 200 (2400 B)
# 4000     shmem_Gfn     Float32 × 3 × 200 (2400 B)
# 6400     shmem_SRn     Float32 × 3 × 3 × 200 (7200 B)
# 13600    shmem_Sshell  Float32 × 3 × 200 (2400 B)
# 16000    shmem_S2shell Float32 × 3 × 200 (2400 B)
# 18400    shmem_n       Int32   × 200    (800 B)
# 19200    shmem_lapd    Float32 × 200    (800 B)
# 20000    sdata_f       Float32 × 32     (128 B — block-reduce scratch)
# 20128    sdata_i       Int32   × 32     (128 B — block-reduce scratch)
# 20256    END (~20 KB)
#
# Fits twice in the 48 KB/block shmem budget on Ada / Ampere / Hopper.

# Shared-memory offsets (bytes)
const _FSA_OFF_RAD     = 0
const _FSA_OFF_FBAR    = _FSA_OFF_RAD     + sizeof(Float32) * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_GN      = _FSA_OFF_FBAR    + sizeof(Float32) * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_GFN     = _FSA_OFF_GN      + sizeof(Float32) * 3 * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_SRN     = _FSA_OFF_GFN     + sizeof(Float32) * 3 * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_SSH     = _FSA_OFF_SRN     + sizeof(Float32) * 9 * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_S2SH    = _FSA_OFF_SSH     + sizeof(Float32) * 3 * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_N       = _FSA_OFF_S2SH    + sizeof(Float32) * 3 * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_LAPD    = _FSA_OFF_N       + sizeof(Int32)   * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_SDATA_F = _FSA_OFF_LAPD    + sizeof(Float32) * Int(_MAX_SHELLS_GPU)
const _FSA_OFF_SDATA_I = _FSA_OFF_SDATA_F + sizeof(Float32) * 32
const _FSA_SHMEM_BYTES = _FSA_OFF_SDATA_I + sizeof(Int32)   * 32

function _fused_shell_analysis_kernel!(
    # Outputs (indexed by peak)
    RTHL_out, Fbarx_out, e_v_out, p_v_out,
    strain_final_out, eigs_out,
    Srb_out, Sbar_out, Sbar2_out,
    gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
    zvir_half_out, mask_out,
    # Field inputs
    delta, etax, etay, etaz, eta2x, eta2y, eta2z, lapd,
    # Peak centers
    peaks_i, peaks_j, peaks_k,
    # Shell offset tables
    off_di, off_dj, off_dk, shell_start, shell_count, shell_r2,
    # Collapse table
    ct_table,
    ct_X1::Float32, ct_Y1::Float32, ct_Z1::Float32,
    ct_dxi::Float32, ct_dyi::Float32, ct_dzi::Float32,
    ct_nx::Int32, ct_ny::Int32, ct_nz::Int32, ct_out_val::Float32,
    # SPH kernel LUT
    akk_tab,
    # Per-peak collapse state
    ZZon_pp, fcrit_pp, Rfclvi_r2::Float32,
    # Scaling constants
    wRnor::Float32, aRnor::Float32, hlatt_1::Float32, hlatt_2::Float32,
    ir2min::Int32,
    nshells::Int32, nshells_max::Int32,
    n1::Int32, n2::Int32, n3::Int32, nbuff::Int32, update_mask::Bool,
)
    peak_id = blockIdx().x
    tid     = threadIdx().x
    bdim    = blockDim().x

    # ---------- Shared-memory views (see layout comment above) ----------
    shmem_rad     = CuDynamicSharedArray(Float32, _MAX_SHELLS_GPU, _FSA_OFF_RAD)
    shmem_Fbar    = CuDynamicSharedArray(Float32, _MAX_SHELLS_GPU, _FSA_OFF_FBAR)
    shmem_Gn      = CuDynamicSharedArray(Float32, (Int32(3), _MAX_SHELLS_GPU), _FSA_OFF_GN)
    shmem_Gfn     = CuDynamicSharedArray(Float32, (Int32(3), _MAX_SHELLS_GPU), _FSA_OFF_GFN)
    shmem_SRn     = CuDynamicSharedArray(Float32, (Int32(3), Int32(3), _MAX_SHELLS_GPU), _FSA_OFF_SRN)
    shmem_Sshell  = CuDynamicSharedArray(Float32, (Int32(3), _MAX_SHELLS_GPU), _FSA_OFF_SSH)
    shmem_S2shell = CuDynamicSharedArray(Float32, (Int32(3), _MAX_SHELLS_GPU), _FSA_OFF_S2SH)
    shmem_n       = CuDynamicSharedArray(Int32,   _MAX_SHELLS_GPU, _FSA_OFF_N)
    shmem_lapd    = CuDynamicSharedArray(Float32, _MAX_SHELLS_GPU, _FSA_OFF_LAPD)
    sdata_f       = CuDynamicSharedArray(Float32, 32, _FSA_OFF_SDATA_F)
    sdata_i       = CuDynamicSharedArray(Int32,   32, _FSA_OFF_SDATA_I)

    @inbounds ci = peaks_i[peak_id]
    @inbounds cj = peaks_j[peak_id]
    @inbounds ck = peaks_k[peak_id]

    # =====================================================================
    # Phase A: Gather + inline normalisation.
    # Per-shell: threads cooperatively accumulate 23 Float32 + 1 Int32 sums,
    # block-reduce each, thread 0 stages normalised values into shmem.
    # dFbar_prev is kept in thread 0's registers across the shell loop.
    # =====================================================================
    dFbar_prev = 0.0f0
    s = _I1
    while s <= nshells
        @inbounds s0 = shell_start[s]
        @inbounds nc = shell_count[s]

        local_F   = 0.0f0
        local_Sx  = 0.0f0; local_Sy  = 0.0f0; local_Sz  = 0.0f0
        local_S2x = 0.0f0; local_S2y = 0.0f0; local_S2z = 0.0f0
        local_Gx  = 0.0f0; local_Gy  = 0.0f0; local_Gz  = 0.0f0
        local_Gfx = 0.0f0; local_Gfy = 0.0f0; local_Gfz = 0.0f0
        local_SRxx = 0.0f0; local_SRxy = 0.0f0; local_SRxz = 0.0f0
        local_SRyx = 0.0f0; local_SRyy = 0.0f0; local_SRyz = 0.0f0
        local_SRzx = 0.0f0; local_SRzy = 0.0f0; local_SRzz = 0.0f0
        local_lapd = 0.0f0
        local_n = _I0

        c = tid
        while c <= nc
            idx = s0 + c - _I1
            @inbounds di = off_di[idx]
            @inbounds dj = off_dj[idx]
            @inbounds dk = off_dk[idx]

            iv1 = ci + di
            iv2 = cj + dj
            iv3 = ck + dk

            if (_I1 <= iv1) & (iv1 <= n1) & (_I1 <= iv2) & (iv2 <= n2) & (_I1 <= iv3) & (iv3 <= n3)
                @inbounds d_val  = delta[iv1, iv2, iv3]
                @inbounds df_val = delta[iv3, iv1, iv2]  # transpose for Gf
                @inbounds ex  = etax[iv1, iv2, iv3]
                @inbounds ey  = etay[iv1, iv2, iv3]
                @inbounds ez  = etaz[iv1, iv2, iv3]
                @inbounds e2x = eta2x[iv1, iv2, iv3]
                @inbounds e2y = eta2y[iv1, iv2, iv3]
                @inbounds e2z = eta2z[iv1, iv2, iv3]

                local_F += d_val
                local_Sx  += ex;  local_Sy  += ey;  local_Sz  += ez
                local_S2x += e2x; local_S2y += e2y; local_S2z += e2z

                fdi = Float32(di); fdj = Float32(dj); fdk = Float32(dk)
                local_Gx  += d_val  * fdi; local_Gy  += d_val  * fdj; local_Gz  += d_val  * fdk
                local_Gfx += df_val * fdi; local_Gfy += df_val * fdj; local_Gfz += df_val * fdk

                local_SRxx += ex * fdi; local_SRxy += ex * fdj; local_SRxz += ex * fdk
                local_SRyx += ey * fdi; local_SRyy += ey * fdj; local_SRyz += ey * fdk
                local_SRzx += ez * fdi; local_SRzy += ez * fdj; local_SRzz += ez * fdk

                @inbounds local_lapd += lapd[iv1, iv2, iv3]
                local_n += _I1
            end
            c += bdim
        end

        # 22 Float32 + 1 Int32 block reductions per shell.
        F_tot   = block_reduce_sum_f32(local_F,   sdata_f); sync_threads()
        Sx_tot  = block_reduce_sum_f32(local_Sx,  sdata_f); sync_threads()
        Sy_tot  = block_reduce_sum_f32(local_Sy,  sdata_f); sync_threads()
        Sz_tot  = block_reduce_sum_f32(local_Sz,  sdata_f); sync_threads()
        S2x_tot = block_reduce_sum_f32(local_S2x, sdata_f); sync_threads()
        S2y_tot = block_reduce_sum_f32(local_S2y, sdata_f); sync_threads()
        S2z_tot = block_reduce_sum_f32(local_S2z, sdata_f); sync_threads()
        Gx_tot  = block_reduce_sum_f32(local_Gx,  sdata_f); sync_threads()
        Gy_tot  = block_reduce_sum_f32(local_Gy,  sdata_f); sync_threads()
        Gz_tot  = block_reduce_sum_f32(local_Gz,  sdata_f); sync_threads()
        Gfx_tot = block_reduce_sum_f32(local_Gfx, sdata_f); sync_threads()
        Gfy_tot = block_reduce_sum_f32(local_Gfy, sdata_f); sync_threads()
        Gfz_tot = block_reduce_sum_f32(local_Gfz, sdata_f); sync_threads()
        SRxx_tot = block_reduce_sum_f32(local_SRxx, sdata_f); sync_threads()
        SRxy_tot = block_reduce_sum_f32(local_SRxy, sdata_f); sync_threads()
        SRxz_tot = block_reduce_sum_f32(local_SRxz, sdata_f); sync_threads()
        SRyx_tot = block_reduce_sum_f32(local_SRyx, sdata_f); sync_threads()
        SRyy_tot = block_reduce_sum_f32(local_SRyy, sdata_f); sync_threads()
        SRyz_tot = block_reduce_sum_f32(local_SRyz, sdata_f); sync_threads()
        SRzx_tot = block_reduce_sum_f32(local_SRzx, sdata_f); sync_threads()
        SRzy_tot = block_reduce_sum_f32(local_SRzy, sdata_f); sync_threads()
        SRzz_tot = block_reduce_sum_f32(local_SRzz, sdata_f); sync_threads()
        lapd_tot = block_reduce_sum_f32(local_lapd, sdata_f); sync_threads()
        n_tot    = block_reduce_sum_i32(local_n,   sdata_i); sync_threads()

        if tid == _I1
            @inbounds r2_val = shell_r2[s]
            if s == _I1
                # Shell 1 (center cell, rad=0): special-case. G/Gf/SR are
                # zero, Fbar[1] = dFbar_1 = F/n.
                dFbar_1 = n_tot > _I0 ? F_tot / Float32(n_tot) : 0.0f0
                @inbounds shmem_rad[1]  = 0.0f0
                @inbounds shmem_Fbar[1] = dFbar_1
                @inbounds shmem_n[1]    = n_tot
                for L in _I1:Int32(3)
                    @inbounds shmem_Gn[L, 1]  = 0.0f0
                    @inbounds shmem_Gfn[L, 1] = 0.0f0
                    for K in _I1:Int32(3)
                        @inbounds shmem_SRn[L, K, 1] = 0.0f0
                    end
                end
                @inbounds shmem_Sshell[_I1, 1] = Sx_tot
                @inbounds shmem_Sshell[Int32(2), 1] = Sy_tot
                @inbounds shmem_Sshell[Int32(3), 1] = Sz_tot
                @inbounds shmem_S2shell[_I1, 1] = S2x_tot
                @inbounds shmem_S2shell[Int32(2), 1] = S2y_tot
                @inbounds shmem_S2shell[Int32(3), 1] = S2z_tot
                @inbounds shmem_lapd[1] = lapd_tot
                dFbar_prev = dFbar_1
            else
                rad_s = sqrt(Float32(r2_val))
                rad_s_inv = 1.0f0 / rad_s
                dFbar_s = n_tot > _I0 ? F_tot / Float32(n_tot) : 0.0f0
                # Trapezoidal integration for Fbar(r) (matches post-process)
                @inbounds rad_sp = shmem_rad[s - _I1]
                rad3p = rad_sp * rad_sp * rad_sp
                rad3  = rad_s * rad_s * rad_s
                @inbounds Fbar_prev = shmem_Fbar[s - _I1]
                Fbar_s = (rad3p * Fbar_prev + 0.5f0 * (dFbar_prev + dFbar_s) * (rad3 - rad3p)) / rad3

                @inbounds shmem_rad[s]  = rad_s
                @inbounds shmem_Fbar[s] = Fbar_s
                @inbounds shmem_n[s]    = n_tot
                # Normalise G/Gf/SR by rad_s inline
                @inbounds shmem_Gn[_I1, s]       = Gx_tot  * rad_s_inv
                @inbounds shmem_Gn[Int32(2), s]  = Gy_tot  * rad_s_inv
                @inbounds shmem_Gn[Int32(3), s]  = Gz_tot  * rad_s_inv
                @inbounds shmem_Gfn[_I1, s]      = Gfx_tot * rad_s_inv
                @inbounds shmem_Gfn[Int32(2), s] = Gfy_tot * rad_s_inv
                @inbounds shmem_Gfn[Int32(3), s] = Gfz_tot * rad_s_inv
                @inbounds shmem_SRn[_I1, _I1, s]      = SRxx_tot * rad_s_inv
                @inbounds shmem_SRn[_I1, Int32(2), s] = SRxy_tot * rad_s_inv
                @inbounds shmem_SRn[_I1, Int32(3), s] = SRxz_tot * rad_s_inv
                @inbounds shmem_SRn[Int32(2), _I1, s]      = SRyx_tot * rad_s_inv
                @inbounds shmem_SRn[Int32(2), Int32(2), s] = SRyy_tot * rad_s_inv
                @inbounds shmem_SRn[Int32(2), Int32(3), s] = SRyz_tot * rad_s_inv
                @inbounds shmem_SRn[Int32(3), _I1, s]      = SRzx_tot * rad_s_inv
                @inbounds shmem_SRn[Int32(3), Int32(2), s] = SRzy_tot * rad_s_inv
                @inbounds shmem_SRn[Int32(3), Int32(3), s] = SRzz_tot * rad_s_inv
                @inbounds shmem_Sshell[_I1, s]       = Sx_tot
                @inbounds shmem_Sshell[Int32(2), s]  = Sy_tot
                @inbounds shmem_Sshell[Int32(3), s]  = Sz_tot
                @inbounds shmem_S2shell[_I1, s]      = S2x_tot
                @inbounds shmem_S2shell[Int32(2), s] = S2y_tot
                @inbounds shmem_S2shell[Int32(3), s] = S2z_tot
                @inbounds shmem_lapd[s] = lapd_tot

                dFbar_prev = dFbar_s
            end
        end
        sync_threads()
        s += _I1
    end

    # =====================================================================
    # Phase B: sequential post-process on thread 0 (reads from shmem).
    # Matches `_post_process_kernel!` phases 3–8 verbatim, with global
    # reads replaced by shmem_* reads. Inlining here avoids the duplicate
    # shmem-staging loop that kernel would otherwise perform.
    # =====================================================================
    if tid != _I1
        return
    end

    @inbounds ZZon  = ZZon_pp[peak_id]
    @inbounds fcrit = fcrit_pp[peak_id]
    m = nshells_max

    # ----- Phase 3: gradient at filter scale (gradpkrf) -----
    mrf = _I1
    rmrf = 10.0f0
    mrfi = _I1
    s_mrf = Int32(2)
    while s_mrf <= nshells_max
        @inbounds r2_s = shell_r2[s_mrf]
        dist_mrf = abs(Float32(r2_s) - Rfclvi_r2)
        if dist_mrf < rmrf
            mrf = mrfi
            rmrf = dist_mrf
        end
        if s_mrf < nshells_max
            mrfi += _I1
        end
        s_mrf += _I1
    end
    @inbounds rlow_rf = max(shmem_rad[mrf] - 2.0f0, 0.0f0)
    mlow_rf = mrf
    if mrf > _I1
        m1_rf = mrf - _I1
        while m1_rf >= _I1
            @inbounds if shmem_rad[m1_rf] > rlow_rf
                mlow_rf = m1_rf
            end
            m1_rf -= _I1
        end
    end
    @inbounds rupp_rf = shmem_rad[mrf] + 2.0f0
    mupp_rf = mrf
    m1_rf = mrf + _I1
    while m1_rf <= m
        @inbounds if shmem_rad[m1_rf] < rupp_rf
            mupp_rf = m1_rf
        end
        m1_rf += _I1
    end
    (_, _, _, _, _, _,
     gradpkrf_x, gradpkrf_y, gradpkrf_z, _, _, _) = _kernel_strain_f32(
        shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
        mrf, mlow_rf, mupp_rf, wRnor, aRnor, hlatt_1, hlatt_2)

    # ----- Phase 4a: m0/mupp initialisation -----
    m0 = _I1
    ifcrit = _I1
    rupp_init = sqrt(Float32(ir2min)) + 2.0f0
    mupp = _I1
    s = Int32(2)
    while s <= nshells_max
        @inbounds r2 = shell_r2[s]
        if ifcrit == _I1
            m0 = s
            if Float32(r2) > Float32(ir2min)
                @inbounds rupp_init = shmem_rad[m0] + 2.0f0
                ifcrit = _I0
                mupp = s
            end
        else
            @inbounds if shmem_rad[s] < rupp_init
                mupp = s
            end
        end
        s += _I1
    end

    # ----- Phase 4b: fcrit crossing -----
    @inbounds Fbar_m0 = shmem_Fbar[m0]
    found_cross = false
    if Fbar_m0 >= fcrit
        mm = m0
        while mm <= m
            @inbounds if shmem_Fbar[mm] < fcrit
                m0 = mm
                found_cross = true
                break
            end
            mm += _I1
        end
        if !found_cross
            mupp = m
        end
        if found_cross
            @inbounds rupp_m0 = shmem_rad[m0] + 2.0f0
            mupp = m0
            mm = m0 + _I1
            while mm <= m
                @inbounds if shmem_rad[mm] < rupp_m0
                    mupp = mm
                end
                mm += _I1
            end
        end
    else
        mstart = m0 > _I1 ? (m0 - _I1) : _I1
        mp = mstart
        while mp >= _I1
            @inbounds if shmem_Fbar[mp] >= fcrit
                found_cross = true
                break
            end
            m0 = mp
            mp -= _I1
        end
        if !found_cross
            _write_no_collapse_outputs!(
                RTHL_out, Fbarx_out, e_v_out, p_v_out,
                strain_final_out, eigs_out,
                Srb_out, Sbar_out, Sbar2_out,
                gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
                zvir_half_out, peak_id)
            return
        end
        @inbounds rupp_m0 = shmem_rad[m0] + 2.0f0
        mupp_new = m0
        mm = m0 + _I1
        while mm <= m
            @inbounds if shmem_rad[mm] < rupp_m0
                mupp_new = mm
            end
            mm += _I1
        end
        mupp = mupp_new
    end

    # ----- Phase 5: strain at m0 + inward walk -----
    mlow = m0
    @inbounds rlow = max(shmem_rad[m0] - 2.0f0, 0.0f0)
    if m0 > _I1
        m1 = m0 - _I1
        while m1 >= _I1
            @inbounds if shmem_rad[m1] > rlow
                mlow = m1
            end
            m1 -= _I1
        end
    end

    (E11_m0, E12_m0, E13_m0, E22_m0, E23_m0, E33_m0,
     gx_m0, gy_m0, gz_m0, gfx_m0, gfy_m0, gfz_m0) = _kernel_strain_f32(
        shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
        m0, mlow, mupp, wRnor, aRnor, hlatt_1, hlatt_2)
    (lam1_m0, lam2_m0, lam3_m0) = _eig3_symmetric_f32(
        E11_m0, E22_m0, E33_m0, E12_m0, E13_m0, E23_m0)
    Frho_m0_raw = lam1_m0 + lam2_m0 + lam3_m0
    e_v_m0 = Frho_m0_raw > 0.0f0 ? 0.5f0 * (lam3_m0 - lam1_m0) / Frho_m0_raw : 0.0f0
    p_v_m0 = Frho_m0_raw > 0.0f0 ? 0.5f0 * (lam3_m0 + lam1_m0 - 2.0f0 * lam2_m0) / Frho_m0_raw : 0.0f0
    Frhoh_m0 = Frho_m0_raw
    @inbounds Frho_m0_fbar = shmem_Fbar[m0]
    zvir1p_m0 = -1.0f0

    Frhoh = Frhoh_m0
    Frho_val = Frho_m0_fbar
    e_v_curr = e_v_m0
    p_v_curr = p_v_m0
    E11_last = E11_m0; E12_last = E12_m0; E13_last = E13_m0
    E22_last = E22_m0; E23_last = E23_m0; E33_last = E33_m0
    E11_prev = E11_m0; E12_prev = E12_m0; E13_prev = E13_m0
    E22_prev = E22_m0; E23_prev = E23_m0; E33_prev = E33_m0
    gx_last = gx_m0;  gy_last = gy_m0;  gz_last = gz_m0
    gfx_last = gfx_m0; gfy_last = gfy_m0; gfz_last = gfz_m0
    gx_prev = gx_m0;  gy_prev = gy_m0;  gz_prev = gz_m0
    gfx_prev = gfx_m0; gfy_prev = gfy_m0; gfz_prev = gfz_m0
    Frhpk = Frhoh_m0; Fnupk = Frho_m0_fbar; Fevpk = e_v_m0; Fpvpk = p_v_m0

    collapsed = false
    zvir1 = -1.0f0
    zvir1p = zvir1p_m0

    if zvir1p_m0 >= ZZon
        collapsed = true
    else
        zvir1 = zvir1p_m0
        if m0 == _I1
            _write_no_collapse_outputs!(
                RTHL_out, Fbarx_out, e_v_out, p_v_out,
                strain_final_out, eigs_out,
                Srb_out, Sbar_out, Sbar2_out,
                gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
                zvir_half_out, peak_id)
            return
        end

        mp = m0 - _I1
        while mp >= _I1
            @inbounds rupp_mp = shmem_rad[mp] + 2.0f0
            @inbounds rlow_mp = max(shmem_rad[mp] - 2.0f0, 0.0f0)

            @inbounds if shmem_rad[mupp] > rupp_mp
                mupp_new = mp
                m1 = mp + _I1
                while m1 <= mupp
                    @inbounds if shmem_rad[m1] < rupp_mp
                        mupp_new = m1
                    end
                    m1 += _I1
                end
                mupp = mupp_new
            end
            @inbounds if mlow > _I1 && shmem_rad[mlow - _I1] > rlow_mp
                mlownew = mlow
                m1 = mlow - _I1
                while m1 >= _I1
                    @inbounds if shmem_rad[m1] > rlow_mp
                        mlownew = m1
                    end
                    m1 -= _I1
                end
                mlow = mlownew
            end

            (E11_mp, E12_mp, E13_mp, E22_mp, E23_mp, E33_mp,
             gx_mp, gy_mp, gz_mp, gfx_mp, gfy_mp, gfz_mp) = _kernel_strain_f32(
                shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
                mp, mlow, mupp, wRnor, aRnor, hlatt_1, hlatt_2)
            (lam1_mp, lam2_mp, lam3_mp) = _eig3_symmetric_f32(
                E11_mp, E22_mp, E33_mp, E12_mp, E13_mp, E23_mp)
            Frho_mp_raw = lam1_mp + lam2_mp + lam3_mp
            e_v_mp = Frho_mp_raw > 0.0f0 ? 0.5f0 * (lam3_mp - lam1_mp) / Frho_mp_raw : 0.0f0
            p_v_mp = Frho_mp_raw > 0.0f0 ? 0.5f0 * (lam3_mp + lam1_mp - 2.0f0 * lam2_mp) / Frho_mp_raw : 0.0f0

            Frhoh_mp = Frho_mp_raw
            @inbounds Frho_mp_fbar = shmem_Fbar[mp]

            zvir1p_mp = -1.0f0
            if Frho_mp_fbar > 0.0f0
                poe_mp = e_v_mp < 1.0f-5 ? 0.0f0 : p_v_mp / e_v_mp
                zvir1p_mp = _interp3_trilinear_f32(
                    ct_table, log10(Frho_mp_fbar), e_v_mp, poe_mp,
                    ct_X1, ct_Y1, ct_Z1, ct_dxi, ct_dyi, ct_dzi,
                    ct_nx, ct_ny, ct_nz, ct_out_val)
            end

            if zvir1p_mp >= ZZon
                m0 = mp + _I1
                zvir1p = zvir1p_mp
                E11_last = E11_mp; E12_last = E12_mp; E13_last = E13_mp
                E22_last = E22_mp; E23_last = E23_mp; E33_last = E33_mp
                gx_last = gx_mp;   gy_last = gy_mp;   gz_last = gz_mp
                gfx_last = gfx_mp; gfy_last = gfy_mp; gfz_last = gfz_mp
                Frhoh = Frhoh_mp
                Frho_val = Frho_mp_fbar
                e_v_curr = e_v_mp
                p_v_curr = p_v_mp
                collapsed = true
                break
            end

            zvir1 = zvir1p_mp
            Frhpk = Frhoh_mp; Fnupk = Frho_mp_fbar; Fevpk = e_v_mp; Fpvpk = p_v_mp
            E11_prev = E11_mp; E12_prev = E12_mp; E13_prev = E13_mp
            E22_prev = E22_mp; E23_prev = E23_mp; E33_prev = E33_mp
            gx_prev = gx_mp;   gy_prev = gy_mp;   gz_prev = gz_mp
            gfx_prev = gfx_mp; gfy_prev = gfy_mp; gfz_prev = gfz_mp
            Frhoh = Frhoh_mp
            Frho_val = Frho_mp_fbar
            e_v_curr = e_v_mp
            p_v_curr = p_v_mp
            E11_last = E11_mp; E12_last = E12_mp; E13_last = E13_mp
            E22_last = E22_mp; E23_last = E23_mp; E33_last = E33_mp
            gx_last = gx_mp;   gy_last = gy_mp;   gz_last = gz_mp
            gfx_last = gfx_mp; gfy_last = gfy_mp; gfz_last = gfz_mp
            mp -= _I1
        end
    end

    if !collapsed
        _write_no_collapse_outputs!(
            RTHL_out, Fbarx_out, e_v_out, p_v_out,
            strain_final_out, eigs_out,
            Srb_out, Sbar_out, Sbar2_out,
            gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
            zvir_half_out, peak_id)
        return
    end

    # ----- Phase 6: RTHL interpolation -----
    dZvir = zvir1p - zvir1
    RTHL = 0.0f0
    Fbarx = 0.0f0
    frac = 0.0f0
    if zvir1 > 0.0f0 && dZvir != 0.0f0
        @inbounds radmp = shmem_rad[m0 - _I1]
        @inbounds radm  = shmem_rad[m0]
        RTHL3 = radmp * radmp * radmp +
                (radm * radm * radm - radmp * radmp * radmp) * (zvir1p - ZZon) / dZvir
        RTHL = cbrt(RTHL3)
    else
        @inbounds RTHL = shmem_rad[m0 - _I1]
    end

    if RTHL <= 0.0f0
        _write_no_collapse_outputs!(
            RTHL_out, Fbarx_out, e_v_out, p_v_out,
            strain_final_out, eigs_out,
            Srb_out, Sbar_out, Sbar2_out,
            gradpk_out, gradpkf_out, gradpkrf_out, d2F_out,
            zvir_half_out, peak_id)
        return
    end

    @inbounds rad3p = shmem_rad[m0 - _I1] * shmem_rad[m0 - _I1] * shmem_rad[m0 - _I1]
    @inbounds rad3  = shmem_rad[m0] * shmem_rad[m0] * shmem_rad[m0]
    drad3 = rad3 - rad3p
    RTHL3 = RTHL * RTHL * RTHL
    @inbounds dFbar = shmem_Fbar[m0] - shmem_Fbar[m0 - _I1]

    if zvir1 > 0.0f0 && drad3 != 0.0f0
        frac = (RTHL3 - rad3p) / drad3
        @inbounds Fbarx = shmem_Fbar[m0 - _I1] + frac * dFbar
        Frhpk = Frhoh + frac * (Frhpk - Frhoh)
        Fnupk = Frho_val + frac * (Fnupk - Frho_val)
        Fevpk = e_v_curr + frac * (Fevpk - e_v_curr)
        Fpvpk = p_v_curr + frac * (Fpvpk - p_v_curr)
        E11_mat = E11_last + frac * (E11_prev - E11_last)
        E12_mat = E12_last + frac * (E12_prev - E12_last)
        E13_mat = E13_last + frac * (E13_prev - E13_last)
        E22_mat = E22_last + frac * (E22_prev - E22_last)
        E23_mat = E23_last + frac * (E23_prev - E23_last)
        E33_mat = E33_last + frac * (E33_prev - E33_last)
        gx_mat  = gx_last  + frac * (gx_prev  - gx_last)
        gy_mat  = gy_last  + frac * (gy_prev  - gy_last)
        gz_mat  = gz_last  + frac * (gz_prev  - gz_last)
        gfx_mat = gfx_last + frac * (gfx_prev - gfx_last)
        gfy_mat = gfy_last + frac * (gfy_prev - gfy_last)
        gfz_mat = gfz_last + frac * (gfz_prev - gfz_last)
    else
        @inbounds Fbarx = shmem_Fbar[m0 - _I1]
        Fnupk = Frho_val; Fevpk = e_v_curr; Fpvpk = p_v_curr
        E11_mat = E11_last; E12_mat = E12_last; E13_mat = E13_last
        E22_mat = E22_last; E23_mat = E23_last; E33_mat = E33_last
        gx_mat  = gx_prev;  gy_mat  = gy_prev;  gz_mat  = gz_prev
        gfx_mat = gfx_prev; gfy_mat = gfy_prev; gfz_mat = gfz_prev
    end

    tr_E = E11_mat + E22_mat + E33_mat
    if tr_E != 0.0f0
        sc = Fbarx / tr_E
        E11n = E11_mat * sc; E12n = E12_mat * sc; E13n = E13_mat * sc
        E22n = E22_mat * sc; E23n = E23_mat * sc; E33n = E33_mat * sc
    else
        E11n = E11_mat; E12n = E12_mat; E13n = E13_mat
        E22n = E22_mat; E23n = E23_mat; E33n = E33_mat
    end
    (lamf1, lamf2, lamf3) = _eig3_symmetric_f32(E11n, E22n, E33n, E12n, E13n, E23n)
    Frhoc = lamf1 + lamf2 + lamf3
    e_v_f = Frhoc > 0.0f0 ? 0.5f0 * (lamf3 - lamf1) / Frhoc : 0.0f0
    p_v_f = Frhoc > 0.0f0 ? 0.5f0 * (lamf3 + lamf1 - 2.0f0 * lamf2) / Frhoc : 0.0f0

    if Frhoh != 0.0f0
        sc2 = Fbarx / Frhoh
        E11f = E11_last * sc2; E12f = E12_last * sc2; E13f = E13_last * sc2
        E22f = E22_last * sc2; E23f = E23_last * sc2; E33f = E33_last * sc2
    else
        E11f = E11_last; E12f = E12_last; E13f = E13_last
        E22f = E22_last; E23f = E23_last; E33f = E33_last
    end

    # ----- Phase 7: Sbar, Sbar2, Srb -----
    Sbar_x  = 0.0f0; Sbar_y  = 0.0f0; Sbar_z  = 0.0f0
    Sbar2_x = 0.0f0; Sbar2_y = 0.0f0; Sbar2_z = 0.0f0
    nSbar = _I0
    m1 = _I1
    while m1 <= (m0 - _I1)
        @inbounds n_m1 = shmem_n[m1]
        if n_m1 > _I0
            @inbounds Sbar_x += shmem_Sshell[_I1, m1]
            @inbounds Sbar_y += shmem_Sshell[Int32(2), m1]
            @inbounds Sbar_z += shmem_Sshell[Int32(3), m1]
            @inbounds Sbar2_x += shmem_S2shell[_I1, m1]
            @inbounds Sbar2_y += shmem_S2shell[Int32(2), m1]
            @inbounds Sbar2_z += shmem_S2shell[Int32(3), m1]
            nSbar += n_m1
        end
        m1 += _I1
    end
    if nSbar > _I0
        invN = 1.0f0 / Float32(nSbar)
        Sbar_x *= invN; Sbar_y *= invN; Sbar_z *= invN
        Sbar2_x *= invN; Sbar2_y *= invN; Sbar2_z *= invN
    end

    Srb = 0.0f0
    rad5p = 0.0f0
    if m0 > Int32(2)
        mp = Int32(2)
        while mp <= (m0 - _I1)
            @inbounds rad_mp = shmem_rad[mp]
            rad5 = rad_mp * rad_mp * rad_mp * rad_mp * rad_mp
            @inbounds Srb += 0.5f0 * (shmem_Fbar[mp - _I1] + shmem_Fbar[mp]) * (rad5 - rad5p)
            rad5p = rad5
            mp += _I1
        end
    end
    RTHL5 = RTHL3 * RTHL * RTHL
    if zvir1 > 0.0f0 && dZvir != 0.0f0
        @inbounds Srb += 0.5f0 * (shmem_Fbar[m0 - _I1] + Fbarx) * (RTHL5 - rad5p)
    end
    denom_Srb = Fbarx * RTHL5
    if denom_Srb > 0.0f0
        Srb /= denom_Srb
    end

    # ----- Phase 7b: d2F Laplacian average within RTHL -----
    d2F = 0.0f0
    nd2 = _I0
    s_lap = _I1
    while s_lap <= nshells
        @inbounds r2_s = shell_r2[s_lap]
        sqrt_r2 = sqrt(Float32(r2_s))
        if sqrt_r2 > RTHL
            break
        end
        @inbounds d2F += shmem_lapd[s_lap]
        @inbounds nd2 += shmem_n[s_lap]
        s_lap += _I1
    end
    if nd2 > _I0
        d2F /= Float32(nd2)
    end

    # ----- Phase 8: zvir_half -----
    zvir_half = -1.0f0
    if RTHL >= 3.0f0
        jj_int = Int32(1)
        while jj_int <= Int32(10)
            tfrac = (Float32(jj_int - _I1)) * (1.0f0 / 9.0f0) * 0.5f0 + 0.5f0
            rcur = RTHL * cbrt(tfrac)

            rupp_rc = rcur + 2.0f0
            mupp_rc = _I1
            m1_u = Int32(2)
            while m1_u <= nshells_max - _I1
                @inbounds r_u = shmem_rad[m1_u]
                if r_u > rupp_rc
                    mupp_rc = m1_u
                    break
                elseif r_u == 0.0f0
                    mupp_rc = m1_u - _I1
                    break
                end
                mupp_rc = m1_u
                m1_u += _I1
            end

            rlow_rc = rcur - 2.0f0
            if rlow_rc < 0.0f0
                rlow_rc = 0.0f0
            end
            mlow_rc = _I1
            m1_l = _I1
            while m1_l <= nshells_max - _I1
                @inbounds if shmem_rad[m1_l] > rlow_rc
                    mlow_rc = m1_l
                    break
                end
                m1_l += _I1
            end

            m0_rc = _I1
            m1_0 = _I1
            while m1_0 <= nshells_max - _I1
                @inbounds if shmem_rad[m1_0] > rcur
                    m0_rc = m1_0
                    break
                end
                m0_rc = m1_0
                m1_0 += _I1
            end

            (E11_rc, E12_rc, E13_rc, E22_rc, E23_rc, E33_rc,
             _, _, _, _, _, _) = _kernel_strain_f32(
                shmem_rad, shmem_Gn, shmem_Gfn, shmem_SRn, akk_tab,
                m0_rc, mlow_rc, mupp_rc, wRnor, aRnor, 1.0f0, 1.0f0)

            (lam1_rc, lam2_rc, lam3_rc) = _eig3_symmetric_f32(
                E11_rc, E22_rc, E33_rc, E12_rc, E13_rc, E23_rc)
            Frho_rc = lam1_rc + lam2_rc + lam3_rc
            e_v_rc = 0.0f0; p_v_rc = 0.0f0
            if Frho_rc > 0.0f0
                e_v_rc = 0.5f0 * (lam3_rc - lam1_rc) / Frho_rc
                p_v_rc = 0.5f0 * (lam3_rc + lam1_rc - 2.0f0 * lam2_rc) / Frho_rc
            end

            @inbounds Frhoc_rc = shmem_Fbar[m0_rc]
            poe_rc = e_v_rc < 1.0f-5 ? 0.0f0 : p_v_rc / e_v_rc
            zvir_rc = -1.0f0
            if Frhoc_rc > 0.0f0
                zvir_rc = _interp3_trilinear_f32(
                    ct_table, log10(Frhoc_rc), e_v_rc, poe_rc,
                    ct_X1, ct_Y1, ct_Z1, ct_dxi, ct_dyi, ct_dzi,
                    ct_nx, ct_ny, ct_nz, ct_out_val)
            end
            zrc_shifted = zvir_rc - 1.0f0
            if zrc_shifted > zvir_half
                zvir_half = zrc_shifted
            end

            jj_int += _I1
        end
    end

    # ----- Write outputs -----
    @inbounds RTHL_out[peak_id]  = RTHL
    @inbounds Fbarx_out[peak_id] = Fbarx
    @inbounds e_v_out[peak_id]   = e_v_f
    @inbounds p_v_out[peak_id]   = p_v_f
    @inbounds Srb_out[peak_id]   = Srb
    @inbounds d2F_out[peak_id]   = d2F
    @inbounds zvir_half_out[peak_id] = zvir_half
    @inbounds Sbar_out[_I1, peak_id]  = Sbar_x
    @inbounds Sbar_out[Int32(2), peak_id] = Sbar_y
    @inbounds Sbar_out[Int32(3), peak_id] = Sbar_z
    @inbounds Sbar2_out[_I1, peak_id] = Sbar2_x
    @inbounds Sbar2_out[Int32(2), peak_id] = Sbar2_y
    @inbounds Sbar2_out[Int32(3), peak_id] = Sbar2_z
    @inbounds gradpk_out[_I1, peak_id]  = gx_mat
    @inbounds gradpk_out[Int32(2), peak_id] = gy_mat
    @inbounds gradpk_out[Int32(3), peak_id] = gz_mat
    @inbounds gradpkf_out[_I1, peak_id]  = gfx_mat
    @inbounds gradpkf_out[Int32(2), peak_id] = gfy_mat
    @inbounds gradpkf_out[Int32(3), peak_id] = gfz_mat
    @inbounds gradpkrf_out[_I1, peak_id]  = gradpkrf_x
    @inbounds gradpkrf_out[Int32(2), peak_id] = gradpkrf_y
    @inbounds gradpkrf_out[Int32(3), peak_id] = gradpkrf_z
    @inbounds strain_final_out[1, 1, peak_id] = E11f
    @inbounds strain_final_out[1, 2, peak_id] = E12f
    @inbounds strain_final_out[1, 3, peak_id] = E13f
    @inbounds strain_final_out[2, 1, peak_id] = E12f
    @inbounds strain_final_out[2, 2, peak_id] = E22f
    @inbounds strain_final_out[2, 3, peak_id] = E23f
    @inbounds strain_final_out[3, 1, peak_id] = E13f
    @inbounds strain_final_out[3, 2, peak_id] = E23f
    @inbounds strain_final_out[3, 3, peak_id] = E33f
    @inbounds eigs_out[1, peak_id] = lamf1
    @inbounds eigs_out[2, peak_id] = lamf2
    @inbounds eigs_out[3, peak_id] = lamf3

    # ----- Mask side-effect -----
    if update_mask
        s_mask = _I1
        while s_mask <= nshells
            @inbounds r2_s = shell_r2[s_mask]
            if sqrt(Float32(r2_s)) > RTHL
                break
            end
            @inbounds s0 = shell_start[s_mask]
            @inbounds nc = shell_count[s_mask]
            c_mask = _I1
            while c_mask <= nc
                idx = s0 + c_mask - _I1
                @inbounds di = off_di[idx]
                @inbounds dj = off_dj[idx]
                @inbounds dk = off_dk[idx]
                iv1 = ci + di; iv2 = cj + dj; iv3 = ck + dk
                if (_I1 <= iv1) & (iv1 <= n1) & (_I1 <= iv2) & (iv2 <= n2) & (_I1 <= iv3) & (iv3 <= n3)
                    if (iv1 > nbuff) & (iv1 <= n1 - nbuff) &
                       (iv2 > nbuff) & (iv2 <= n2 - nbuff) &
                       (iv3 > nbuff) & (iv3 <= n3 - nbuff)
                        @inbounds mask_out[iv1, iv2, iv3] = Int8(1)
                    end
                end
                c_mask += _I1
            end
            s_mask += _I1
        end
    end
    return
end

# ============================================================
# Host entry point: analyse_peak_gpu_cuda
# ============================================================

function PeakPatch.analyse_peak_gpu_cuda(
        delta_h, etax_h, etay_h, etaz_h, eta2x_h, eta2y_h, eta2z_h,
        peaks_i, peaks_j, peaks_k,
        stab::ShellTables,
        ct_table_h::AbstractArray{<:Real,3},
        X1::Real, X2::Real, Y1::Real, Y2::Real, Z1::Real, Z2::Real,
        alatt::Real, ir2min::Integer, ZZon::Real, Rfclvi::Real;
        fcrit_override::Union{Nothing, Real}=nothing,
        growth_tables=nothing,
        rmax2rs::Real=0.0, threads::Int=128, ct_out_val::Real=-1.0,
        lapd::Union{Nothing, AbstractArray{<:Real,3}}=nothing,
        mask::Union{Nothing, AbstractArray{<:Integer,3}}=nothing,
        nbuff::Integer=0)
    # Resolve fcrit using the same precedence as CPU analyse_peak_gpu:
    #   explicit override → growth_tables → hard default 1.686
    fcrit = if fcrit_override !== nothing
        Float64(fcrit_override)
    elseif growth_tables !== nothing
        PeakPatch.RadialShell.fsc_of_z(Float64(ZZon) - 1.0, growth_tables)
    else
        1.686
    end
    npeaks = length(peaks_i)
    nshells = stab.nshells

    # Upload all host-side inputs to GPU once.
    delta_d  = CuArray{Float32}(delta_h)
    etax_d   = CuArray{Float32}(etax_h)
    etay_d   = CuArray{Float32}(etay_h)
    etaz_d   = CuArray{Float32}(etaz_h)
    eta2x_d  = CuArray{Float32}(eta2x_h)
    eta2y_d  = CuArray{Float32}(eta2y_h)
    eta2z_d  = CuArray{Float32}(eta2z_h)
    lapd_d   = lapd === nothing ? CUDA.zeros(Float32, size(delta_h)...) :
                                   (@assert size(lapd) == size(delta_h); CuArray{Float32}(lapd))
    pi_d2 = CuArray{Int32}(peaks_i)
    pj_d2 = CuArray{Int32}(peaks_j)
    pk_d2 = CuArray{Int32}(peaks_k)
    stab_d = ShellTablesGPU(stab)

    # Allocate gather output buffers on GPU (no host round-trip)
    Fshell_d    = CUDA.zeros(Float32, nshells, npeaks)
    nshell_d    = CUDA.zeros(Int32,   nshells, npeaks)
    Sshell_d    = CUDA.zeros(Float32, 3, nshells, npeaks)
    S2shell_d   = CUDA.zeros(Float32, 3, nshells, npeaks)
    Gshell_d    = CUDA.zeros(Float32, 3, nshells, npeaks)
    Gfshell_d   = CUDA.zeros(Float32, 3, nshells, npeaks)
    SRshell_d   = CUDA.zeros(Float32, 3, 3, nshells, npeaks)
    lapdshell_d = CUDA.zeros(Float32, nshells, npeaks)

    # Launch gather kernel (device-side only, no host traffic)
    _launch_shell_gather_full!(
        Fshell_d, nshell_d, Sshell_d, S2shell_d, Gshell_d, Gfshell_d, SRshell_d,
        lapdshell_d,
        delta_d, etax_d, etay_d, etaz_d, eta2x_d, eta2y_d, eta2z_d, lapd_d,
        pi_d2, pj_d2, pk_d2, stab_d; threads=threads)

    # After the gather kernel, the field inputs (delta_d, etas...) are no
    # longer needed by the post-process kernel. Free their pinned pool slots
    # so the mask grid and output arrays have room.
    CUDA.unsafe_free!(delta_d)
    CUDA.unsafe_free!(etax_d); CUDA.unsafe_free!(etay_d); CUDA.unsafe_free!(etaz_d)
    CUDA.unsafe_free!(eta2x_d); CUDA.unsafe_free!(eta2y_d); CUDA.unsafe_free!(eta2z_d)
    CUDA.unsafe_free!(lapd_d)

    shell_r2_d = CuArray{Int32}(stab.shell_r2)
    ct_table_d = CuArray{Float32}(ct_table_h)
    nx, ny, nz = size(ct_table_h)
    dxi = Float32((nx - 1) / (X2 - X1))
    dyi = Float32((ny - 1) / (Y2 - Y1))
    dzi = Float32((nz - 1) / (Z2 - Z1))

    # Determine nshells_max based on rmax2rs.
    #
    # Match RadialShell.analyse_peak which processes nshells - 1 shells (its
    # per-cell loop never finalises the last shell). Without this cap the GPU
    # path is more permissive than the real CPU at the outermost boundary.
    nshells_max = nshells - 1
    if rmax2rs > 0.0
        r2_limit = (Float64(rmax2rs) * Float64(Rfclvi) / Float64(alatt))^2
        count = 0
        for s in 1:nshells - 1
            stab.shell_r2[s] <= r2_limit || break
            count = s
        end
        nshells_max = count
    end

    # Scaling constants from CPU analyse_peak_gpu
    hlatt = 1.0
    anor = 1.0 / hlatt^3
    aRnor = 3.0 * anor / alatt
    wnor = anor / (4.0 * pi)
    wRnor = 3.0 * wnor / alatt
    hlatt_1 = 1.0  # corresponds to CPU default
    hlatt_2 = 1.0

    # Allocate outputs
    RTHL_d          = CUDA.zeros(Float32, npeaks)
    Fbarx_d         = CUDA.zeros(Float32, npeaks)
    e_v_d           = CUDA.zeros(Float32, npeaks)
    p_v_d           = CUDA.zeros(Float32, npeaks)
    strain_final_d  = CUDA.zeros(Float32, 3, 3, npeaks)
    eigs_d          = CUDA.zeros(Float32, 3, npeaks)
    Srb_d           = CUDA.zeros(Float32, npeaks)
    Sbar_d          = CUDA.zeros(Float32, 3, npeaks)
    Sbar2_d         = CUDA.zeros(Float32, 3, npeaks)
    gradpk_d        = CUDA.zeros(Float32, 3, npeaks)
    gradpkf_d       = CUDA.zeros(Float32, 3, npeaks)
    gradpkrf_d      = CUDA.zeros(Float32, 3, npeaks)
    d2F_d           = CUDA.zeros(Float32, npeaks)
    zvir_half_d     = CUDA.zeros(Float32, npeaks)

    # Mask: if provided, upload as CuArray{Int8,3}. Otherwise allocate a
    # 1×1×1 dummy so the kernel argument is always a valid CuArray;
    # update_mask flag decides whether to actually write to it.
    update_mask = mask !== nothing
    n1_delta, n2_delta, n3_delta = size(delta_h)
    mask_d = if update_mask
        @assert size(mask) == size(delta_h)
        CuArray{Int8}(mask)
    else
        CuArray{Int8}(zeros(Int8, 1, 1, 1))
    end

    # Shared memory: rad + Fbar + 3*Gn + 3*Gfn + 9*SRn = 17 Float32 × MAX_SHELLS
    shmem_bytes = 17 * Int(_MAX_SHELLS_GPU) * sizeof(Float32)

    Rfclvi_r2 = Float32((Rfclvi / alatt)^2)

    # Per-peak ZZon / fcrit arrays. For ievol==0 we broadcast the scalar
    # values; for ievol==1 the caller provides per-peak arrays in the
    # multirf entry point. This single-call entry always broadcasts.
    ZZon_pp_d  = CUDA.fill(Float32(ZZon),  npeaks)
    fcrit_pp_d = CUDA.fill(Float32(fcrit), npeaks)

    @cuda threads=threads blocks=npeaks shmem=shmem_bytes _post_process_kernel!(
        RTHL_d, Fbarx_d, e_v_d, p_v_d, strain_final_d, eigs_d,
        Srb_d, Sbar_d, Sbar2_d, gradpk_d, gradpkf_d, gradpkrf_d, d2F_d,
        zvir_half_d, mask_d,
        Fshell_d, nshell_d, Sshell_d, S2shell_d, Gshell_d, Gfshell_d, SRshell_d,
        lapdshell_d,
        pi_d2, pj_d2, pk_d2,
        stab_d.offsets_di, stab_d.offsets_dj, stab_d.offsets_dk,
        stab_d.shell_start, stab_d.shell_count,
        shell_r2_d,
        ct_table_d,
        Float32(X1), Float32(Y1), Float32(Z1),
        dxi, dyi, dzi,
        Int32(nx), Int32(ny), Int32(nz), Float32(ct_out_val),
        _AKK_GPU[],
        ZZon_pp_d, fcrit_pp_d, Rfclvi_r2,
        Float32(wRnor), Float32(aRnor), Float32(hlatt_1), Float32(hlatt_2),
        Int32(ir2min),
        Int32(nshells), Int32(nshells_max),
        Int32(n1_delta), Int32(n2_delta), Int32(n3_delta), Int32(nbuff), update_mask,
    )

    return (
        RTHL         = Array(RTHL_d),
        Fbarx        = Array(Fbarx_d),
        e_v          = Array(e_v_d),
        p_v          = Array(p_v_d),
        strain_final = Array(strain_final_d),
        eigs         = Array(eigs_d),
        Srb          = Array(Srb_d),
        Sbar         = Array(Sbar_d),
        Sbar2        = Array(Sbar2_d),
        gradpk       = Array(gradpk_d),
        gradpkf      = Array(gradpkf_d),
        gradpkrf     = Array(gradpkrf_d),
        d2F          = Array(d2F_d),
        zvir_half    = Array(zvir_half_d),
        mask         = update_mask ? Array(mask_d) : nothing,
    )
end

# ============================================================
# Host entry point: analyse_peaks_gpu_cuda_multirf
# ============================================================
#
# Batched-per-Rf variant. Uploads the seven field arrays + lapd + mask ONCE,
# then runs (gather + post-process) for each entry in `batches`. Used by
# `run_multitile_split(use_gpu=true)` which analyses each tile at 3 filter
# scales — the single-call variant re-uploads 200+ MB of fields per Rf,
# totalling several GB of redundant PCIe traffic per run.
#
# `batches` entries: `(peaks_i, peaks_j, peaks_k, Rf::Real, ir2min::Integer)`.

function PeakPatch.analyse_peaks_gpu_cuda_multirf(
        delta_h, etax_h, etay_h, etaz_h, eta2x_h, eta2y_h, eta2z_h,
        batches::AbstractVector,
        stab::ShellTables,
        ct_table_h::AbstractArray{<:Real,3},
        X1::Real, X2::Real, Y1::Real, Y2::Real, Z1::Real, Z2::Real,
        alatt::Real, ZZon::Real;
        fcrit_override::Union{Nothing, Real}=nothing,
        growth_tables=nothing,
        rmax2rs::Real=0.0, threads::Int=128, ct_out_val::Real=-1.0,
        lapd::Union{Nothing, AbstractArray{<:Real,3}}=nothing,
        mask::Union{Nothing, AbstractArray{<:Integer,3}}=nothing,
        nbuff::Integer=0)

    fcrit = if fcrit_override !== nothing
        Float64(fcrit_override)
    elseif growth_tables !== nothing
        PeakPatch.RadialShell.fsc_of_z(Float64(ZZon) - 1.0, growth_tables)
    else
        1.686
    end
    nshells = stab.nshells

    # ---------------- One-time uploads ----------------
    # Any of (delta_h, eta*_h, eta2*_h, lapd) may already be a device array —
    # in that case we reuse it directly and avoid a redundant upload of the
    # same O(nmesh³ × Float32) volume (~250 MB at nmesh=399). Caller owns
    # inputs, so we only free the ones we just uploaded.
    _upload(x) = x isa CuArray{Float32,3} ? (x, false) : (CuArray{Float32}(x), true)
    delta_d,  _free_delta  = _upload(delta_h)
    etax_d,   _free_etax   = _upload(etax_h)
    etay_d,   _free_etay   = _upload(etay_h)
    etaz_d,   _free_etaz   = _upload(etaz_h)
    eta2x_d,  _free_eta2x  = _upload(eta2x_h)
    eta2y_d,  _free_eta2y  = _upload(eta2y_h)
    eta2z_d,  _free_eta2z  = _upload(eta2z_h)
    _free_lapd = false
    if lapd === nothing
        lapd_d = CUDA.zeros(Float32, size(delta_d)...)
        _free_lapd = true
    elseif lapd isa CuArray{Float32,3}
        @assert size(lapd) == size(delta_d)
        lapd_d = lapd
    else
        @assert size(lapd) == size(delta_d)
        lapd_d = CuArray{Float32}(lapd)
        _free_lapd = true
    end
    stab_d     = ShellTablesGPU(stab)
    shell_r2_d = CuArray{Int32}(stab.shell_r2)
    ct_table_d = CuArray{Float32}(ct_table_h)

    update_mask = mask !== nothing
    n1_delta, n2_delta, n3_delta = size(delta_d)
    mask_d = if update_mask
        @assert size(mask) == size(delta_d)
        CuArray{Int8}(mask)
    else
        CuArray{Int8}(zeros(Int8, 1, 1, 1))
    end

    nx, ny, nz = size(ct_table_h)
    dxi = Float32((nx - 1) / (X2 - X1))
    dyi = Float32((ny - 1) / (Y2 - Y1))
    dzi = Float32((nz - 1) / (Z2 - Z1))

    # Scaling constants
    hlatt = 1.0
    anor = 1.0 / hlatt^3
    aRnor = 3.0 * anor / alatt
    wnor = anor / (4.0 * pi)
    wRnor = 3.0 * wnor / alatt
    hlatt_1 = 1.0
    hlatt_2 = 1.0

    # Post-process kernel shmem budget (rad, Fbar, Gn, Gfn, SRn staging).
    shmem_bytes = 17 * Int(_MAX_SHELLS_GPU) * sizeof(Float32)

    results = Vector{NamedTuple}(undef, length(batches))

    # Persistent per-call scratch pool — sized to the largest batch so we
    # can reuse every buffer across every Rf without per-batch
    # CUDA.zeros + unsafe_free! churn. npeaks varies per batch (one batch
    # per unique Rf) but all batches in a single multirf call share the
    # same nshells / field geometry.
    #
    # Empirically the fused single-kernel variant (see
    # `_fused_shell_analysis_kernel!`, kept in the file for future reuse)
    # costs more than it saves on SM_86: combining gather + post-process
    # raises per-block shmem from 256 B to ~20 KB, which drops occupancy
    # from 12 blocks/SM to 2 blocks/SM during the memory-bound gather
    # phase. The lost latency-hiding outweighs the ~100 MB/batch of
    # global traffic saved. We keep the two-kernel pipeline here but
    # reuse the persistent scratch below to cut per-batch allocation
    # overhead, which is pure win.
    max_np = 0
    for b in batches
        max_np = max(max_np, length(b[1]))
    end
    pi_scratch       = CUDA.zeros(Int32,   max_np)
    pj_scratch       = CUDA.zeros(Int32,   max_np)
    pk_scratch       = CUDA.zeros(Int32,   max_np)
    ZZon_pp_scratch  = CUDA.zeros(Float32, max_np)
    fcrit_pp_scratch = CUDA.zeros(Float32, max_np)
    # Per-shell intermediates between gather and post-process (sized to max_np).
    Fshell_scratch    = CUDA.zeros(Float32, nshells, max_np)
    nshell_scratch    = CUDA.zeros(Int32,   nshells, max_np)
    Sshell_scratch    = CUDA.zeros(Float32, 3, nshells, max_np)
    S2shell_scratch   = CUDA.zeros(Float32, 3, nshells, max_np)
    Gshell_scratch    = CUDA.zeros(Float32, 3, nshells, max_np)
    Gfshell_scratch   = CUDA.zeros(Float32, 3, nshells, max_np)
    SRshell_scratch   = CUDA.zeros(Float32, 3, 3, nshells, max_np)
    lapdshell_scratch = CUDA.zeros(Float32, nshells, max_np)
    # Per-peak outputs.
    RTHL_scratch         = CUDA.zeros(Float32, max_np)
    Fbarx_scratch        = CUDA.zeros(Float32, max_np)
    e_v_scratch          = CUDA.zeros(Float32, max_np)
    p_v_scratch          = CUDA.zeros(Float32, max_np)
    strain_final_scratch = CUDA.zeros(Float32, 3, 3, max_np)
    eigs_scratch         = CUDA.zeros(Float32, 3, max_np)
    Srb_scratch          = CUDA.zeros(Float32, max_np)
    Sbar_scratch         = CUDA.zeros(Float32, 3, max_np)
    Sbar2_scratch        = CUDA.zeros(Float32, 3, max_np)
    gradpk_scratch       = CUDA.zeros(Float32, 3, max_np)
    gradpkf_scratch      = CUDA.zeros(Float32, 3, max_np)
    gradpkrf_scratch     = CUDA.zeros(Float32, 3, max_np)
    d2F_scratch          = CUDA.zeros(Float32, max_np)
    zvir_half_scratch    = CUDA.zeros(Float32, max_np)

    # ---------------- Per-Rf batch work ----------------
    for (bi, batch) in enumerate(batches)
        # Batch tuple is either (peaks_i, peaks_j, peaks_k, Rf, ir2min) for
        # the ievol==0 scalar-broadcast path, or extended to length 7 with
        # trailing (ZZon_pp_b, fcrit_pp_b) per-peak arrays for ievol==1.
        peaks_i = batch[1]; peaks_j = batch[2]; peaks_k = batch[3]
        Rf      = batch[4]; ir2min  = batch[5]
        has_pp_b = length(batch) >= 7
        ZZon_pp_b_h  = has_pp_b ? batch[6] : nothing
        fcrit_pp_b_h = has_pp_b ? batch[7] : nothing
        npeaks = length(peaks_i)

        # Slice views into persistent scratch (npeaks ≤ max_np).
        pi_d = view(pi_scratch, 1:npeaks)
        pj_d = view(pj_scratch, 1:npeaks)
        pk_d = view(pk_scratch, 1:npeaks)
        copyto!(pi_d, Int32.(peaks_i))
        copyto!(pj_d, Int32.(peaks_j))
        copyto!(pk_d, Int32.(peaks_k))

        # Per-shell intermediate views.
        Fshell_d    = view(Fshell_scratch,    :, 1:npeaks)
        nshell_d    = view(nshell_scratch,    :, 1:npeaks)
        Sshell_d    = view(Sshell_scratch,    :, :, 1:npeaks)
        S2shell_d   = view(S2shell_scratch,   :, :, 1:npeaks)
        Gshell_d    = view(Gshell_scratch,    :, :, 1:npeaks)
        Gfshell_d   = view(Gfshell_scratch,   :, :, 1:npeaks)
        SRshell_d   = view(SRshell_scratch,   :, :, :, 1:npeaks)
        lapdshell_d = view(lapdshell_scratch, :, 1:npeaks)

        _launch_shell_gather_full!(
            Fshell_d, nshell_d, Sshell_d, S2shell_d, Gshell_d, Gfshell_d, SRshell_d,
            lapdshell_d,
            delta_d, etax_d, etay_d, etaz_d, eta2x_d, eta2y_d, eta2z_d, lapd_d,
            pi_d, pj_d, pk_d, stab_d; threads=threads)

        # Per-Rf nshells_max (match CPU: excludes last shell)
        nshells_max = nshells - 1
        if rmax2rs > 0.0
            r2_limit = (Float64(rmax2rs) * Float64(Rf) / Float64(alatt))^2
            count = 0
            for s in 1:nshells - 1
                stab.shell_r2[s] <= r2_limit || break
                count = s
            end
            nshells_max = count
        end

        RTHL_d         = view(RTHL_scratch,         1:npeaks)
        Fbarx_d        = view(Fbarx_scratch,        1:npeaks)
        e_v_d          = view(e_v_scratch,          1:npeaks)
        p_v_d          = view(p_v_scratch,          1:npeaks)
        strain_final_d = view(strain_final_scratch, :, :, 1:npeaks)
        eigs_d         = view(eigs_scratch,         :, 1:npeaks)
        Srb_d          = view(Srb_scratch,          1:npeaks)
        Sbar_d         = view(Sbar_scratch,         :, 1:npeaks)
        Sbar2_d        = view(Sbar2_scratch,        :, 1:npeaks)
        gradpk_d       = view(gradpk_scratch,       :, 1:npeaks)
        gradpkf_d      = view(gradpkf_scratch,      :, 1:npeaks)
        gradpkrf_d     = view(gradpkrf_scratch,     :, 1:npeaks)
        d2F_d          = view(d2F_scratch,          1:npeaks)
        zvir_half_d    = view(zvir_half_scratch,    1:npeaks)

        Rfclvi_r2 = Float32((Rf / alatt)^2)

        # Per-peak ZZon / fcrit vectors (sliced views into persistent scratch).
        ZZon_pp_d  = view(ZZon_pp_scratch,  1:npeaks)
        fcrit_pp_d = view(fcrit_pp_scratch, 1:npeaks)
        if has_pp_b
            copyto!(ZZon_pp_d,  Float32.(ZZon_pp_b_h))
            copyto!(fcrit_pp_d, Float32.(fcrit_pp_b_h))
        else
            fill!(ZZon_pp_d,  Float32(ZZon))
            fill!(fcrit_pp_d, Float32(fcrit))
        end

        @cuda threads=threads blocks=npeaks shmem=shmem_bytes _post_process_kernel!(
            RTHL_d, Fbarx_d, e_v_d, p_v_d, strain_final_d, eigs_d,
            Srb_d, Sbar_d, Sbar2_d, gradpk_d, gradpkf_d, gradpkrf_d, d2F_d,
            zvir_half_d, mask_d,
            Fshell_d, nshell_d, Sshell_d, S2shell_d, Gshell_d, Gfshell_d, SRshell_d,
            lapdshell_d,
            pi_d, pj_d, pk_d,
            stab_d.offsets_di, stab_d.offsets_dj, stab_d.offsets_dk,
            stab_d.shell_start, stab_d.shell_count,
            shell_r2_d,
            ct_table_d,
            Float32(X1), Float32(Y1), Float32(Z1),
            dxi, dyi, dzi,
            Int32(nx), Int32(ny), Int32(nz), Float32(ct_out_val),
            _AKK_GPU[],
            ZZon_pp_d, fcrit_pp_d, Rfclvi_r2,
            Float32(wRnor), Float32(aRnor), Float32(hlatt_1), Float32(hlatt_2),
            Int32(ir2min),
            Int32(nshells), Int32(nshells_max),
            Int32(n1_delta), Int32(n2_delta), Int32(n3_delta), Int32(nbuff), update_mask,
        )

        results[bi] = (
            RTHL         = Array(RTHL_d),
            Fbarx        = Array(Fbarx_d),
            e_v          = Array(e_v_d),
            p_v          = Array(p_v_d),
            strain_final = Array(strain_final_d),
            eigs         = Array(eigs_d),
            Srb          = Array(Srb_d),
            Sbar         = Array(Sbar_d),
            Sbar2        = Array(Sbar2_d),
            gradpk       = Array(gradpk_d),
            gradpkf      = Array(gradpkf_d),
            gradpkrf     = Array(gradpkrf_d),
            d2F          = Array(d2F_d),
            zvir_half    = Array(zvir_half_d),
        )
        # Nothing to free per-batch — all buffers are persistent scratch.
    end

    # Release the per-call persistent scratch pool.
    CUDA.unsafe_free!(pi_scratch); CUDA.unsafe_free!(pj_scratch); CUDA.unsafe_free!(pk_scratch)
    CUDA.unsafe_free!(ZZon_pp_scratch); CUDA.unsafe_free!(fcrit_pp_scratch)
    CUDA.unsafe_free!(Fshell_scratch); CUDA.unsafe_free!(nshell_scratch)
    CUDA.unsafe_free!(Sshell_scratch); CUDA.unsafe_free!(S2shell_scratch)
    CUDA.unsafe_free!(Gshell_scratch); CUDA.unsafe_free!(Gfshell_scratch)
    CUDA.unsafe_free!(SRshell_scratch); CUDA.unsafe_free!(lapdshell_scratch)
    CUDA.unsafe_free!(RTHL_scratch); CUDA.unsafe_free!(Fbarx_scratch)
    CUDA.unsafe_free!(e_v_scratch); CUDA.unsafe_free!(p_v_scratch)
    CUDA.unsafe_free!(strain_final_scratch); CUDA.unsafe_free!(eigs_scratch)
    CUDA.unsafe_free!(Srb_scratch); CUDA.unsafe_free!(Sbar_scratch); CUDA.unsafe_free!(Sbar2_scratch)
    CUDA.unsafe_free!(gradpk_scratch); CUDA.unsafe_free!(gradpkf_scratch); CUDA.unsafe_free!(gradpkrf_scratch)
    CUDA.unsafe_free!(d2F_scratch); CUDA.unsafe_free!(zvir_half_scratch)

    # Free the one-time uploads (only the buffers we allocated — caller-owned
    # CuArrays that were passed in are left alone).
    _free_delta && CUDA.unsafe_free!(delta_d)
    _free_etax  && CUDA.unsafe_free!(etax_d)
    _free_etay  && CUDA.unsafe_free!(etay_d)
    _free_etaz  && CUDA.unsafe_free!(etaz_d)
    _free_eta2x && CUDA.unsafe_free!(eta2x_d)
    _free_eta2y && CUDA.unsafe_free!(eta2y_d)
    _free_eta2z && CUDA.unsafe_free!(eta2z_d)
    _free_lapd  && CUDA.unsafe_free!(lapd_d)
    CUDA.unsafe_free!(shell_r2_d); CUDA.unsafe_free!(ct_table_d)

    mask_out = update_mask ? Array(mask_d) : nothing
    CUDA.unsafe_free!(mask_d)
    return (results=results, mask=mask_out)
end

# ============================================================
# GPU isolated FFT convolution (Phase 1 of GPU plan)
# ============================================================
#
# Replaces CPU `MultiResolution._isolated_convolve`:
#   noise (n_input³) → zero-pad to (2*n_input)³ → rfft → multiply by
#   √P(k) [optional kernel_fn(kx,ky,kz,k²)] → irfft → extract tile.
#
# GPU pipeline:
# 1. Build P(k) lookup (sqrt-amplitude in log-k space) on host, upload once
# 2. Upload noise, zero-pad
# 3. cuFFT plan_rfft / mul!(plan_fwd) into complex buffer
# 4. Element-wise CUDA kernel applies √P(k) and selected transfer kernel
# 5. cuFFT plan_irfft / mul!(plan_inv) back to real
# 6. Extract central tile region

@inline function _interp_pk_table(pk_table, log_k::Float32, log_k_min::Float32,
                                  d_log_k::Float32, ntable::Int32)
    fpos = (log_k - log_k_min) / d_log_k
    i0 = Int32(floor(fpos)) + _I1
    if i0 < _I1
        @inbounds return pk_table[1]
    elseif i0 >= ntable
        @inbounds return pk_table[ntable]
    else
        frac = fpos - Float32(i0 - _I1)
        @inbounds return pk_table[i0] + frac * (pk_table[i0 + _I1] - pk_table[i0])
    end
end

# Transfer kernel: each thread handles one (ix, iy, iz) of the rfft buffer.
# `padded_k` is a CuArray{ComplexF32,3} of shape (n2÷2+1, n2, n2).
function _transfer_kernel!(padded_k, pk_table, log_k_min::Float32,
                            d_log_k::Float32, ntable::Int32,
                            dk::Float64, n2::Int32,
                            kernel_fn_id::Int32, dim1::Int32, dim2::Int32)
    nk = (n2 >> _I1) + _I1   # rfft first-dim length
    total = Int64(nk) * Int64(n2) * Int64(n2)
    idx0 = Int64(blockIdx().x - _I1) * Int64(blockDim().x) + Int64(threadIdx().x)
    if idx0 > total
        return
    end
    # Decompose linear idx (1-based) into (ix, iy, iz) — column-major
    li = idx0 - Int64(1)
    ix = Int32(li % Int64(nk)) + _I1
    rem = li ÷ Int64(nk)
    iy = Int32(rem % Int64(n2)) + _I1
    iz = Int32(rem ÷ Int64(n2)) + _I1

    # Frequencies in Float64 — matches CPU `_isolated_convolve` which casts
    # kx_arr[ix] to Float64 before computing transfer function. Without this
    # the per-mode error compounds to ~1e-5 at n≥64, beyond Float32 noise.
    half = (n2 >> _I1)  # n2 / 2
    kx = Float64(ix - _I1) * dk
    iy_signed = iy <= half ? (iy - _I1) : (iy - n2 - _I1)
    iz_signed = iz <= half ? (iz - _I1) : (iz - n2 - _I1)
    ky = Float64(iy_signed) * dk
    kz = Float64(iz_signed) * dk
    k2 = kx * kx + ky * ky + kz * kz

    if k2 == 0.0
        @inbounds padded_k[ix, iy, iz] = ComplexF32(0)
        return
    end
    sqrt_pk = Float64(_interp_pk_table(pk_table, Float32(log(k2) * 0.5),
                                        log_k_min, d_log_k, ntable))

    @inbounds val = padded_k[ix, iy, iz]

    if kernel_fn_id == _I0
        # δ: just √P(k)
        @inbounds padded_k[ix, iy, iz] = ComplexF32(val * sqrt_pk)
    elseif kernel_fn_id == _I1
        # 1LPT: i * ki / k² × √P(k)
        ki = (dim1 == _I1) ? kx : ((dim1 == Int32(2)) ? ky : kz)
        coeff = sqrt_pk * ki / k2
        a = Float64(real(val)); b = Float64(imag(val))
        @inbounds padded_k[ix, iy, iz] = ComplexF32(-b * coeff, a * coeff)
    elseif kernel_fn_id == Int32(2)
        # 2LPT: -i * ki / k² × √P(k)
        ki = (dim1 == _I1) ? kx : ((dim1 == Int32(2)) ? ky : kz)
        coeff = sqrt_pk * ki / k2
        a = Float64(real(val)); b = Float64(imag(val))
        @inbounds padded_k[ix, iy, iz] = ComplexF32(b * coeff, -a * coeff)
    elseif kernel_fn_id == Int32(3)
        # φ_ij: -ki*kj/k² × √P(k)  (real coefficient)
        ki = (dim1 == _I1) ? kx : ((dim1 == Int32(2)) ? ky : kz)
        kj = (dim2 == _I1) ? kx : ((dim2 == Int32(2)) ? ky : kz)
        coeff = Float32(-sqrt_pk * ki * kj / k2)
        @inbounds padded_k[ix, iy, iz] = val * coeff
    elseif kernel_fn_id == Int32(4)
        # Laplacian: k² × √P(k)
        coeff = Float32(sqrt_pk * k2)
        @inbounds padded_k[ix, iy, iz] = val * coeff
    end
    return
end

# Enforce Hermitian symmetry on the DC (ix=1) and Nyquist (ix=n2/2+1) planes
# of an rfft buffer. cuFFT's irfft (CUFFT_C2R) gives results that diverge
# from FFTW's irfft by ~5% when the input has Hermitian-violating content on
# these special planes (e.g., from the φ_ij transfer for off-diagonal ij at
# n≥64). FFTW evidently averages the conjugate pairs implicitly; cuFFT does
# something different. We match FFTW by symmetrising explicitly.
#
# For each (iy, iz) in {DC, Nyquist} plane, the conjugate-pair index is
# (iy_c, iz_c) where iy_c = iy==1 ? 1 : (n2 - iy + 2), similarly iz_c.
# Self-conjugate cells (iy in {1, n2/2+1} AND iz in {1, n2/2+1}) → imag := 0.
# Otherwise: avg = 0.5*(F[iy,iz] + conj(F[iy_c,iz_c])); F[iy,iz] = avg;
# the partner cell will independently compute its own avg = conj(this).
function _hermitize_planes_kernel!(padded_k, n2::Int32, plane_ix1::Int32,
                                    plane_ix2::Int32)
    idx0 = (blockIdx().x - _I1) * blockDim().x + threadIdx().x
    total = Int64(n2) * Int64(n2)
    if idx0 > total
        return
    end
    li = idx0 - _I1
    iy = Int32(li % Int64(n2)) + _I1
    iz = Int32(li ÷ Int64(n2)) + _I1
    iy_c = iy == _I1 ? _I1 : (n2 - iy + Int32(2))
    iz_c = iz == _I1 ? _I1 : (n2 - iz + Int32(2))

    # Avoid read/write races: only the "lower" of each conjugate pair does
    # the work, writing BOTH slots. Lower defined lexicographically by
    # (iy, iz) < (iy_c, iz_c).
    is_lower = (iy < iy_c) | ((iy == iy_c) & (iz < iz_c))
    is_self  = (iy == iy_c) & (iz == iz_c)
    if !(is_lower | is_self)
        return
    end

    for ix in (plane_ix1, plane_ix2)
        if is_self
            # Self-conjugate cell: must be real
            @inbounds val = padded_k[ix, iy, iz]
            @inbounds padded_k[ix, iy, iz] = ComplexF32(real(val), 0.0f0)
        else
            # Average with conjugate partner; write both slots
            @inbounds vme = padded_k[ix, iy, iz]
            @inbounds vco = padded_k[ix, iy_c, iz_c]
            avg = (vme + conj(vco)) * 0.5f0
            @inbounds padded_k[ix, iy, iz]     = avg
            @inbounds padded_k[ix, iy_c, iz_c] = conj(avg)
        end
    end
    return
end

function PeakPatch.isolated_convolve_gpu(noise::AbstractArray{<:Real,3},
                                          pk, boxsize_local::Real, n::Integer;
                                          kernel_fn_id::Integer=0,
                                          dim1::Integer=0, dim2::Integer=0,
                                          nshell::Integer=0,
                                          pk_table_n::Integer=8192,
                                          threads::Integer=256)
    n_input = size(noise, 1)
    @assert size(noise) == (n_input, n_input, n_input)
    n2 = 2 * n_input
    dx = Float64(boxsize_local) / n
    dk = 2π / (n2 * dx)

    # ----- Pre-compute √(P(k) · dk³ · n2³) lookup in log-k space -----
    # Range: k_min just above 0 to k_max ≈ √3 × Nyquist
    k_min_raw = dk
    k_max_raw = sqrt(3.0) * (n2 / 2) * dk
    log_k = range(log(k_min_raw * 0.5), log(k_max_raw * 1.5), length=pk_table_n)
    norm = dk^3 * Float64(n2)^3
    pk_vals = Float32[sqrt(Float64(pk(exp(lk))) * norm) for lk in log_k]
    pk_table_d = CuArray(pk_vals)
    log_k_min_f = Float32(log_k[1])
    d_log_k_f   = Float32(step(log_k))

    # ----- Upload noise + zero-pad on GPU -----
    padded_d = CUDA.zeros(Float32, n2, n2, n2)
    noise_d = CuArray{Float32}(noise)
    @inbounds padded_d[1:n_input, 1:n_input, 1:n_input] .= noise_d
    CUDA.unsafe_free!(noise_d)

    # ----- Forward FFT -----
    padded_k_d = _cufft_rfft(padded_d)
    CUDA.unsafe_free!(padded_d)

    # ----- Apply transfer function -----
    nk = size(padded_k_d, 1)
    total = nk * n2 * n2
    blocks = cld(total, threads)
    @cuda threads=threads blocks=blocks _transfer_kernel!(
        padded_k_d, pk_table_d,
        log_k_min_f, d_log_k_f, Int32(pk_table_n),
        Float64(dk), Int32(n2),
        Int32(kernel_fn_id), Int32(dim1), Int32(dim2),
    )

    # ----- Enforce Hermitian symmetry on DC and Nyquist planes -----
    # Necessary so that cuFFT.irfft matches FFTW.irfft. cuFFT and FFTW
    # diverge on Hermitian-violating inputs at these planes (~5% relerr
    # for off-diagonal φ_ij at n≥64). FFTW averages conjugate pairs
    # implicitly; we do it explicitly.
    plane_total = n2 * n2
    plane_blocks = cld(plane_total, threads)
    @cuda threads=threads blocks=plane_blocks _hermitize_planes_kernel!(
        padded_k_d, Int32(n2), _I1, Int32(n2 ÷ 2 + 1))

    # ----- Inverse FFT -----
    result_d = _cufft_irfft(padded_k_d, n2)
    CUDA.unsafe_free!(padded_k_d)

    # ----- Extract tile region -----
    s1 = nshell + 1
    s2 = nshell + n
    tile_d = result_d[s1:s2, s1:s2, s1:s2]
    CUDA.unsafe_free!(result_d)

    return Array(tile_d)
end

# -------------------------------------------------------------------
# Multi-output isolated convolution: shares one (upload + zero-pad +
# forward rFFT) of the noise across N output kernels. Used in the tile
# loop to compute (δ, ψ₁_x, ψ₁_y, ψ₁_z) from the same residual noise.
# Each entry in `kernels` is a 3-tuple (kernel_fn_id, dim1, dim2)
# matching the IDs accepted by isolated_convolve_gpu.
#
# Returns Vector{Array{Float32,3}} of length(kernels), each (n,n,n).
# Hermitisation must be applied per-kernel because some transfer
# kernels (1LPT i·kx/k²) break Hermitian symmetry on the Nyquist
# planes; we re-hermitise on the per-kernel scratch buffer.
# -------------------------------------------------------------------
function PeakPatch.isolated_convolve_gpu_multi(noise::AbstractArray{<:Real,3},
                                                pk, boxsize_local::Real, n::Integer;
                                                kernels::AbstractVector,
                                                nshell::Integer=0,
                                                pk_table_n::Integer=8192,
                                                threads::Integer=256,
                                                return_device::Bool=false)
    n_input = size(noise, 1)
    @assert size(noise) == (n_input, n_input, n_input)
    n2 = 2 * n_input
    dx = Float64(boxsize_local) / n
    dk = 2π / (n2 * dx)

    # Shared P(k) lookup
    k_min_raw = dk
    k_max_raw = sqrt(3.0) * (n2 / 2) * dk
    log_k = range(log(k_min_raw * 0.5), log(k_max_raw * 1.5), length=pk_table_n)
    norm = dk^3 * Float64(n2)^3
    pk_vals = Float32[sqrt(Float64(pk(exp(lk))) * norm) for lk in log_k]
    pk_table_d = CuArray(pk_vals)
    log_k_min_f = Float32(log_k[1])
    d_log_k_f   = Float32(step(log_k))

    # Shared upload + zero-pad + forward rFFT of the noise
    padded_d = CUDA.zeros(Float32, n2, n2, n2)
    noise_d  = CuArray{Float32}(noise)
    @inbounds padded_d[1:n_input, 1:n_input, 1:n_input] .= noise_d
    CUDA.unsafe_free!(noise_d)
    base_k_d = _cufft_rfft(padded_d)
    CUDA.unsafe_free!(padded_d)

    nk = size(base_k_d, 1)
    total = nk * n2 * n2
    blocks = cld(total, threads)
    plane_total = n2 * n2
    plane_blocks = cld(plane_total, threads)
    s1 = nshell + 1
    s2 = nshell + n

    if return_device
        out_dev = Vector{CuArray{Float32,3}}(undef, length(kernels))
        for (ki, ktup) in enumerate(kernels)
            kernel_fn_id, dim1, dim2 = ktup
            scratch_k_d = copy(base_k_d)
            @cuda threads=threads blocks=blocks _transfer_kernel!(
                scratch_k_d, pk_table_d,
                log_k_min_f, d_log_k_f, Int32(pk_table_n),
                Float64(dk), Int32(n2),
                Int32(kernel_fn_id), Int32(dim1), Int32(dim2),
            )
            @cuda threads=threads blocks=plane_blocks _hermitize_planes_kernel!(
                scratch_k_d, Int32(n2), _I1, Int32(n2 ÷ 2 + 1))
            result_d = _cufft_irfft(scratch_k_d, n2)
            CUDA.unsafe_free!(scratch_k_d)
            # Materialise the (s1:s2)³ slice into an owning CuArray so the
            # caller can free the result directly (slicing returns a view
            # backed by result_d, which we must not free while the caller
            # still holds it).
            out_dev[ki] = CuArray{Float32}(result_d[s1:s2, s1:s2, s1:s2])
            CUDA.unsafe_free!(result_d)
        end
        CUDA.unsafe_free!(base_k_d)
        CUDA.unsafe_free!(pk_table_d)
        return out_dev
    else
        out = Vector{Array{Float32,3}}(undef, length(kernels))
        for (ki, ktup) in enumerate(kernels)
            kernel_fn_id, dim1, dim2 = ktup
            scratch_k_d = copy(base_k_d)
            @cuda threads=threads blocks=blocks _transfer_kernel!(
                scratch_k_d, pk_table_d,
                log_k_min_f, d_log_k_f, Int32(pk_table_n),
                Float64(dk), Int32(n2),
                Int32(kernel_fn_id), Int32(dim1), Int32(dim2),
            )
            @cuda threads=threads blocks=plane_blocks _hermitize_planes_kernel!(
                scratch_k_d, Int32(n2), _I1, Int32(n2 ÷ 2 + 1))
            result_d = _cufft_irfft(scratch_k_d, n2)
            CUDA.unsafe_free!(scratch_k_d)
            tile_d = result_d[s1:s2, s1:s2, s1:s2]
            CUDA.unsafe_free!(result_d)
            out[ki] = Array(tile_d)
            CUDA.unsafe_free!(tile_d)
        end
        CUDA.unsafe_free!(base_k_d)
        CUDA.unsafe_free!(pk_table_d)
        return out
    end
end

# ===================================================================
# Tile-local periodic k-space kernel (2LPT + Laplacian)
# ===================================================================
# Matches `MultiResolution._apply_kernel_inplace!` — applies a kernel
# multiplier to each mode of a rfft buffer, and zeros the DC origin
# plus the Nyquist row/plane (ix=nk, iy=nyq+1, iz=nyq+1). No P(k).
#
# kernel_fn_id (subset of isolated-convolve IDs, no √P(k) factor):
#   2 = -i * k_{dim1} / k²               (2LPT displacement)
#   3 = -k_{dim1} * k_{dim2} / k²        (φ_ij)
#   4 = k²                               (Laplacian)
function _periodic_kernel!(arr_k, dk::Float64, n::Int32,
                            kernel_fn_id::Int32, dim1::Int32, dim2::Int32)
    nk = (n >> _I1) + _I1      # n/2 + 1 (rfft first-dim length)
    nyq = n >> _I1             # n/2
    total = Int64(nk) * Int64(n) * Int64(n)
    idx0 = Int64(blockIdx().x - _I1) * Int64(blockDim().x) + Int64(threadIdx().x)
    if idx0 > total
        return
    end
    li = idx0 - Int64(1)
    ix = Int32(li % Int64(nk)) + _I1
    rem = li ÷ Int64(nk)
    iy = Int32(rem % Int64(n)) + _I1
    iz = Int32(rem ÷ Int64(n)) + _I1

    # Signed frequency indices (fftfreq convention)
    iy_signed = iy <= (nyq + _I1) ? (iy - _I1) : (iy - n - _I1)
    iz_signed = iz <= (nyq + _I1) ? (iz - _I1) : (iz - n - _I1)
    kx = Float64(ix - _I1) * dk
    ky = Float64(iy_signed) * dk
    kz = Float64(iz_signed) * dk
    k2 = kx * kx + ky * ky + kz * kz

    # CPU zeroing: DC origin, Nyquist rows/plane
    if k2 == 0.0 || ix == nk || iy == nyq + _I1 || iz == nyq + _I1
        @inbounds arr_k[ix, iy, iz] = ComplexF32(0)
        return
    end

    @inbounds val = arr_k[ix, iy, iz]
    if kernel_fn_id == Int32(2)
        # 2LPT: -i * ki / k²
        ki = (dim1 == _I1) ? kx : ((dim1 == Int32(2)) ? ky : kz)
        coeff = ki / k2
        a = Float64(real(val)); b = Float64(imag(val))
        @inbounds arr_k[ix, iy, iz] = ComplexF32(b * coeff, -a * coeff)
    elseif kernel_fn_id == Int32(3)
        # φ_ij: -ki * kj / k² (real)
        ki = (dim1 == _I1) ? kx : ((dim1 == Int32(2)) ? ky : kz)
        kj = (dim2 == _I1) ? kx : ((dim2 == Int32(2)) ? ky : kz)
        coeff = Float32(-ki * kj / k2)
        @inbounds arr_k[ix, iy, iz] = val * coeff
    elseif kernel_fn_id == Int32(4)
        # Laplacian: k²
        coeff = Float32(k2)
        @inbounds arr_k[ix, iy, iz] = val * coeff
    end
    return
end

function PeakPatch.compute_2lpt_gpu(delta_tile::AbstractArray{<:Real,3},
                                      nmesh::Integer, boxsize_local::Real;
                                      threads::Integer=256,
                                      return_device::Bool=false)
    n = Int(nmesh)
    @assert size(delta_tile) == (n, n, n)
    dk = 2π / Float64(boxsize_local)

    # Accept either host Array{Float32,3} or device CuArray{Float32,3}. If we
    # already have a CuArray, reuse it (don't free — caller owns it); otherwise
    # upload a fresh copy that we free after FFT.
    delta_d_is_input = delta_tile isa CuArray{Float32,3}
    delta_d = delta_d_is_input ? delta_tile : CuArray{Float32}(delta_tile)
    delta_k_d = _cufft_rfft(delta_d)

    # src2 = 0.5 * δ²
    src2_d = delta_d .^ 2 .* 0.5f0
    if !delta_d_is_input
        CUDA.unsafe_free!(delta_d)
    end

    nk = size(delta_k_d, 1)
    total = nk * n * n
    blocks = cld(total, threads)

    # Six φ_ij terms: diagonal (d,d) with coef 0.5, off-diag with coef 1.0.
    for (di, dj, coef) in ((1, 1, 0.5f0), (2, 2, 0.5f0), (3, 3, 0.5f0),
                            (1, 2, 1.0f0), (1, 3, 1.0f0), (2, 3, 1.0f0))
        phi_k_d = copy(delta_k_d)
        @cuda threads=threads blocks=blocks _periodic_kernel!(
            phi_k_d, Float64(dk), Int32(n),
            Int32(3), Int32(di), Int32(dj))
        phi_r_d = _cufft_irfft(phi_k_d, n)
        CUDA.unsafe_free!(phi_k_d)
        src2_d .-= phi_r_d .^ 2 .* coef
        CUDA.unsafe_free!(phi_r_d)
    end
    CUDA.unsafe_free!(delta_k_d)

    # Forward FFT of src2, then ψ₂ for each dimension
    src2_k_d = _cufft_rfft(src2_d)
    CUDA.unsafe_free!(src2_d)

    if return_device
        psi2_dev = Vector{CuArray{Float32,3}}(undef, 3)
        for dim in 1:3
            psi2_k_d = copy(src2_k_d)
            @cuda threads=threads blocks=blocks _periodic_kernel!(
                psi2_k_d, Float64(dk), Int32(n),
                Int32(2), Int32(dim), Int32(0))
            psi2_dev[dim] = _cufft_irfft(psi2_k_d, n)
            CUDA.unsafe_free!(psi2_k_d)
        end
        CUDA.unsafe_free!(src2_k_d)
        return psi2_dev
    else
        psi2_host = Vector{Array{Float32,3}}(undef, 3)
        for dim in 1:3
            psi2_k_d = copy(src2_k_d)
            @cuda threads=threads blocks=blocks _periodic_kernel!(
                psi2_k_d, Float64(dk), Int32(n),
                Int32(2), Int32(dim), Int32(0))
            psi2_d = _cufft_irfft(psi2_k_d, n)
            CUDA.unsafe_free!(psi2_k_d)
            psi2_host[dim] = Array(psi2_d)
            CUDA.unsafe_free!(psi2_d)
        end
        CUDA.unsafe_free!(src2_k_d)
        return psi2_host
    end
end

function PeakPatch.compute_laplacian_gpu(delta_tile::AbstractArray{<:Real,3},
                                          nmesh::Integer, boxsize_local::Real;
                                          threads::Integer=256,
                                          return_device::Bool=false)
    n = Int(nmesh)
    @assert size(delta_tile) == (n, n, n)
    dk = 2π / Float64(boxsize_local)

    delta_d_is_input = delta_tile isa CuArray{Float32,3}
    delta_d = delta_d_is_input ? delta_tile : CuArray{Float32}(delta_tile)
    delta_k_d = _cufft_rfft(delta_d)
    if !delta_d_is_input
        CUDA.unsafe_free!(delta_d)
    end

    nk = size(delta_k_d, 1)
    total = nk * n * n
    blocks = cld(total, threads)
    @cuda threads=threads blocks=blocks _periodic_kernel!(
        delta_k_d, Float64(dk), Int32(n),
        Int32(4), Int32(0), Int32(0))

    lapd_d = _cufft_irfft(delta_k_d, n)
    CUDA.unsafe_free!(delta_k_d)
    if return_device
        return lapd_d
    else
        out = Array(lapd_d)
        CUDA.unsafe_free!(lapd_d)
        return out
    end
end

# ===================================================================
# Peak finding: smooth_field + find_peaks on GPU
# ===================================================================
# Matches `Filters.smooth_field` + `PeakFind.find_peaks` for the periodic
# tile-local rfft. CPU `smooth_field` does NOT zero Nyquist planes here
# (unlike `_apply_kernel_inplace!`) because the window is real-valued and
# Hermitian symmetry is preserved automatically.
#
# wsmooth: 0 = Gaussian (Fortran-convention, kR/2), 1 = top-hat,
#          3 = k²-weighted Gaussian (lapd diagnostic).
function _smooth_window_kernel!(arr_k, dk::Float64, Rf::Float64, n::Int32,
                                  wsmooth::Int32)
    nk = (n >> _I1) + _I1
    nyq = n >> _I1
    total = Int64(nk) * Int64(n) * Int64(n)
    idx0 = Int64(blockIdx().x - _I1) * Int64(blockDim().x) + Int64(threadIdx().x)
    if idx0 > total
        return
    end
    li = idx0 - Int64(1)
    ix = Int32(li % Int64(nk)) + _I1
    rem = li ÷ Int64(nk)
    iy = Int32(rem % Int64(n)) + _I1
    iz = Int32(rem ÷ Int64(n)) + _I1

    iy_signed = iy <= (nyq + _I1) ? (iy - _I1) : (iy - n - _I1)
    iz_signed = iz <= (nyq + _I1) ? (iz - _I1) : (iz - n - _I1)
    kx = Float64(ix - _I1) * dk
    ky = Float64(iy_signed) * dk
    kz = Float64(iz_signed) * dk
    k2 = kx * kx + ky * ky + kz * kz
    # k=0 mode: leave unchanged (matches CPU `fkR == 0 && continue`)
    if k2 == 0.0
        return
    end
    k = sqrt(k2)
    kR = k * Rf
    w = 0.0
    if wsmooth == _I0
        # Gaussian (Fortran convention): W = exp(-kR²/8)  [= exp(-(kR/2)²/2)]
        w = exp(-kR * kR * 0.125)
    elseif wsmooth == Int32(3)
        # k²-weighted Gaussian (lapd diagnostic)
        w = exp(-kR * kR * 0.125) * k2
    else
        # wsmooth == 1 → top-hat
        w = 3.0 * (sin(kR) - kR * cos(kR)) / (kR * kR * kR)
    end
    @inbounds val = arr_k[ix, iy, iz]
    @inbounds arr_k[ix, iy, iz] = val * Float32(w)
    return
end

# Local-maximum finder. Thread per cell in the inner cube
# (nbuff+1 : n-nbuff)³. Atomic counter reserves output slots.
function _find_peaks_kernel!(peak_count, peak_i, peak_j, peak_k,
                               peak_delta, peak_lapd, peak_filter_idx,
                               mask_arr, delta_s, lapd_s,
                               n::Int32, nbuff::Int32, fcrit::Float32,
                               has_lapd::Int32, filter_idx::Int32,
                               max_peaks::Int32)
    inner_n = n - Int32(2) * nbuff
    total = Int64(inner_n) * Int64(inner_n) * Int64(inner_n)
    idx0 = Int64(blockIdx().x - _I1) * Int64(blockDim().x) + Int64(threadIdx().x)
    if idx0 > total
        return
    end
    li = idx0 - Int64(1)
    i_inner = Int32(li % Int64(inner_n)) + _I1
    rem = li ÷ Int64(inner_n)
    j_inner = Int32(rem % Int64(inner_n)) + _I1
    k_inner = Int32(rem ÷ Int64(inner_n)) + _I1
    i = i_inner + nbuff
    j = j_inner + nbuff
    k = k_inner + nbuff

    @inbounds if mask_arr[i, j, k] == Int8(1)
        return
    end
    @inbounds ff = delta_s[i, j, k]
    if ff < fcrit
        return
    end

    # Strict local maximum over the 3×3×3 stencil. Because `ff > delta[i,j,k]`
    # is never true (they're the same element), checking `<` over the full
    # 27 cells — including the centre — is equivalent to the CPU version.
    is_max = true
    for kk in (k - _I1):(k + _I1)
        for jj in (j - _I1):(j + _I1)
            for ii in (i - _I1):(i + _I1)
                @inbounds if ff < delta_s[ii, jj, kk]
                    is_max = false
                end
            end
        end
    end
    if !is_max
        return
    end

    # Atomically claim a slot (1-based). peak_count[1] holds the running total.
    old = CUDA.atomic_add!(CUDA.pointer(peak_count, 1), Int32(1))
    slot = old + _I1
    if slot > max_peaks
        return
    end
    @inbounds mask_arr[i, j, k] = Int8(1)
    @inbounds peak_i[slot] = i
    @inbounds peak_j[slot] = j
    @inbounds peak_k[slot] = k
    @inbounds peak_delta[slot] = ff
    @inbounds peak_filter_idx[slot] = filter_idx
    if has_lapd == _I1
        @inbounds peak_lapd[slot] = lapd_s[i, j, k]
    end
    return
end

function PeakPatch.peak_find_tile_gpu(delta_tile::AbstractArray{<:Real,3},
                                        filters::AbstractVector,
                                        fcrits::AbstractVector{<:Real},
                                        tile_masks_in::AbstractArray{<:Integer,3},
                                        xbx::Real, ybx::Real, zbx::Real,
                                        alatt::Real, nbuff::Integer,
                                        wsmooth::Integer, ioutshear::Integer;
                                        threads::Integer=256,
                                        max_peaks::Integer=2_000_000)
    n = Int(size(delta_tile, 1))
    @assert size(delta_tile) == (n, n, n)
    @assert size(tile_masks_in) == (n, n, n)
    @assert length(filters) == length(fcrits)
    dk = 2π / (n * Float64(alatt))    # periodic: boxsize_local = n * alatt
    has_lapd = ioutshear >= 1 ? 1 : 0

    # Upload δ if host, reuse if already on device. Caller owns the input
    # CuArray (we never free it). FFT once, keep k-space buffer on device
    # across the filter loop.
    delta_d_is_input = delta_tile isa CuArray{Float32,3}
    delta_d = delta_d_is_input ? delta_tile : CuArray{Float32}(delta_tile)
    delta_k_d = _cufft_rfft(delta_d)
    if !delta_d_is_input
        CUDA.unsafe_free!(delta_d)
    end

    # Mask on device (updated across filters)
    mask_d = CuArray(Int8.(tile_masks_in))

    # Pre-allocate output buffers sized for worst case
    peak_count_d = CUDA.zeros(Int32, 1)
    peak_i_d     = CUDA.zeros(Int32, max_peaks)
    peak_j_d     = CUDA.zeros(Int32, max_peaks)
    peak_k_d     = CUDA.zeros(Int32, max_peaks)
    peak_delta_d = CUDA.zeros(Float32, max_peaks)
    peak_filter_d = CUDA.zeros(Int32, max_peaks)
    peak_lapd_d  = has_lapd == 1 ? CUDA.zeros(Float32, max_peaks) :
                                    CUDA.zeros(Float32, 1)   # dummy, unused

    nk_freq = size(delta_k_d, 1)
    smooth_total = nk_freq * n * n
    smooth_blocks = cld(smooth_total, threads)

    inner_n = n - 2 * Int(nbuff)
    find_total = inner_n * inner_n * inner_n
    find_blocks = cld(find_total, threads)

    smoothed_k_d = similar(delta_k_d)    # scratch for k-space window
    lapd_k_d     = has_lapd == 1 ? similar(delta_k_d) : delta_k_d  # alias if unused

    for (fi, (filter, fcrit)) in enumerate(zip(filters, fcrits))
        Rf = Float64(filter[3])
        # ---- smooth_field: window → irfft ----
        copyto!(smoothed_k_d, delta_k_d)
        @cuda threads=threads blocks=smooth_blocks _smooth_window_kernel!(
            smoothed_k_d, dk, Rf, Int32(n), Int32(wsmooth))
        delta_s_d = _cufft_irfft(smoothed_k_d, n)

        # ---- optional wsmooth=3 field for lapd diagnostic ----
        lapd_s_d = if has_lapd == 1
            copyto!(lapd_k_d, delta_k_d)
            @cuda threads=threads blocks=smooth_blocks _smooth_window_kernel!(
                lapd_k_d, dk, Rf, Int32(n), Int32(3))
            _cufft_irfft(lapd_k_d, n)
        else
            delta_s_d  # dummy, kernel won't read
        end

        # ---- find_peaks ----
        @cuda threads=threads blocks=find_blocks _find_peaks_kernel!(
            peak_count_d, peak_i_d, peak_j_d, peak_k_d,
            peak_delta_d, peak_lapd_d, peak_filter_d,
            mask_d, delta_s_d, lapd_s_d,
            Int32(n), Int32(nbuff), Float32(fcrit),
            Int32(has_lapd), Int32(fi), Int32(max_peaks))

        CUDA.unsafe_free!(delta_s_d)
        if has_lapd == 1
            CUDA.unsafe_free!(lapd_s_d)
        end
    end

    CUDA.unsafe_free!(smoothed_k_d)
    if has_lapd == 1
        CUDA.unsafe_free!(lapd_k_d)
    end
    CUDA.unsafe_free!(delta_k_d)

    # Download results
    npeaks = Int(Array(peak_count_d)[1])
    npeaks = min(npeaks, max_peaks)
    peak_i_h     = Array(peak_i_d[1:npeaks])
    peak_j_h     = Array(peak_j_d[1:npeaks])
    peak_k_h     = Array(peak_k_d[1:npeaks])
    peak_delta_h = Array(peak_delta_d[1:npeaks])
    peak_filter_h = Array(peak_filter_d[1:npeaks])
    peak_lapd_h  = has_lapd == 1 ? Array(peak_lapd_d[1:npeaks]) : Float32[]
    mask_out = Array(mask_d)

    CUDA.unsafe_free!(peak_count_d)
    CUDA.unsafe_free!(peak_i_d)
    CUDA.unsafe_free!(peak_j_d)
    CUDA.unsafe_free!(peak_k_d)
    CUDA.unsafe_free!(peak_delta_d)
    CUDA.unsafe_free!(peak_filter_d)
    CUDA.unsafe_free!(peak_lapd_d)
    CUDA.unsafe_free!(mask_d)

    # Assemble PeakCandidate vector on host (match CPU coordinate computation)
    cen1 = 0.5 * (n + 1)
    peaks = Vector{PeakPatch.PeakCandidate}(undef, npeaks)
    Rf_out  = Vector{Float64}(undef, npeaks)
    FcRf_out = Vector{Float32}(undef, npeaks)
    d2Rf_out = Vector{Float32}(undef, npeaks)
    for p in 1:npeaks
        i = Int(peak_i_h[p]); j = Int(peak_j_h[p]); k = Int(peak_k_h[p])
        xc = Float64(xbx) + Float64(alatt) * (i - cen1)
        yc = Float64(ybx) + Float64(alatt) * (j - cen1)
        zc = Float64(zbx) + Float64(alatt) * (k - cen1)
        ipp = i + (j - 1) * n + (k - 1) * n * n
        fi_ = Int(peak_filter_h[p])
        Rf_out[p] = Float64(filters[fi_][3])
        FcRf_out[p] = peak_delta_h[p]
        d2Rf_out[p] = has_lapd == 1 ? peak_lapd_h[p] : 0.0f0
        peaks[p] = PeakPatch.PeakCandidate{Float32}(
            i, j, k, ipp, xc, yc, zc, peak_delta_h[p], Rf_out[p])
    end

    return (peaks=peaks, Rf=Rf_out, FcRf=FcRf_out, d2Rf=d2Rf_out, mask=mask_out)
end

# ===================================================================
# Coarse → fine tile tri-cubic (Catmull-Rom) interpolation
# ===================================================================
# GPU port of `MultiResolution._interpolate_to_tile`. Matches the CPU
# output element-for-element at Float32 precision.

@inline function _catmull_rom_f32(t::Float32)
    t2 = t * t; t3 = t2 * t
    w0 = -0.5f0 * t3 + t2 - 0.5f0 * t
    w1 =  1.5f0 * t3 - 2.5f0 * t2 + 1.0f0
    w2 = -1.5f0 * t3 + 2.0f0 * t2 + 0.5f0 * t
    w3 =  0.5f0 * t3 - 0.5f0 * t2
    return (w0, w1, w2, w3)
end

@inline function _mod1_i32(i::Int32, M::Int32)
    # Equivalent to Julia's mod1(i, M): result in 1..M for any integer i.
    r = (i - _I1) % M
    if r < _I0
        r += M
    end
    return r + _I1
end

function _tricubic_kernel!(tile_out, coarse, i0::Int32, j0::Int32, k0::Int32,
                             block::Int32, nmesh::Int32, N::Int32, M::Int32)
    total = Int64(nmesh) * Int64(nmesh) * Int64(nmesh)
    idx0 = Int64(blockIdx().x - _I1) * Int64(blockDim().x) + Int64(threadIdx().x)
    if idx0 > total
        return
    end
    li_ = idx0 - Int64(1)
    li = Int32(li_ % Int64(nmesh)) + _I1
    rem = li_ ÷ Int64(nmesh)
    lj = Int32(rem % Int64(nmesh)) + _I1
    lk = Int32(rem ÷ Int64(nmesh)) + _I1

    # Global fine-grid index (wrapped) → fractional coarse-grid coordinate.
    # CPU: gi = mod1(i0 + li - 1, N); fx = (gi - 0.5)/block + 0.5
    gi = _mod1_i32(i0 + li - _I1, N)
    gj = _mod1_i32(j0 + lj - _I1, N)
    gk = _mod1_i32(k0 + lk - _I1, N)

    # Use Float64 for fractional coord (matches CPU where fx is Float64)
    fx = (Float64(gi) - 0.5) / Float64(block) + 0.5
    fy = (Float64(gj) - 0.5) / Float64(block) + 0.5
    fz = (Float64(gk) - 0.5) / Float64(block) + 0.5

    ix = Int32(floor(fx)); dx_ = Float32(fx - Float64(ix))
    iy = Int32(floor(fy)); dy_ = Float32(fy - Float64(iy))
    iz = Int32(floor(fz)); dz_ = Float32(fz - Float64(iz))

    (wx0, wx1, wx2, wx3) = _catmull_rom_f32(dx_)
    (wy0, wy1, wy2, wy3) = _catmull_rom_f32(dy_)
    (wz0, wz1, wz2, wz3) = _catmull_rom_f32(dz_)

    # Precompute wrapped indices for the 4-wide stencils along each axis
    # (avoids 4 calls inside the 64-op inner loop).
    ix_arr = (_mod1_i32(ix - _I1, M), _mod1_i32(ix, M),
              _mod1_i32(ix + _I1, M), _mod1_i32(ix + Int32(2), M))
    iy_arr = (_mod1_i32(iy - _I1, M), _mod1_i32(iy, M),
              _mod1_i32(iy + _I1, M), _mod1_i32(iy + Int32(2), M))
    iz_arr = (_mod1_i32(iz - _I1, M), _mod1_i32(iz, M),
              _mod1_i32(iz + _I1, M), _mod1_i32(iz + Int32(2), M))
    wx = (wx0, wx1, wx2, wx3)
    wy = (wy0, wy1, wy2, wy3)
    wz = (wz0, wz1, wz2, wz3)

    val = 0.0f0
    @inbounds for kk_idx in 1:4
        kk = iz_arr[kk_idx]
        wk = wz[kk_idx]
        for jj_idx in 1:4
            jj = iy_arr[jj_idx]
            wj = wy[jj_idx]
            for ii_idx in 1:4
                ii = ix_arr[ii_idx]
                wi = wx[ii_idx]
                val += wi * wj * wk * coarse[ii, jj, kk]
            end
        end
    end
    @inbounds tile_out[li, lj, lk] = val
    return
end

function PeakPatch.interpolate_to_tile_gpu(coarse::AbstractArray{<:Real,3},
                                             it::Integer, jt::Integer, kt::Integer,
                                             nsub::Integer, nmesh::Integer,
                                             N::Integer, M::Integer;
                                             threads::Integer=256,
                                             return_device::Bool=false)
    @assert size(coarse) == (M, M, M)
    block = Int(N) ÷ Int(M)
    nbuff = (Int(nmesh) - Int(nsub)) ÷ 2
    i0 = (Int(it) - 1) * Int(nsub) + 1 - nbuff
    j0 = (Int(jt) - 1) * Int(nsub) + 1 - nbuff
    k0 = (Int(kt) - 1) * Int(nsub) + 1 - nbuff

    coarse_d_is_input = coarse isa CuArray{Float32,3}
    coarse_d = coarse_d_is_input ? coarse : CuArray{Float32}(coarse)
    tile_d = CUDA.zeros(Float32, Int(nmesh), Int(nmesh), Int(nmesh))

    total = Int(nmesh)^3
    blocks = cld(total, threads)
    @cuda threads=threads blocks=blocks _tricubic_kernel!(
        tile_d, coarse_d, Int32(i0), Int32(j0), Int32(k0),
        Int32(block), Int32(nmesh), Int32(N), Int32(M))

    if !coarse_d_is_input
        CUDA.unsafe_free!(coarse_d)
    end
    if return_device
        return tile_d
    else
        out = Array(tile_d)
        CUDA.unsafe_free!(tile_d)
        return out
    end
end

# Multi-GPU dispatcher helper: bind the current task to device `id` (0-based).
# Implementation of the stub declared in `src/PeakPatch.jl`.
function PeakPatch.set_cuda_device!(id::Int)
    CUDA.device!(id)
    return nothing
end

end # module CUDAExt
