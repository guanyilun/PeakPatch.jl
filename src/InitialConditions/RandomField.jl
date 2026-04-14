module RandomField

using FFTW
using Random123: threefry
import ..LCG

"""
    generate_grf_lcg(n, pk, boxsize, seed; T=Float32)

Same as `generate_grf` but uses the Fortran 48-bit multiplicative LCG
for white-noise generation, producing bit-for-bit identical ICs to the
Fortran hpkvd code when given the same seed.

The Fortran seeds via `rans(1, seed)` which for single-task gives
`seed_state = (abs(seed), 0, 0, 0)`, then fills `delta(i,j,k)` with
`gaussdev()` in column-major (i-innermost) order.
"""
function generate_grf_lcg(n::Int, pk, boxsize::Real, seed::Integer; T::Type=Float32,
                          fortran_compat::Bool=false)
    # --- White noise via Fortran LCG ---
    # rans(1, seed): for N=1 (odd), nn=1, returns seed_state=(abs(seed),0,0,0)
    seeds = LCG.rans(1, seed)
    lcg_seed = seeds[1]
    state = LCG.GaussState()

    noise = Array{T,3}(undef, n, n, n)
    for idx in eachindex(noise)  # column-major: i varies fastest
        g, lcg_seed = LCG.gaussdev(lcg_seed, state)
        noise[idx] = g
    end

    return _convolve_noise!(noise, pk, boxsize, n; fortran_compat=fortran_compat)
end

"""
    _threefry_gaussian(seed, pair_idx) -> (g1, g2)

Generate a pair of standard-normal deviates from a single Threefry2x call.
Uses the Box-Muller transform on two uniform deviates derived from the
counter-based RNG output.

`seed` is the global RNG seed; `pair_idx` is a 0-based counter that uniquely
identifies which pair of Gaussian deviates to produce.  The output is fully
deterministic given (seed, pair_idx) with no state dependencies.
"""
@inline function _threefry_gaussian(seed::UInt64, pair_idx::UInt64)
    u1_raw, u2_raw = threefry((seed, UInt64(0)), (pair_idx, UInt64(0)), Val(20))
    # Convert UInt64 -> Float64 in (0, 1] — exclude exact 0 for log safety
    u1 = (Float64(u1_raw >>> 11) + 1.0) * 1.1102230246251565e-16  # (x+1) * 2^-53
    u2 = (Float64(u2_raw >>> 11) + 1.0) * 1.1102230246251565e-16
    # Box-Muller transform
    r = sqrt(-2.0 * log(u1))
    θ = 2π * u2
    return r * cos(θ), r * sin(θ)
end

"""
    fill_noise_threefry!(noise, n, seed)

Fill a 3D array (or subregion) with Gaussian white noise using the Threefry2x
counter-based RNG.  Each cell's value depends only on `(seed, linear_index)`,
making this trivially parallelizable — any rank can generate any cell
independently with no communication.

Linear index mapping (column-major, 1-based): `idx = i + (j-1)*n + (k-1)*n*n`

For a full `n×n×n` array, call as `fill_noise_threefry!(noise, n, seed)`.
For a distributed sub-region (e.g. PencilArray), use `fill_noise_threefry_region!`.
"""
function fill_noise_threefry!(noise::AbstractArray{T,3}, n::Int, seed::Integer) where T
    s = UInt64(seed)
    n3 = n * n * n
    Threads.@threads for k in 1:n
        @inbounds for j in 1:n, i in 1:n
            idx = i + (j - 1) * n + (k - 1) * n * n  # 1-based linear index
            pair_idx = UInt64((idx - 1) >> 1)          # 0-based pair index
            g1, g2 = _threefry_gaussian(s, pair_idx)
            noise[i, j, k] = T(iseven(idx - 1) ? g1 : g2)
        end
    end
    # If n³ is odd, the last cell uses g1 from the last pair (g2 is discarded)
    return noise
end

"""
    fill_noise_threefry_region!(arr, global_ranges, n, seed)

Fill a sub-region of a distributed 3D field with Threefry Gaussian noise.
`global_ranges` is a tuple `(irange, jrange, krange)` of the global indices
this region covers.  `n` is the full grid dimension (for linear index computation).

Each cell's value is identical to what `fill_noise_threefry!` would produce
for the same global index — decomposition-independent and deterministic.
"""
function fill_noise_threefry_region!(arr, global_ranges, n::Int, seed::Integer)
    s = UInt64(seed)
    irange, jrange, krange = global_ranges
    Threads.@threads for (lk, gk) in collect(enumerate(krange))
        @inbounds for (lj, gj) in enumerate(jrange)
            for (li, gi) in enumerate(irange)
                idx = gi + (gj - 1) * n + (gk - 1) * n * n
                pair_idx = UInt64((idx - 1) >> 1)
                g1, g2 = _threefry_gaussian(s, pair_idx)
                arr[li, lj, lk] = iseven(idx - 1) ? g1 : g2
            end
        end
    end
    return arr
end

"""
    generate_grf(n, pk, boxsize, seed; T=Float32, fortran_compat=false)

Generate a 3D Gaussian random field with power spectrum `pk(k)` on an
`n×n×n` periodic grid of side length `boxsize` (Mpc/h).

Uses the Threefry2x counter-based RNG: each cell's noise value depends only
on `(seed, linear_index)`, making this suitable for distributed generation
where each rank fills its own portion independently.

Returns a 3D `Array{T,3}` with zero mean and variance set by the input P(k).
"""
function generate_grf(n::Int, pk, boxsize::Real, seed::Integer; T::Type=Float32,
                      fortran_compat::Bool=false)
    noise = Array{T,3}(undef, n, n, n)
    fill_noise_threefry!(noise, n, seed)
    return _convolve_noise!(noise, pk, boxsize, n; fortran_compat=fortran_compat)
end

"""
Internal: take white-noise array, FFT → convolve with √P(k) → iFFT.

When `fortran_compat=true`, compute k-values and P(k) amplitudes in Float32
to match Fortran's single-precision arithmetic in `convolve_noise`.
"""
function _convolve_noise!(noise, pk, boxsize, n; fortran_compat::Bool=false)

    # --- Forward FFT ---
    # Plan on a dummy array to avoid FFTW.MEASURE destroying the input
    dummy = similar(noise)
    plan = plan_rfft(dummy; flags=FFTW.MEASURE)
    noise_k = plan * noise

    # --- Convolve with √P(k) in Fourier space ---
    Tk = fortran_compat ? Float32 : Float64  # k-space computation precision
    dk = Tk(2π / boxsize)

    kx_arr = Tk.(FFTW.rfftfreq(n, n * dk))   # [0, dk, 2dk, ..., n/2*dk]
    ky_arr = Tk.(FFTW.fftfreq(n, n * dk))     # [0, dk, ..., n/2*dk, ..., -dk]
    kz_arr = Tk.(FFTW.fftfreq(n, n * dk))

    nk = size(noise_k, 1)  # n÷2 + 1

    Threads.@threads for iz in 1:n
        for iy in 1:n, ix in 1:nk
            kx = kx_arr[ix]
            ky = ky_arr[iy]
            kz = kz_arr[iz]
            k2 = kx^2 + ky^2 + kz^2
            if k2 == zero(Tk)
                noise_k[ix, iy, iz] = 0  # zero DC mode
                continue
            end
            k = sqrt(k2)
            amp = sqrt(Tk(pk(Float64(k))) * dk^3 * Tk(n)^3)
            noise_k[ix, iy, iz] *= amp
        end
    end

    # --- Inverse FFT ---
    dummy_k = similar(noise_k)
    invplan = plan_irfft(dummy_k, n; flags=FFTW.MEASURE)
    field = invplan * noise_k

    return field
end

end # module RandomField
