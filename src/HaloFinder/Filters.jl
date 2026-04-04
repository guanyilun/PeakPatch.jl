module Filters

# Window functions in Fourier space. All normalized to W(0) = 1.

"""
    gaussian_window(kR)

Gaussian window in k-space: W(kR) = exp(-k²R²/2).

# Fortran R/2 convention
The Fortran `smooth_field` passes `kR/2` to this function, not `kR`.  
This means the effective real-space Gaussian has σ = R/2, not R.  
Callers must pass `k * R/2` (or wrap with `gaussian_window_fortran`) to
match the Fortran output exactly:

    gaussian_window_fortran(kR) = gaussian_window(kR / 2)   # = exp(-k²R²/8)
"""
gaussian_window(kR::Real) = exp(-kR^2 / 2)

"""
    gaussian_window_fortran(kR)

Gaussian window matching the Fortran `smooth_field` convention, which
passes `kR/2` internally: W(kR) = exp(-k²R²/8).
Use this to reproduce Fortran-smoothed fields exactly.
"""
gaussian_window_fortran(kR::Real) = gaussian_window(kR / 2)

"""Top-hat (spherical) window: W(k) = 3(sin(kR) - kR cos(kR)) / (kR)³."""
function tophat_window(kR::Real)
    if kR == 0.0
        return 1.0
    end
    return 3.0 * (sin(kR) - kR * cos(kR)) / kR^3
end

"""Read a filter bank file (ASCII, Fortran formatted).

Format:
  Line 1: nic  (number of filter scales)
  Lines 2..nic+1: id fcrit Rf

Returns a vector of (id::Int, fcrit::Float64, Rf::Float64) tuples.
"""
function read_filterbank(path::String)
    lines = readlines(path)
    nic = parse(Int, strip(lines[1]))
    result = Vector{Tuple{Int, Float64, Float64}}(undef, nic)
    for i in 1:nic
        parts = split(strip(lines[i + 1]))
        id_ = parse(Int, parts[1])
        fcrit = parse(Float64, parts[2])
        Rf = parse(Float64, parts[3])
        result[i] = (id_, fcrit, Rf)
    end
    return result
end

using FFTW

"""Smooth a field in Fourier space.

Window types (`wsmooth`):
- `0`: Gaussian (Fortran R/2 convention)
- `1`: top-hat (spherical)
- `3`: SIGMA_2 — k²-weighted Gaussian (for Laplacian diagnostic)

Arguments:
- `delta_k`: rfft of the density field (ComplexF32 or ComplexF64, not mutated)
- `n`: grid size (assumed cubic)
- `boxsize`: side length in Mpc/h
- `Rf`: smoothing radius
- `wsmooth`: window type

Returns an `n×n×n` Float32 real-space array.
"""
function smooth_field(delta_k, n::Int, boxsize::Float64, Rf::Float64, wsmooth::Int;
                      fortran_compat::Bool=false)
    T = fortran_compat ? Float32 : Float64
    dk = T(2π / boxsize)
    kx = dk .* T.(FFTW.rfftfreq(n, n))
    ky = dk .* T.(FFTW.fftfreq(n, n))
    kz = dk .* T.(FFTW.fftfreq(n, n))
    Rf_use = T(Rf)

    smoothed_k = copy(delta_k)
    for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2)
        fkR = k * Rf_use
        fkR == zero(T) && continue
        if wsmooth == 0
            w = gaussian_window_fortran(fkR)
        elseif wsmooth == 3
            w = gaussian_window_fortran(fkR) * k^2
        else
            w = tophat_window(fkR)
        end
        smoothed_k[ix, iy, iz] *= Float32(w)
    end
    return irfft(smoothed_k, n)
end

end # module Filters
