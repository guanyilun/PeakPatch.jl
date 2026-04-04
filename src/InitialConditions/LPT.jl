module LPT

using FFTW

"""
    displacements_1lpt(delta_k, n, boxsize)

Compute 1LPT (Zel'dovich) displacement fields from the density field in
Fourier space.

# Arguments
- `delta_k`: rfft-transformed density field, shape `(n÷2+1, n, n)`
- `n`: grid size (assumed cubic)
- `boxsize`: side length in Mpc/h

# Returns
- `(psi_x, psi_y, psi_z)`: three `n×n×n` arrays of displacements (same precision as input)

Formula: ψᵢ(k) = im * kᵢ / k² * δ̃(k)
"""
function displacements_1lpt(delta_k, n::Int, boxsize::Real)
    dk = 2π / Float64(boxsize)
    nk = size(delta_k, 1)  # n÷2+1
    CT = eltype(delta_k)   # Complex{Float32} or Complex{Float64}

    kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
    ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
    kz_arr = Float64.(FFTW.fftfreq(n, n * dk))

    psi_x_k = similar(delta_k, CT)
    psi_y_k = similar(delta_k, CT)
    psi_z_k = similar(delta_k, CT)

    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = kx_arr[ix]
        ky = ky_arr[iy]
        kz = kz_arr[iz]
        k2 = kx^2 + ky^2 + kz^2

        if k2 == 0.0
            psi_x_k[ix, iy, iz] = 0
            psi_y_k[ix, iy, iz] = 0
            psi_z_k[ix, iy, iz] = 0
            continue
        end

        d = delta_k[ix, iy, iz]

        # Zero Nyquist modes on each axis to match Fortran behaviour:
        # ψ at k_Nyq cannot be properly encoded by the rfft (ambiguous ±k_Nyq).
        nyq = n ÷ 2 + 1
        if ix == nk || iy == nyq || iz == nyq
            psi_x_k[ix, iy, iz] = 0
            psi_y_k[ix, iy, iz] = 0
            psi_z_k[ix, iy, iz] = 0
            continue
        end

        psi_x_k[ix, iy, iz] = im * kx / k2 * d
        psi_y_k[ix, iy, iz] = im * ky / k2 * d
        psi_z_k[ix, iy, iz] = im * kz / k2 * d
    end

    dummy_k = similar(psi_x_k)
    invplan = plan_irfft(dummy_k, n; flags=FFTW.MEASURE)
    return (
        invplan * psi_x_k,
        invplan * psi_y_k,
        invplan * psi_z_k,
    )
end

"""
    displacements_2lpt(delta_k, n, boxsize)

Compute 2LPT displacement fields from the density field in Fourier space.

# Returns
- `(psi2_x, psi2_y, psi2_z)`: second-order displacement fields

Algorithm:
1. Compute φᵢⱼ(k) = -(kᵢ kⱼ / k²) δ̃(k) → irfft → real-space φᵢⱼ
2. Form source: s = φ₁₁ φ₂₂ - φ₁₂² + φ₁₁ φ₃₃ - φ₁₃² + φ₂₂ φ₃₃ - φ₂₃²
3. Gradient of s in k-space gives ψ⁽²⁾
"""
function displacements_2lpt(delta_k, n::Int, boxsize::Real)
    dk = 2π / Float64(boxsize)
    nk = size(delta_k, 1)
    CT = eltype(delta_k)

    kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
    ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
    kz_arr = Float64.(FFTW.fftfreq(n, n * dk))

    # Compute phi_ij in k-space and transform to real space
    phi11_k = similar(delta_k, CT)
    phi22_k = similar(delta_k, CT)
    phi33_k = similar(delta_k, CT)
    phi12_k = similar(delta_k, CT)
    phi13_k = similar(delta_k, CT)
    phi23_k = similar(delta_k, CT)

    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = kx_arr[ix]
        ky = ky_arr[iy]
        kz = kz_arr[iz]
        k2 = kx^2 + ky^2 + kz^2

        if k2 == 0.0
            phi11_k[ix, iy, iz] = 0
            phi22_k[ix, iy, iz] = 0
            phi33_k[ix, iy, iz] = 0
            phi12_k[ix, iy, iz] = 0
            phi13_k[ix, iy, iz] = 0
            phi23_k[ix, iy, iz] = 0
            continue
        end

        d = delta_k[ix, iy, iz]
        phi11_k[ix, iy, iz] = -kx * kx / k2 * d
        phi22_k[ix, iy, iz] = -ky * ky / k2 * d
        phi33_k[ix, iy, iz] = -kz * kz / k2 * d
        phi12_k[ix, iy, iz] = -kx * ky / k2 * d
        phi13_k[ix, iy, iz] = -kx * kz / k2 * d
        phi23_k[ix, iy, iz] = -ky * kz / k2 * d
    end

    dummy_k = similar(phi11_k)
    invplan = plan_irfft(dummy_k, n; flags=FFTW.MEASURE)

    phi11 = invplan * phi11_k
    phi22 = invplan * phi22_k
    phi33 = invplan * phi33_k
    phi12 = invplan * phi12_k
    phi13 = invplan * phi13_k
    phi23 = invplan * phi23_k

    # Source term: s = det(φ_{ij}) in 3D (all 3 eigenvalue products)
    src2 = phi11 .* phi22 .- phi12 .^ 2 .+
           phi11 .* phi33 .- phi13 .^ 2 .+
           phi22 .* phi33 .- phi23 .^ 2

    # FFT source to k-space
    dummy_src = similar(src2)
    fwdplan = plan_rfft(dummy_src; flags=FFTW.MEASURE)
    src2_k = fwdplan * src2

    # Gradient of source → 2LPT displacement
    psi2_x_k = similar(src2_k, CT)
    psi2_y_k = similar(src2_k, CT)
    psi2_z_k = similar(src2_k, CT)

    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = kx_arr[ix]
        ky = ky_arr[iy]
        kz = kz_arr[iz]
        k2 = kx^2 + ky^2 + kz^2

        if k2 == 0.0
            psi2_x_k[ix, iy, iz] = 0
            psi2_y_k[ix, iy, iz] = 0
            psi2_z_k[ix, iy, iz] = 0
            continue
        end

        # Zero Nyquist modes (same as 1LPT, matching Fortran)
        nyq = n ÷ 2 + 1
        if ix == nk || iy == nyq || iz == nyq
            psi2_x_k[ix, iy, iz] = 0
            psi2_y_k[ix, iy, iz] = 0
            psi2_z_k[ix, iy, iz] = 0
            continue
        end

        s = src2_k[ix, iy, iz]
        # Sign is NEGATIVE: Fortran uses -im for 2LPT gradients.
        # The downstream factor (-3/7 Ω_m^{-1/143}) adds another sign flip,
        # giving net +3/7 — which would be wrong with +im here.
        psi2_x_k[ix, iy, iz] = -im * kx / k2 * s
        psi2_y_k[ix, iy, iz] = -im * ky / k2 * s
        psi2_z_k[ix, iy, iz] = -im * kz / k2 * s
    end

    dummy_k2 = similar(psi2_x_k)
    invplan2 = plan_irfft(dummy_k2, n; flags=FFTW.MEASURE)
    return (
        invplan2 * psi2_x_k,
        invplan2 * psi2_y_k,
        invplan2 * psi2_z_k,
    )
end

end # module LPT
