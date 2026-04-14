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

    Threads.@threads for iz in 1:n
        for iy in 1:n, ix in 1:nk
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

    Threads.@threads for iz in 1:n
        for iy in 1:n, ix in 1:nk
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

    Threads.@threads for iz in 1:n
        for iy in 1:n, ix in 1:nk
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

            nyq = n ÷ 2 + 1
            if ix == nk || iy == nyq || iz == nyq
                psi2_x_k[ix, iy, iz] = 0
                psi2_y_k[ix, iy, iz] = 0
                psi2_z_k[ix, iy, iz] = 0
                continue
            end

            s = src2_k[ix, iy, iz]
            psi2_x_k[ix, iy, iz] = -im * kx / k2 * s
            psi2_y_k[ix, iy, iz] = -im * ky / k2 * s
            psi2_z_k[ix, iy, iz] = -im * kz / k2 * s
        end
    end

    dummy_k2 = similar(psi2_x_k)
    invplan2 = plan_irfft(dummy_k2, n; flags=FFTW.MEASURE)
    return (
        invplan2 * psi2_x_k,
        invplan2 * psi2_y_k,
        invplan2 * psi2_z_k,
    )
end

"""
    displacement_1lpt_component(delta_k, n, boxsize, dim) -> Array

Compute a single 1LPT displacement component (dim=1,2,3 for x,y,z).
Returns the real-space field ψ_dim(x).

Kernel: ψ_dim(k) = im * k_dim / k² * δ̃(k)
"""
function displacement_1lpt_component(delta_k, n::Int, boxsize::Real, dim::Int)
    dk = 2π / Float64(boxsize)
    nk = size(delta_k, 1)
    CT = eltype(delta_k)

    kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
    ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
    kz_arr = Float64.(FFTW.fftfreq(n, n * dk))
    k_arrs = (kx_arr, ky_arr, kz_arr)

    out_k = similar(delta_k, CT)
    nyq = n ÷ 2 + 1

    Threads.@threads for iz in 1:n
        for iy in 1:n, ix in 1:nk
            kx = kx_arr[ix]; ky = ky_arr[iy]; kz = kz_arr[iz]
            k2 = kx^2 + ky^2 + kz^2

            if k2 == 0.0 || ix == nk || iy == nyq || iz == nyq
                out_k[ix, iy, iz] = 0
                continue
            end
            out_k[ix, iy, iz] = im * k_arrs[dim][dim == 1 ? ix : (dim == 2 ? iy : iz)] / k2 * delta_k[ix, iy, iz]
        end
    end

    dummy_k = similar(out_k)
    invplan = plan_irfft(dummy_k, n; flags=FFTW.MEASURE)
    return invplan * out_k
end

"""
    displacement_2lpt_component(src2_k, n, boxsize, dim) -> Array

Compute a single 2LPT displacement component from the source term in k-space.
Returns the real-space field ψ²_dim(x).

Kernel: ψ²_dim(k) = -im * k_dim / k² * s̃(k)
"""
function displacement_2lpt_component(src2_k, n::Int, boxsize::Real, dim::Int)
    dk = 2π / Float64(boxsize)
    nk = size(src2_k, 1)
    CT = eltype(src2_k)

    kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
    ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
    kz_arr = Float64.(FFTW.fftfreq(n, n * dk))
    k_arrs = (kx_arr, ky_arr, kz_arr)

    out_k = similar(src2_k, CT)
    nyq = n ÷ 2 + 1

    Threads.@threads for iz in 1:n
        for iy in 1:n, ix in 1:nk
            kx = kx_arr[ix]; ky = ky_arr[iy]; kz = kz_arr[iz]
            k2 = kx^2 + ky^2 + kz^2

            if k2 == 0.0 || ix == nk || iy == nyq || iz == nyq
                out_k[ix, iy, iz] = 0
                continue
            end
            out_k[ix, iy, iz] = -im * k_arrs[dim][dim == 1 ? ix : (dim == 2 ? iy : iz)] / k2 * src2_k[ix, iy, iz]
        end
    end

    dummy_k = similar(out_k)
    invplan = plan_irfft(dummy_k, n; flags=FFTW.MEASURE)
    return invplan * out_k
end

"""
    compute_src2_k(delta_k, n, boxsize) -> src2_k

Compute the 2LPT source term in k-space with minimal peak memory.

Instead of computing all 6 phi_ij simultaneously (6 full-grid arrays),
computes them incrementally: off-diagonals one at a time (free immediately),
then diagonals in a memory-efficient order.

Peak memory: delta_k + src2 + 3 diagonal phi arrays (vs 6 simultaneous in `displacements_2lpt`).
"""
function compute_src2_k(delta_k, n::Int, boxsize::Real)
    dk = 2π / Float64(boxsize)
    nk = size(delta_k, 1)
    CT = eltype(delta_k)
    Tf = real(CT)

    kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
    ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
    kz_arr = Float64.(FFTW.fftfreq(n, n * dk))

    # Accumulator for the source term
    src2 = zeros(Tf, n, n, n)

    # Reusable k-space buffer
    phi_k = similar(delta_k, CT)

    # FFTW plans (reuse for all phi_ij computations)
    dummy_k = similar(phi_k)
    invplan = plan_irfft(dummy_k, n; flags=FFTW.MEASURE)

    # Helper: compute phi_ij in k-space
    function _fill_phi_k!(phi_k, di::Int, dj::Int)
        k_arrs = (kx_arr, ky_arr, kz_arr)
        Threads.@threads for iz in 1:n
            for iy in 1:n, ix in 1:nk
                kx = kx_arr[ix]; ky = ky_arr[iy]; kz = kz_arr[iz]
                k2 = kx^2 + ky^2 + kz^2
                if k2 == 0.0
                    phi_k[ix, iy, iz] = 0
                    continue
                end
                ki = k_arrs[di][di == 1 ? ix : (di == 2 ? iy : iz)]
                kj = k_arrs[dj][dj == 1 ? ix : (dj == 2 ? iy : iz)]
                phi_k[ix, iy, iz] = -ki * kj / k2 * delta_k[ix, iy, iz]
            end
        end
    end

    # --- Off-diagonals (each used once, free immediately) ---
    for (di, dj) in ((1,2), (1,3), (2,3))
        _fill_phi_k!(phi_k, di, dj)
        phi_ij = invplan * phi_k
        src2 .-= phi_ij .^ 2
    end

    # --- Diagonals (need pairwise products) ---
    _fill_phi_k!(phi_k, 1, 1)
    phi_11 = invplan * phi_k

    _fill_phi_k!(phi_k, 2, 2)
    phi_22 = invplan * phi_k

    _fill_phi_k!(phi_k, 3, 3)
    phi_33 = invplan * phi_k

    src2 .+= phi_11 .* phi_22 .+ phi_11 .* phi_33 .+ phi_22 .* phi_33

    # Free diagonals
    phi_11 = phi_22 = phi_33 = nothing

    # FFT source to k-space
    fwdplan = plan_rfft(src2; flags=FFTW.MEASURE)
    src2_k = fwdplan * src2

    return src2_k
end

"""
    compute_laplacian_field(delta_k, n, boxsize) -> Array

Compute k²-weighted field: irfft(k² × δ_k).
Matches the Fortran convention (positive k², no sign flip).
"""
function compute_laplacian_field(delta_k, n::Int, boxsize::Real)
    Tf = real(eltype(delta_k))
    dk = Tf(2π / boxsize)
    nk = n ÷ 2 + 1
    nyq = n ÷ 2

    lapd_k = similar(delta_k)
    Threads.@threads for iz in 1:n
        for iy in 1:n, ix in 1:nk
            kx = Tf(ix - 1) * dk
            ky = Tf(iy <= nyq + 1 ? iy - 1 : iy - 1 - n) * dk
            kz = Tf(iz <= nyq + 1 ? iz - 1 : iz - 1 - n) * dk
            k2 = kx^2 + ky^2 + kz^2
            lapd_k[ix, iy, iz] = k2 * delta_k[ix, iy, iz]
        end
    end
    lapd = irfft(lapd_k, n) ./ Tf(n)^3
    return Tf.(lapd)
end

end # module LPT
