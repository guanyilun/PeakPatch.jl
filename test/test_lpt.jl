@testset "LPT" begin
    using FFTW
    using Statistics: mean, std

    @testset "1LPT div(ψ) = -δ (k-space identity)" begin
        # In k-space: ψᵢ(k) = im * kᵢ/k² * δ̃(k)
        # div(ψ)(k) = Σᵢ im*kᵢ * ψᵢ(k) = Σᵢ im²*kᵢ²/k² * δ̃(k) = -δ̃(k)
        #
        # Nyquist modes are now zeroed in displacements_1lpt itself (matching Fortran),
        # so those modes won't satisfy div(ψ) = -δ. Skip them in the check.
        n = 64
        boxsize = 200.0
        pk_uniform(k) = 1.0

        field = generate_grf(n, pk_uniform, boxsize, 42)
        plan = plan_rfft(copy(field))
        delta_k = plan * field

        psi_x, psi_y, psi_z = displacements_1lpt(delta_k, n, boxsize)

        # FFT displacements back to k-space (need separate plans for each real array)
        psi_plan = plan_rfft(copy(psi_x))
        psi_x_k = psi_plan * psi_x
        psi_y_k = psi_plan * psi_y
        psi_z_k = psi_plan * psi_z

        dk = 2π / Float64(boxsize)
        kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
        ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
        kz_arr = Float64.(FFTW.fftfreq(n, n * dk))
        nk = size(delta_k, 1)
        nyq = n ÷ 2 + 1

        max_err = 0.0
        n_checked = 0
        for iz in 1:n, iy in 1:n, ix in 1:nk
            # Skip Nyquist modes (zeroed in the implementation)
            (ix == nk || iy == nyq || iz == nyq) && continue
            kx, ky, kz = kx_arr[ix], ky_arr[iy], kz_arr[iz]
            k2 = kx^2 + ky^2 + kz^2
            k2 == 0 && continue
            div_k = im * kx * psi_x_k[ix, iy, iz] +
                    im * ky * psi_y_k[ix, iy, iz] +
                    im * kz * psi_z_k[ix, iy, iz]
            err = abs(div_k + delta_k[ix, iy, iz])
            max_err = max(max_err, err)
            n_checked += 1
        end

        @test n_checked > 0  # sanity: we actually checked something
        # Float32 round-trip: error scales with field amplitude (~√n³ after n³ fix)
        # Tolerance 0.01 is appropriate for Float32 with n=64 grid
        @test max_err < 0.01
    end

    @testset "1LPT displacement physical properties" begin
        n = 64
        boxsize = 200.0
        pk_uniform(k) = 1.0

        field = generate_grf(n, pk_uniform, boxsize, 42)
        plan = plan_rfft(field)
        delta_k = plan * field

        psi_x, psi_y, psi_z = displacements_1lpt(delta_k, n, boxsize)

        # Displacements should be real and finite
        @test all(isfinite, psi_x)
        @test all(isfinite, psi_y)
        @test all(isfinite, psi_z)

        # Mean displacement should be approximately zero (DC mode zeroed)
        @test abs(mean(psi_x)) < 0.1
        @test abs(mean(psi_y)) < 0.1
        @test abs(mean(psi_z)) < 0.1
    end

    @testset "2LPT displacements are finite and nonzero" begin
        n = 32
        boxsize = 100.0
        pk_uniform(k) = 1.0

        field = generate_grf(n, pk_uniform, boxsize, 42)
        plan = plan_rfft(field)
        delta_k = plan * field

        psi2_x, psi2_y, psi2_z = displacements_2lpt(delta_k, n, boxsize)

        @test all(isfinite, psi2_x)
        @test all(isfinite, psi2_y)
        @test all(isfinite, psi2_z)

        # 2LPT should produce nonzero displacements
        @test std(psi2_x) > 0
        @test std(psi2_y) > 0
        @test std(psi2_z) > 0
    end

    @testset "2LPT mean displacement ≈ 0" begin
        n = 64
        boxsize = 200.0
        pk_uniform(k) = 1.0

        field = generate_grf(n, pk_uniform, boxsize, 42)
        plan = plan_rfft(field)
        delta_k = plan * field

        psi2_x, psi2_y, psi2_z = displacements_2lpt(delta_k, n, boxsize)

        @test abs(mean(psi2_x)) < 0.01
        @test abs(mean(psi2_y)) < 0.01
        @test abs(mean(psi2_z)) < 0.01
    end

    @testset "1LPT + 2LPT cross-correlation" begin
        # 1LPT and 2LPT should be correlated but not identical
        n = 64
        boxsize = 200.0
        pk_uniform(k) = 1.0

        field = generate_grf(n, pk_uniform, boxsize, 42)
        plan = plan_rfft(field)
        delta_k = plan * field

        psi_x, _, _ = displacements_1lpt(delta_k, n, boxsize)
        psi2_x, _, _ = displacements_2lpt(delta_k, n, boxsize)

        # They should be different (not trivially related)
        @test !(psi_x ≈ psi2_x)
    end

    @testset "2LPT sign: div(ψ²)(k) = +s(k) with -im gradient" begin
        # With the Fortran -im convention:
        #   ψ²ᵢ(k) = -im * kᵢ/k² * s(k)
        # The k-space divergence is:
        #   div(ψ²)(k) = Σᵢ im*kᵢ * ψ²ᵢ(k)
        #              = Σᵢ im*kᵢ * (-im*kᵢ/k²) * s(k)
        #              = Σᵢ (kᵢ²/k²) * s(k) = s(k)
        # So div(ψ²)(k) = +s(k).  We verify this directly.
        n = 32
        boxsize = 100.0
        pk_uniform(k) = 1.0

        field = generate_grf(n, pk_uniform, boxsize, 7)
        fwd_plan = plan_rfft(copy(field))
        delta_k = fwd_plan * field

        psi2_x, psi2_y, psi2_z = displacements_2lpt(delta_k, n, boxsize)

        # FFT the displacements back to k-space
        fwd_plan2 = plan_rfft(copy(psi2_x))
        psi2_x_k = fwd_plan2 * psi2_x
        psi2_y_k = fwd_plan2 * psi2_y
        psi2_z_k = fwd_plan2 * psi2_z

        # Independently recompute the 2LPT source s in k-space
        # (copy the phi_ij computation from LPT.jl)
        dk = 2π / Float64(boxsize)
        nk = size(delta_k, 1)
        kx_arr = Float64.(FFTW.rfftfreq(n, n * dk))
        ky_arr = Float64.(FFTW.fftfreq(n, n * dk))
        kz_arr = Float64.(FFTW.fftfreq(n, n * dk))

        phi11_k = similar(delta_k, ComplexF32)
        phi22_k = similar(delta_k, ComplexF32)
        phi33_k = similar(delta_k, ComplexF32)
        phi12_k = similar(delta_k, ComplexF32)
        phi13_k = similar(delta_k, ComplexF32)
        phi23_k = similar(delta_k, ComplexF32)
        for iz in 1:n, iy in 1:n, ix in 1:nk
            kx, ky, kz = kx_arr[ix], ky_arr[iy], kz_arr[iz]
            k2 = kx^2 + ky^2 + kz^2
            d = k2 == 0 ? ComplexF32(0) : delta_k[ix, iy, iz]
            inv_k2 = k2 == 0 ? 0f0 : Float32(1 / k2)
            phi11_k[ix, iy, iz] = -kx * kx * inv_k2 * d
            phi22_k[ix, iy, iz] = -ky * ky * inv_k2 * d
            phi33_k[ix, iy, iz] = -kz * kz * inv_k2 * d
            phi12_k[ix, iy, iz] = -kx * ky * inv_k2 * d
            phi13_k[ix, iy, iz] = -kx * kz * inv_k2 * d
            phi23_k[ix, iy, iz] = -ky * kz * inv_k2 * d
        end
        inv_phi = plan_irfft(phi11_k, n)
        phi11 = inv_phi * phi11_k
        phi22 = inv_phi * phi22_k
        phi33 = inv_phi * phi33_k
        phi12 = inv_phi * phi12_k
        phi13 = inv_phi * phi13_k
        phi23 = inv_phi * phi23_k
        src2 = phi11 .* phi22 .- phi12 .^ 2 .+
               phi11 .* phi33 .- phi13 .^ 2 .+
               phi22 .* phi33 .- phi23 .^ 2
        src2_k = fwd_plan * src2  # s(k)

        # Now check: div(ψ²)(k) ≈ s(k)  for non-Nyquist, non-DC modes
        nyq = n ÷ 2 + 1
        max_rel_err = 0.0
        n_checked = 0
        for iz in 1:n, iy in 1:n, ix in 1:nk
            (ix == nk || iy == nyq || iz == nyq) && continue
            kx, ky, kz = kx_arr[ix], ky_arr[iy], kz_arr[iz]
            k2 = kx^2 + ky^2 + kz^2
            k2 == 0 && continue
            div2 = im * kx * psi2_x_k[ix, iy, iz] +
                   im * ky * psi2_y_k[ix, iy, iz] +
                   im * kz * psi2_z_k[ix, iy, iz]
            s_ref = src2_k[ix, iy, iz]
            abs(s_ref) < 0.1 && continue  # skip near-zero modes
            rel = abs(div2 - s_ref) / abs(s_ref)
            max_rel_err = max(max_rel_err, rel)
            n_checked += 1
            n_checked >= 200 && break
        end
        @test n_checked > 0
        @test max_rel_err < 0.05  # within 5% (Float32 arithmetic + FFT round-trip)
    end
end
