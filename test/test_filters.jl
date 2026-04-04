@testset "Filters" begin
    @testset "Gaussian window" begin
        # W(0) = 1
        @test gaussian_window(0.0) ≈ 1.0
        # W(kR) = exp(-k²R²/2) is positive and monotonically decreasing
        for kR in [0.1, 0.5, 1.0, 2.0, 5.0]
            @test gaussian_window(kR) ≈ exp(-kR^2 / 2) rtol=1e-10
        end
        # Gaussian window decays rapidly
        @test gaussian_window(10.0) < 1e-20
    end

    @testset "Gaussian window: Fortran R/2 convention" begin
        # Fortran smooth_field divides kR by 2 before calling wgauss,
        # giving W(kR) = exp(-k²R²/8), i.e. effective σ = R/2.
        # gaussian_window_fortran must match this.
        for kR in [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]
            @test gaussian_window_fortran(kR) ≈ exp(-kR^2 / 8) rtol=1e-10
        end
        # gaussian_window_fortran is strictly narrower than gaussian_window
        @test gaussian_window_fortran(1.0) > gaussian_window(1.0)
    end

    @testset "Top-hat window" begin
        # W(0) = 1 (by L'Hôpital)
        @test tophat_window(0.0) ≈ 1.0

        # Known values: W(kR) = 3(sin(kR) - kR cos(kR)) / (kR)³
        for kR in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
            expected = 3.0 * (sin(kR) - kR * cos(kR)) / kR^3
            @test tophat_window(kR) ≈ expected rtol=1e-10
        end

        # Top-hat has zeros (first zero near kR ≈ 4.493)
        # W changes sign after first zero
        @test tophat_window(4.0) > 0
        @test tophat_window(5.0) < 0
    end

    @testset "Filter bank file I/O" begin
        # Write a test filter bank file
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "3")
            println(io, "1  1.686  2.0")
            println(io, "2  1.686  4.0")
            println(io, "3  1.686  8.0")
        end

        fb = read_filterbank(tmpfile)
        rm(tmpfile)

        @test length(fb) == 3
        @test fb[1] == (1, 1.686, 2.0)
        @test fb[2] == (2, 1.686, 4.0)
        @test fb[3] == (3, 1.686, 8.0)
    end

    @testset "Window function normalization" begin
        # ∫₀^∞ W²(kR) k² dk = 1/(2π² R³) for top-hat (real space)
        # For Gaussian: ∫ W²(kR) 4π k² dk = (2π)^(3/2) / R³
        # Just check the windows are well-behaved at small kR (Taylor expand)
        # Gaussian: 1 - (kR)²/2 + ...
        @test gaussian_window(0.01) > 0.999
        # Top-hat: 1 - (kR)²/10 + ...
        @test tophat_window(0.01) > 0.999
    end
end
