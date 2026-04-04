@testset "PowerSpectrum" begin
    using Interpolations

    @testset "Load and interpolate test P(k)" begin
        pkfile = joinpath(@__DIR__, "data", "test_pk.dat")
        pk = load_pk(pkfile)

        # P(k) should be positive at tabulated points
        @test pk(0.01)  > 0.0
        @test pk(0.1)   > 0.0
        @test pk(1.0)   > 0.0

        # P(k) should be roughly decreasing at high k
        @test pk(1.0) < pk(0.01)
        @test pk(5.0) < pk(0.1)
    end

    @testset "Log-log interpolation accuracy" begin
        # Create a known power-law: P(k) = A * k^n
        # Write to temp file, load, and verify
        A, n = 1000.0, -1.5
        tmpfile = tempname()
        open(tmpfile, "w") do io
            for logk in -3.0:0.5:1.0
                k = 10.0^logk
                pk = A * k^n
                println(io, "$k  $pk")
            end
        end

        pk = load_pk(tmpfile)
        rm(tmpfile)

        # Check accuracy at tabulated and intermediate points
        for k in [0.005, 0.05, 0.5, 5.0]
            expected = A * k^n
            @test pk(k) ≈ expected rtol=0.01  # log-linear interp on power law
        end
    end

    @testset "Extrapolation" begin
        pkfile = joinpath(@__DIR__, "data", "test_pk.dat")
        pk = load_pk(pkfile)

        # Extrapolation beyond table range should return positive values
        @test pk(1e-6) > 0.0
        @test pk(100.0) > 0.0
    end

    @testset "Integrability (for σ_R computation)" begin
        pkfile = joinpath(@__DIR__, "data", "test_pk.dat")
        pk = load_pk(pkfile)

        # Variance σ²_R = (1/2π²) ∫ P(k) k² W²(kR) dk
        # Should be finite for Gaussian and top-hat windows
        using QuadGK
        R = 8.0  # Mpc/h
        integrand(k) = pk(k) * k^2 * exp(-k^2 * R^2) / (2π^2)
        sigma2, _ = quadgk(integrand, 1e-4, 100.0)
        @test sigma2 > 0.0
        @test isfinite(sigma2)
    end
end
