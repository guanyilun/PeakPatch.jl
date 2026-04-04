@testset "Cosmology" begin
    # Planck 2018 parameters
    c = CosmologyParams(0.2645 + 0.0493, 0.0493, 1.0 - 0.2645 - 0.0493, 0.6735, 0.9649, 0.8111)

    @testset "E²(z)" begin
        # E²(0) = Om + OL = 1 for flat universe
        @test E2(0.0, c) ≈ 1.0 atol=1e-10
        # E²(z) > 1 for z > 0 (matter dominated)
        @test E2(1.0, c) > 1.0
        # E²(z) = Om*(1+z)^3 + OL  (exact for flat LCDM)
        z = 2.5
        @test E2(z, c) ≈ c.Om * (1+z)^3 + c.OL
    end

    @testset "H(z)" begin
        # H(0) = 100 * h km/s/Mpc
        @test H(0.0, c) ≈ 100.0 * 0.6735 rtol=1e-10
        # H(z) = 100*h*sqrt(E2(z))
        z = 1.5
        @test H(z, c) ≈ 100.0 * c.h * sqrt(E2(z, c)) rtol=1e-10
    end

    @testset "Growth factor D(z)" begin
        # D(0) = 1 by construction
        @test growth_factor(0.0, c) ≈ 1.0 rtol=1e-10
        # D(z) ~ 1/(1+z) at high z (matter domination, Ω_m correction)
        @test growth_factor(10.0, c) ≈ 1.0 / 11.0 rtol=0.3
        # D is monotonically decreasing
        @test growth_factor(1.0, c) < growth_factor(0.0, c)
        @test growth_factor(2.0, c) < growth_factor(1.0, c)
        # Known values from CLASS/Eisenstein-Hu for Planck18
        # D(z=0) = 1, D(z=0.5) ≈ 0.78, D(z=1) ≈ 0.61 (CPT92 fitting formula)
        @test 0.75 < growth_factor(0.5, c) < 0.82
        @test 0.58 < growth_factor(1.0, c) < 0.65
    end

    @testset "Comoving distance χ(z)" begin
        # χ(0) = 0
        @test chi(0.0, c) ≈ 0.0 atol=1e-10
        # χ(z) > 0 for z > 0
        @test chi(1.0, c) > 0.0
        # χ is monotonically increasing
        @test chi(2.0, c) > chi(1.0, c)
        # Known value: χ(z=1) ≈ 1550 Mpc (≈ 2300 Mpc/h) for Planck18
        # Note: chi() returns Mpc (includes h factor)
        @test 1500 < chi(1.0, c) < 1600
    end

    @testset "Growth rate f(z)" begin
        # f(z) ≈ Ω_m(z)^0.55 (Linder approximation)
        # At high z, Ω_m(z) → 1, so f → 1
        f_highz = growth_rate(100.0, c)
        @test abs(f_highz - 1.0) < 0.05
        # f(0) ≈ 0.53 for Planck18
        @test 0.48 < growth_rate(0.0, c) < 0.58
    end

    @testset "Critical overdensity δ_c(z) — Nakamura & Suto 1997" begin
        # Baseline: (3/20)(12π)^{2/3} ≈ 1.6865
        dc_base = (3.0/20.0) * (12π)^(2/3)
        @test dc_base ≈ 1.6865 rtol=1e-4

        # δ_c(0) should be close to 1.686 but slightly adjusted by N&S97 correction
        dc0 = delta_c(0.0, c)
        @test 1.67 < dc0 < 1.69

        # Formula: (3/20)(12π)^{2/3} * (1 + 0.0123 * log10(Ω_m(z)))
        # At z=0, Ω_m < 1 for ΛCDM, so log10(Ω_m) < 0, giving δ_c slightly < 1.686
        @test dc0 < dc_base

        # δ_c is weakly dependent on z
        dc1 = delta_c(1.0, c)
        @test abs(dc1 - dc0) / dc0 < 0.02

        # At high z (matter dominated), Ω_m(z) → 1, log10(1) = 0, δ_c → dc_base
        dc_highz = delta_c(100.0, c)
        @test abs(dc_highz - dc_base) / dc_base < 0.001

        # Spot-check exact formula against hand-computed value
        # For Om=0.314, at z=0: Omz = 0.314, log10(0.314) = -0.503
        c_test = CosmologyParams(0.314, 0.049, 0.686, 0.67, 0.96, 0.82)
        Omz0 = c_test.Om / (c_test.Om + c_test.OL)
        expected = dc_base * (1.0 + 0.0123 * log10(Omz0))
        @test delta_c(0.0, c_test) ≈ expected rtol=1e-10
    end

    @testset "Self-consistency: D(z) and f(z)" begin
        # Numerical derivative of D(z) should match growth rate
        # f(z) = -d ln D / d ln(1+z) = -(1+z)/D * dD/dz
        z = 1.0
        dz = 0.001
        D_plus  = growth_factor(z + dz, c)
        D_minus = growth_factor(z - dz, c)
        dDdz = (D_plus - D_minus) / (2 * dz)
        D_z = growth_factor(z, c)
        f_numerical = -(1 + z) / D_z * dDdz
        @test f_numerical ≈ growth_rate(z, c) rtol=0.05
    end
end
