@testset "EllipsoidalCollapse" begin

using .PeakPatch

# Planck 2018 cosmology (matches test_cosmology.jl)
cosmo = CosmologyParams(0.315, 0.049, 0.685, 0.674, 0.965, 0.808)

@testset "Carlson RD" begin
    # RD(x,x,x) = x^{-3/2}
    for x in [1.0, 2.0, 0.5, 10.0]
        rd = PeakPatch._elliptic_rd(x, x, x)
        @test rd ≈ x^(-1.5) rtol=1e-10
    end

    # Just check non-trivial case is positive and finite
    @test PeakPatch._elliptic_rd(1.0, 2.0, 3.0) > 0
    @test isfinite(PeakPatch._elliptic_rd(1.0, 2.0, 3.0))
end

@testset "Shape integrals get_b_2" begin
    # Spherical case: (1,1,1) -> (1/3, 1/3, 1/3)
    b = get_b_2(1.0, 1.0, 1.0, 1)
    @test b[1] ≈ 1/3 rtol=1e-10
    @test b[2] ≈ 1/3 rtol=1e-10
    @test b[3] ≈ 1/3 rtol=1e-10

    # Sum should always be 1
    for (a1, a2, a3) in [(0.5, 0.8, 1.2), (0.3, 0.7, 1.5), (0.9, 0.95, 1.0)]
        b = get_b_2(a1, a2, a3, 1)
        @test sum(b) ≈ 1.0 rtol=1e-10
    end

    # For a1 <= a2 <= a3, the r_3 >= 0.999 check always triggers (since r_3 = a3/a1 >= 1)
    # so b = (1/3, 1/3, 1/3) for normal ellipsoid evolution with a1 <= a2 <= a3.
    # The triaxial branch (Carlson RD) is only reached for unusual axis orderings.
    b = get_b_2(0.5, 0.8, 1.2, 1)
    @test b[1] ≈ 1/3 rtol=1e-10
end

@testset "DlinearTables" begin
    tables = Dlinear_tables(cosmo)

    # D(a=1) should be 1.0
    dlin, HD_Ha, D_a = Dlinear_ab(1.0, tables)
    @test dlin ≈ 1.0 rtol=1e-3

    # D should be positive for all a > 0
    for a in [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5]
        dlin, HD_Ha, D_a = Dlinear_ab(a, tables)
        @test dlin > 0
        @test isfinite(dlin)
    end

    # Dfnofa at a=1 should give 1.0
    @test Dfnofa(1.0, tables) ≈ 1.0 rtol=1e-3

    # D should be proportional to a at early times
    d1 = Dfnofa(0.01, tables)
    d2 = Dfnofa(0.02, tables)
    @test d2 / d1 ≈ 2.0 rtol=0.05

    # Compare with CPT92 approximation
    for z in [0.0, 0.5, 1.0, 2.0]
        a = 1.0 / (1.0 + z)
        d_table = Dfnofa(a, tables)
        d_cpt = growth_factor(z, cosmo)
        @test d_table ≈ d_cpt rtol=5e-3
    end
end

@testset "Single ellipsoid: spherical case" begin
    ep = EllipsoidParams(cosmo)

    # Spherical collapse (e=0, p=0) with Frho = 1.686 should collapse near z=0
    Frho = 1.686
    zvir, Dvir, fcvir = evolve_ellipse_full(Frho, 0.0, 0.0, ep)
    # z_vir should be close to 0 (positive for overdensity)
    @test zvir > -0.5
    @test zvir > 0  # for this overdensity, should collapse
    @test isfinite(zvir)
end

@testset "Single ellipsoid: known values" begin
    ep = EllipsoidParams(cosmo)

    # High Frho should give high z_vir
    zvir_high, _, _ = evolve_ellipse_full(5.0, 0.0, 0.0, ep)
    zvir_low, _, _ = evolve_ellipse_full(1.5, 0.0, 0.0, ep)
    @test zvir_high > zvir_low

    # Increasing eccentricity should give different z_vir
    zvir_e0, _, _ = evolve_ellipse_full(3.0, 0.0, 0.0, ep)
    zvir_e3, _, _ = evolve_ellipse_full(3.0, 0.3, 0.0, ep)
    @test zvir_e0 > 0
    @test zvir_e3 > 0
    @test isfinite(zvir_e0)
    @test isfinite(zvir_e3)
end

@testset "CollapseTable round-trip" begin
    ep = EllipsoidParams(cosmo)

    # Small table for fast tests
    tp = CollapseTableParams(Nx=5, Ny=3, Nz=3)
    table = make_table(ep, tp; verbose=false)

    # Check dimensions
    @test size(table) == (5, 3, 3)

    # Round-trip I/O
    tmppath = tempname()
    write_homeltab(tmppath, table, tp)
    table2, tp2 = read_homeltab(tmppath)
    rm(tmppath)

    @test size(table2) == size(table)
    @test table2 ≈ table

    # Verify header
    @test tp2.Nx == tp.Nx
    @test tp2.Ny == tp.Ny
    @test tp2.Nz == tp.Nz
    @test tp2.X1 ≈ tp.X1 rtol=1e-6
    @test tp2.X2 ≈ tp.X2 rtol=1e-6
end

@testset "Trilinear interpolation" begin
    # Create a simple linear table to test interpolation
    tp = CollapseTableParams(Nx=11, Ny=11, Nz=11,
        X1=0.0, X2=1.0, Y1=0.0, Y2=1.0, Z1=0.0, Z2=1.0)

    # Table = x + y + z (linear function should be interpolated exactly)
    dX = 0.1; dY = 0.1; dZ = 0.1
    table = Float32[
        Float32((i-1)*dX + (j-1)*dY + (k-1)*dZ)
        for i in 1:11, j in 1:11, k in 1:11
    ]

    ct = CollapseTableInterp(table, tp)

    # Test exact interpolation at grid points
    @test PeakPatch.interpolate(ct, 0.0, 0.0, 0.0) ≈ 0.0 atol=1e-5
    @test PeakPatch.interpolate(ct, 0.5, 0.0, 0.0) ≈ 0.5 atol=1e-5
    @test PeakPatch.interpolate(ct, 0.0, 0.5, 0.0) ≈ 0.5 atol=1e-5

    # Test interpolation at intermediate points (should be exact for linear)
    @test PeakPatch.interpolate(ct, 0.25, 0.25, 0.25) ≈ 0.75 atol=1e-4

    # Out of bounds should return -1
    @test PeakPatch.interpolate(ct, -0.1, 0.0, 0.0) == -1.0
    @test PeakPatch.interpolate(ct, 1.1, 0.0, 0.0) == -1.0
end

end # EllipsoidalCollapse testset
