@testset "RadialShell" begin

@testset "hRinteg boundary values" begin
    # u0=0: full integral from 0 to 2 — matches Fortran return value 1.4*u12-31/35
    r0 = hRinteg(0.0, 0.0)
    @test r0 ≈ 1.4 * 0.0 - 31.0 / 35.0 atol = 1e-12

    # u0=0 with nonzero u12: verify Fortran formula
    r0b = hRinteg(0.0, 1.5)
    @test r0b ≈ 1.4 * 1.5 - 31.0 / 35.0 atol = 1e-12

    # u0 >= 2: integral is zero
    @test hRinteg(2.0, 0.0) == 0.0
    @test hRinteg(3.0, 1.0) == 0.0

    # u0=1: on the boundary between the two kernel regions
    r1 = hRinteg(1.0, 0.0)
    # Fortran: hRinteg=1.6*0-6.4/7, then x0=1, x02=1, x04=1
    # h1 = 1*(2-1*(2.4-1*(1-1/7))) = 2 - 2.4 + 1 - 1/7 = 0.6 - 1/7 = 3.2/7
    # h2 = 1*(4-1*(4-1*(1.5-0.2))) = 4 - 4 + 1.5 - 0.2 = 1.3
    # result = -6.4/7 - 1.3*0 + 3.2/7 = -3.2/7
    @test r1 ≈ -3.2 / 7.0 atol = 1e-12

    # symmetry: negative u0 should give same result
    @test hRinteg(0.5, 2.0) ≈ hRinteg(-0.5, 2.0)
end

@testset "hRinteg u12 dependence" begin
    # With u12 > 0, the integral includes extra u12 terms
    r_a = hRinteg(0.5, 1.0)
    r_b = hRinteg(0.5, 2.0)
    # Should differ
    @test r_a != r_b
end

@testset "atab4" begin
    akk = atab4()
    @test length(akk) == 1101
    # at i=1 (u=0, ru=0): akk = ca*(2/3 - 0 + 0) = ca * 2/3
    ca = 3.0 / (2.0 * pi)
    @test akk[1] ≈ ca * (2.0 / 3.0) atol = 1e-14
    # at ru=2 (u=4, i=401): akk should be 0
    @test akk[401] == 0.0
    # beyond ru=2: all zero
    @test akk[500] == 0.0
    @test akk[1100] == 0.0
    # continuity at ru=1 (u=1, i=101)
    ca_val = ca * (2.0/3.0 - 1.0 + 0.5 * 1.0 * 1.0)  # ca*(2/3 - 1 + 0.5)
    cb = 1.0 / (4.0 * pi)
    cb_val = cb * (2.0 - 1.0)^3  # cb * 1
    @test ca_val ≈ cb_val atol = 1e-14
end

@testset "precompute_shells" begin
    # rmax=1: center + 6 face neighbors = 7 cells
    cells = precompute_shells(1)
    @test length(cells) == 7
    @test cells[1].r2 == 0  # center cell first
    # All r2 <= 1
    for c in cells
        @test c.r2 <= 1
    end

    # rmax=2: should have all cells with r^2 <= 4
    cells2 = precompute_shells(2)
    @test length(cells2) == 33  # 1 + 6 + 12 + 8 + 6
    # sorted by r2
    for i in 2:length(cells2)
        @test cells2[i].r2 >= cells2[i-1].r2
    end
end

@testset "get_evals" begin
    # Identity matrix → eigenvalues all 1
    I3 = Float64[1 0 0; 0 1 0; 0 0 1]
    lam, vecs, iflag = get_evals(I3)
    @test iflag == 0
    @test lam ≈ [1.0, 1.0, 1.0] atol = 1e-12

    # Diagonal matrix
    D = Float64[1 0 0; 0 2 0; 0 0 3]
    lam2, _, _ = get_evals(D)
    @test lam2 ≈ [1.0, 2.0, 3.0] atol = 1e-12
end

@testset "normalize_strain!" begin
    S = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
    normalize_strain!(S, 12.0)
    @test S[1,1] + S[2,2] + S[3,3] ≈ 12.0 atol = 1e-12
    @test S[1,1] ≈ 2.0 atol = 1e-12
    @test S[2,2] ≈ 4.0 atol = 1e-12
    @test S[3,3] ≈ 6.0 atol = 1e-12
end

@testset "analyse_peak synthetic" begin
    # 20^3 grid, peak at center
    n = 20
    delta = zeros(Float32, n, n, n)
    etax = zeros(Float32, n, n, n)
    etay = zeros(Float32, n, n, n)
    etaz = zeros(Float32, n, n, n)

    # Put a strong overdensity at the center
    i0, j0, k0 = 10, 10, 10
    delta[i0, j0, k0] = 10.0f0
    # Add some surrounding overdensity
    for di in -1:1, dj in -1:1, dk in -1:1
        delta[i0+di, j0+dj, k0+dk] = 5.0f0
    end

    pg = PeakGrid(delta, etax, etay, etaz, nothing, nothing, nothing, nothing,
                  (n, n, n), nothing)

    # Create a minimal collapse table (read from file or synthetic)
    # Use a small table with all values = 1.0 (zvir = 1 for all entries)
    tp = PeakPatch.CollapseTableParams(Nx=5, Ny=3, Nz=3)
    table = fill(Float32(1.0), 5, 3, 3)
    ct = PeakPatch.CollapseTableInterp(table, tp)

    shells = precompute_shells(8)

    # ipp = center cell in column-major 3D grid
    # ipp = i + (j-1)*n1 + (k-1)*n1*n2
    ipp = i0 + (j0 - 1) * n + (k0 - 1) * n * n

    result = analyse_peak(pg, ipp, 1.0, 4, 0.0, 4.0, ct, shells)
    # With a strong overdensity and zvir table = 1.0, should collapse
    @test result.RTHL > 0.0
    @test result.Srb > 0.0
end

@testset "analyse_peak no collapse" begin
    n = 20
    delta = zeros(Float32, n, n, n)
    etax = zeros(Float32, n, n, n)
    etay = zeros(Float32, n, n, n)
    etaz = zeros(Float32, n, n, n)

    # Low overdensity everywhere
    delta .= 0.5f0

    pg = PeakGrid(delta, etax, etay, etaz, nothing, nothing, nothing, nothing,
                  (n, n, n), nothing)

    # Table with very high zvir values — nothing collapses
    tp = PeakPatch.CollapseTableParams(Nx=5, Ny=3, Nz=3)
    table = fill(Float32(100.0), 5, 3, 3)
    ct = PeakPatch.CollapseTableInterp(table, tp)

    shells = precompute_shells(8)
    i0, j0, k0 = 10, 10, 10
    ipp = i0 + (j0 - 1) * n + (k0 - 1) * n * n

    result = analyse_peak(pg, ipp, 1.0, 4, 0.0, 4.0, ct, shells)
    @test result.RTHL == -1.0
end

end # testset RadialShell
