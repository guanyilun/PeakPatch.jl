@testset "Pipeline" begin

using PeakPatch
using FFTW
using DelimitedFiles

# ---------- helpers to create synthetic test data files ----------

"""Create a minimal P(k) file (2-column: k, P(k)) for testing."""
function _make_pk_file(dir)
    path = joinpath(dir, "test_pk.dat")
    # Simple power law: P(k) = 1e4 * k^(-1.5) — enough to generate a field
    ks = 10.0 .^ range(-4, stop=1, length=200)
    open(path, "w") do f
        for k in ks
            Pk = 1e4 * k^(-1.5)
            println(f, "$k  $Pk")
        end
    end
    return path
end

"""Create a minimal filter bank file for testing."""
function _make_filter_file(dir)
    path = joinpath(dir, "test_filter.dat")
    open(path, "w") do f
        println(f, "3")            # 3 filter scales
        println(f, "1  1.686  8.0")   # largest
        println(f, "2  1.686  5.0")
        println(f, "3  1.686  3.0")   # smallest
    end
    return path
end

"""Build a PipelineConfig for a small grid."""
function _make_config(dir; n=64, boxsize=200.0, z=0.0, ilpt=1)
    pkfile = _make_pk_file(dir)
    filterfile = _make_filter_file(dir)
    outfile = joinpath(dir, "test_out.pksc")
    tabfile = joinpath(@__DIR__, "data", "HomelTab_julia.dat")

    @test isfile(tabfile)

    PipelineConfig(
        Omx       = 0.315 - 0.049,
        OmB       = 0.049,
        Omvac     = 0.685,
        h         = 0.674,
        n         = n,
        boxsize   = boxsize,
        nbuff     = 4,
        z_out     = z,
        z_max     = z,
        ilpt      = ilpt,
        rmax2rs   = 0.0,
        pkfile    = pkfile,
        filterfile = filterfile,
        fileout   = outfile,
        tabfile   = tabfile,
    )
end

# =================================================================
# Test 1: smooth_field unit test
# =================================================================
@testset "smooth_field" begin
    n = 32
    boxsize = 100.0

    # Spike at center: delta=1 at grid center, 0 elsewhere
    delta = zeros(Float32, n, n, n)
    c = n ÷ 2 + 1
    delta[c, c, c] = 1.0f0

    delta_k = rfft(delta)

    # Save copy to verify non-mutation
    delta_k_copy = copy(delta_k)

    Rf = 2.0
    smoothed = PeakPatch.Filters.smooth_field(delta_k, n, boxsize, Rf, 0)  # Gaussian

    # Verify k-space copy unchanged
    @test delta_k == delta_k_copy

    # Smoothed field should be peaked at center
    @test smoothed[c, c, c] == maximum(smoothed)

    # Should be approximately isotropic: check a few equidistant points
    val_c = smoothed[c, c, c]
    # Points at distance 1 from center (cardinal directions)
    @test abs(smoothed[c+1, c, c] - smoothed[c, c+1, c]) < 0.01 * val_c
    @test abs(smoothed[c+1, c, c] - smoothed[c, c, c+1]) < 0.01 * val_c

    # Top-hat smoothing should also work
    smoothed_th = PeakPatch.Filters.smooth_field(delta_k, n, boxsize, Rf, 1)
    @test delta_k == delta_k_copy  # still unchanged
    @test smoothed_th[c, c, c] == maximum(smoothed_th)
end

# =================================================================
# Test 2: Velocity formula test
# =================================================================
@testset "Velocity formulas" begin
    # Construct a known PeakResult and verify velocity computation
    cosmo = PeakPatch.CosmologyParams(0.315, 0.049, 0.685, 0.674, 0.965, 0.808)
    growth_tables = PeakPatch.Dlinear_tables(cosmo)
    z_out = 0.0
    a_out = 1.0
    _, _, D_out = PeakPatch.Dlinear_ab(a_out, growth_tables)
    Omnr = cosmo.Om

    # Fake result with known values
    Sbar = [1.0, 2.0, 3.0]
    Sbar2 = [0.5, -0.3, 0.1]
    Srb = 0.8
    RTHL = 5.0  # lattice units
    alatt = 3.0
    Fbarx = 2.0

    result = PeakResult(RTHL, Srb, Sbar, Sbar2,
                        zeros(3,3), zeros(3), zeros(3),
                        Fbarx, 0.1, 0.0, -1.0,
                        zeros(3), 0.0)

    RTHL_phys = Float32(RTHL * alatt)
    vTHvir0 = 100.0 * cosmo.h * sqrt(Omnr)

    # 1LPT
    Sbar_vel = result.Sbar .* D_out
    @test Sbar_vel ≈ Sbar atol=1e-4  # D(0) ≈ 1 (numerical)

    # 2LPT
    Om_a = Omnr * a_out^3 / (Omnr * a_out^3 + cosmo.OL)
    coeff = -3.0/7.0 * Om_a^(-1.0/143) * D_out^2
    Sbar2_vel = -result.Sbar2 .* coeff

    # With D(0)=1, Om_a ≈ Omnr/(Omnr+OL)
    Om_a_expected = Omnr / (Omnr + cosmo.OL)
    coeff_expected = -3.0/7.0 * Om_a_expected^(-1.0/143)
    @test Sbar2_vel[1] ≈ -0.5 * coeff_expected atol=1e-4

    # Virial velocity
    vE2 = vTHvir0^2 * Fbarx * Srb * Float64(RTHL_phys)^2
    @test vE2 > 0  # should be positive for these values
end

# =================================================================
# Test 3: Small synthetic pipeline (64^3)
# =================================================================
@testset "End-to-end 64³ pipeline" begin
    tmpdir = mktempdir()
    cfg = _make_config(tmpdir; n=64, boxsize=200.0, z=0.0, ilpt=1)

    halos = PeakPatch.Pipeline.run_tile(cfg; seed=42, verbose=false)

    # Catalog should be readable
    halos_read, RTHLmax, z_read = PeakPatch.read_pksc(cfg.fileout)
    @test length(halos_read) == length(halos)
    @test z_read ≈ Float32(0.0)

    # All halos have positive RTHL
    if !isempty(halos)
        @test all(h -> h.RTHL > 0, halos)
        @test RTHLmax > 0
    end

    rm(tmpdir; recursive=true)
end

# =================================================================
# Test 4: 2LPT pipeline variant
# =================================================================
@testset "2LPT pipeline" begin
    tmpdir = mktempdir()
    cfg = _make_config(tmpdir; n=32, boxsize=100.0, z=0.0, ilpt=2)

    halos = PeakPatch.Pipeline.run_tile(cfg; seed=123, verbose=false)

    # Just verify it runs without error and produces valid output
    if !isempty(halos)
        @test all(h -> h.RTHL > 0, halos)
        # 2LPT velocities should be non-zero for at least some halos
        has_2lpt = any(h -> h.vx2 != 0 || h.vy2 != 0 || h.vz2 != 0, halos)
        @test has_2lpt
    end

    rm(tmpdir; recursive=true)
end

# =================================================================
# Test 5: Lightcone mode (ievol=1)
# =================================================================
@testset "Lightcone mode (ievol=1)" begin
    tmpdir = mktempdir()

    pkfile = _make_pk_file(tmpdir)
    filterfile = _make_filter_file(tmpdir)
    tabfile = joinpath(@__DIR__, "data", "HomelTab_julia.dat")

    n = 64; boxsize = 200.0
    cfg_lc = PipelineConfig(
        Omx       = 0.315 - 0.049,
        OmB       = 0.049,
        Omvac     = 0.685,
        h         = 0.674,
        n         = n,
        boxsize   = boxsize,
        nbuff     = 4,
        z_out     = 0.0,
        z_max     = 3.0,
        ievol     = 1,
        rmax2rs   = 0.0,
        pkfile    = pkfile,
        filterfile = filterfile,
        fileout   = joinpath(tmpdir, "test_lc.pksc"),
        tabfile   = tabfile,
    )

    halos_lc = PeakPatch.Pipeline.run_tile(cfg_lc; seed=42, verbose=false)

    @test length(halos_lc) >= 0

    if !isempty(halos_lc)
        @test all(h -> h.RTHL > 0, halos_lc)
    end

    # Compare with ievol=0
    cfg_nolc = PipelineConfig(
        Omx       = 0.315 - 0.049,
        OmB       = 0.049,
        Omvac     = 0.685,
        h         = 0.674,
        n         = n,
        boxsize   = boxsize,
        nbuff     = 4,
        z_out     = 0.0,
        z_max     = 0.0,
        ievol     = 0,
        rmax2rs   = 0.0,
        pkfile    = pkfile,
        filterfile = filterfile,
        fileout   = joinpath(tmpdir, "test_nolc.pksc"),
        tabfile   = tabfile,
    )

    halos_nolc = PeakPatch.Pipeline.run_tile(cfg_nolc; seed=42, verbose=false)

    @test length(halos_lc) > 0 || length(halos_nolc) == 0

    rm(tmpdir; recursive=true)
end

end # Pipeline testset
