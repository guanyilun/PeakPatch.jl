@testset "RadialShell vs Fortran get_homel" begin

    # =========================================================
    # Fortran reference values from test_get_homel (single peak)
    # Synthetic 142^3 grid: delta(71,71,71)=10, neighbours=5
    # Non-zero displacements at (71,71,71) and face neighbours
    # Params: alatt=1, ir2min=4, ZZon=1.0, Rfclvi=4.0
    # =========================================================

    # Fortran fcrit (computed from collapse table via fsc_of_z)
    fortran_fcrit = 1.7161525487899780

    fortran_RTHL      = 1.7320507764816284
    fortran_Srb       = 1.9020472764968872
    fortran_Fbarx     = 2.3230779170989990
    fortran_Sbar      = (2.5925926864147186E-02,
                         5.1851853728294373E-02,
                         7.7777773141860962E-02)
    fortran_Sbar2     = (0.0, 0.0, 0.0)
    # Fortran strain_mat: diagonal with entries [0.387, 0.774, 1.162]
    fortran_strain    = [ 3.8717970252037048E-01  -0.0                    -0.0;
                         -0.0                     7.7435940504074097E-01  -0.0;
                         -0.0                    -0.0                     1.1615388393402100E+00]
    fortran_gradpk    = (0.0, 0.0, 0.0)
    fortran_gradpkrf  = (0.0, 0.0, 0.0)
    fortran_e_v       = 1.6666662693023682E-01
    fortran_p_v       = -5.1315236504478889E-08
    fortran_zvir_half = -1.0

    # =========================================================
    # Build identical synthetic grid in Julia
    # =========================================================
    n = 142
    delta = zeros(Float32, n, n, n)
    etax  = zeros(Float32, n, n, n)
    etay  = zeros(Float32, n, n, n)
    etaz  = zeros(Float32, n, n, n)

    # Overdensity
    delta[71,71,71] = 10.0f0
    delta[70,71,71] = 5.0f0; delta[72,71,71] = 5.0f0
    delta[71,70,71] = 5.0f0; delta[71,72,71] = 5.0f0
    delta[71,71,70] = 5.0f0; delta[71,71,72] = 5.0f0

    # 1LPT displacements (same as Fortran test)
    etax[71,71,71] = 0.1f0
    etax[70,71,71] = 0.05f0; etax[72,71,71] = 0.15f0
    etax[71,70,71] = 0.1f0;  etax[71,72,71] = 0.1f0
    etax[71,71,70] = 0.1f0;  etax[71,71,72] = 0.1f0

    etay[71,71,71] = 0.2f0
    etay[70,71,71] = 0.2f0; etay[72,71,71] = 0.2f0
    etay[71,70,71] = 0.1f0; etay[71,72,71] = 0.3f0
    etay[71,71,70] = 0.2f0; etay[71,71,72] = 0.2f0

    etaz[71,71,71] = 0.3f0
    etaz[70,71,71] = 0.3f0; etaz[72,71,71] = 0.3f0
    etaz[71,70,71] = 0.3f0; etaz[71,72,71] = 0.3f0
    etaz[71,71,70] = 0.15f0; etaz[71,71,72] = 0.45f0

    pg = PeakGrid(delta, etax, etay, etaz, nothing, nothing, nothing,
                  nothing, (n, n, n), nothing)

    # Precompute shells (same as icloud with nhunt=8)
    shells = precompute_shells(8)

    # Load collapse table (same as Fortran)
    table, tp = read_homeltab(joinpath(@__DIR__, "data", "HomelTab_fortran.dat"))
    ct = CollapseTableInterp(table, tp)

    # Peak index: column-major 1-based (71,71,71)
    ipp = 71 + 70*142 + 70*142*142

    # Run analyse_peak with Fortran's fcrit
    result = analyse_peak(pg, ipp, 1.0, 4, 1.0, 4.0, ct, shells;
                          growth_tables = nothing,
                          fcrit_override = fortran_fcrit)

    # =========================================================
    # Compare results
    # =========================================================

    # RTHL — should match to high precision
    @test result.RTHL > 0
    @test abs(result.RTHL - fortran_RTHL) / fortran_RTHL < 0.01

    # Srb
    @test abs(result.Srb - fortran_Srb) / abs(fortran_Srb) < 0.01

    # Fbarx
    @test abs(result.Fbarx - fortran_Fbarx) / abs(fortran_Fbarx) < 0.01

    # Sbar (displacements)
    for i in 1:3
        @test abs(result.Sbar[i] - fortran_Sbar[i]) < 0.01
    end

    # Sbar2 (2LPT displacements — all zero)
    for i in 1:3
        @test abs(result.Sbar2[i] - fortran_Sbar2[i]) < 0.01
    end

    # strain_mat
    for i in 1:3, j in 1:3
        val = fortran_strain[i, j]
        if abs(val) > 1e-10
            @test abs(result.strain_mat[i,j] - val) / abs(val) < 0.01
        else
            @test abs(result.strain_mat[i,j]) < 0.01
        end
    end

    # gradpk
    for i in 1:3
        @test abs(result.gradpk[i] - fortran_gradpk[i]) < 0.01
    end

    # gradpkrf
    for i in 1:3
        @test abs(result.gradpkrf[i] - fortran_gradpkrf[i]) < 0.01
    end

    # e_v, p_v
    @test abs(result.e_v - fortran_e_v) / abs(fortran_e_v) < 0.01
    @test abs(result.p_v - fortran_p_v) < 0.01

    # zvir_half
    @test abs(result.zvir_half - fortran_zvir_half) < 0.01
end
