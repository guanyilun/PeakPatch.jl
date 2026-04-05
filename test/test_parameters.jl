@testset "Parameters" begin
    @testset "Round-trip write/read params" begin
        p = FortranParams(
            Int32(0),    Int32(1),    0.5f0,       5.0f0,       Int32(1),
            0.31f0,      0.049f0,     0.69f0,      0.67f0,
            Int32(4),    Int32(4),    Int32(4),
            0.0f0,       1000.0f0,    0.0f0,       0.0f0,       0.0f0,
            Int32(10),   Int32(0),    Int32(0),
            Int32(1),    0.0f0,       0.0f0,       0.0f0,       200.0f0,
            Int32(0),
            Int32(100),  Int32(50),   Int32(50),
            -2.0f0,      2.0f0,       0.0f0,       0.5f0,       -0.5f0,      0.5f0,
            Int32(0),    3.0f0,       Int32(0),
            Int32(0),    0.0f0,       0.0f0,       0.0f0,       0.0f0,
            Int32(2),    Int32(0),    Int32(0),
            "fields",    "dens.dat",  "hpkvd_params.bin",
            "pk.dat",    "filters.dat", "output",  "HomelTab.dat"
        )

        tmppath = tempname()
        write_params_bin(tmppath, p)
        p2 = read_params_bin(tmppath)
        rm(tmppath)

        # Check all fields match
        @test p2.ireadfield == p.ireadfield
        @test p2.ioutshear == p.ioutshear
        @test p2.global_redshift == p.global_redshift
        @test p2.maximum_redshift == p.maximum_redshift
        @test p2.num_redshifts == p.num_redshifts
        @test p2.Omx == p.Omx
        @test p2.OmB == p.OmB
        @test p2.Omvac == p.Omvac
        @test p2.h == p.h
        @test p2.nlx == p.nlx
        @test p2.nly == p.nly
        @test p2.nlz == p.nlz
        @test p2.dL_box == p.dL_box
        @test p2.nbuff == p.nbuff
        @test p2.dcrit == p.dcrit
        @test p2.ilpt == p.ilpt
        @test p2.wsmooth == p.wsmooth
        @test p2.NonGauss == p.NonGauss
        @test p2.fNL == p.fNL
        @test p2.fielddir == p.fielddir
        @test p2.pkfile == p.pkfile
        @test p2.fileout == p.fileout
        @test p2.TabInterpFile == p.TabInterpFile
    end

    @testset "Binary format compatibility" begin
        # Write params and verify byte-level structure
        p = FortranParams(
            Int32(0),    Int32(1),    0.5f0,       5.0f0,       Int32(1),
            0.31f0,      0.049f0,     0.69f0,      0.67f0,
            Int32(4),    Int32(4),    Int32(4),
            0.0f0,       1000.0f0,    0.0f0,       0.0f0,       0.0f0,
            Int32(10),   Int32(0),    Int32(0),
            Int32(1),    0.0f0,       0.0f0,       0.0f0,       200.0f0,
            Int32(0),
            Int32(100),  Int32(50),   Int32(50),
            -2.0f0,      2.0f0,       0.0f0,       0.5f0,       -0.5f0,      0.5f0,
            Int32(0),    3.0f0,       Int32(0),
            Int32(0),    0.0f0,       0.0f0,       0.0f0,       0.0f0,
            Int32(2),    Int32(0),    Int32(0),
            "fields",    "dens.dat",  "hpkvd_params.bin",
            "pk.dat",    "filters.dat", "output",  "HomelTab.dat"
        )

        tmppath = tempname()
        write_params_bin(tmppath, p)

        # Check file size: 46 * 4 bytes (Float32) + 7 * (4 bytes + string length)
        expected_size = 46 * 4 + sum(s -> 4 + length(s),
            ["fields", "dens.dat", "hpkvd_params.bin", "pk.dat", "filters.dat", "output", "HomelTab.dat"])
        actual_size = filesize(tmppath)
        @test actual_size == expected_size
        rm(tmppath)
    end

    @testset "PipelineConfig from TOML" begin
        config = Dict{String,Any}(
            "cosmology" => Dict{String,Any}(
                "Om" => 0.315, "OB" => 0.049, "OL" => 0.685, "h" => 0.674
            ),
            "grid" => Dict{String,Any}(
                "n" => 128, "boxsize" => 300.0, "nbuff" => 6
            ),
            "run" => Dict{String,Any}(
                "z_out" => 1.0, "ilpt" => 2, "ioutshear" => 1,
                "wsmooth" => 0, "rmax2rs" => 1.5, "NonGauss" => 0, "fNL" => 0.0
            ),
            "files" => Dict{String,Any}(
                "pk" => "my_pk.dat",
                "filterbank" => "my_filters.dat",
                "homeltab" => "my_HomelTab.dat",
                "output" => "my_catalog.pksc"
            )
        )
        cfg = PipelineConfig(config)

        # Cosmology: Omx = Om - OB
        @test cfg.Omx ≈ 0.315 - 0.049
        @test cfg.OmB ≈ 0.049
        @test cfg.Omvac ≈ 0.685
        @test cfg.h ≈ 0.674

        # Grid
        @test cfg.n == 128
        @test cfg.boxsize ≈ 300.0
        @test cfg.nbuff == 6

        # Run
        @test cfg.z_out ≈ 1.0
        @test cfg.ilpt == 2
        @test cfg.ioutshear == 1
        @test cfg.wsmooth == 0
        @test cfg.rmax2rs ≈ 1.5

        # Files
        @test cfg.pkfile == "my_pk.dat"
        @test cfg.filterfile == "my_filters.dat"
        @test cfg.tabfile == "my_HomelTab.dat"
        @test cfg.fileout == "my_catalog.pksc"
    end

    @testset "PipelineConfig from TOML — defaults" begin
        cfg = PipelineConfig(Dict{String,Any}())
        @test cfg.Omx ≈ 0.315 - 0.049
        @test cfg.n == 142
        @test cfg.z_out ≈ 0.0
        @test cfg.ilpt == 2
        @test cfg.pkfile == "pk.dat"
    end

    @testset "PipelineConfig from TOML — dL_box backward compat" begin
        config = Dict{String,Any}(
            "grid" => Dict{String,Any}("dL_box" => 500.0)
        )
        cfg = PipelineConfig(config)
        @test cfg.boxsize ≈ 500.0
    end

    @testset "PipelineConfig from FortranParams" begin
        sp = FortranParams(
            Int32(0), Int32(1), 0.5f0, 5.0f0, Int32(1),
            0.31f0, 0.049f0, 0.69f0, 0.67f0,
            Int32(64), Int32(64), Int32(64),
            0.0f0, 200.0f0, 1.0f0, 2.0f0, 3.0f0,
            Int32(4), Int32(0), Int32(1),
            Int32(1), 0.0f0, 0.0f0, 0.0f0, 200.0f0, Int32(0),
            Int32(50), Int32(20), Int32(20),
            -2.0f0, 2.0f0, 0.0f0, 0.5f0, -0.5f0, 0.5f0,
            Int32(0), 1.5f0, Int32(0),
            Int32(0), 0.0f0, 0.0f0, 0.0f0, 0.0f0,
            Int32(2), Int32(0), Int32(0),
            "", "", "", "pk.dat", "filters.dat", "out.pksc", "HomelTab.dat"
        )
        cfg = PipelineConfig(sp)

        @test cfg.Omx ≈ Float64(sp.Omx)
        @test cfg.n == 64
        @test cfg.boxsize ≈ 200.0
        @test cfg.z_out ≈ Float64(sp.global_redshift)
        @test cfg.z_max ≈ Float64(sp.maximum_redshift)
        @test cfg.ievol == 1
        @test cfg.cenx ≈ 1.0
        @test cfg.ceny ≈ 2.0
        @test cfg.cenz ≈ 3.0
        @test cfg.ilpt == 2
        @test cfg.rmax2rs ≈ Float64(sp.rmax2rs)
        @test cfg.pkfile == "pk.dat"
        @test cfg.tabfile == "HomelTab.dat"
    end

    @testset "Read Fortran-produced hpkvd_params.bin" begin
        # Read the actual params file produced by peak-patch.py for the 142³ test
        params_path = joinpath(@__DIR__, "..", "fulltile_test_42", "hpkvd_params.bin")
        if isfile(params_path)
            p = read_params_bin(params_path)
            # These integer fields were the ones broken by the Float32 reinterpretation bug
            @test p.nlx == Int32(1)    # single tile: nlx=1 (was read as 0 before fix)
            @test p.nly == Int32(1)
            @test p.nlz == Int32(1)
            @test p.ireadfield == Int32(0) || p.ireadfield == Int32(1)
            @test p.num_redshifts >= Int32(1)
            @test p.nbuff > Int32(0)
            @test p.ilpt >= Int32(1)
            # Float fields should be physically reasonable
            @test 0.0f0 < p.Omx < 1.0f0
            @test 0.0f0 < p.OmB < 1.0f0
            @test 0.0f0 < p.h < 2.0f0
            @test p.dL_box > 0.0f0
        end
    end
end
