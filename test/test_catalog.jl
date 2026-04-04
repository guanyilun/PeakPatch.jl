@testset "Catalog" begin
    @testset "Round-trip write/read" begin
        halos_in = [
            HaloRecord(1.0f0, 2.0f0, 3.0f0, 0.1f0, 0.2f0, 0.3f0,
                       5.0f0, -0.1f0, -0.2f0, -0.3f0, 1.686f0),
            HaloRecord(10.0f0, 20.0f0, 30.0f0, 1.0f0, 2.0f0, 3.0f0,
                       8.0f0, -1.0f0, -2.0f0, -3.0f0, 2.5f0),
        ]
        RTHLmax = 8.0f0
        z_out = 0.5f0

        tmppath = tempname()
        write_pksc(tmppath, halos_in, RTHLmax, z_out)
        halos_out, RTHLmax_out, z_out_out = read_pksc(tmppath)
        rm(tmppath)

        @test length(halos_out) == length(halos_in)
        @test RTHLmax_out == RTHLmax
        @test z_out_out == z_out

        for i in eachindex(halos_in)
            @test halos_out[i].x  == halos_in[i].x
            @test halos_out[i].y  == halos_in[i].y
            @test halos_out[i].z  == halos_in[i].z
            @test halos_out[i].vx == halos_in[i].vx
            @test halos_out[i].vy == halos_in[i].vy
            @test halos_out[i].vz == halos_in[i].vz
            @test halos_out[i].RTHL == halos_in[i].RTHL
            @test halos_out[i].vx2 == halos_in[i].vx2
            @test halos_out[i].vy2 == halos_in[i].vy2
            @test halos_out[i].vz2 == halos_in[i].vz2
            @test halos_out[i].overdensity == halos_in[i].overdensity
        end
    end

    @testset "Empty catalog" begin
        halos_in = HaloRecord[]
        tmppath = tempname()
        write_pksc(tmppath, halos_in, 0.0f0, 0.0f0)
        halos_out, RTHLmax, z_out = read_pksc(tmppath)
        rm(tmppath)

        @test length(halos_out) == 0
        @test RTHLmax == 0.0f0
        @test z_out == 0.0f0
    end

    @testset "Large catalog" begin
        # Write 1000 halos to verify no off-by-one errors
        N = 1000
        halos_in = [HaloRecord(
            Float32(i), Float32(2i), Float32(3i),
            Float32(0.1i), Float32(0.2i), Float32(0.3i),
            Float32(sqrt(i)),
            Float32(-0.1i), Float32(-0.2i), Float32(-0.3i),
            Float32(1.0 + 0.001i)
        ) for i in 1:N]

        tmppath = tempname()
        write_pksc(tmppath, halos_in, Float32(sqrt(N)), 1.0f0)
        halos_out, _, _ = read_pksc(tmppath)
        rm(tmppath)

        @test length(halos_out) == N
        # Spot-check first and last
        @test halos_out[1].x == 1.0f0
        @test halos_out[N].x == Float32(N)
        @test halos_out[N].overdensity ≈ Float32(1.0 + 0.001N)
    end

    @testset "Extended catalog round-trip (33 fields)" begin
        halos_in = [
            ExtHaloRecord(1.0f0, 2.0f0, 3.0f0, 0.1f0, 0.2f0, 0.3f0,
                          5.0f0, -0.1f0, -0.2f0, -0.3f0, 1.686f0,
                          0.15f0, 0.05f0,     # e_v, p_v
                          0.1f0, 0.2f0, 0.3f0, 0.01f0, 0.02f0, 0.03f0,  # strain
                          0.5f0, 1.2f0,       # d2F, zform
                          0.01f0, 0.02f0, 0.03f0,  # grad
                          0.04f0, 0.05f0, 0.06f0,  # gradf
                          4.0f0,              # Rf
                          1.5f0, 0.3f0,       # FcollvRf, d2FRf
                          0.07f0, 0.08f0, 0.09f0)  # gradrf
        ]
        tmppath = tempname()
        write_pksc(tmppath, halos_in, 5.0f0, 0.0f0)

        # File size: header (12) + 1 × 33 × 4 = 144
        @test filesize(tmppath) == 12 + 33 * 4

        halos_out, RTHLmax, z_out = read_pksc(tmppath)
        rm(tmppath)

        @test length(halos_out) == 1
        @test halos_out[1] isa ExtHaloRecord
        h = halos_out[1]
        @test h.x == 1.0f0
        @test h.RTHL == 5.0f0
        @test h.e_v == 0.15f0
        @test h.p_v == 0.05f0
        @test h.strain_11 == 0.1f0
        @test h.d2F == 0.5f0
        @test h.zform == 1.2f0
        @test h.Rf == 4.0f0
        @test h.FcollvRf == 1.5f0
        @test h.gradrf_z == 0.09f0
    end

    @testset "Binary layout" begin
        # Verify the binary format matches Fortran expectation:
        # Int32 (N) + Float32 (RTHLmax) + Float32 (z_out) + N * (11 * Float32)
        N = 5
        halos = [HaloRecord(1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0, 6.0f0,
                            7.0f0, 8.0f0, 9.0f0, 10.0f0, 11.0f0)
                 for _ in 1:N]
        tmppath = tempname()
        write_pksc(tmppath, halos, 7.0f0, 0.5f0)

        expected_size = 4 + 4 + 4 + N * 11 * 4  # header + data
        @test filesize(tmppath) == expected_size
        rm(tmppath)
    end
end
