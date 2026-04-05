#!/usr/bin/env julia
# Run with: mpiexecjl -np N julia --project=. test/test_mpi_multitile.jl
#
# Tests that run_multitile_mpi produces results consistent with serial run_multitile.
# Note: distributed FFT uses Float64, serial uses Float32, so results are approximate.

using MPI
MPI.Init()

using PeakPatch
using PencilFFTs
using PencilArrays
using Test

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nranks = MPI.Comm_size(comm)

# ---------- helpers ----------

function _make_pk_file(dir)
    path = joinpath(dir, "test_pk.dat")
    ks = 10.0 .^ range(-4, stop=1, length=200)
    open(path, "w") do f
        for k in ks
            Pk = 1e4 * k^(-1.5)
            println(f, "$k  $Pk")
        end
    end
    return path
end

function _make_filter_file(dir)
    path = joinpath(dir, "test_filter.dat")
    open(path, "w") do f
        println(f, "3")
        println(f, "1  1.686  8.0")
        println(f, "2  1.686  5.0")
        println(f, "3  1.686  3.0")
    end
    return path
end

function _make_config(dir; n=64, boxsize=200.0, z=0.0, ilpt=1, nbuff=4,
                      ioutshear=0, rmax2rs=0.0)
    pkfile = _make_pk_file(dir)
    filterfile = _make_filter_file(dir)
    outfile = joinpath(dir, "test_out.pksc")
    tabfile = joinpath(@__DIR__, "data", "HomelTab_julia.dat")

    PipelineConfig(
        Omx       = 0.315 - 0.049,
        OmB       = 0.049,
        Omvac     = 0.685,
        h         = 0.674,
        n         = n,
        boxsize   = boxsize,
        nbuff     = nbuff,
        z_out     = z,
        z_max     = z,
        ilpt      = ilpt,
        ioutshear = ioutshear,
        rmax2rs   = rmax2rs,
        pkfile    = pkfile,
        filterfile = filterfile,
        fileout   = outfile,
        tabfile   = tabfile,
    )
end

# ============================================================
# Test 1: ntile=1, np=1 should approximately match serial
# ============================================================
if nranks == 1
    @testset "MPI np=1, ntile=1 ≈ serial" begin
        tmpdir = mktempdir()

        cfg_serial = _make_config(tmpdir; n=32, boxsize=100.0, z=0.0, ilpt=1, nbuff=3)
        halos_serial = PeakPatch.MultiTile.run_multitile(cfg_serial; ntile=1, seed=42,
                                                          verbose=false)

        cfg_mpi = PipelineConfig(
            Omx = cfg_serial.Omx, OmB = cfg_serial.OmB,
            Omvac = cfg_serial.Omvac, h = cfg_serial.h,
            n = cfg_serial.n, boxsize = cfg_serial.boxsize, nbuff = cfg_serial.nbuff,
            z_out = cfg_serial.z_out, z_max = cfg_serial.z_max, ilpt = cfg_serial.ilpt,
            rmax2rs = cfg_serial.rmax2rs, pkfile = cfg_serial.pkfile,
            filterfile = cfg_serial.filterfile,
            fileout = joinpath(tmpdir, "test_out_mpi.pksc"),
            tabfile = cfg_serial.tabfile,
        )

        halos_mpi = PeakPatch.run_multitile_mpi(cfg_mpi; ntile=1, seed=42,
                                                  verbose=false, comm=comm)

        @test length(halos_serial) == length(halos_mpi)

        if !isempty(halos_serial) && !isempty(halos_mpi)
            sort!(halos_serial; by=h -> (h.x, h.y, h.z))
            sort!(halos_mpi; by=h -> (h.x, h.y, h.z))

            # Approximate match: Float64 distributed FFT vs Float32 serial FFT
            for i in 1:length(halos_serial)
                @test halos_serial[i].x ≈ halos_mpi[i].x atol=0.5
                @test halos_serial[i].y ≈ halos_mpi[i].y atol=0.5
                @test halos_serial[i].z ≈ halos_mpi[i].z atol=0.5
                @test halos_serial[i].RTHL ≈ halos_mpi[i].RTHL rtol=0.1
            end
        end

        rm(tmpdir; recursive=true)
    end
end

# ============================================================
# Test 2: ntile=2 runs without error and produces halos
# ============================================================
@testset "MPI ntile=2 end-to-end (np=$nranks)" begin
    tmpdir = mktempdir()
    cfg = _make_config(tmpdir; n=20, boxsize=60.0, z=0.0, ilpt=1, nbuff=3)

    halos = PeakPatch.run_multitile_mpi(cfg; ntile=2, seed=42, verbose=rank==0, comm=comm)

    @test length(halos) >= 0  # should not error

    if !isempty(halos)
        @test all(h -> h.RTHL > 0, halos)
    end

    rm(tmpdir; recursive=true)
end

# ============================================================
# Test 3: ntile=2 multi-rank matches serial (if nranks >= 2)
# ============================================================
if nranks >= 2
    @testset "MPI ntile=2 np=$nranks ≈ serial" begin
        tmpdir = mktempdir()
        cfg = _make_config(tmpdir; n=20, boxsize=60.0, z=0.0, ilpt=1, nbuff=3)

        # Serial reference (rank 0 only)
        halos_serial = nothing
        if rank == 0
            halos_serial = PeakPatch.MultiTile.run_multitile(cfg; ntile=2, seed=42,
                                                              verbose=false)
        end

        # MPI run
        halos_mpi = PeakPatch.run_multitile_mpi(cfg; ntile=2, seed=42,
                                                  verbose=false, comm=comm)

        # run_multitile_mpi returns all halos on rank 0, empty on others
        if rank == 0
            @test length(halos_serial) == length(halos_mpi)

            if !isempty(halos_serial) && !isempty(halos_mpi)
                sort!(halos_serial; by=h -> (h.x, h.y, h.z))
                sort!(halos_mpi; by=h -> (h.x, h.y, h.z))

                for i in 1:min(length(halos_serial), length(halos_mpi))
                    @test halos_serial[i].x ≈ halos_mpi[i].x atol=0.5
                    @test halos_serial[i].y ≈ halos_mpi[i].y atol=0.5
                    @test halos_serial[i].z ≈ halos_mpi[i].z atol=0.5
                end
            end
        else
            @test isempty(halos_mpi)
        end

        rm(tmpdir; recursive=true)
    end
end

rank == 0 && println("\nAll MPI tests passed (nranks=$nranks)")

MPI.Finalize()
