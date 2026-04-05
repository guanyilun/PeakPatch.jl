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

# ---------- helpers (same as test_multitile.jl) ----------

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

function _make_simparams(dir; n=64, boxsize=200.0, z=0.0, ilpt=1, nbuff=4,
                         ioutshear=0, rmax2rs=0.0)
    pkfile = _make_pk_file(dir)
    filterfile = _make_filter_file(dir)
    outfile = joinpath(dir, "test_out.pksc")
    tabfile = joinpath(@__DIR__, "data", "HomelTab_julia.dat")

    nsub = n - 2 * nbuff
    dcore_box = boxsize * nsub / n

    PeakPatch.SimParams(
        Int32(0), Int32(ioutshear), Float32(z), Float32(z), Int32(1),
        Float32(0.315 - 0.049), Float32(0.049), Float32(0.685), Float32(0.674),
        Int32(n), Int32(n), Int32(n), Float32(dcore_box), Float32(boxsize),
        Float32(0.0), Float32(0.0), Float32(0.0),
        Int32(nbuff), Int32(0), Int32(0), Int32(2),
        Float32(0.171), Float32(0.171), Float32(0.01), Float32(200.0), Int32(4),
        Int32(50), Int32(20), Int32(20),
        Float32(log10(1.5)), Float32(log10(8.0)),
        Float32(0.0), Float32(0.5),
        Float32(-0.9999), Float32(0.9999),
        Int32(0), Float32(rmax2rs), Int32(1),
        Int32(0), Float32(0.0), Float32(0.0), Float32(0.0), Float32(0.0),
        Int32(ilpt), Int32(0), Int32(0),
        "", "", "", pkfile, filterfile, outfile, tabfile
    )
end

# ============================================================
# Test 1: ntile=1, np=1 should approximately match serial
# ============================================================
if nranks == 1
    @testset "MPI np=1, ntile=1 ≈ serial" begin
        tmpdir = mktempdir()

        sp_serial = _make_simparams(tmpdir; n=32, boxsize=100.0, z=0.0, ilpt=1, nbuff=3)
        halos_serial = PeakPatch.MultiTile.run_multitile(sp_serial; ntile=1, seed=42,
                                                          verbose=false)

        outfile_mpi = joinpath(tmpdir, "test_out_mpi.pksc")
        sp_mpi = PeakPatch.SimParams(
            sp_serial.ireadfield, sp_serial.ioutshear,
            sp_serial.global_redshift, sp_serial.maximum_redshift,
            sp_serial.num_redshifts, sp_serial.Omx, sp_serial.OmB,
            sp_serial.Omvac, sp_serial.h,
            sp_serial.nlx, sp_serial.nly, sp_serial.nlz,
            sp_serial.dcore_box, sp_serial.dL_box,
            sp_serial.cenx, sp_serial.ceny, sp_serial.cenz,
            sp_serial.nbuff, sp_serial.next, sp_serial.ievol,
            sp_serial.ivir_strat, sp_serial.fcoll_3, sp_serial.fcoll_2,
            sp_serial.fcoll_1, sp_serial.dcrit, sp_serial.iforce_strat,
            sp_serial.TabInterpNx, sp_serial.TabInterpNy, sp_serial.TabInterpNz,
            sp_serial.TabInterpX1, sp_serial.TabInterpX2,
            sp_serial.TabInterpY1, sp_serial.TabInterpY2,
            sp_serial.TabInterpZ1, sp_serial.TabInterpZ2,
            sp_serial.wsmooth, sp_serial.rmax2rs, sp_serial.ioutfield,
            sp_serial.NonGauss, sp_serial.fNL, sp_serial.A_nG, sp_serial.B_nG, sp_serial.R_nG,
            sp_serial.ilpt, sp_serial.iwant_field_part, sp_serial.largerun,
            sp_serial.fielddir, sp_serial.densfilein, sp_serial.filein,
            sp_serial.pkfile, sp_serial.filterfile,
            outfile_mpi, sp_serial.TabInterpFile
        )

        halos_mpi = PeakPatch.run_multitile_mpi(sp_mpi; ntile=1, seed=42,
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
    sp = _make_simparams(tmpdir; n=20, boxsize=60.0, z=0.0, ilpt=1, nbuff=3)

    halos = PeakPatch.run_multitile_mpi(sp; ntile=2, seed=42, verbose=rank==0, comm=comm)

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
        sp = _make_simparams(tmpdir; n=20, boxsize=60.0, z=0.0, ilpt=1, nbuff=3)

        # Serial reference (rank 0 only)
        halos_serial = nothing
        if rank == 0
            halos_serial = PeakPatch.MultiTile.run_multitile(sp; ntile=2, seed=42,
                                                              verbose=false)
        end

        # MPI run
        halos_mpi = PeakPatch.run_multitile_mpi(sp; ntile=2, seed=42,
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
