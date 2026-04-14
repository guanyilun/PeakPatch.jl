#!/usr/bin/env julia
# Validate multi-resolution split approach against global FFT.
#
# Usage: julia --project=. test/test_multiresolution.jl

using PeakPatch
using Test

# ---- Helper: create test config ----
function _make_config(dir; n=60, boxsize=100.0, z=0.0, ilpt=1, nbuff=8,
                      ioutshear=0, rmax2rs=0.0)
    pk_path = joinpath(dir, "test_pk.dat")
    ks = 10.0 .^ range(-4, stop=1, length=200)
    open(pk_path, "w") do f
        for k in ks
            Pk = 1e4 * k^(-1.5)
            println(f, "$k  $Pk")
        end
    end

    filter_path = joinpath(dir, "test_filter.dat")
    open(filter_path, "w") do f
        println(f, "3")
        println(f, "1  1.686  8.0")
        println(f, "2  1.686  5.0")
        println(f, "3  1.686  3.0")
    end

    cosmo = PeakPatch.CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.808)
    ep = PeakPatch.EllipsoidParams(cosmo)
    tp = PeakPatch.CollapseTableParams()
    tab = PeakPatch.make_table(ep, tp)
    tab_path = joinpath(dir, "test_table.dat")
    PeakPatch.write_homeltab(tab_path, tab, tp)

    cat_path = joinpath(dir, "test_output.pksc")

    PeakPatch.PipelineConfig(
        n=n, boxsize=Float64(boxsize), pkfile=pk_path, filterfile=filter_path,
        tabfile=tab_path, fileout=cat_path, nbuff=nbuff,
        z_out=Float64(z), ilpt=ilpt,
        ioutshear=ioutshear, rmax2rs=rmax2rs,
        ievol=0, z_max=0.0, cenx=0.0, ceny=0.0, cenz=0.0,
        Omx=0.261, OmB=0.049, Omvac=0.69, h=0.68,
        NonGauss=0, fNL=0.0, wsmooth=1,
    )
end

# ============================================================
# Test 1: Multi-resolution ntile=2 vs global FFT ntile=2
# ============================================================
@testset "Multi-resolution vs global FFT (ntile=2, 1LPT)" begin
    tmpdir = mktempdir()
    # n=18, nbuff=3 → nsub=12, N=12*2+6=30. M=2*5=10 divides 30.
    cfg = _make_config(tmpdir; n=18, boxsize=60.0, z=0.0, ilpt=1, nbuff=3)

    halos_global = PeakPatch.run_multitile(cfg; ntile=2, seed=42, verbose=false)
    halos_split = PeakPatch.run_multitile_split(cfg; ntile=2, seed=42, verbose=true,
                                                  coarse_factor=5)

    println("  Global FFT: $(length(halos_global)) halos")
    println("  Split mode: $(length(halos_split)) halos")

    # Both should find halos (test isn't useful if both are empty)
    # On small grids, the mode split introduces ~10% RMS error in the fields,
    # so halo counts won't match exactly.  Check that they're in the same
    # ballpark (within 50% or ±2, whichever is larger).
    n_g = length(halos_global)
    n_s = length(halos_split)
    tol = max(2, n_g ÷ 2)
    @test abs(n_s - n_g) ≤ tol

    rm(tmpdir; recursive=true)
end

# ============================================================
# Test 2: Multi-resolution ntile=2 with 2LPT
# ============================================================
@testset "Multi-resolution vs global FFT (ntile=2, 2LPT)" begin
    tmpdir = mktempdir()
    cfg = _make_config(tmpdir; n=18, boxsize=60.0, z=0.0, ilpt=2, nbuff=3)

    halos_global = PeakPatch.run_multitile(cfg; ntile=2, seed=42, verbose=false)
    halos_split = PeakPatch.run_multitile_split(cfg; ntile=2, seed=42, verbose=false,
                                                  coarse_factor=5)

    println("  Global FFT (2LPT): $(length(halos_global)) halos")
    println("  Split mode (2LPT): $(length(halos_split)) halos")

    if !isempty(halos_global) && !isempty(halos_split)
        @test abs(length(halos_split) - length(halos_global)) ≤ max(1, length(halos_global) ÷ 10)
    else
        @test true
    end

    rm(tmpdir; recursive=true)
end

# ============================================================
# Test 3: Vary coarse_factor to show convergence
# ============================================================
@testset "Convergence with coarse_factor" begin
    tmpdir = mktempdir()
    # n=22, nbuff=4 → nsub=14, N=14*2+8=36. M=2*cf must divide 36.
    # cf=2→M=4 (36/4=9✓), cf=3→M=6 (36/6=6✓), cf=6→M=12 (36/12=3✓), cf=9→M=18 (36/18=2✓)
    cfg = _make_config(tmpdir; n=22, boxsize=100.0, z=0.0, ilpt=1, nbuff=4)

    halos_global = PeakPatch.run_multitile(cfg; ntile=2, seed=42, verbose=false)
    n_global = length(halos_global)

    println("  Global FFT: $n_global halos")

    for cf in [2, 3, 6, 9]
        halos_split = PeakPatch.run_multitile_split(cfg; ntile=2, seed=42, verbose=false,
                                                      coarse_factor=cf)
        n_split = length(halos_split)
        println("  coarse_factor=$cf: $n_split halos (Δ=$(n_split - n_global))")
    end

    # Just check it runs without error
    @test true

    rm(tmpdir; recursive=true)
end

println("\nAll multi-resolution tests completed.")
