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

# Simple P(k) interpolator for field-level tests
function _make_pk()
    tmpfile = tempname() * ".dat"
    ks = 10.0 .^ range(-4, stop=1, length=200)
    open(tmpfile, "w") do f
        for k in ks
            Pk = 1e4 * k^(-1.5)
            println(f, "$k  $Pk")
        end
    end
    pk = PeakPatch.load_pk(tmpfile)
    rm(tmpfile)
    return pk
end

# ============================================================
# Test 1: Field-level comparison (main accuracy test)
# ============================================================
@testset "Field-level accuracy (ntile=2, N=120)" begin
    # N=120, ntile=2, nbuff=8 → nsub=52, nmesh=68
    # coarse_factor=5 → M=10, block=12
    N = 120; ntile = 2; nbuff = 8; cf = 5
    boxsize = 400.0; seed = 42
    pk = _make_pk()

    results = compare_fields_split(pk, N, boxsize, seed, ntile, cf;
                                    nbuff=nbuff, verbose=true)

    for r in results
        println("  Tile $(r.tile): corr_δ=$(round(r.corr_delta; digits=5)), " *
                "rms_δ=$(round(r.rms_delta*100; digits=2))%, " *
                "corr_ψ₁=$(round(r.corr_psi1; digits=5)), " *
                "rms_ψ₁=$(round(r.rms_psi1*100; digits=2))%")
    end

    for r in results
        @test r.corr_delta > 0.999
        @test r.rms_delta < 0.05
        @test r.corr_psi1 > 0.93
        @test r.rms_psi1 < 0.40
    end
end

# ============================================================
# Test 2: Coarse factor sweep
# ============================================================
@testset "Coarse factor sweep (N=120)" begin
    N = 120; ntile = 2; nbuff = 8
    boxsize = 400.0; seed = 42
    pk = _make_pk()

    # Accuracy is NOT monotonic in cf. Small cf → better delta (residual is
    # more short-range) but worse ψ₁ (coarse grid too coarse for 1/k² modes).
    # Large cf → worse delta (residual has long modes) but better ψ₁.
    # cf=3-5 is optimal for halo counts (delta-dominated).
    for cf in [3, 5, 6, 10, 20]
        results = compare_fields_split(pk, N, boxsize, seed, ntile, cf;
                                        nbuff=nbuff, verbose=false)
        mean_corr_d = sum(r.corr_delta for r in results) / length(results)
        mean_rms_d = sum(r.rms_delta for r in results) / length(results)
        mean_corr_p = sum(r.corr_psi1 for r in results) / length(results)
        mean_rms_p = sum(r.rms_psi1 for r in results) / length(results)
        println("  cf=$cf (M=$(2*cf), block=$(N÷(2*cf))): " *
                "δ: corr=$(round(mean_corr_d; digits=5)) rms=$(round(mean_rms_d*100; digits=2))% | " *
                "ψ₁: corr=$(round(mean_corr_p; digits=5)) rms=$(round(mean_rms_p*100; digits=2))%")

        @test mean_corr_d > 0.99
        @test mean_rms_d < 0.15
    end
end

# ============================================================
# Test 3: Halo count comparison (end-to-end, 1LPT)
# ============================================================
@testset "Halo count comparison (ntile=2, 1LPT)" begin
    tmpdir = mktempdir()
    # n=68, nbuff=8 → nsub=52, N=120. cf=5→M=10 divides 120.
    cfg = _make_config(tmpdir; n=68, boxsize=226.67, z=0.0, ilpt=1, nbuff=8)

    halos_global = PeakPatch.run_multitile(cfg; ntile=2, seed=42, verbose=false)
    halos_split = PeakPatch.run_multitile_split(cfg; ntile=2, seed=42, verbose=true,
                                                  coarse_factor=5)

    n_g = length(halos_global)
    n_s = length(halos_split)
    println("  Global FFT: $n_g halos")
    println("  Split mode: $n_s halos")

    if n_g > 0
        frac_diff = abs(n_s - n_g) / n_g
        println("  Fractional difference: $(round(frac_diff*100; digits=1))%")
        @test frac_diff < 0.05 || abs(n_s - n_g) <= 3  # <5% or ≤3 halos
    end

    rm(tmpdir; recursive=true)
end

# ============================================================
# Test 4: Halo count comparison with 2LPT
# ============================================================
@testset "Halo count comparison (ntile=2, 2LPT)" begin
    tmpdir = mktempdir()
    cfg = _make_config(tmpdir; n=68, boxsize=226.67, z=0.0, ilpt=2, nbuff=8)

    halos_global = PeakPatch.run_multitile(cfg; ntile=2, seed=42, verbose=false)
    halos_split = PeakPatch.run_multitile_split(cfg; ntile=2, seed=42, verbose=false,
                                                  coarse_factor=5)

    n_g = length(halos_global)
    n_s = length(halos_split)
    println("  Global FFT (2LPT): $n_g halos")
    println("  Split mode (2LPT): $n_s halos")

    if n_g > 0
        frac_diff = abs(n_s - n_g) / n_g
        println("  Fractional difference: $(round(frac_diff*100; digits=1))%")
        @test frac_diff < 0.05 || abs(n_s - n_g) <= 3
    end

    rm(tmpdir; recursive=true)
end

println("\nAll multi-resolution tests completed.")
