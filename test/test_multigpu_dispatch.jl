#!/usr/bin/env julia
# Multi-GPU tile dispatcher parity test.
#
# Verifies that `run_multitile_split(...; use_gpu=true, devices=[0,1])`
# produces the same halo catalog as the sequential GPU path.
# Requires 2 visible CUDA devices and Julia launched with --threads>=2.

using Test, PeakPatch, CUDA, Random

@assert CUDA.functional()
ndev = length(CUDA.devices())
println("Visible CUDA devices: $ndev")
@assert ndev >= 2 "multi-GPU test needs at least 2 devices visible"
@assert Threads.nthreads() >= 2 "start Julia with --threads>=2"

function make_config(dir; n, boxsize, ilpt, nbuff=8, ioutshear=0)
    pk_path = joinpath(dir, "pk.dat")
    open(pk_path, "w") do f
        for k in 10.0 .^ range(-4, stop=1, length=200)
            println(f, "$k  $(1e4 * k^(-1.5))")
        end
    end
    filter_path = joinpath(dir, "filter.dat")
    open(filter_path, "w") do f
        println(f, "3"); println(f, "1  1.686  8.0")
        println(f, "2  1.686  5.0"); println(f, "3  1.686  3.0")
    end
    cosmo = PeakPatch.CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.808)
    tp = PeakPatch.CollapseTableParams()
    tab = PeakPatch.make_table(PeakPatch.EllipsoidParams(cosmo), tp)
    tab_path = joinpath(dir, "tab.dat")
    PeakPatch.write_homeltab(tab_path, tab, tp)
    PeakPatch.PipelineConfig(
        n=n, boxsize=Float64(boxsize), pkfile=pk_path, filterfile=filter_path,
        tabfile=tab_path, fileout=joinpath(dir, "out.pksc"), nbuff=nbuff,
        z_out=0.0, ilpt=ilpt, ioutshear=ioutshear, rmax2rs=0.0,
        ievol=0, z_max=0.0, cenx=0.0, ceny=0.0, cenz=0.0,
        Omx=0.261, OmB=0.049, Omvac=0.69, h=0.68,
        NonGauss=0, fNL=0.0, wsmooth=1,
    )
end

# Small config: ntile=2 → 8 tiles, distributed 4 per GPU.
# nmesh=68 as in the existing gpu parity test.
const NTILE = 2; const NMESH = 68; const NBUFF = 8
const NSUB = NMESH - 2 * NBUFF
const N = NSUB * NTILE + 2 * NBUFF
const BOXSIZE = Float64(N) * 226.67 / 68.0

@testset "multi-GPU dispatcher vs single-GPU (ilpt=$ilpt)" for ilpt in (1, 2)
    tmp = mktempdir()
    cfg = make_config(tmp; n=NMESH, boxsize=BOXSIZE, ilpt=ilpt, nbuff=NBUFF)

    # Single-GPU baseline
    println(">> single-GPU (ilpt=$ilpt)...")
    t1 = @elapsed h1 = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=42,
                                                      verbose=false, use_gpu=true)
    # 2-GPU parallel
    println(">> 2-GPU dispatch (ilpt=$ilpt)...")
    t2 = @elapsed h2 = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=42,
                                                      verbose=false, use_gpu=true,
                                                      devices=[0, 1])

    @test length(h1) == length(h2)
    println("   single=$(length(h1))  multi=$(length(h2))   " *
            "times: 1×=$(round(t1,digits=2))s  2×=$(round(t2,digits=2))s  " *
            "speedup=$(round(t1/t2, digits=2))x")

    # Sort by coordinates to compare (tile ordering differs between paths).
    key(h) = (round(h.x, digits=4), round(h.y, digits=4), round(h.z, digits=4))
    s1 = sort(h1; by=key); s2 = sort(h2; by=key)
    # Check top-10 by radius
    top10(h) = sort(h; by=r -> -r.RTHL)[1:min(10, length(h))]
    for (a, b) in zip(top10(s1), top10(s2))
        @test isapprox(a.x, b.x; atol=1e-4)
        @test isapprox(a.y, b.y; atol=1e-4)
        @test isapprox(a.z, b.z; atol=1e-4)
        @test isapprox(a.RTHL, b.RTHL; rtol=1e-4)
    end
    rm(tmp; recursive=true)
end
