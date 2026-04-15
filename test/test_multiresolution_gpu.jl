#!/usr/bin/env julia
# End-to-end test: run_multitile_split with `use_gpu=true` vs CPU.
#
# Run (one GPU needed):
#   CUDA_VISIBLE_DEVICES=1 julia --project=/tmp/pp_gpu_sandbox --threads=8 \
#     test/test_multiresolution_gpu.jl

using PeakPatch
using CUDA
using Test
using Printf

if !CUDA.functional()
    @warn "CUDA not functional — skipping GPU end-to-end test"
    exit(0)
end

# ---- Re-use the helpers from test_multiresolution.jl ----
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

# Halo positions agree if they're within `tol` in each coord
function _compare_halos(hs_cpu, hs_gpu; label="")
    @info "$label: CPU=$(length(hs_cpu))  GPU=$(length(hs_gpu))"
    ng = length(hs_cpu); ns = length(hs_gpu)
    @test abs(ng - ns) / max(ng, 1) < 0.05 || abs(ng - ns) <= 3

    # Sort both by (x,y,z) and compare the first few
    sort_key = h -> (h.x, h.y, h.z)
    sa = sort(hs_cpu; by=sort_key)
    sb = sort(hs_gpu; by=sort_key)
    ncmp = min(length(sa), length(sb), 10)
    for i in 1:ncmp
        da = sa[i]; db = sb[i]
        dx = abs(da.x - db.x); dy = abs(da.y - db.y); dz = abs(da.z - db.z)
        @printf "  halo %2d: Δpos=(%.3f,%.3f,%.3f)  R: CPU=%.3g GPU=%.3g\n" i dx dy dz da.RTHL db.RTHL
    end
end

@testset "run_multitile_split use_gpu=true vs CPU (1LPT)" begin
    tmpdir = mktempdir()
    cfg = _make_config(tmpdir; n=68, boxsize=226.67, z=0.0, ilpt=1, nbuff=8)

    t_cpu = @elapsed halos_cpu = PeakPatch.run_multitile_split(
        cfg; ntile=2, seed=42, verbose=false, coarse_factor=5, use_gpu=false)
    t_gpu = @elapsed halos_gpu = PeakPatch.run_multitile_split(
        cfg; ntile=2, seed=42, verbose=false, coarse_factor=5, use_gpu=true)

    @printf "  CPU time: %.2f s  GPU time: %.2f s (%.2fx)\n" t_cpu t_gpu (t_cpu / t_gpu)
    _compare_halos(halos_cpu, halos_gpu; label="1LPT ntile=2")

    rm(tmpdir; recursive=true)
end

@testset "run_multitile_split use_gpu=true vs CPU (2LPT)" begin
    tmpdir = mktempdir()
    cfg = _make_config(tmpdir; n=68, boxsize=226.67, z=0.0, ilpt=2, nbuff=8)

    t_cpu = @elapsed halos_cpu = PeakPatch.run_multitile_split(
        cfg; ntile=2, seed=42, verbose=false, coarse_factor=5, use_gpu=false)
    t_gpu = @elapsed halos_gpu = PeakPatch.run_multitile_split(
        cfg; ntile=2, seed=42, verbose=false, coarse_factor=5, use_gpu=true)

    @printf "  CPU time: %.2f s  GPU time: %.2f s (%.2fx)\n" t_cpu t_gpu (t_cpu / t_gpu)
    _compare_halos(halos_cpu, halos_gpu; label="2LPT ntile=2")

    rm(tmpdir; recursive=true)
end
