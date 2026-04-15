#!/usr/bin/env julia
# Per-stage profile for run_multitile_split. Runs CPU then GPU path and
# prints stage accounting so we can pick the next GPU port target.
#
# Run:
#   CUDA_VISIBLE_DEVICES=1 julia --project=/tmp/pp_gpu_sandbox --threads=64 \
#       bench/profile_multitile.jl
#
# Knobs (env):
#   NTILE (default 2), NMESH (default 128), ILPT (default 2),
#   IOUTSHEAR (default 0), NBUFF (default 8)

using PeakPatch
using CUDA
using Printf

const NTILE = parse(Int, get(ENV, "NTILE", "2"))
const NMESH = parse(Int, get(ENV, "NMESH", "128"))
const ILPT  = parse(Int, get(ENV, "ILPT",  "2"))
const IOUTSHEAR = parse(Int, get(ENV, "IOUTSHEAR", "0"))
const NBUFF = parse(Int, get(ENV, "NBUFF", "8"))
# Lightcone mode (ievol=1) evaluates each tile at the redshift implied by
# its observer distance. Set IEVOL=1 to exercise the per-peak-ZZon path.
const IEVOL = parse(Int, get(ENV, "IEVOL", "0"))
const Z_MAX = parse(Float64, get(ENV, "ZMAX", "3.0"))

function make_config(dir; n, boxsize, z=0.0, ilpt=2, nbuff=8, ioutshear=0,
                      ievol=0, z_max=0.0,
                      cenx=0.0, ceny=0.0, cenz=0.0)
    pk_path = joinpath(dir, "test_pk.dat")
    open(pk_path, "w") do f
        for k in 10.0 .^ range(-4, stop=1, length=200)
            println(f, "$k  $(1e4 * k^(-1.5))")
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
        z_out=Float64(z), ilpt=ilpt, ioutshear=ioutshear, rmax2rs=0.0,
        ievol=ievol, z_max=Float64(z_max),
        cenx=Float64(cenx), ceny=Float64(ceny), cenz=Float64(cenz),
        Omx=0.261, OmB=0.049, Omvac=0.69, h=0.68,
        NonGauss=0, fNL=0.0, wsmooth=1,
    )
end

println("\n" * "="^60)
println("Per-stage profile: CPU vs GPU")
println("="^60)
nsub = NMESH - 2 * NBUFF
N = nsub * NTILE + 2 * NBUFF
@printf "ntile=%d nmesh=%d nsub=%d N=%d  ilpt=%d ioutshear=%d  ievol=%d  z_max=%.2f\n" NTILE NMESH nsub N ILPT IOUTSHEAR IEVOL Z_MAX
@printf "threads=%d GPU=%s\n" Threads.nthreads() (CUDA.functional() ? CUDA.name(CUDA.device()) : "n/a")
println()

tmpdir = mktempdir()
boxsize = Float64(N) * 226.67 / 68.0
cfg = make_config(tmpdir; n=NMESH, boxsize=boxsize, ilpt=ILPT, nbuff=NBUFF,
                  ioutshear=IOUTSHEAR, ievol=IEVOL, z_max=Z_MAX,
                  cenx=(IEVOL == 1 ? boxsize/2 : 0.0),
                  ceny=(IEVOL == 1 ? boxsize/2 : 0.0),
                  cenz=(IEVOL == 1 ? boxsize/2 : 0.0))

# Warm both paths first (JIT + plan build)
println(">> warmup CPU..."); flush(stdout)
_ = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=1, verbose=false, use_gpu=false)
println(">> warmup GPU..."); flush(stdout)
_ = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=1, verbose=false, use_gpu=true)

println("\n>> CPU profile..."); flush(stdout)
t_cpu = @elapsed h_cpu = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=42,
                                                        verbose=false, use_gpu=false,
                                                        profile=true)
@printf "CPU total wall: %.2f s  (%d halos)\n" t_cpu length(h_cpu)

println("\n>> GPU profile..."); flush(stdout)
t_gpu = @elapsed h_gpu = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=42,
                                                        verbose=false, use_gpu=true,
                                                        profile=true)
@printf "GPU total wall: %.2f s  (%d halos)\n" t_gpu length(h_gpu)
@printf "Speedup: %.2fx   ΔN=%d\n" (t_cpu / t_gpu) (length(h_gpu) - length(h_cpu))

rm(tmpdir; recursive=true)
