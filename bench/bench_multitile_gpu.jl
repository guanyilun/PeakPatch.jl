#!/usr/bin/env julia
# End-to-end benchmark: run_multitile_split CPU vs GPU.
#
# Run (one GPU needed):
#   CUDA_VISIBLE_DEVICES=1 julia --project=/tmp/pp_gpu_sandbox --threads=64 \
#     bench/bench_multitile_gpu.jl
#
# Env knobs:
#   NTILE (default 2), NMESH (default 192), ILPT (default 2), IOUTSHEAR (default 0)

using PeakPatch
using CUDA
using Printf

const NTILE = parse(Int, get(ENV, "NTILE", "2"))
const NMESH = parse(Int, get(ENV, "NMESH", "192"))
const ILPT  = parse(Int, get(ENV, "ILPT",  "2"))
const IOUTSHEAR = parse(Int, get(ENV, "IOUTSHEAR", "0"))
const NBUFF = parse(Int, get(ENV, "NBUFF", "8"))
const COARSE_FACTOR = parse(Int, get(ENV, "COARSE_FACTOR", "5"))

function make_config(dir; n, boxsize, z=0.0, ilpt=2, nbuff=8,
                     ioutshear=0, rmax2rs=0.0)
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
        z_out=Float64(z), ilpt=ilpt, ioutshear=ioutshear, rmax2rs=rmax2rs,
        ievol=0, z_max=0.0, cenx=0.0, ceny=0.0, cenz=0.0,
        Omx=0.261, OmB=0.049, Omvac=0.69, h=0.68,
        NonGauss=0, fNL=0.0, wsmooth=1,
    )
end

println("\n" * "="^60)
println("End-to-end multi-tile: CPU vs GPU")
println("="^60)
nsub = NMESH - 2 * NBUFF
N = nsub * NTILE + 2 * NBUFF
@printf "ntile=%d nmesh=%d nsub=%d N=%d  ilpt=%d ioutshear=%d\n" NTILE NMESH nsub N ILPT IOUTSHEAR
@printf "padded isolated FFT: %d³ (%.2f GB Float32)\n" (2*NMESH) (Float64(2*NMESH)^3 * 4 / 1e9)
@printf "Threads (Julia/FFTW): %d\n" Threads.nthreads()
@printf "GPU: %s\n" CUDA.functional() ? CUDA.name(CUDA.device()) : "n/a"
println()

tmpdir = mktempdir()
# boxsize scaled to keep comparable fluctuations
boxsize = Float64(N) * 226.67 / 68.0
cfg = make_config(tmpdir; n=NMESH, boxsize=boxsize, ilpt=ILPT,
                    nbuff=NBUFF, ioutshear=IOUTSHEAR)

# Warm up (cuFFT plan build, JIT for both paths). Reuse a tiny ntile=1 config
# just to compile the GPU path once.
println(">> warming up GPU path..."); flush(stdout)
_ = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=1, verbose=true,
                                   use_gpu=true)
println(">> warming up CPU path..."); flush(stdout)
_ = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=1, verbose=true,
                                   use_gpu=false)

println("\n>> timing CPU run (seed=42)..."); flush(stdout)
t_cpu = @elapsed halos_cpu = PeakPatch.run_multitile_split(
    cfg; ntile=NTILE, seed=42, verbose=true, use_gpu=false)

println("\n>> timing GPU run (seed=42)..."); flush(stdout)
t_gpu = @elapsed halos_gpu = PeakPatch.run_multitile_split(
    cfg; ntile=NTILE, seed=42, verbose=true, use_gpu=true)

println()
@printf "CPU: %.2f s  (%d halos)\n" t_cpu length(halos_cpu)
@printf "GPU: %.2f s  (%d halos)\n" t_gpu length(halos_gpu)
@printf "Speedup: %.2fx\n" (t_cpu / t_gpu)
@printf "ΔN_halos: %d (%.2f%%)\n" (length(halos_gpu) - length(halos_cpu)) (
    100 * abs(length(halos_gpu) - length(halos_cpu)) / max(length(halos_cpu), 1))

# Compare top halos by sorted position
sort_key = h -> (h.x, h.y, h.z)
sa = sort(halos_cpu; by=sort_key)
sb = sort(halos_gpu; by=sort_key)
ncmp = min(length(sa), length(sb), 5)
println("\nTop $ncmp halos (sorted by position):")
for i in 1:ncmp
    da = sa[i]; db = sb[i]
    @printf "  %2d: CPU=(%.3f,%.3f,%.3f) R=%.3f  |  GPU=(%.3f,%.3f,%.3f) R=%.3f\n" i da.x da.y da.z da.RTHL db.x db.y db.z db.RTHL
end
println("="^60)

rm(tmpdir; recursive=true)
