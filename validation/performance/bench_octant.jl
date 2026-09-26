#!/usr/bin/env julia
# Benchmark driver: same pipeline as validation/websky_6144/run_gpu_octant.jl (split GPU
# multi-resolution → merge → finalize → write), with each stage timed separately and an
# optional warm-up run first so JIT compilation is excluded from the timed run.
#
#   julia --project=validation -t T validation/performance/bench_octant.jl \
#         <config.toml> <ndev> <outdir> [warmup_config.toml]
#
# Prints one machine-readable line:  BENCH label=... ndev=... t_pipeline=... (seconds)
# GPU memory is recorded outside Julia (nvidia-smi logging in the SLURM script) because the
# CUDA.jl pool caches allocations — device-level usage is the real footprint.
using TOML, Printf
using PeakPatch
using CUDA

function run_once(config_path, devices, outdir; verbose=true)
    config = TOML.parsefile(config_path)
    cfg = PipelineConfig(config)
    rc = config["run"]
    ntile = rc["ntile"]; cf = rc["coarse_factor"]; seed = get(rc, "seed", 12345)
    t = Dict{String,Float64}()

    t0 = time()
    halos = run_multitile_split(cfg; ntile=ntile, seed=seed, coarse_factor=cf, use_gpu=true,
                                devices=devices, verbose=verbose, profile=verbose)
    t["pipeline"] = time() - t0
    npre = length(halos)

    t0 = time()
    halos = length(halos) > 1 ? merge_catalog(halos; verbose=verbose) : halos
    t["merge"] = time() - t0

    t0 = time()
    cosmo = CosmologyParams(cfg.Omx + cfg.OmB, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)
    obs = (Float64(cfg.cenx), Float64(cfg.ceny), Float64(cfg.cenz))
    isempty(halos) || (halos = finalize_eulerian(halos, cosmo, obs; ievol=cfg.ievol, z_out=Float64(cfg.z_out)))
    t["finalize"] = time() - t0

    t0 = time()
    path = joinpath(outdir, basename(cfg.fileout))
    RTHLmax = isempty(halos) ? 0f0 : maximum(h.RTHL for h in halos)
    write_pksc(path, halos, RTHLmax, Float32(cfg.z_out))
    t["write"] = time() - t0
    return t, npre, length(halos), path
end

function main()
    length(ARGS) >= 3 || error("usage: bench_octant.jl <config> <ndev> <outdir> [warmup_config]")
    config_path, ndev, outdir = ARGS[1], parse(Int, ARGS[2]), ARGS[3]
    mkpath(outdir)
    @assert CUDA.functional()
    length(CUDA.devices()) >= ndev || error("need $ndev GPUs, have $(length(CUDA.devices()))")
    devices = collect(0:ndev-1)
    gpu = CUDA.name(CUDA.device())
    @info "bench" config=config_path ndev threads=Threads.nthreads() gpu

    if length(ARGS) >= 4
        tw = time(); run_once(ARGS[4], devices, outdir; verbose=false)
        @info "warm-up done" seconds=round(time() - tw; digits=1)
        GC.gc(true); CUDA.reclaim()
    end

    t, npre, npost, path = run_once(config_path, devices, outdir)
    label = splitext(basename(config_path))[1]
    @printf("BENCH label=%s gpu=\"%s\" ndev=%d threads=%d npre=%d npost=%d t_pipeline=%.1f t_merge=%.1f t_finalize=%.1f t_write=%.1f bytes=%d\n",
            label, gpu, ndev, Threads.nthreads(), npre, npost, t["pipeline"], t["merge"],
            t["finalize"], t["write"], filesize(path))
    for d in devices
        CUDA.device!(d)
        @printf("CUDA_POOL dev=%d used_GB=%.2f\n", d, CUDA.used_memory() / 2^30)
    end
end

main()
