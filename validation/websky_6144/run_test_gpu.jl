#!/usr/bin/env julia
#
# Small-scale GPU pipeline test.
# Runs a 120^3 box (ntile=2, nmesh=68) on a single GPU.
#
# Usage:
#   julia --project=validation -t 8 validation/websky_6144/run_test_gpu.jl

using TOML
using PeakPatch
using CUDA

function main()
    config_path = joinpath(@__DIR__, "config_test_gpu_small.toml")
    config = TOML.parsefile(config_path)
    cfg = PipelineConfig(config)

    run_cfg = get(config, "run", Dict{String,Any}())
    seed = get(run_cfg, "seed", 42)
    ntile = get(run_cfg, "ntile", 2)
    coarse_factor = get(run_cfg, "coarse_factor", 5)

    out_cfg = get(config, "output", Dict{String,Any}())
    outdir = get(out_cfg, "path", ".")
    isdir(outdir) || mkpath(outdir)

    # ---- Init GPU ----
    @assert CUDA.functional() "CUDA not functional"
    ndev = length(CUDA.devices())
    @info "CUDA OK — $ndev device(s): $(CUDA.device!(0) |> CUDA.name)"

    # Use just 1 GPU for the test
    devices = [0]

    # ---- Geometry ----
    nmesh = cfg.n
    nbuff = cfg.nbuff
    nsub = nmesh - 2 * nbuff
    N = nsub * ntile + 2 * nbuff
    @info "Test config" N=N ntile=ntile nmesh=nmesh nbuff=nbuff nsub=nsub

    # ---- Run ----
    t0 = time()
    halos = run_multitile_split(cfg;
        ntile=ntile,
        seed=seed,
        coarse_factor=coarse_factor,
        use_gpu=true,
        devices=devices,
        verbose=true,
        profile=true)
    elapsed = time() - t0

    @info "Pipeline done" halos=length(halos) elapsed_s=round(elapsed; digits=1)

    # ---- Merge ----
    if length(halos) > 1
        halos = merge_catalog(halos; verbose=true)
        @info "After merge: $(length(halos)) halos"
    end

    # ---- Write output ----
    z_out = Float32(cfg.z_out)
    RTHLmax = isempty(halos) ? Float32(0) : maximum(h.RTHL for h in halos)
    pksc_path = joinpath(outdir, basename(cfg.fileout))
    write_pksc(pksc_path, halos, RTHLmax, z_out)

    @info "TEST PASSED — $(length(halos)) halos, $(round(elapsed; digits=1))s"
    @info "Output: $pksc_path"
end

main()
