#!/usr/bin/env julia
#
# GPU driver for Websky 6144^3 multi-resolution on Killarney.
# Reads a TOML config and calls run_multitile_split with GPU acceleration.
#
# Usage:
#   julia --project=validation -t 32 validation/websky_6144/run_gpu_octant.jl \
#       validation/websky_6144/config_websky_6144_oct000.toml
#
# On SLURM: invoked by run_websky_6144_killarney.slurm

using TOML
using PeakPatch
using CUDA

function main()
    if isempty(ARGS)
        println(stderr, "Usage: julia --project=validation run_gpu_octant.jl <config.toml>")
        exit(1)
    end

    config_path = ARGS[1]
    if !isfile(config_path)
        println(stderr, "Error: config not found: $config_path")
        exit(1)
    end

    config = TOML.parsefile(config_path)
    cfg = PipelineConfig(config)
    # Reproduction metadata, captured first so a dirty code tree fails before any GPU time is spent
    # (PEAKPATCH_ALLOW_DIRTY=1 embeds the diff instead). Written next to the catalog as <catalog>.prov.toml.
    prov = capture_provenance(config_path; scripts=[@__FILE__])

    # Run parameters
    run_cfg = get(config, "run", Dict{String,Any}())
    seed       = get(run_cfg, "seed", 42)
    ntile      = get(run_cfg, "ntile", 16)
    coarse_factor = get(run_cfg, "coarse_factor", 4)
    gen_table  = get(run_cfg, "generate_table", false)
    ode_solver = Symbol(get(run_cfg, "ode_solver", "rk4"))
    # Fortran merge_pkvd volume reduction after exclusion (off in every catalog up to v3fs;
    # validation/paper/CLUSTERING_EXCESS_2026-09-28.md)
    vol_red    = Bool(get(run_cfg, "volume_reduction", false))

    # Output parameters
    out_cfg = get(config, "output", Dict{String,Any}())
    outdir  = get(out_cfg, "path", ".")
    isdir(outdir) || mkpath(outdir)

    # ---- Optional: generate collapse table ----
    if gen_table
        Om_total = cfg.Omx + cfg.OmB
        cosmo_tab = CosmologyParams(Om_total, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)
        ep = EllipsoidParams(cosmo_tab; solver=ode_solver)
        tp = CollapseTableParams()
        @info "Generating collapse table (solver=$ode_solver)..."
        table = make_table_threaded(ep, tp; verbose=true)
        write_homeltab(cfg.tabfile, table, tp)
        @info "Wrote collapse table: $(cfg.tabfile)"
    end

    # ---- Initialize GPU ----
    @assert CUDA.functional() "CUDA is not functional on this node"
    ndev = length(CUDA.devices())
    @info "CUDA devices available: $ndev"

    if ndev < 1
        error("No CUDA devices found")
    end

    devices = collect(0:ndev-1)
    @info "Using $ndev GPUs: $devices"
    prov["env.gpu"] = join(unique(CUDA.name(d) for d in CUDA.devices()), ", ") * " × $ndev"
    prov["env.cuda_runtime"] = string(CUDA.runtime_version())
    prov["env.cuda_driver"] = string(CUDA.driver_version())

    # ---- Geometry summary ----
    nmesh = cfg.n
    nbuff = cfg.nbuff
    nsub, N = grid_layout(cfg, ntile)
    alatt = cfg.boxsize / nmesh
    boxsize_full = N * alatt

    @info "Websky 6144^3 GPU multi-resolution" N=N ntile=ntile nmesh=nmesh nbuff=nbuff nsub=nsub boxsize_full=round(boxsize_full; digits=1) periodic_cores=cfg.periodic_cores seed=seed coarse_factor=coarse_factor ievol=cfg.ievol z_max=cfg.z_max ilpt=cfg.ilpt ioutshear=cfg.ioutshear peak_threshold_per_tile=cfg.peak_threshold_per_tile volume_reduction=vol_red

    # ---- Run pipeline ----
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

    @info "Pipeline complete" halos=length(halos) elapsed_min=round(elapsed/60; digits=1)

    # ---- Merge (Lagrangian-space exclusion) ----
    if length(halos) > 1
        @info "Running merger" exclusion=true volume_reduction=vol_red
        halos = merge_catalog(halos; verbose=true, volume_reduction=vol_red)
        @info "After merge: $(length(halos)) halos"
    end

    # ---- Finalize: Eulerian positions + km/s peculiar velocities (merge_pkvd conversion) ----
    if !isempty(halos)
        @info "Finalizing (Eulerian position + km/s velocity)..."
        cosmo_fin = CosmologyParams(cfg.Omx + cfg.OmB, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)
        obs_fin = (Float64(cfg.cenx), Float64(cfg.ceny), Float64(cfg.cenz))
        halos = finalize_eulerian(halos, cosmo_fin, obs_fin; ievol=cfg.ievol, z_out=Float64(cfg.z_out))
    end

    # ---- Write output ----
    z_out = Float32(cfg.z_out)
    RTHLmax = isempty(halos) ? Float32(0) : maximum(h.RTHL for h in halos)

    pksc_path = joinpath(outdir, basename(cfg.fileout))
    write_pksc(pksc_path, halos, RTHLmax, z_out)
    @info "Wrote catalog: $pksc_path ($(length(halos)) halos)"
    prov["run.n_halos"] = string(length(halos))
    prov["run.pipeline_minutes"] = string(round(elapsed / 60; digits=2))
    prov["run.product"] = "raw catalog (merged + finalized; Eulerian positions, km/s velocities; masses before AM)"
    side = write_provenance_sidecar(pksc_path, prov)
    @info "Wrote provenance: $side"

    @info "Done: $(length(halos)) halos in $(round(elapsed/60; digits=1)) minutes"
end

main()
