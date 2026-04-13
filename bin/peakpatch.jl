#!/usr/bin/env julia
#
# PeakPatch driver — reads a TOML config file and runs the halo-finding pipeline.
#
# Usage:
#   julia --project=. bin/peakpatch.jl config.toml
#   julia --project=. bin/peakpatch.jl config.toml --verbose
#
# For MPI runs:
#   mpiexec -np N julia --project=. bin/peakpatch.jl config.toml
#

using TOML
using PeakPatch

function main()
    if isempty(ARGS)
        println(stderr, "Usage: julia bin/peakpatch.jl <config.toml> [--verbose]")
        exit(1)
    end

    config_path = ARGS[1]
    verbose = "--verbose" in ARGS || "-v" in ARGS

    if !isfile(config_path)
        println(stderr, "Error: config file not found: $config_path")
        exit(1)
    end

    config = TOML.parsefile(config_path)

    # Build PipelineConfig from TOML
    cfg = PipelineConfig(config)

    # Run parameters (not part of PipelineConfig)
    run_cfg = get(config, "run", Dict{String,Any}())
    seed    = get(run_cfg, "seed", 42)
    ntile   = get(run_cfg, "ntile", 1)
    use_mpi = get(run_cfg, "use_mpi", false)
    merge   = get(run_cfg, "merge", true)
    gen_table  = get(run_cfg, "generate_table", false)
    ode_solver = Symbol(get(run_cfg, "ode_solver", "rk4"))
    use_lcg       = get(run_cfg, "use_lcg", false)
    fortran_compat = get(run_cfg, "fortran_compat", false)

    # Output parameters
    out_cfg = get(config, "output", Dict{String,Any}())
    format  = get(out_cfg, "format", "pksc")
    outdir  = get(out_cfg, "path", ".")
    isdir(outdir) || mkpath(outdir)

    verbose && @info "PeakPatch driver" config=config_path ntile seed format

    # ---- Optional: generate collapse table on-the-fly ----
    if gen_table
        Om_total = cfg.Omx + cfg.OmB
        cosmo_tab = CosmologyParams(Om_total, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)
        ep = EllipsoidParams(cosmo_tab; solver=ode_solver)
        tp = CollapseTableParams()
        verbose && @info "Generating collapse table (solver=$ode_solver)..."
        table = make_table_threaded(ep, tp; verbose=verbose)
        write_homeltab(cfg.tabfile, table, tp)
        verbose && @info "Wrote collapse table: $(cfg.tabfile)"
    end

    # ---- Run pipeline ----
    is_rank0 = true
    halos = if use_mpi
        # MPI path: requires MPI.jl loaded at top level before calling main()
        @eval using MPI
        MPI.Initialized() || MPI.Init()
        is_rank0 = MPI.Comm_rank(MPI.COMM_WORLD) == 0
        is_rank0 && @info "Running MPI multi-tile pipeline (ntile=$ntile)"
        run_multitile_mpi(cfg; ntile=ntile, seed=seed, verbose=verbose)
    elseif ntile > 1
        run_multitile(cfg; ntile=ntile, seed=seed, verbose=verbose,
                      use_lcg=use_lcg, fortran_compat=fortran_compat)
    else
        run_tile(cfg; seed=seed, verbose=verbose,
                 use_lcg=use_lcg, fortran_compat=fortran_compat)
    end

    # Post-processing and output only on rank 0 (all halos gathered there)
    if !is_rank0
        MPI.Barrier(MPI.COMM_WORLD)
        return
    end

    verbose && @info "Pipeline complete: $(length(halos)) halos"

    # ---- Merge ----
    if merge && length(halos) > 1
        verbose && @info "Running merger (exclusion + volume reduction)..."
        halos = merge_catalog(halos; verbose=verbose)
        verbose && @info "After merge: $(length(halos)) halos"
    end

    # ---- Write output ----
    z_out = Float32(cfg.z_out)
    RTHLmax = isempty(halos) ? Float32(0) : maximum(h.RTHL for h in halos)

    if format in ("pksc", "both")
        pksc_path = joinpath(outdir, basename(cfg.fileout))
        write_pksc(pksc_path, halos, RTHLmax, z_out)
        verbose && @info "Wrote pksc: $pksc_path ($(length(halos)) halos)"
    end

    if format in ("hdf5", "both")
        hdf5_path = joinpath(outdir, replace(basename(cfg.fileout), r"\.\w+$" => ".h5"))
        Om_total = cfg.Omx + cfg.OmB
        cosmo = CosmologyParams(Om_total, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)
        write_catalog_hdf5(hdf5_path, halos, cosmo)
        verbose && @info "Wrote HDF5: $hdf5_path ($(length(halos)) halos)"
    end

    @info "Done: $(length(halos)) halos written (format=$format)"

    use_mpi && MPI.Barrier(MPI.COMM_WORLD)
end

main()
