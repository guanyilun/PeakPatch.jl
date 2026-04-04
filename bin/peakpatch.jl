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

    # Build SimParams from TOML
    sp = SimParams(config)

    # Run parameters
    run_cfg = get(config, "run", Dict{String,Any}())
    seed    = get(run_cfg, "seed", 42)
    ntile   = get(run_cfg, "ntile", 1)
    use_mpi = get(run_cfg, "use_mpi", false)
    merge   = get(run_cfg, "merge", true)

    # Output parameters
    out_cfg = get(config, "output", Dict{String,Any}())
    format  = get(out_cfg, "format", "pksc")
    outdir  = get(out_cfg, "path", ".")
    isdir(outdir) || mkpath(outdir)

    verbose && @info "PeakPatch driver" config=config_path ntile seed format

    # ---- Run pipeline ----
    halos = if use_mpi
        # MPI path: requires MPI.jl + PencilFFTs + PencilArrays loaded
        @info "Running MPI multi-tile pipeline (ntile=$ntile)"
        run_multitile_mpi(sp; ntile=ntile, seed=seed)
    elseif ntile > 1
        run_multitile(sp; ntile=ntile, seed=seed, verbose=verbose)
    else
        run_tile(sp; seed=seed, verbose=verbose)
    end

    verbose && @info "Pipeline complete: $(length(halos)) halos"

    # ---- Merge ----
    if merge && length(halos) > 1
        verbose && @info "Running merger (exclusion + volume reduction)..."
        halos = merge_catalog(halos; verbose=verbose)
        verbose && @info "After merge: $(length(halos)) halos"
    end

    # ---- Write output ----
    z_out = Float32(sp.global_redshift)
    RTHLmax = isempty(halos) ? Float32(0) : maximum(h.RTHL for h in halos)

    if format in ("pksc", "both")
        pksc_path = joinpath(outdir, basename(sp.fileout))
        write_pksc(pksc_path, halos, RTHLmax, z_out)
        verbose && @info "Wrote pksc: $pksc_path ($(length(halos)) halos)"
    end

    if format in ("hdf5", "both")
        hdf5_path = joinpath(outdir, replace(basename(sp.fileout), r"\.\w+$" => ".h5"))
        # Build cosmology for coordinate conversion
        Om_total = Float64(sp.Omx) + Float64(sp.OmB)
        cosmo = CosmologyParams(Om_total, Float64(sp.OmB), Float64(sp.Omvac),
                                Float64(sp.h), 0.965, 0.808)
        write_catalog_hdf5(hdf5_path, halos, cosmo)
        verbose && @info "Wrote HDF5: $hdf5_path ($(length(halos)) halos)"
    end

    @info "Done: $(length(halos)) halos written (format=$format)"
end

main()
