#!/usr/bin/env julia
#
# PeakPatch driver — reads a TOML config file and runs the halo-finding pipeline:
#   field + peaks (serial / multi-tile / MPI / GPU multi-resolution)
#   → merge (Lagrangian exclusion) → finalize (Eulerian positions + km/s velocities)
#   → optional abundance matching → pksc / HDF5 output.
#
# Usage:
#   julia --project=. bin/peakpatch.jl config.toml [--verbose]
#   julia --project=validation -t 32 bin/peakpatch.jl config.toml     # [run] use_gpu = true
#   mpiexec -np N julia --project=. bin/peakpatch.jl config.toml     # [run] use_mpi = true
#
# GPU runs need an environment that provides CUDA.jl (CUDA is a weak dependency of
# PeakPatch): use the script environment `validation/Project.toml`.
#

using TOML
using PeakPatch

"""Load CUDA (triggering ext/CUDAExt.jl) and return the device ids to use."""
function _init_gpu(run_cfg)
    try
        @eval Main using CUDA
    catch err
        error("[run] use_gpu = true needs CUDA.jl in the active environment " *
              "(it is a weak dependency). Run with --project=validation. ($err)")
    end
    Base.invokelatest(() -> Main.CUDA.functional()) || error("CUDA is not functional on this node")
    ndev = Base.invokelatest(() -> length(Main.CUDA.devices()))
    devices = get(run_cfg, "devices", "all")
    devices = devices == "all" ? collect(0:ndev-1) : Int.(devices)
    isempty(devices) && error("no CUDA devices found")
    length(devices) > Threads.nthreads() &&
        @warn "$(length(devices)) GPUs but $(Threads.nthreads()) threads; start julia with -t >= $(length(devices))"
    return devices
end

"""Sky fraction covered by a lightcone catalog: an observer at a grid corner sees one
octant (1/8), at the grid centre the full sky (1). Anything else must be explicit."""
function _infer_fsky(cfg, ntile)
    nsub = cfg.n - 2 * cfg.nbuff
    half = (nsub * ntile + 2 * cfg.nbuff) * (cfg.boxsize / cfg.n) / 2   # boxsize_full / 2
    obs = (cfg.cenx, cfg.ceny, cfg.cenz)
    all(o -> abs(o) < 1e-3 * half, obs) && return 1.0
    all(o -> isapprox(abs(o), half; rtol=0.02), obs) && return 1 / 8
    error("cannot infer the sky fraction for observer $obs (box half-size $(round(half; digits=1)) " *
          "Mpc/h); set [abundance_match] fsky explicitly")
end

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
    use_gpu = get(run_cfg, "use_gpu", false)
    multires = get(run_cfg, "multires", use_gpu)   # GPU implies the multi-resolution path
    coarse_factor = get(run_cfg, "coarse_factor", 0) # 0 = auto; ≈32 for velocity/potential products
    merge   = get(run_cfg, "merge", true)
    finalize = get(run_cfg, "finalize", true)       # Eulerian positions + km/s velocities
    gen_table  = get(run_cfg, "generate_table", false)
    ode_solver = Symbol(get(run_cfg, "ode_solver", "rk4"))
    use_lcg       = get(run_cfg, "use_lcg", false)
    fortran_compat = get(run_cfg, "fortran_compat", false)
    lowmem        = get(run_cfg, "lowmem", false)

    use_mpi && (use_gpu || multires) && error("[run] use_mpi cannot be combined with use_gpu/multires")
    use_gpu && !multires && error("[run] use_gpu = true requires the multi-resolution path (multires = true)")
    multires && (use_lcg || fortran_compat) &&
        error("[run] multires/use_gpu uses the Threefry RNG; use_lcg/fortran_compat are not supported")

    # Output parameters
    out_cfg = get(config, "output", Dict{String,Any}())
    format  = get(out_cfg, "format", "pksc")
    outdir  = get(out_cfg, "path", ".")
    isdir(outdir) || mkpath(outdir)

    # Abundance matching (optional; lightcone catalogs only)
    am_cfg = get(config, "abundance_match", Dict{String,Any}())
    do_am  = get(am_cfg, "enabled", false)
    do_am && cfg.ievol != 1 && error("[abundance_match] needs a lightcone run (ievol = 1): halos are matched in redshift shells")
    do_am && !finalize && error("[abundance_match] needs finalize = true (it bins halos by their Eulerian distance)")

    Om_total = cfg.Omx + cfg.OmB
    cosmo = CosmologyParams(Om_total, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)

    verbose && @info "PeakPatch driver" config=config_path ntile seed format use_gpu multires use_mpi

    # ---- Optional: generate collapse table on-the-fly ----
    if gen_table
        ep = EllipsoidParams(cosmo; solver=ode_solver)
        tp = CollapseTableParams()
        verbose && @info "Generating collapse table (solver=$ode_solver)..."
        table = make_table_threaded(ep, tp; verbose=verbose)
        write_homeltab(cfg.tabfile, table, tp)
        verbose && @info "Wrote collapse table: $(cfg.tabfile)"
    end

    # ---- Run pipeline ----
    is_rank0 = true
    t0 = time()
    halos = if use_mpi
        # MPI path: load MPI + PencilFFTs to trigger MPIExt package extension,
        # then @invokelatest to call methods defined in the new world age.
        @eval using MPI
        @eval using PencilFFTs
        Base.invokelatest(MPI.Init)
        is_rank0 = Base.invokelatest(MPI.Comm_rank, MPI.COMM_WORLD) == 0
        is_rank0 && @info "Running MPI multi-tile pipeline (ntile=$ntile)"
        Base.invokelatest(run_multitile_mpi, cfg; ntile=ntile, seed=seed, verbose=verbose, lowmem=lowmem)
    elseif multires
        devices = use_gpu ? _init_gpu(run_cfg) : nothing
        @info "Running multi-resolution pipeline" ntile coarse_factor use_gpu devices
        Base.invokelatest(run_multitile_split, cfg; ntile=ntile, seed=seed, verbose=verbose,
                          coarse_factor=coarse_factor, use_gpu=use_gpu, devices=devices)
    elseif ntile > 1 && lowmem
        run_multitile_lowmem(cfg; ntile=ntile, seed=seed, verbose=verbose,
                              use_lcg=use_lcg, fortran_compat=fortran_compat)
    elseif ntile > 1
        run_multitile(cfg; ntile=ntile, seed=seed, verbose=verbose,
                      use_lcg=use_lcg, fortran_compat=fortran_compat)
    else
        run_tile(cfg; seed=seed, verbose=verbose,
                 use_lcg=use_lcg, fortran_compat=fortran_compat)
    end

    # Post-processing and output only on rank 0 (all halos gathered there)
    if !is_rank0
        Base.invokelatest(MPI.Barrier, MPI.COMM_WORLD)
        return
    end

    verbose && @info "Pipeline complete: $(length(halos)) halos in $(round((time() - t0) / 60; digits=1)) min"

    # ---- Merge ----
    if merge && length(halos) > 1
        verbose && @info "Running merger (Lagrangian exclusion)..."
        halos = merge_catalog(halos; verbose=verbose)
        verbose && @info "After merge: $(length(halos)) halos"
    end

    # ---- Finalize: Lagrangian q + displacement → Eulerian position + km/s velocity ----
    # (the Fortran merge_pkvd conversion). Without it the catalog holds peak positions
    # and raw displacements, NOT what readers of a pksc catalog expect.
    if finalize && !isempty(halos)
        verbose && @info "Finalizing (Eulerian position + km/s velocity)..."
        halos = finalize_eulerian(halos, cosmo, (Float64(cfg.cenx), Float64(cfg.ceny), Float64(cfg.cenz));
                                  ievol=cfg.ievol, z_out=Float64(cfg.z_out))
    elseif !finalize
        @warn "[run] finalize = false: catalog holds Lagrangian positions and displacements, not Eulerian positions/km/s velocities"
    end

    catalogs = [("", halos)]

    # ---- Abundance matching ----
    if do_am && !isempty(halos)
        # The field-generation pkfile is P(k)/(2π)³; σ(M) for the target HMF needs the
        # physical P(k), which is exactly (2π)³ × that (checked on pk_websky.dat).
        pkfield = PeakPatch.PowerSpectrum.load_pk(get(am_cfg, "pkfile", cfg.pkfile))
        pk_phys = haskey(am_cfg, "pkfile") ? pkfield : (k -> (2π)^3 * pkfield(k))
        z_max_am = Float64(get(am_cfg, "z_max", cfg.z_max))
        fsky = Float64(get(am_cfg, "fsky", 0.0))
        fsky > 0 || (fsky = _infer_fsky(cfg, ntile))
        hmf = Symbol(get(am_cfg, "hmf", "tinker"))
        obs = (Float64(cfg.cenx), Float64(cfg.ceny), Float64(cfg.cenz))
        verbose && @info "Abundance matching" hmf z_max_am fsky obs
        table = build_abundance_table(halos, cosmo, pk_phys; hmf=hmf, z_max=z_max_am,
                                      nzbins=max(2, round(Int, z_max_am / 0.1)), obs=obs,   # ≥2: linear lookup in z
                                      fsky=fsky, verbose=verbose)
        halos_am = abundance_match(halos, table, cosmo; obs=obs)
        catalogs = get(am_cfg, "keep_raw", true) ? [("", halos), ("_AM", halos_am)] : [("", halos_am)]
    end

    # ---- Write output ----
    z_out = Float32(cfg.z_out)
    for (suffix, hs) in catalogs
        base = basename(cfg.fileout)
        ext  = something(match(r"\.\w+$", base), (match = ".pksc",)).match
        stem = replace(base, r"\.\w+$" => "") * suffix
        RTHLmax = isempty(hs) ? Float32(0) : maximum(h.RTHL for h in hs)

        if format in ("pksc", "both")
            pksc_path = joinpath(outdir, stem * ext)
            write_pksc(pksc_path, hs, RTHLmax, z_out)
            verbose && @info "Wrote pksc: $pksc_path ($(length(hs)) halos)"
        end

        if format in ("hdf5", "both")
            @eval using HDF5
            hdf5_path = joinpath(outdir, stem * ".h5")
            Base.invokelatest(write_catalog_hdf5, hdf5_path, hs, cosmo)
            verbose && @info "Wrote HDF5: $hdf5_path ($(length(hs)) halos)"
        end
    end

    @info "Done: $(length(halos)) halos written (format=$format$(do_am ? ", + abundance-matched" : ""))"

    use_mpi && Base.invokelatest(MPI.Barrier, MPI.COMM_WORLD)
end

main()
