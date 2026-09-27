#!/usr/bin/env julia
# GPU driver for the field-matter lightcone κ/mass maps (Phase A of
# docs/field_lightcone_plan.md). Paints ALL lattice cells (full matter — halo + field in
# one map) of one octant onto HEALPix maps via run_multitile_fieldmap.
#
# Usage: julia --project=validation -t 32 run_fieldmap_octant.jl <config.toml> [nside] [outdir] [gpu_paint|-] [exclude=<catalog.pksc>]
# Optional 4th arg "gpu_paint" switches to Phase-B on-device pixelization (own RING
# ang2pix + device atomics; enables Nside 4096 / subdiv 5 at low cost); "-" skips.
# Optional 5th arg "exclude=<path>" masks the Lagrangian spheres of the catalog halos
# (Websky field component). The catalog must be OLD-convention (x,y,z = Lagrangian);
# output files then get a "_field" tag suffix.
using TOML, PeakPatch, CUDA, Healpix, Printf, Statistics

# Bulk-read Lagrangian centers + RTHL from a 33-field ext pksc (26 GB → 4 columns).
function read_exclusion_catalog(path::String)
    open(path, "r") do io
        nhalo = Int(read(io, Int32)); read(io, Float32); read(io, Float32)
        nfields = (filesize(path) - 12) ÷ (nhalo * 4)
        @assert nfields in (11, 33) "unexpected pksc record size: $nfields fields"
        hx = Vector{Float32}(undef, nhalo); hy = similar(hx); hz = similar(hx)
        hR = similar(hx)
        chunk = 4_000_000
        buf = Matrix{Float32}(undef, nfields, chunk)
        done = 0
        while done < nhalo
            n = min(chunk, nhalo - done)
            n < chunk && (buf = Matrix{Float32}(undef, nfields, n))
            read!(io, buf)
            @views begin
                hx[done+1:done+n] = buf[1, 1:n]
                hy[done+1:done+n] = buf[2, 1:n]
                hz[done+1:done+n] = buf[3, 1:n]
                hR[done+1:done+n] = buf[7, 1:n]
            end
            done += n
        end
        return (x=hx, y=hy, z=hz, R=hR)
    end
end

function main()
    config_path = ARGS[1]
    nside = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2048
    outdir = length(ARGS) >= 3 ? ARGS[3] : "/home/yguan/scratch/websky_6144/fieldmaps"
    gpu_paint = length(ARGS) >= 4 && ARGS[4] == "gpu_paint"
    excl_path = ""
    for a in ARGS[4:end]
        startswith(a, "exclude=") && (excl_path = a[9:end])
    end
    isdir(outdir) || mkpath(outdir)

    config = TOML.parsefile(config_path)
    cfg = PipelineConfig(config)
    run_cfg = get(config, "run", Dict{String,Any}())
    seed = get(run_cfg, "seed", 42)
    ntile = get(run_cfg, "ntile", 16)
    coarse_factor = get(run_cfg, "coarse_factor", 4)
    tag = replace(basename(config_path), "config_" => "", ".toml" => "")

    exclude_halos = nothing
    if !isempty(excl_path)
        @info "reading exclusion catalog (Lagrangian centers + RTHL)..." excl_path
        exclude_halos = read_exclusion_catalog(excl_path)
        @info "exclusion catalog" nhalo=length(exclude_halos.x) Rmax=maximum(exclude_halos.R)
        tag *= "_field"
    end

    @assert CUDA.functional() "CUDA is not functional on this node"
    devices = collect(0:length(CUDA.devices())-1)
    res = Resolution(nside)
    npix = nside2npix(nside)
    v2p = (x, y, z) -> Healpix.vec2pixRing(res, x, y, z)

    # chi_star: Websky/pks2map hardwire 14.2 Gpc for the CMB source plane -> Mpc/h
    chi_star = 14200.0 * cfg.h

    @info "fieldmap octant" tag nside npix devices seed ntile coarse_factor chi_star gpu_paint
    t0 = time()
    maps = run_multitile_fieldmap(cfg; ntile=ntile, seed=seed, coarse_factor=coarse_factor,
                                  npix=npix, vec2pix=v2p,
                                  kernels=[:kappa, :mass, :tau, :ksz, :isw],
                                  chi_star=chi_star, subdiv_max=(gpu_paint ? 5 : 3),
                                  exclude_halos=exclude_halos,
                                  use_gpu=true, devices=devices,
                                  gpu_paint=gpu_paint, nside=(gpu_paint ? nside : 0),
                                  verbose=true)
    elapsed = time() - t0
    @info "painting complete" elapsed_min = round(elapsed / 60; digits=1)

    # ---- diagnostics ----
    mtot = sum(maps[:mass])
    kmean = mean(maps[:kappa])
    kon = count(>(0), maps[:mass]) / npix
    @printf("total painted mass: %.6e Msun/h\n", mtot)
    @printf("mean kappa (full sky, unsubtracted): %.6e\n", kmean)
    @printf("sky fraction covered: %.4f (octant = 0.125)\n", kon)

    # ---- save ----
    for (kern, m) in maps
        hm = HealpixMap{Float64,RingOrder}(m)
        path = joinpath(outdir, "$(kern)_$(tag)_nside$(nside).fits")
        Healpix.saveToFITS(hm, "!" * path, typechar="E")
        @info "wrote $path"
    end
end

main()
