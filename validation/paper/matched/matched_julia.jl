#!/usr/bin/env julia
# Same-field Julia-vs-Fortran comparison (MATCHED_FORTRAN_2026-10.md), Julia side.
# One snapshot box at z = 0.7, production cell (0.852 Mpc/h), nbuff 26, ntile 4, NON-periodic cores
# (the Fortran tiling): N = nsub*ntile + 2*nbuff = 1056.
#   mode "field" : write the exact global linear field δ(z=0) that run_multitile generates for this seed,
#                  in Fortran Fvec format (N³ Float32, x fastest) for hpkvd ireadfield = 1
#   mode "exact" : CPU run_multitile on that same global field (no multires split)
#   mode "split" : GPU run_multitile_split, same noise (Threefry counter → same white noise), cf 22 →
#                  coarse block 12 as in production (6144 / (16·32))
# Raw and merged (exclusion + volume reduction, Lagrangian + displacements) catalogs are written as pksc.
#   usage: julia --project=validation -t N matched_julia.jl <config.toml> <field|exact|split> <outdir>
using TOML, Printf
using PeakPatch
import PeakPatch.PowerSpectrum: load_pk
const MODE = ARGS[2]
MODE == "split" && @eval using CUDA

function main()
    cfgd = TOML.parsefile(ARGS[1]); cfg = PipelineConfig(cfgd); rc = cfgd["run"]
    outdir = ARGS[3]; mkpath(outdir)
    ntile = rc["ntile"]; seed = rc["seed"]
    nsub, N = grid_layout(cfg, ntile); alatt = cfg.boxsize / cfg.n; L = N * alatt
    @info "matched box" MODE N nsub alatt L
    if MODE == "field"
        δ = generate_grf(N, load_pk(cfg.pkfile), L, seed)
        fn = joinpath(outdir, "Fvec_jmatched")
        open(io -> write(io, δ), fn, "w")
        m = sum(Float64, δ) / length(δ); s = sqrt(sum(x -> (x - m)^2, δ) / length(δ))
        @printf("wrote %s (%d bytes, expect %d); mean %.3e sigma %.5f\n", fn, filesize(fn), 4N^3, m, s)
        return
    end
    raw = if MODE == "exact"
        run_multitile(cfg; ntile=ntile, seed=seed, verbose=true)
    else
        run_multitile_split(cfg; ntile=ntile, seed=seed, coarse_factor=rc["coarse_factor"], use_gpu=true,
                            devices=collect(0:length(CUDA.devices())-1), verbose=true)
    end
    raw = filter(h -> h.RTHL > 0, raw)
    zf = Float32(cfg.z_out)
    write_pksc(joinpath(outdir, "julia_$(MODE)_raw.pksc"), raw, maximum(h.RTHL for h in raw), zf)
    m = merge_catalog(raw; verbose=true, volume_reduction=true)
    write_pksc(joinpath(outdir, "julia_$(MODE)_merged.pksc"), m, maximum(h.RTHL for h in m), zf)
    @info "done" MODE raw = length(raw) merged = length(m)
end
main()
