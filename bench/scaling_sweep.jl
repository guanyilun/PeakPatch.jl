#!/usr/bin/env julia
# Scaling sweep: run GPU `run_multitile_split` at several (nmesh, ntile) points
# to measure per-tile cost and extrapolate to production volumes (e.g. N=6144).
#
# Run:
#   CUDA_VISIBLE_DEVICES=0 julia --project=/tmp/pp_gpu_sandbox --threads=32 \
#       bench/scaling_sweep.jl
#
# Knobs (env):
#   NMESH_LIST (default "128,192,256,320,384")
#   NTILE (default 2)
#   ILPT  (default 2)
#   IOUTSHEAR (default 0)
#   NBUFF (default 8)
#   TARGET_N (default 6144; used for extrapolation)

using PeakPatch
using CUDA
using Printf

const NMESH_LIST = [parse(Int, s) for s in split(get(ENV, "NMESH_LIST", "128,192,256,320,384"), ",")]
const NTILE      = parse(Int, get(ENV, "NTILE", "2"))
const ILPT       = parse(Int, get(ENV, "ILPT",  "2"))
const IOUTSHEAR  = parse(Int, get(ENV, "IOUTSHEAR", "0"))
const NBUFF      = parse(Int, get(ENV, "NBUFF", "8"))
const TARGET_N   = parse(Int, get(ENV, "TARGET_N", "6144"))
const ALATT_REF  = 226.67 / 68.0  # Mpc/h per cell (matches other benches)

function make_config(dir; n, boxsize, z=0.0, ilpt=2, nbuff=8, ioutshear=0)
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
        ievol=0, z_max=0.0, cenx=0.0, ceny=0.0, cenz=0.0,
        Omx=0.261, OmB=0.049, Omvac=0.69, h=0.68,
        NonGauss=0, fNL=0.0, wsmooth=1,
    )
end

@inline _mem_used_mb() = (CUDA.total_memory() - CUDA.free_memory()) / (1024^2)
@inline _mem_total_mb() = CUDA.total_memory() / (1024^2)

println("\n" * "="^72)
println("GPU scaling sweep — run_multitile_split")
println("="^72)
if CUDA.functional()
    dev = CUDA.device()
    @printf "GPU: %s  (%.1f GB)   threads=%d\n" CUDA.name(dev) (CUDA.totalmem(dev)/1024^3) Threads.nthreads()
else
    error("CUDA not functional — this benchmark requires a GPU.")
end
@printf "ntile=%d  ilpt=%d  ioutshear=%d  nbuff=%d   target N=%d\n" NTILE ILPT IOUTSHEAR NBUFF TARGET_N
@printf "nmesh sweep: %s\n\n" join(NMESH_LIST, ", ")

# Warmup: smallest size only, both CPU (for comparison at small N) and GPU.
nmesh_warm = minimum(NMESH_LIST)
nsub_warm  = nmesh_warm - 2 * NBUFF
N_warm     = nsub_warm * NTILE + 2 * NBUFF
box_warm   = Float64(N_warm) * ALATT_REF
tmpdir_warm = mktempdir()
cfg_warm = make_config(tmpdir_warm; n=nmesh_warm, boxsize=box_warm, ilpt=ILPT, nbuff=NBUFF, ioutshear=IOUTSHEAR)
println(">> GPU warmup (nmesh=$nmesh_warm)..."); flush(stdout)
_ = PeakPatch.run_multitile_split(cfg_warm; ntile=NTILE, seed=1, verbose=false, use_gpu=true)
rm(tmpdir_warm; recursive=true)

results = NamedTuple[]

for nmesh in NMESH_LIST
    nsub = nmesh - 2 * NBUFF
    N    = nsub * NTILE + 2 * NBUFF
    boxsize = Float64(N) * ALATT_REF
    @printf "\n---- nmesh=%d  nsub=%d  N=%d  ntiles=%d  box=%.1f Mpc/h ----\n" nmesh nsub N NTILE^3 boxsize
    flush(stdout)

    # Pre-bench: drop any cached state, measure baseline memory
    GC.gc(); CUDA.reclaim()
    baseline_used = _mem_used_mb()

    tmpdir = mktempdir()
    cfg = make_config(tmpdir; n=nmesh, boxsize=boxsize, ilpt=ILPT, nbuff=NBUFF, ioutshear=IOUTSHEAR)

    # Poll GPU memory in the background to capture peak usage mid-run.
    peak_tracker = Ref(baseline_used)
    stop_poll = Ref(false)
    poller = Threads.@spawn begin
        while !stop_poll[]
            u = _mem_used_mb()
            if u > peak_tracker[]; peak_tracker[] = u; end
            sleep(0.25)
        end
    end

    # Run with profile=true to get per-stage breakdown
    t_gpu = @elapsed h_gpu = PeakPatch.run_multitile_split(cfg; ntile=NTILE, seed=42,
                                                          verbose=false, use_gpu=true,
                                                          profile=true)
    stop_poll[] = true; wait(poller)
    peak_used = peak_tracker[]

    n_tiles = NTILE^3
    per_tile = t_gpu / n_tiles

    @printf "   wall=%.2fs  per_tile=%.3fs  halos=%d  mem_used~%.0f MB (baseline %.0f MB)\n" t_gpu per_tile length(h_gpu) peak_used baseline_used
    flush(stdout)

    push!(results, (nmesh=nmesh, nsub=nsub, N=N, boxsize=boxsize, ntiles=n_tiles,
                    t_total=t_gpu, t_per_tile=per_tile, nhalos=length(h_gpu),
                    mem_mb=peak_used))

    rm(tmpdir; recursive=true)
    GC.gc(); CUDA.reclaim()
end

println("\n" * "="^72)
println("Summary table")
println("="^72)
@printf "%-6s %-6s %-6s %-10s %-10s %-9s %-8s %-10s\n" "nmesh" "nsub" "N" "wall(s)" "per-tile" "halos" "mem(MB)" "cell/μs/GPU"
for r in results
    cells_per_us = (r.N^3) / (r.t_total * 1e6)
    @printf "%-6d %-6d %-6d %-10.2f %-10.3f %-9d %-8.0f %-10.2f\n" r.nmesh r.nsub r.N r.t_total r.t_per_tile r.nhalos r.mem_mb cells_per_us
end

println("\n" * "="^72)
println("Extrapolation to N=$TARGET_N")
println("="^72)
# Choose the largest nmesh that fit (best amortised per-tile cost).
best = results[end]
@printf "Reference point: nmesh=%d, per-tile=%.3fs, mem=%.0f MB\n" best.nmesh best.t_per_tile best.mem_mb

println("\nProjected full-run cost on a SINGLE GPU using this nmesh:")
@printf "%-10s %-10s %-12s %-14s %-14s\n" "nmesh_t" "ntile_t" "tiles_tot" "wall_1GPU(h)" "wall_8GPU(h)"
for nmesh_t in (best.nmesh, 256, 192)
    nsub_t = nmesh_t - 2 * NBUFF
    nsub_t <= 0 && continue
    ntile_t = cld(TARGET_N - 2 * NBUFF, nsub_t)
    # Adjust so N matches exactly if possible (approximate).
    tiles_tot = ntile_t^3
    # Per-tile cost for this nmesh_t (take closest measurement, else scale from `best` by nmesh³).
    idx = findfirst(r -> r.nmesh == nmesh_t, results)
    per_tile = idx === nothing ? best.t_per_tile * (nmesh_t / best.nmesh)^3 : results[idx].t_per_tile
    wall_1 = per_tile * tiles_tot / 3600
    wall_8 = wall_1 / 8
    @printf "%-10d %-10d %-12d %-14.2f %-14.2f\n" nmesh_t ntile_t tiles_tot wall_1 wall_8
end

println("""

Notes:
 - Wall-times assume tile-level work dominates. Coarse-grid setup (phase 1a) is
   shared across all tiles; at target scales it is negligible relative to tile loop.
 - Multi-GPU numbers assume embarrassingly parallel tile dispatch (not yet
   implemented in the pipeline — requires minor refactor to distribute the
   outer `for tid in tile_ids` loop across CUDA devices / MPI ranks).
 - Memory figure is post-run CUDA usage; peak during execution is somewhat
   higher (cuFFT plans + intermediate buffers). Keep ~2× headroom.
""")
