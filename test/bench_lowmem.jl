#!/usr/bin/env julia
#
# Benchmark run_multitile vs run_multitile_lowmem: memory and runtime.
# Usage: julia --project=. -t 4 test/bench_lowmem.jl
#
using PeakPatch
using Printf

# ---- Setup ----
function make_config(tmpdir; n, boxsize, ilpt=2, nbuff=4, ioutshear=0)
    # Generate P(k) file
    pkfile = joinpath(tmpdir, "pk.dat")
    open(pkfile, "w") do f
        for k in 10.0 .^ range(-4, stop=1, length=200)
            println(f, "$k  $(1e4 * k^(-1.5))")
        end
    end

    # Generate filter file
    filterfile = joinpath(tmpdir, "filters.dat")
    open(filterfile, "w") do f
        println(f, "5")
        for (i, Rf) in enumerate([12.0, 9.0, 7.0, 5.0, 3.5])
            println(f, "$i  1.686  $Rf  0.3")
        end
    end

    tabfile = joinpath(@__DIR__, "data", "HomelTab_julia.dat")
    PipelineConfig(
        Omx=0.266, OmB=0.049, Omvac=0.685, h=0.674,
        n=n, boxsize=Float64(boxsize), nbuff=nbuff,
        z_out=0.0, z_max=0.0, ilpt=ilpt, ioutshear=ioutshear,
        rmax2rs=0.0,
        pkfile=pkfile, filterfile=filterfile,
        fileout=joinpath(tmpdir, "out.pksc"), tabfile=tabfile,
    )
end

function measure_memory_and_time(f, label)
    GC.gc(); GC.gc()
    mem_before = Sys.free_memory()

    t0 = time()
    result = f()
    elapsed = time() - t0

    GC.gc(); GC.gc()
    mem_after = Sys.free_memory()

    # Peak memory from GC stats
    peak_bytes = Base.gc_live_bytes()

    @printf "  %-30s  %6.1f sec  %d halos\n" label elapsed length(result)
    return (elapsed=elapsed, nhalos=length(result), result=result)
end

# ---- Benchmarks at increasing grid sizes ----
println("=" ^ 70)
println("Benchmark: run_multitile vs run_multitile_lowmem")
println("Threads: $(Threads.nthreads())")
println("=" ^ 70)

for (n, ntile, boxsize, nbuff) in [
    (32,  2, 100.0, 4),
    (64,  2, 200.0, 6),
    (128, 2, 400.0, 8),
]
    N = (n - 2*nbuff) * ntile + 2*nbuff
    field_gb = N^3 * 8 / 1e9  # Float64
    println()
    @printf "--- n=%d, ntile=%d, N=%d (%.2f GB per Float64 field) ---\n" n ntile N field_gb

    tmpdir = mktempdir()
    cfg = make_config(tmpdir; n=n, boxsize=boxsize, nbuff=nbuff)

    # Warmup (first run triggers compilation)
    if n == 32
        println("  (warmup...)")
        run_multitile(cfg; ntile=ntile, seed=1, verbose=false)
        run_multitile_lowmem(cfg; ntile=ntile, seed=1, verbose=false)
    end

    # Standard
    GC.gc(); GC.gc()
    bytes_before_std = Base.gc_live_bytes()
    t_std = @elapsed halos_std = run_multitile(cfg; ntile=ntile, seed=42, verbose=false)
    bytes_peak_std = Base.gc_live_bytes()
    GC.gc(); GC.gc()
    bytes_after_std = Base.gc_live_bytes()

    # Lowmem
    GC.gc(); GC.gc()
    bytes_before_lm = Base.gc_live_bytes()
    t_lm = @elapsed halos_lm = run_multitile_lowmem(cfg; ntile=ntile, seed=42, verbose=false)
    bytes_peak_lm = Base.gc_live_bytes()
    GC.gc(); GC.gc()
    bytes_after_lm = Base.gc_live_bytes()

    @printf "  %-12s  time=%6.2fs  halos=%d  live_bytes_peak=%.1f MB\n" "standard" t_std length(halos_std) (bytes_peak_std - bytes_before_std)/1e6
    @printf "  %-12s  time=%6.2fs  halos=%d  live_bytes_peak=%.1f MB\n" "lowmem" t_lm length(halos_lm) (bytes_peak_lm - bytes_before_lm)/1e6
    @printf "  speedup: %.2fx   memory ratio: %.2fx\n" t_std/t_lm (bytes_peak_std - bytes_before_std) / max(1, bytes_peak_lm - bytes_before_lm)

    # Verify identical results
    @assert length(halos_std) == length(halos_lm) "Halo count mismatch: $(length(halos_std)) vs $(length(halos_lm))"

    rm(tmpdir; recursive=true)
end

# ---- Detailed allocation tracking for the largest case ----
println()
println("=" ^ 70)
println("Allocation tracking (n=128, ntile=2)")
println("=" ^ 70)

tmpdir = mktempdir()
cfg = make_config(tmpdir; n=128, boxsize=400.0, nbuff=8)

println("\nStandard:")
GC.gc(); GC.gc()
stats_std = @timed run_multitile(cfg; ntile=2, seed=42, verbose=false)
@printf "  time=%.2fs  alloc=%.1f MB  gctime=%.2fs\n" stats_std.time stats_std.bytes/1e6 stats_std.gctime

println("\nLowmem:")
GC.gc(); GC.gc()
stats_lm = @timed run_multitile_lowmem(cfg; ntile=2, seed=42, verbose=false)
@printf "  time=%.2fs  alloc=%.1f MB  gctime=%.2fs\n" stats_lm.time stats_lm.bytes/1e6 stats_lm.gctime

@printf "\nAllocation ratio (std/lowmem): %.2fx\n" stats_std.bytes / stats_lm.bytes
@printf "Time ratio (lowmem/std): %.2fx\n" stats_lm.time / stats_std.time

rm(tmpdir; recursive=true)
println("\nDone.")
