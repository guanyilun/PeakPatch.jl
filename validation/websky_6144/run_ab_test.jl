#!/usr/bin/env julia
# A/B field-construction test (CPU, z=0):
#   A = run_multitile_split  (coarse+residual MUSIC, GPU production path, here on CPU)
#   B = run_multitile        (global FFT, Fortran-equivalent)
# Same cfg/seed/filters. Prints mass function + smoothed-peak-amplitude vs Rf
# for both, to see whether the multi-res field produces the fake-giant tail.
#
# Usage: julia --project=. -t 8 validation/websky_6144/run_ab_test.jl

using TOML, PeakPatch, Printf

cfg_path = isempty(ARGS) ? joinpath(@__DIR__, "config_ab_test.toml") : ARGS[1]
config = TOML.parsefile(cfg_path)
cfg = PipelineConfig(config)
rc = config["run"]
seed = get(rc, "seed", 13579); ntile = get(rc, "ntile", 4); cf = get(rc, "coarse_factor", 4)

rho_m = 2.775e11 * (cfg.Omx + cfg.OmB)   # Msun/h per (Mpc/h)^3

function report(label, halos)
    n = length(halos)
    @printf("\n==== %s : %d raw halos ====\n", label, n)
    isempty(halos) && return
    # mass = 4/3 pi RTHL^3 rho_m ; RTHL is physical Mpc/h
    massbin = Dict{Float64,Int}()
    rfbin   = Dict{Float64,Vector{Float32}}()
    rmax = 0.0
    for h in halos
        R = Float64(h.RTHL); rmax = max(rmax, R)
        M = 4/3*pi*R^3*rho_m
        M > 0 && (massbin[floor(log10(M)*4)/4] = get(massbin, floor(log10(M)*4)/4, 0) + 1)
        rf = round(Float64(h.Rf), digits=2)
        push!(get!(rfbin, rf, Float32[]), h.FcollvRf)
    end
    @printf("  max RTHL = %.1f Mpc/h\n", rmax)
    println("  log10(M)  N   (high-mass tail flagged with * if it RISES)")
    ks = sort(collect(keys(massbin)))
    prev = 0
    for k in ks
        k < 13.0 && continue
        flag = massbin[k] > prev ? " *" : ""
        @printf("    %5.2f  %8d%s\n", k, massbin[k], flag)
        prev = massbin[k]
    end
    println("  Rf [Mpc/h] | N | median smoothed-delta-at-peak (should FALL with Rf):")
    for rf in sort(collect(keys(rfbin)))
        v = sort(rfbin[rf])
        @printf("    Rf=%6.2f  %7d   median=%.3f\n", rf, length(v), v[length(v)÷2+1])
    end
end

@info "Running A: run_multitile_split (coarse+residual)..."
halos_split = run_multitile_split(cfg; ntile=ntile, seed=seed, coarse_factor=cf,
                                  use_gpu=false, verbose=false)
@info "Running B: run_multitile (global FFT, Fortran-equivalent)..."
halos_global = run_multitile(cfg; ntile=ntile, seed=seed, verbose=false)

report("A  split (coarse+residual)", halos_split)
report("B  global FFT (Fortran-equiv)", halos_global)

@printf("\nSUMMARY: split=%d vs global=%d raw halos (ratio %.2f)\n",
        length(halos_split), length(halos_global),
        length(halos_split)/max(1,length(halos_global)))
