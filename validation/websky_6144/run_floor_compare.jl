#!/usr/bin/env julia
# Julia side of the floor-completeness test. Runs the SAME config as the Fortran
# hpkvd floor run (next=256, box=320.83, prod cellsize 1.25326, websky, z=0
# snapshot, pk_websky /(2pi)^3, filters_websky) and prints N(>M) at the SAME mass
# thresholds as the Fortran catalog, for both the global-FFT (Fortran-equivalent)
# and multi-res-split paths. Direct apples-to-apples vs Fortran's:
#   N(>5e11)=220889 1e12=218002 1.69e12=212806 3e12=199604 1e13=128920
#   3e13=50252 1e14=11294 3e14=1969   (total 224181)
#
# Usage: julia --project=. -t 8 validation/websky_6144/run_floor_compare.jl
using TOML, PeakPatch, Printf

cfg_path = joinpath(@__DIR__, "config_ab_test.toml")
config = TOML.parsefile(cfg_path)
cfg = PipelineConfig(config)
rc = config["run"]
seed = get(rc, "seed", 13579); ntile = get(rc, "ntile", 4); cf = get(rc, "coarse_factor", 4)
rho_m = 2.775e11 * (cfg.Omx + cfg.OmB)

THRESH = [(5e11,"5e11"),(1e12,"1e12"),(1.69e12,"1.69e12"),(3e12,"3e12"),
          (1e13,"1e13"),(3e13,"3e13"),(1e14,"1e14"),(3e14,"3e14")]
FORT = Dict("5e11"=>220889,"1e12"=>218002,"1.69e12"=>212806,"3e12"=>199604,
            "1e13"=>128920,"3e13"=>50252,"1e14"=>11294,"3e14"=>1969)

function report(label, halos)
    M = [4/3*pi*Float64(h.RTHL)^3*rho_m for h in halos]
    sort!(M)
    @printf("\n==== %s : %d total halos (max RTHL=%.2f) ====\n",
            label, length(halos), isempty(halos) ? 0.0 : maximum(h.RTHL for h in halos))
    @printf("  %-10s  %10s  %10s  %7s\n","M>","Julia","Fortran","J/F")
    for (m0,lbl) in THRESH
        njl = count(>(m0), M)
        nf = FORT[lbl]
        @printf("  %-10s  %10d  %10d  %6.2f\n", lbl, njl, nf, njl/max(1,nf))
    end
end

@info "Running GLOBAL (run_multitile, Fortran-equivalent)..."
halos_global = run_multitile(cfg; ntile=ntile, seed=seed, verbose=false)
@info "Running SPLIT (run_multitile_split, multi-res)..."
halos_split = run_multitile_split(cfg; ntile=ntile, seed=seed, coarse_factor=cf,
                                  use_gpu=false, verbose=false)
report("GLOBAL (Fortran-equiv)", halos_global)
report("SPLIT (multi-res)", halos_split)
@printf("\nFortran total = 224181\n")
