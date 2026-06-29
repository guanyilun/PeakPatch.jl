#!/usr/bin/env julia
# Capstone: confirm Julia finder == Fortran by applying merge_catalog to a Julia
# global run and overlaying its full N(>M) on Fortran-merged
# (total 71864; >1e13=31075; >1e14=2211).
using TOML, PeakPatch, Printf
import PeakPatch.Merger: merge_catalog

cfg = PipelineConfig(TOML.parsefile(joinpath(@__DIR__, "config_ab_test.toml")))
rho_m = 2.775e11 * (cfg.Omx + cfg.OmB)
NgtM(hs,M0) = count(h -> (4/3*pi*Float64(h.RTHL)^3*rho_m) > M0, hs)
rep(l,hs) = @printf("%-26s N=%8d | >1.69e12=%8d >1e13=%7d >1e14=%6d\n", l,
    length(hs), NgtM(hs,1.69e12), NgtM(hs,1e13), NgtM(hs,1e14))

@info "Running global..."
halos = run_multitile(cfg; ntile=4, seed=13579, verbose=false)
rep("JULIA raw (masked-find)", halos)
merged = merge_catalog(halos; verbose=true)
rep("JULIA + merge_catalog", merged)
@printf("\nFORTRAN+merge reference:   N=   71864 | >1.69e12=   66206 >1e13=  31075 >1e14=  2211\n")
