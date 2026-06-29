#!/usr/bin/env julia
# Test the find-masking hypothesis. PP_NO_FINDMASK disables Julia's cross-scale
# peak masking. If nomask+merge jumps to ~49643 (= Fortran-peaks+merge) while
# masked+merge=25875, the premature find-masking is confirmed as the under-finding bug.
ENV["PP_NO_FINDMASK"] = "1"
using TOML, PeakPatch, Printf
import PeakPatch.Merger: merge_catalog

cfg = PipelineConfig(TOML.parsefile(joinpath(@__DIR__,"config_ab_test.toml")))
rho_m = 2.775e11*(cfg.Omx+cfg.OmB)
NgtM(hs,M0) = count(h -> (4/3*pi*Float64(h.RTHL)^3*rho_m) > M0, hs)
rep(l,hs) = @printf("%-26s N=%8d | >1.69e12=%8d >1e13=%7d >1e14=%6d\n", l,
    length(hs), NgtM(hs,1.69e12), NgtM(hs,1e13), NgtM(hs,1e14))

@info "Running global with PP_NO_FINDMASK..."
halos = run_multitile(cfg; ntile=4, seed=13579, verbose=false)
rep("NOMASK raw", halos)
merged = merge_catalog(halos; verbose=true)
rep("NOMASK + merge_catalog", merged)
@printf("\nReference: masked+merge=25875 ; Fortran-peaks+merge=49643 (>1e13=15561,>1e14=1341)\n")
