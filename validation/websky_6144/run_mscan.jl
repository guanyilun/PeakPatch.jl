#!/usr/bin/env julia
# Does the multi-res coarse grid M affect the floor count? Production uses M=64; my
# verification used M=16. Run run_multitile_split at M=16/32/64 (coarse_factor=4/8/16,
# ntile=4, N=256) + run_multitile (global, exact) on the SAME z=0 field, and compare
# floor densities. If the floor count drops as M grows, the multi-res reconstruction
# loses small-scale (floor) power at high M — a production-deficit suspect.
using TOML, PeakPatch, Printf
import PeakPatch.Merger: merge_catalog

cfg = PipelineConfig(TOML.parsefile(joinpath(@__DIR__, "config_ab_test.toml")))
rho = 2.775e11*(cfg.Omx+cfg.OmB)
NgtM(hs,M0)=count(h->(4/3*pi*Float64(h.RTHL)^3*rho)>M0, hs)
rep(l,hs)=@printf("%-22s N=%8d | >1.69e12=%8d >1e13=%7d >1e14=%6d\n", l,
    length(hs), NgtM(hs,1.69e12), NgtM(hs,1e13), NgtM(hs,1e14))

@info "global (exact)..."
rep("GLOBAL", run_multitile(cfg; ntile=4, seed=13579))
for cf in (4, 8, 16)
    M = 4*cf
    @info "split M=$M (block=$(256÷M))..."
    rep("SPLIT M=$M", run_multitile_split(cfg; ntile=4, seed=13579, coarse_factor=cf, use_gpu=false))
end
@printf("\nProduction M=64 (block=96). If floor stable across M here, multi-res is not the loss.\n")
