#!/usr/bin/env julia
# Isolate the production-specific GPU shell-analysis path. Run run_multitile_split
# with use_gpu=true vs use_gpu=false on the IDENTICAL config (z=0 snapshot), apply
# the SAME merge_catalog to both, and compare floor densities. Same field, same M,
# same merge — the ONLY difference is GPU vs CPU shell analysis. If GPU < CPU at the
# floor, the GPU kernel rejects floor halos differently (the production-deficit suspect).
using CUDA
using TOML, PeakPatch, Printf
import PeakPatch.Merger: merge_catalog

cfg = PipelineConfig(TOML.parsefile(joinpath(@__DIR__, "config_ab_test.toml")))
rho = 2.775e11*(cfg.Omx+cfg.OmB)
NgtM(hs,M0)=count(h->(4/3*pi*Float64(h.RTHL)^3*rho)>M0, hs)
function rep(l,hs)
    @printf("%-24s N=%8d | >1.69e12=%8d >1e13=%7d >1e14=%6d  maxRTHL=%.2f\n",
        l, length(hs), NgtM(hs,1.69e12), NgtM(hs,1e13), NgtM(hs,1e14),
        isempty(hs) ? 0.0 : maximum(h.RTHL for h in hs))
end

@info "CUDA" functional=CUDA.functional() dev=(CUDA.functional() ? name(CUDA.device()) : "none")

@info "CPU split..."
hc = run_multitile_split(cfg; ntile=4, seed=13579, coarse_factor=4, use_gpu=false)
@info "GPU split..."
hg = run_multitile_split(cfg; ntile=4, seed=13579, coarse_factor=4, use_gpu=true, devices=[0])

rep("CPU raw", hc); rep("GPU raw", hg)
rep("CPU + merge", merge_catalog(hc)); rep("GPU + merge", merge_catalog(hg))
@printf("\nGPU/CPU raw floor ratio = %.4f ; if <1 the GPU shell path under-finds at the floor\n",
        NgtM(hg,1.69e12)/max(1,NgtM(hc,1.69e12)))
