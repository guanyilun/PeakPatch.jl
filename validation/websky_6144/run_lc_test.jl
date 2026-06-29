#!/usr/bin/env julia
# Test the lightcone machinery (ievol=1) in a small box where z spans only 0-0.19
# (observer at corner). At z~0 the lightcone count MUST match the ievol=0 snapshot.
# If ievol=1 gives far fewer, the lightcone path under-produces (production bug).
using TOML, PeakPatch, Printf

rho(cfg) = 2.775e11*(cfg.Omx+cfg.OmB)
NgtM(hs,M0,rm) = count(h -> (4/3*pi*Float64(h.RTHL)^3*rm) > M0, hs)
function rep(l,hs,rm)
    @printf("%-22s N=%8d | >1.69e12=%8d >1e13=%7d >1e14=%6d  maxRTHL=%.2f\n", l,
        length(hs), NgtM(hs,1.69e12,rm), NgtM(hs,1e13,rm), NgtM(hs,1e14,rm),
        isempty(hs) ? 0.0 : maximum(h.RTHL for h in hs))
end

cfg0 = PipelineConfig(TOML.parsefile(joinpath(@__DIR__,"config_ab_test.toml")))  # ievol=0
cfg1 = PipelineConfig(TOML.parsefile(joinpath(@__DIR__,"config_lc_test.toml")))  # ievol=1, obs corner
rm = rho(cfg0)

@info "ievol=0 snapshot (split)..."
h0 = run_multitile_split(cfg0; ntile=4, seed=13579, coarse_factor=4, use_gpu=false)
@info "ievol=1 lightcone, observer corner, z<0.19 (split)..."
h1 = run_multitile_split(cfg1; ntile=4, seed=13579, coarse_factor=4, use_gpu=false)
rep("ievol=0 snapshot", h0, rm)
rep("ievol=1 lightcone", h1, rm)
@printf("\nAt z~0 these should MATCH. ratio lc/snap = %.3f\n", length(h1)/max(1,length(h0)))
