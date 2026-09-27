using TOML, PeakPatch, Printf
import PeakPatch.Merger: merge_catalog
NgtM(hs,M0,rm)=count(h->(4/3*pi*Float64(h.RTHL)^3*rm)>M0, hs)
function rep(l,hs,rm)
    @printf("%-30s N=%8d | >1.23e12=%8d >1.69e12=%8d >3e12=%8d >1e13=%7d >1e14=%6d\n",
        l, length(hs), NgtM(hs,1.23e12,rm),NgtM(hs,1.69e12,rm),NgtM(hs,3e12,rm),NgtM(hs,1e13,rm),NgtM(hs,1e14,rm))
end
for (tag,cf) in [("STANDARD 21filt rmincell1.65","config_ab_test.toml"),
                 ("DENSE 35filt rmincell1.2","config_ab_dense.toml")]
    cfg=PipelineConfig(TOML.parsefile(joinpath(@__DIR__,cf)))
    rm=2.775e11*(cfg.Omx+cfg.OmB)
    @info "running $tag..."
    h=run_multitile_split(cfg; ntile=4, seed=13579, coarse_factor=4, use_gpu=false)
    rep(tag*" RAW", h, rm)
    rep(tag*" MERGED", merge_catalog(h), rm)
end
@printf("\nApples-to-apples = MERGED rows (production exclusion/merge applied). RAW inflates dense via cross-scale dupes.\n")
@printf("Websky/deg² ratios (ours/websky): 0.41@1.23e12 ... 1.0@3e13. Filter bank confirmed if MERGED dense > MERGED standard at low mass.\n")
