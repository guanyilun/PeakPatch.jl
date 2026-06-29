using TOML, PeakPatch, Printf
import PeakPatch.Merger: merge_catalog
cfg=PipelineConfig(TOML.parsefile(joinpath(@__DIR__,"config_cellfix.toml")))
rho=2.775e11*(cfg.Omx+cfg.OmB)
NgtM(hs,M0)=count(h->(4/3*pi*Float64(h.RTHL)^3*rho)>M0,hs)
rep(l,hs)=@printf("%-24s N=%8d | >1.69e12=%8d >1e13=%7d >1e14=%6d\n",l,length(hs),NgtM(hs,1.69e12),NgtM(hs,1e13),NgtM(hs,1e14))
@info "Julia run_multitile N=468 cellsize 1.2533 (matches corrected Fortran)..."
h=run_multitile(cfg; ntile=2, seed=13579, verbose=false)
rep("JULIA raw (N=468)",h)
rep("JULIA + merge",merge_catalog(h))
@printf("\nCompare to corrected Fortran (same n=468, cellsize 1.2533, box 586.5)\n")
