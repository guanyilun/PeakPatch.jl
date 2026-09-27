#!/usr/bin/env julia
# Fast (CPU, small-box) confirmation of the factor-h cellsize MECHANISM, in parallel with the
# full GPU octant. Same N=256 tiling, z=0 snapshot, two cellsizes:
#   COARSE: cellsize 1.2533 Mpc/h, filters_websky.dat (Rf_min 2.068 -> M_min 3.18e12) = our old run
#   FINE:   cellsize 0.8522 Mpc/h, filters_websky_finecell.dat (Rf_min 1.406 -> M_min 1.00e12) = Websky
# CLAIM: the fine run's raw mass function extends down to ~1e12; the coarse one turns over ~3e12.
using TOML, PeakPatch, Printf

base = TOML.parsefile(joinpath(@__DIR__, "config_ab_test.toml"))
ntile = 4; seed = 13579

function run_one(boxsize_pertile, filterbank)
    cfg_d = deepcopy(base)
    cfg_d["grid"]["boxsize"] = boxsize_pertile
    cfg_d["files"]["filterbank"] = filterbank
    cfg = PipelineConfig(cfg_d)
    halos = run_multitile_split(cfg; ntile=ntile, seed=seed, coarse_factor=4, use_gpu=false)
    nbuff = cfg.nbuff; nmesh = cfg.n; nsub = nmesh - 2*nbuff
    N = nsub*ntile + 2*nbuff
    alatt = cfg.boxsize / cfg.n
    boxfull = N * alatt
    rho = 2.775e11*(cfg.Omx+cfg.OmB)
    M = [ (4/3*pi*rho*Float64(h.RTHL)^3) for h in halos if h.RTHL > 0 ]
    return M, boxfull, alatt
end

# differential MF dn/dlnM per (Mpc/h)^3
edges = 10.0 .^ (11.6:0.15:14.0)
ctr = [sqrt(edges[i]*edges[i+1]) for i in 1:length(edges)-1]
function mf(M, boxfull)
    V = boxfull^3; h = zeros(length(ctr))
    for m in M; i=searchsortedfirst(edges,m)-1; (1<=i<=length(ctr)) && (h[i]+=1.0); end
    [h[i]/V/(log(edges[i+1])-log(edges[i])) for i in 1:length(ctr)]
end

@info "COARSE run (cell 1.2533, filters_websky)..."
Mc, boxc, ac = run_one(125.326, "validation/websky_6144/data/filters_websky.dat")
@info "FINE run (cell 0.8522, filters_websky_finecell)..."
Mf, boxf, af = run_one(85.221, "validation/websky_6144/data/filters_websky_finecell.dat")

dc = mf(Mc, boxc); df = mf(Mf, boxf)
rho = 2.775e11*0.359
Mmin_c = 4/3*pi*rho*(1.65*ac)^3; Mmin_f = 4/3*pi*rho*(1.65*af)^3
@printf("\nCOARSE: cell=%.4f Mpc/h box=%.1f  N(halos)=%d  smallest-filter M=%.2e\n", ac, boxc, length(Mc), Mmin_c)
@printf("FINE:   cell=%.4f Mpc/h box=%.1f  N(halos)=%d  smallest-filter M=%.2e\n", af, boxf, length(Mf), Mmin_f)
@printf("\n%-10s %-14s %-14s   (dn/dlnM per (Mpc/h)^3)\n","M_center","COARSE","FINE")
for i in 1:length(ctr)
    mk = ""
    (edges[i]<=Mmin_c<edges[i+1]) && (mk*=" <coarse Rf_min")
    (edges[i]<=Mmin_f<edges[i+1]) && (mk*=" <fine Rf_min")
    @printf("%-10.2e %-14.3e %-14.3e%s\n", ctr[i], dc[i], df[i], mk)
end
@printf("\nMECHANISM CONFIRMED if FINE has nonzero dn/dlnM down to ~1e12 while COARSE rolls off ~3e12.\n")
