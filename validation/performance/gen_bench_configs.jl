#!/usr/bin/env julia
# Generate benchmark configs (validation/performance/configs/) at the production cell size
# (0.85221 Mpc/h = Websky), production physics (2LPT, ioutshear=1, finecell filters, ievol=1),
# observer at a grid corner, and z_max chosen so the lightcone horizon covers the whole box
# (every tile is active → work is fixed and load-balanced).
#
#   julia --project=validation validation/performance/gen_bench_configs.jl
#
# Sets:
#   scale   n=414 (= production tile), ntile=4, N=1560, cf=30 (coarse cell 11.1 Mpc/h ~ prod 10.2)
#           → strong scaling 1/2/4 L40S + 1 H100
#   warmup  n=96, ntile=2 — tiny run executed first in each process to take JIT out of timings
#   mem_nXXX  ntile=2, tile size n ∈ {256,320,384,448,512,576} → GPU memory vs tile size
import PeakPatch.Cosmology: CosmologyParams, chi
using Printf

const ALATT = 5236.0 / 6144      # Mpc/h, production cell size
const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.808)
const DATA  = "validation/websky_6144/data"
outdir = joinpath(@__DIR__, "configs"); mkpath(outdir)

function z_of_chi(target)
    lo, hi = 0.0, 20.0
    for _ in 1:100
        mid = (lo + hi) / 2
        chi(mid, COSMO) < target ? (lo = mid) : (hi = mid)
    end
    return hi
end

function write_config(name; n, ntile, nbuff=16, cf)
    nsub = n - 2nbuff; N = nsub * ntile + 2nbuff; M = ntile * cf
    N % M == 0 || error("$name: N=$N not divisible by M=$M")
    L = N * ALATT; cen = -L / 2
    zmax = round(z_of_chi(1.01 * sqrt(3) * L); digits=3)     # horizon beyond the far corner
    open(joinpath(outdir, "$name.toml"), "w") do io
        println(io, "# Benchmark config '$name' (gen_bench_configs.jl): n=$n ntile=$ntile nbuff=$nbuff N=$N ",
                "M=$M block=$(N ÷ M); full box $(round(L; digits=1)) Mpc/h; z_max=$zmax covers the box")
        println(io, """
        [cosmology]
        Om = 0.31
        OB = 0.049
        OL = 0.69
        h  = 0.68

        [grid]
        n       = $n
        boxsize = $(round(n * ALATT; digits=4))
        nbuff   = $nbuff
        cenx    = $(round(cen; digits=3))
        ceny    = $(round(cen; digits=3))
        cenz    = $(round(cen; digits=3))

        [run]
        z_out     = 0.0
        z_max     = $zmax
        ievol     = 1
        ilpt      = 2
        ioutshear = 1
        wsmooth   = 1
        rmax2rs   = 0.0
        seed      = 12345
        ntile     = $ntile
        coarse_factor = $cf

        [files]
        pk         = "$DATA/pk_websky.dat"
        filterbank = "$DATA/filters_websky_finecell.dat"
        homeltab   = "$DATA/HomelTab_websky.dat"
        output     = "bench_$name.pksc"
        """)
    end
    @printf("%-10s n=%3d ntile=%d N=%4d M=%3d block=%2d  L=%7.1f Mpc/h  z_max=%.3f\n",
            name, n, ntile, N, M, N ÷ M, L, zmax)
end

write_config("scale";  n=414, ntile=4, cf=30)
write_config("warmup"; n=96,  ntile=2, cf=10)
for (n, cf) in ((256, 24), (320, 16), (384, 16), (448, 24), (512, 16), (576, 28))
    write_config(@sprintf("mem_n%03d", n); n=n, ntile=2, cf=cf)
end
