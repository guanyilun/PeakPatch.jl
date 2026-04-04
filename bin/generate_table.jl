#!/usr/bin/env julia
#
# Generate ellipsoidal collapse table for PeakPatch.
#
# Usage:
#   julia --project=. bin/generate_table.jl [options]
#
# Options:
#   -o, --output PATH    Output file path (default: HomelTab.dat)
#   --solver SOLVER       ODE solver: "diffeq" (default, adaptive Tsit5) or "rk4"
#   --Om FLOAT            Total Omega_matter (default: 0.315)
#   --OB FLOAT            Omega_baryon (default: 0.049)
#   --OL FLOAT            Omega_Lambda (default: 0.685)
#   --h FLOAT             H0 / (100 km/s/Mpc) (default: 0.674)
#   -v, --verbose         Print progress info
#

using PeakPatch

function parse_args(args)
    output  = "HomelTab.dat"
    solver  = :diffeq
    Om      = 0.315
    OB      = 0.049
    OL      = 0.685
    h       = 0.674
    verbose = false

    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-o", "--output")
            i += 1; output = args[i]
        elseif a == "--solver"
            i += 1; solver = Symbol(args[i])
        elseif a == "--Om"
            i += 1; Om = parse(Float64, args[i])
        elseif a == "--OB"
            i += 1; OB = parse(Float64, args[i])
        elseif a == "--OL"
            i += 1; OL = parse(Float64, args[i])
        elseif a == "--h"
            i += 1; h = parse(Float64, args[i])
        elseif a in ("-v", "--verbose")
            verbose = true
        else
            println(stderr, "Unknown argument: $a")
            println(stderr, "Usage: julia bin/generate_table.jl [-o PATH] [--solver diffeq|rk4] [--Om F] [--OB F] [--OL F] [--h F] [-v]")
            exit(1)
        end
        i += 1
    end

    return (; output, solver, Om, OB, OL, h, verbose)
end

function main()
    opts = parse_args(ARGS)

    cosmo = CosmologyParams(opts.Om, opts.OB, opts.OL, opts.h, 0.965, 0.808)
    ep = EllipsoidParams(cosmo; solver=opts.solver)
    tp = CollapseTableParams()

    opts.verbose && @info "Generating collapse table" solver=opts.solver output=opts.output Om=opts.Om OB=opts.OB OL=opts.OL h=opts.h

    table = make_table_threaded(ep, tp; verbose=opts.verbose)
    write_homeltab(opts.output, table, tp)

    @info "Wrote collapse table: $(opts.output) ($(tp.Nx)×$(tp.Ny)×$(tp.Nz))"
end

main()
