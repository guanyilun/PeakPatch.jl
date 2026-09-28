#!/usr/bin/env julia
# Full-sky abundance matching, following the reference Websky procedure
# (~/work/peakpatch/python/catalogue_tools/abundance_match/make_abundancematch_table.py; Stein+2020 §4.1).
# ONE table from the N(>M|z) counts of all 8 octants against the full-sky Tinker volume
# (fsky=1), instead of a per-octant table (fsky=1/8), and optionally a frozen fractional
# correction above `tail_N` halos per z-bin (AbundanceMatch.build_abundance_table docs;
# validation/paper/V3_RESULTS_2026-09-28.md for why).
#
# Usage:
#   julia ... apply_abundance_match_fullsky.jl table <table.txt> <tail_N> <z_max> <cat_oct000> ... <cat_oct111>
#       (catalog file names must contain "octZYX"; the observer comes from those bits)
#   julia ... apply_abundance_match_fullsky.jl apply <table.txt> <in.pksc> <out.pksc> <ZYX>
using PeakPatch, Printf

const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)   # Om total 0.31 (see apply_abundance_match.jl)
const RHO_M = 2.775e11 * 0.31
const NMBINS, NZBINS, MMIN, MMAX = 10000, 46, 5e11, 1e16                # = apply_abundance_match.jl defaults

# octZYX bit = 1 -> observer at +2618 Mpc/h on that axis
obs_of(zyx) = (zyx[3] == '1' ? 2618.0 : -2618.0, zyx[2] == '1' ? 2618.0 : -2618.0, zyx[1] == '1' ? 2618.0 : -2618.0)
octof(path) = (m = match(r"oct([01]{3})", basename(path)); m === nothing && error("no octZYX in $path"); m.captures[1])

NgtM(hs, M0) = count(h -> (4π / 3) * RHO_M * Float64(h.RTHL)^3 > M0, hs)
report(label, hs) = @printf("%-12s N(>1.7e12)=%.4e  N(>1e14)=%.4e  N(>1e15)=%d  N(>2e15)=%d\n", label,
                            NgtM(hs, 1.7e12), NgtM(hs, 1e14), NgtM(hs, 1e15), NgtM(hs, 2e15))

function make_table(out, tail_N, z_max, cats)
    length(cats) == 8 || error("need the 8 octant catalogs, got $(length(cats))")
    Medge, zedge = am_grid(; nMbins=NMBINS, z_min=0.0, z_max=z_max, nzbins=NZBINS, Mmin=MMIN, Mmax=MMAX)
    Npp = zeros(NMBINS, NZBINS)
    for c in cats
        oct = octof(c); @info "counting" c oct obs = obs_of(oct)
        halos, _, _ = read_pksc(c)
        am_counts!(Npp, halos, COSMO; obs=obs_of(oct), Medge, zedge)
        halos = nothing; GC.gc()
    end
    @printf("full-sky counts: %d halos in [%.1e, %.1e), z<%.2f\n", sum(Npp), MMIN, MMAX, z_max)
    pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "data", "pk_websky_RAW_unnormalized.dat"))
    t = build_abundance_table(Npp, COSMO, pk; z_min=0.0, z_max=z_max, Mmin=MMIN, Mmax=MMAX,
                              hmf=:tinker, fsky=1.0, tail_N=tail_N, verbose=true)
    save_abundance_table(out, t)
    @info "wrote table" out tail_N
end

function apply_table(tpath, inp, out, zyx)
    t = load_abundance_table(tpath)
    halos, _, z_out = read_pksc(inp)
    report("RAW", halos)
    am = abundance_match(halos, t, COSMO; obs=obs_of(zyx))
    report("AM fullsky", am)
    write_pksc(out, am, Float32(maximum(h.RTHL for h in am)), Float32(z_out))
    @info "wrote" out
end

mode = ARGS[1]
if mode == "table"
    make_table(ARGS[2], parse(Float64, ARGS[3]), parse(Float64, ARGS[4]), ARGS[5:end])
elseif mode == "apply"
    apply_table(ARGS[2], ARGS[3], ARGS[4], ARGS[5])
else
    error("mode must be table or apply")
end
