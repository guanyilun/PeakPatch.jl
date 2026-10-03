#!/usr/bin/env julia
# Per-peak harness, Julia side (MATCHED_FORTRAN_2026-10.md step 3): re-run analyse_peak for every Julia raw peak in
# the cube |x+107|,|y+107|,|z+107| < 12 Mpc/h (core of tile (2,2,2)) on the same field and tile arrays as
# run_multitile, dumping the shell-by-shell collapse search through the RadialShell._DBG hook. The instrumented
# Fortran hpkvd (fortran_dbg_patch.py) writes the same quantities for the same cube (run_dbg/dbg_homel_*.txt).
# REQUIRES the uncommitted debug hook: git apply validation/paper/matched/radialshell_dbg_hook.patch (revert after).
#   usage: julia --project=validation -t 8 matched_peak_dbg.jl <config.toml> <rundir>
using TOML, Printf, FFTW
using PeakPatch
import PeakPatch.RadialShell: analyse_peak, precompute_shells, PeakGrid, _DBG
import PeakPatch.CollapseTable: read_homeltab, CollapseTableInterp

function main()
    cfgd = TOML.parsefile(ARGS[1]); cfg = PipelineConfig(cfgd); run = ARGS[2]
    ntile = 4; nmesh = cfg.n; nbuff = cfg.nbuff; nsub, N = grid_layout(cfg, ntile)
    a = cfg.boxsize / nmesh; L = N * a; it = jt = kt = 2
    δ = Array{Float32}(undef, N, N, N); open(io -> read!(io, δ), joinpath(run, "fields", "Fvec_jmatched"))
    δk = rfft(δ)
    px, py, pz = displacements_1lpt(δk, N, L); qx, qy, qz = displacements_2lpt(δk, N, L)
    tl(f) = extract_tile(f, it, jt, kt, nsub, nmesh)
    pg = PeakGrid(tl(δ), tl(px), tl(py), tl(pz), tl(qx), tl(qy), tl(qz), zeros(Int8, nmesh, nmesh, nmesh), (nmesh, nmesh, nmesh), nothing)
    ct = CollapseTableInterp(read_homeltab(cfg.tabfile)...)
    filters = read_filterbank(cfg.filterfile); Rfmax = maximum(f[3] for f in filters)
    shells = precompute_shells(min(nbuff - 1, floor(Int, Rfmax * 1.75 / a)))
    raw = read_pksc(joinpath(run, "julia", "julia_exact_raw.pksc"))[1]
    sel = filter(h -> abs(h.x + 107) < 12 && abs(h.y + 107) < 12 && abs(h.z + 107) < 12, raw)
    @info "peaks in cube" length(sel)
    out = open(joinpath(run, "dbg_julia.txt"), "w"); nbad = 0
    for h in sel
        g = (Float64(h.x), Float64(h.y), Float64(h.z)) ./ a .+ (N + 1) / 2
        loc = round.(Int, g) .- ((it, jt, kt) .- 1) .* nsub
        ipp = loc[1] + (loc[2] - 1) * nmesh + (loc[3] - 1) * nmesh^2
        Rf = Float64(h.Rf); ir2min = min(floor(Int, (1.75Rf / a)^2), floor(Int, (40 / a - 1)^2))
        buf = IOBuffer(); _DBG[] = buf
        r = analyse_peak(pg, ipp, a, ir2min, 1.0 + cfg.z_out, Rf, ct, shells; nbuff=nbuff, rmax2rs=cfg.rmax2rs)
        _DBG[] = nothing
        R = r.RTHL * a
        abs(R - h.RTHL) > 1e-3 * h.RTHL && (nbad += 1)
        @printf(out, "PK %.6f %.6f %.6f %.5f %d\n", h.x, h.y, h.z, Rf, ir2min)
        print(out, String(take!(buf)))
        @printf(out, "R %.7f catalog %.7f\n", R, h.RTHL)
    end
    close(out)
    @info "done" n = length(sel) mismatched_vs_catalog = nbad
end
main()
