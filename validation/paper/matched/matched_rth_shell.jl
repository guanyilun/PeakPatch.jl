#!/usr/bin/env julia
# Signature test for the OTHER-class R_TH disagreements (R_J < R_F, concentrated at small R_f):
# is the Fortran R_TH exactly a lattice shell radius √n·a? That is the footprint of collapse being accepted at the
# first tested shell m0 with no interpolation (RTHL = rad(m0−1)), which the uninitialised Frhoc at
# peakvoidsubs.f90:673 would produce; Julia fixes zvir1p = −1 there and steps inward, interpolating.
#   usage: julia --project=validation -t 16 matched_rth_shell.jl <rundir>
include(joinpath(@__DIR__, "matched_compare.jl"))
onshell(R) = begin q = (R / A)^2; abs(q - round(q)) < 2e-3 * q end          # R = √n·a to float32 precision
function shellsig()
    rd(f) = read_pksc(joinpath(RUN, f))[1]
    Fr = cat_lagr(rd("output/fortran_raw.pksc.12345")); Jr = cat_lagr(rd("julia/julia_exact_raw.pksc"))
    m = match(Fr, Jr; tol=0.5A); ok = findall(>(0), m)
    rat = [Jr.R[m[i]] / Fr.R[i] for i in ok]
    out = open(joinpath(@__DIR__, "..", "results", "matched_rth_shell.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    say("fraction of R_TH values lying exactly on a lattice shell radius √n·a")
    for (lab, sel) in (("J < 0.98 F", k -> rat[k] < 0.98), ("|J/F−1| ≤ 2%", k -> abs(rat[k] - 1) <= 0.02), ("J > 1.02 F", k -> rat[k] > 1.02))
        s = filter(sel, eachindex(ok)); isempty(s) && continue
        say(@sprintf("  %-14s n=%8d   Fortran on-shell %.4f   Julia on-shell %.4f", lab, length(s),
                     count(k -> onshell(Fr.R[ok[k]]), s) / length(s), count(k -> onshell(Jr.R[m[ok[k]]]), s) / length(s)))
    end
    # for J<F: is Julia's R inside the Fortran shell interval, i.e. same crossing shell, different acceptance?
    s = filter(k -> rat[k] < 0.98, eachindex(ok))
    gap = [Fr.R[ok[k]] / A - Jr.R[m[ok[k]]] / A for k in s]
    say(@sprintf("  J<F: (R_F − R_J) in cells: median %.3f  p10 %.3f  p90 %.3f", median(gap), quantile(gap, 0.1), quantile(gap, 0.9)))
    close(out)
end
shellsig()
