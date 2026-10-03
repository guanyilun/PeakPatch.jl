#!/usr/bin/env julia
# Decompose the remaining same-field offset F/JEcap (matched_rth_diag.txt) with hybrids of the Julia raw peaks:
#   JEcap   : Fortran R_TH on the CAP class (R_J beyond the Fortran search limit)
#   JEother : Fortran R_TH on the OTHER class only
#   JEall   : Fortran R_TH on every matched peak (Julia-only peaks kept, Fortran-only peaks absent)
#   JEallU  : JEall + the Fortran-only peaks appended (≈ the Fortran raw catalogue)
# Each is merged (exclusion + volume reduction), rank-matched and measured exactly as in matched_compare.jl.
#   usage: julia --project=validation -t 32 matched_rth_hybrid.jl <rundir>
include(joinpath(@__DIR__, "matched_compare.jl"))
isrep(n) = begin m = n; while m > 0 && m % 4 == 0; m ÷= 4; end; m % 8 != 7 end
nextrep(n) = begin k = n + 1; while !isrep(k); k += 1; end; k end
const NHUNT = min(NB - 1, floor(Int, 30.4358 * 1.75 / A))
rad_m0(Rf) = begin ir2min = min(floor(Int, (1.75Rf / A)^2), floor(Int, (40 / A - 1)^2)); r2 = nextrep(ir2min); r2 > NHUNT^2 ? Inf : sqrt(r2) * A end
withR(h, R) = typeof(h)(ntuple(f -> f == 7 ? Float32(R) : getfield(h, f), fieldcount(typeof(h)))...)

function hybrids()
    rd(f) = read_pksc(joinpath(RUN, f))[1]
    fr = filter(q -> q.RTHL > 0, rd("output/fortran_raw.pksc.12345")); jr = filter(q -> q.RTHL > 0, rd("julia/julia_exact_raw.pksc"))
    Fr = cat_lagr(fr); Jr = cat_lagr(jr); m = match(Fr, Jr; tol=0.5A); ok = findall(>(0), m)
    cap = Dict(m[i] => (Jr.R[m[i]] > rad_m0(Fr.Rf[i]) + 1e-3) for i in ok)
    RF = Dict(m[i] => Fr.R[i] for i in ok)
    mk(sel) = [haskey(RF, j) && sel(j) ? withR(jr[j], RF[j]) : jr[j] for j in eachindex(jr)]
    raws = Dict("JEcap" => mk(j -> cap[j]), "JEother" => mk(j -> !cap[j]), "JEall" => mk(j -> true))
    raws["JEallU"] = vcat(raws["JEall"], fr[m.==0])
    cats = Dict("F" => cat_mpkvd(rd("output/fortran_merge.pksc.12345")), "JE" => cat_lagr(rd("julia/julia_exact_merged.pksc")))
    for (k, v) in raws; cats[k] = cat_lagr(merge_catalog(v; volume_reduction=true)); end
    labs = ["F", "JEallU", "JEall", "JEother", "JEcap", "JE"]
    AM = Dict(l => am_masses(cats[l]) for l in labs)
    dlin = linear_field() .* DZ; Dk = rfft(dlin)
    out = open(joinpath(@__DIR__, "..", "results", "matched_rth_hybrid.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    say("hybrids of the Julia raw peaks with Fortran R_TH (same field); ratios are F / X (jackknife σ)")
    X = Dict(l => (c = AM[l][1]; s = findall(>(5e12), AM[l][2]); xi(c.xE[s], c.yE[s], c.zE[s])) for l in labs)
    say("  ξ(3–15), M>5e12:   " * join([@sprintf("F/%s %.4f  ", l, mean(X["F"][WX]) / mean(X[l][WX])) for l in labs[2:end]]))
    MB = [5e12, 8e12, 1.3e13, 2e13, 3.2e13, 7.9e13, 2.5e14]
    say("  b_E (k<0.1)        " * join([@sprintf("%-17s", "F/" * l) for l in labs[2:end]]))
    for i in 0:length(MB)-1
        B = Dict{String,Tuple{Float64,Vector{Float64}}}()
        for l in labs
            c, M = AM[l]; s = i == 0 ? findall(>(5e12), M) : findall(q -> MB[i] <= q < MB[i+1], M)
            h = deposit(c.xE[s], c.yE[s], c.zE[s]); B[l] = (bx(rfft(h), Dk), bx_jk(h, dlin))
        end
        say(@sprintf("  %-18s", i == 0 ? "M>5e12" : @sprintf("%.1e–%.1e", MB[i], MB[i+1])) *
            join([@sprintf("%.4f±%.4f  ", B["F"][1] / B[l][1], jkerr(B["F"][2] ./ B[l][2])) for l in labs[2:end]]))
    end
    close(out)
end
hybrids()
