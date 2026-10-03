#!/usr/bin/env julia
# Split (production GPU multires path) vs exact global field, same δ (MATCHED_FORTRAN_2026-10.md §8):
#   ξ / b_E ratios, raw-peak cross-match (selection vs R_TH), and hybrids that swap R_TH between the two to tell
#   whether the clustering change comes from WHICH peaks are found or from their collapse radii.
#   usage: MC_N=… MC_NB=… MC_A=… julia --project=validation -t 32 split_vs_exact.jl <rundir with fields/> \
#          <exact dir> <split dir> <label>
include(joinpath(@__DIR__, "matched_compare.jl"))
withR(h, R) = typeof(h)(ntuple(f -> f == 7 ? Float32(R) : getfield(h, f), fieldcount(typeof(h)))...)

function sve_main()
    exd, spd, lab = ARGS[2], ARGS[3], ARGS[4]
    rd(p) = read_pksc(p)[1]
    jer = filter(q -> q.RTHL > 0, rd(joinpath(exd, "julia_exact_raw.pksc"))); jsr = filter(q -> q.RTHL > 0, rd(joinpath(spd, "julia_split_raw.pksc")))
    Er = cat_lagr(jer); Sr = cat_lagr(jsr)
    out = open(joinpath(@__DIR__, "..", "results", "split_vs_exact_$(lab).txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    say("split vs exact, same field: N=$N, region $(round(XHI - XLO; digits=1)) Mpc/h; label $lab")
    ins(c) = findall(i -> inreg(c.xL[i], c.yL[i], c.zL[i]), eachindex(c.R))
    Ei = ins(Er); Si = ins(Sr); Er = sub(Er, Ei); Sr = sub(Sr, Si); jer = jer[Ei]; jsr = jsr[Si]
    m = match(Er, Sr; tol=0.5A); ok = findall(>(0), m)
    sameRf = count(i -> abs(Er.Rf[i] - Sr.Rf[m[i]]) < 1e-3, ok)
    say(@sprintf("raw peaks in region: exact %d, split %d; exact peaks with a split peak in the same cell %.4f (same filter %.4f)",
                 length(Er.R), length(Sr.R), length(ok) / length(Er.R), sameRf / length(Er.R)))
    rr = [Sr.R[m[i]] / Er.R[i] for i in ok]
    say(@sprintf("R_TH split/exact over matched: median %.5f  mean %.5f  p05 %.4f  p95 %.4f  |Δ|>2%% %.4f",
                 median(rr), mean(rr), quantile(rr, 0.05), quantile(rr, 0.95), count(r -> abs(r - 1) > 0.02, rr) / length(rr)))
    for (lo, hi) in ((5e12, 1.3e13), (1.3e13, 5e13), (5e13, 1e16))
        s = filter(i -> lo <= mass(Er.R[i]) < hi, ok); isempty(s) && continue
        r2 = [Sr.R[m[i]] / Er.R[i] for i in s]
        say(@sprintf("  M_exact %.1e–%.1e: n=%d  median %.4f  mean %.4f  frac(S>1.02E) %.4f  frac(S<0.98E) %.4f", lo, hi, length(s),
                     median(r2), mean(r2), count(>(1.02), r2) / length(r2), count(<(0.98), r2) / length(r2)))
    end
    dE = [sqrt((Sr.xE[m[i]] - Er.xE[i])^2 + (Sr.yE[m[i]] - Er.yE[i])^2 + (Sr.zE[m[i]] - Er.zE[i])^2) for i in ok]
    say(@sprintf("Eulerian position difference split−exact: median %.4f  p90 %.4f Mpc/h", median(dE), quantile(dE, 0.9)))
    # hybrids: exact peaks with split R_TH (selection of exact, radii of split) and exact peaks with split displacements
    RS = Dict(i => Sr.R[m[i]] for i in ok)
    hyR = [haskey(RS, i) ? withR(jer[i], RS[i]) : jer[i] for i in eachindex(jer)]
    cats = Dict("E" => cat_lagr(merge_catalog(jer; volume_reduction=true)), "S" => cat_lagr(merge_catalog(jsr; volume_reduction=true)),
                "E_Rsplit" => cat_lagr(merge_catalog(hyR; volume_reduction=true)))
    labs = ["E", "S", "E_Rsplit"]
    AM = Dict(l => am_masses(cats[l]) for l in labs)
    dlin = linear_field() .* DZ; Dk = rfft(dlin)
    X = Dict(l => (c = AM[l][1]; s = findall(>(5e12), AM[l][2]); xi(c.xE[s], c.yE[s], c.zE[s])) for l in labs)
    say(@sprintf("\nξ(3–15), M>5e12:  S/E %.4f   E_Rsplit/E %.4f", mean(X["S"][WX]) / mean(X["E"][WX]), mean(X["E_Rsplit"][WX]) / mean(X["E"][WX])))
    say("  by r (S/E): " * join([@sprintf("%.3f ", X["S"][i] / X["E"][i]) for i in 1:10]))
    MB = [5e12, 8e12, 1.3e13, 3.2e13, 7.9e13, 2.5e14]
    say("b_E (k<0.1)          S/E               E_Rsplit/E")
    for i in 0:length(MB)-1
        B = Dict{String,Tuple{Float64,Vector{Float64}}}()
        for l in labs
            c, M = AM[l]; s = i == 0 ? findall(>(5e12), M) : findall(q -> MB[i] <= q < MB[i+1], M)
            h = deposit(c.xE[s], c.yE[s], c.zE[s]); B[l] = (bx(rfft(h), Dk), bx_jk(h, dlin))
        end
        r(a, b) = @sprintf("%.4f±%.4f   ", B[a][1] / B[b][1], jkerr(B[a][2] ./ B[b][2]))
        say(@sprintf("  %-18s ", i == 0 ? "M>5e12" : @sprintf("%.1e–%.1e", MB[i], MB[i+1])) * r("S", "E") * r("E_Rsplit", "E"))
    end
    close(out)
end
sve_main()
