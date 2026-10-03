#!/usr/bin/env julia
# R_TH disagreement diagnosis on the same-field raw peaks (MATCHED_FORTRAN_2026-10.md, step 2).
# Fortran get_homel stops building shells at ir2upp ≈ (rad(m0)+2)² where m0 is the first shell with r² > ir2min
# (ir2min = min(⌊(1.75 R_f/a)²⌋, ⌊(40/a − 1)²⌋)), and its outward search for the Fbar < fcrit crossing never
# runs (peakvoidsubs.f90:447 reuses jp). So a Fortran R_TH can never exceed rad(m0). Julia builds all shells
# out to nhunt and searches outward. Classes of matched peaks:
#   CAP   : R_J > rad(m0)·a  (Julia's crossing lies beyond where Fortran can look)
#   OTHER : everything else
# Then a hybrid: Julia raw peaks with R_TH replaced by the Fortran value for the CAP class only, merged and
# rank-matched exactly as the others → how much of the F/JE bias and ξ offset the cap accounts for.
#   usage: julia --project=validation -t 32 matched_rth_diag.jl <rundir>
include(joinpath(@__DIR__, "matched_compare.jl"))

isrep(n) = begin m = n; while m > 0 && m % 4 == 0; m ÷= 4; end; m % 8 != 7 end   # sum of three squares
nextrep(n) = begin k = n + 1; while !isrep(k); k += 1; end; k end
const NHUNT = min(NB - 1, floor(Int, 30.4358 * 1.75 / A))
function rad_m0(Rf)                                    # Fortran rad(m0) in Mpc/h, or Inf if beyond the hunt radius
    ir2min = min(floor(Int, (1.75Rf / A)^2), floor(Int, (40 / A - 1)^2))
    r2 = nextrep(ir2min); r2 > NHUNT^2 ? Inf : sqrt(r2) * A
end

function diag()
    rd(f) = read_pksc(joinpath(RUN, f))[1]
    fr = filter(q -> q.RTHL > 0, rd("output/fortran_raw.pksc.12345"))
    jr = filter(q -> q.RTHL > 0, rd("julia/julia_exact_raw.pksc"))
    Fr = cat_lagr(fr); Jr = cat_lagr(jr)
    m = match(Fr, Jr; tol=0.5A); ok = findall(>(0), m)
    out = open(joinpath(@__DIR__, "..", "results", "matched_rth_diag.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    say("same-field raw peaks: Fortran $(length(Fr.R)), Julia $(length(Jr.R)), matched (same cell) $(length(ok)); nhunt = $(NHUNT) cells")
    r0 = [rad_m0(Fr.Rf[i]) for i in ok]
    RF = Fr.R[ok]; RJ = [Jr.R[m[i]] for i in ok]; ratio = RJ ./ RF
    cap = RJ .> r0 .+ 1e-3
    say(@sprintf("  Fortran R_TH above its rad(m0) limit (should be 0): %d", count(RF .> r0 .+ 1e-3)))
    MB = [5e12, 1.3e13, 5e13, 1e16]
    say("\n  by Fortran mass:            n      CAP frac  | CAP: R_J/R_F median  p90 | OTHER: |Δ|>2% frac  R_J>1.02R_F  R_J<0.98R_F  median")
    for i in 1:3
        s = findall(k -> MB[i] <= mass(RF[k]) < MB[i+1], eachindex(ok)); isempty(s) && continue
        c = filter(k -> cap[k], s); o = filter(k -> !cap[k], s)
        say(@sprintf("  %.1e–%.1e  %8d   %.4f    |   %.3f  %.3f   |   %.4f     %.4f     %.4f    %.5f", MB[i], MB[i+1], length(s), length(c) / length(s),
                     isempty(c) ? NaN : median(ratio[c]), isempty(c) ? NaN : quantile(ratio[c], 0.9),
                     count(k -> abs(ratio[k] - 1) > 0.02, o) / length(o), count(k -> ratio[k] > 1.02, o) / length(o),
                     count(k -> ratio[k] < 0.98, o) / length(o), median(ratio[o])))
    end
    # OTHER disagreements: by filter, and how far apart in shells
    o = findall(.!cap)
    say("\n  OTHER class, |R_J/R_F − 1| > 2%, by filter R_f (Mpc/h): count / class count")
    for rf in sort(unique(Fr.Rf[ok[o]]))
        s = filter(k -> abs(Fr.Rf[ok[k]] - rf) < 1e-4, o); d = count(k -> abs(ratio[k] - 1) > 0.02, s)
        say(@sprintf("    R_f %7.4f  n=%8d  disagree %.4f  (J>F %.4f, J<F %.4f)", rf, length(s), d / max(1, length(s)),
                     count(k -> ratio[k] > 1.02, s) / max(1, length(s)), count(k -> ratio[k] < 0.98, s) / max(1, length(s))))
    end
    # unmatched peaks
    mj = match(Jr, Fr; tol=0.5A)
    say(@sprintf("\n  unmatched: Fortran-only %d (%.4f), median M %.2e;  Julia-only %d (%.4f), median M %.2e",
                 length(Fr.R) - length(ok), 1 - length(ok) / length(Fr.R), median(mass.(Fr.R[m.==0])),
                 count(==(0), mj), count(==(0), mj) / length(Jr.R), median(mass.(Jr.R[mj.==0]))))

    # hybrid: Julia raw with Fortran R_TH on the CAP class
    Rnew = Float32.(Jr.R); for (k, i) in enumerate(ok); cap[k] && (Rnew[m[i]] = Float32(RF[k])); end
    hy = [typeof(jr[j])(ntuple(f -> f == 7 ? Rnew[j] : getfield(jr[j], f), fieldcount(typeof(jr[j])))...) for j in eachindex(jr)]
    cats = Dict("F" => cat_mpkvd(rd("output/fortran_merge.pksc.12345")), "JE" => cat_lagr(rd("julia/julia_exact_merged.pksc")),
                "JEcap" => cat_lagr(merge_catalog(hy; volume_reduction=true)))
    labs = ["F", "JEcap", "JE"]
    AM = Dict(l => am_masses(cats[l]) for l in labs)
    dlin = linear_field() .* DZ; Dk = rfft(dlin)
    MB2 = [5e12, 8e12, 1.3e13, 2e13, 3.2e13, 7.9e13, 2.5e14]
    say("\n== hybrid JEcap = Julia raw with Fortran R_TH where Julia's lies beyond the Fortran limit")
    say("  b_E (k<0.1)         F        JEcap    JE       F/JE           F/JEcap        JEcap/JE")
    for i in 0:length(MB2)-1
        B = Dict{String,Tuple{Float64,Vector{Float64}}}()
        for l in labs
            c, M = AM[l]; s = i == 0 ? findall(>(5e12), M) : findall(q -> MB2[i] <= q < MB2[i+1], M)
            h = deposit(c.xE[s], c.yE[s], c.zE[s]); B[l] = (bx(rfft(h), Dk), bx_jk(h, dlin))
        end
        r(a, b) = @sprintf("%.4f±%.4f", B[a][1] / B[b][1], jkerr(B[a][2] ./ B[b][2]))
        say(@sprintf("  %-18s %.4f   %.4f   %.4f   ", i == 0 ? "M>5e12" : @sprintf("%.1e–%.1e", MB2[i], MB2[i+1]),
                     B["F"][1], B["JEcap"][1], B["JE"][1]) * r("F", "JE") * "  " * r("F", "JEcap") * "  " * r("JEcap", "JE"))
    end
    X = Dict(l => (c = AM[l][1]; s = findall(>(5e12), AM[l][2]); xi(c.xE[s], c.yE[s], c.zE[s])) for l in labs)
    say(@sprintf("  ξ(3–15), M>5e12: F/JE %.4f   F/JEcap %.4f   JEcap/JE %.4f", mean(X["F"][WX]) / mean(X["JE"][WX]),
                 mean(X["F"][WX]) / mean(X["JEcap"][WX]), mean(X["JEcap"][WX]) / mean(X["JE"][WX])))
    close(out)
end
diag()
