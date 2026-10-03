#!/usr/bin/env julia
# Code × filter-bank matrix on the same field (MATCHED_FORTRAN_2026-10.md §6):
#   F  : Fortran, production bank (filters_websky_finecell.dat, R_f 1.406→30.44 Mpc/h)   run/
#   F2 : Fortran, paper bank (filters_paper_2cell.dat, 2 a_latt → 36 Mpc = 1.704→24.26)  run_bank2/
#   JE : Julia exact, production bank (= our v4 configuration)                           run/julia
#   JE2: Julia exact, paper bank                                                         run_bank2/julia
#   JSP: Julia production path (GPU multires split, PERIODIC cores, nbuff 25, cf 22 / block 12), same global field
#        (Threefry noise is indexed globally, so N = 1056 gives the identical δ; cell convention identical),
#        production bank — the split-vs-exact leg (nbuff 26 overflows the GPU shell table)  run_split_nb25/
# F2/JE is "Websky as documented, run with the local Fortran" vs "ours"; compare with the full-sky 0.948.
#   usage: julia --project=validation -t 32 matched_bank_compare.jl <fortran_matched dir>
include(joinpath(@__DIR__, "matched_compare.jl"))

function bank_main()
    top = ARGS[1]
    rd(p) = read_pksc(joinpath(top, p))[1]
    cats = Dict("F" => cat_mpkvd(rd("run/output/fortran_merge.pksc.12345")),
                "F2" => cat_mpkvd(rd("run_bank2/output/fortran_merge.pksc.12345")),
                "JE" => cat_lagr(rd("run/julia/julia_exact_merged.pksc")),
                "JE2" => cat_lagr(rd("run_bank2/julia/julia_exact_merged.pksc")))
    labs = ["F", "F2", "JE", "JE2"]
    jsp = joinpath(top, "run_split_nb25", "julia_split_merged.pksc")
    isfile(jsp) && (cats["JSP"] = cat_lagr(rd("run_split_nb25/julia_split_merged.pksc")); push!(labs, "JSP"))
    AM = Dict(l => am_masses(cats[l]) for l in labs)
    dlin = linear_field() .* DZ; Dk = rfft(dlin)
    out = open(joinpath(@__DIR__, "..", "results", "matched_bank.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    for l in labs
        c, M = AM[l]
        say(@sprintf("  %-3s merged %d, in region %d; R_TH of the 5e12 rank %.3f Mpc/h", l, length(cats[l].R), length(c.R), minimum(c.R[M.>5e12])))
    end
    pairs = [("F", "JE"), ("F2", "JE"), ("JE2", "JE"), ("F2", "F"), ("F2", "JE2")]
    "JSP" in labs && append!(pairs, [("JSP", "JE"), ("F2", "JSP")])
    X = Dict(l => (c = AM[l][1]; s = findall(>(5e12), AM[l][2]); xi(c.xE[s], c.yE[s], c.zE[s])) for l in labs)
    say("\n== ξ(3–15 Mpc/h), M>5e12:  " * join([@sprintf("%s/%s %.4f   ", a, b, mean(X[a][WX]) / mean(X[b][WX])) for (a, b) in pairs]))
    say("   by r (F2/JE): " * join([@sprintf("%.3f ", X["F2"][i] / X["JE"][i]) for i in 1:10]))
    MB = [5e12, 8e12, 1.3e13, 2e13, 3.2e13, 7.9e13, 2.5e14]
    say("\n== b_E (k<0.1), ratios with jackknife σ")
    say("  " * " "^18 * join([@sprintf("%-17s", a * "/" * b) for (a, b) in pairs]))
    for i in 0:length(MB)-1
        B = Dict{String,Tuple{Float64,Vector{Float64}}}()
        for l in labs
            c, M = AM[l]; s = i == 0 ? findall(>(5e12), M) : findall(q -> MB[i] <= q < MB[i+1], M)
            h = deposit(c.xE[s], c.yE[s], c.zE[s]); B[l] = (bx(rfft(h), Dk), bx_jk(h, dlin))
        end
        say(@sprintf("  %-18s", i == 0 ? "M>5e12" : @sprintf("%.1e–%.1e", MB[i], MB[i+1])) *
            join([@sprintf("%.4f±%.4f  ", B[a][1] / B[b][1], jkerr(B[a][2] ./ B[b][2])) for (a, b) in pairs]))
    end
    say("  [full-sky Websky/ours v4fs: ξ(3–15) 0.948 ± 0.002; b ≈ 0.946–0.960 below 8e13]")
    close(out)
end
bank_main()
