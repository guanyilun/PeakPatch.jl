#!/usr/bin/env julia
# Split-path block-size scan on the same field (MATCHED_FORTRAN_2026-10.md §7): production GPU multires path,
# periodic cores, nbuff 25, coarse grid M = 4·cf, block = 1056/M, vs Julia exact global field (JE) and the
# Fortran chain with the paper filter bank (F2). If the split raises clustering, JS/JE → 1 as the block grows.
#   usage: julia --project=validation -t 32 matched_split_scan.jl <fortran_matched dir>
include(joinpath(@__DIR__, "matched_compare.jl"))

function scan_main()
    top = ARGS[1]
    rd(p) = read_pksc(joinpath(top, p))[1]
    cats = Dict("JE" => cat_lagr(rd("run/julia/julia_exact_merged.pksc")),
                "F2" => cat_mpkvd(rd("run_bank2/output/fortran_merge.pksc.12345")))
    runs = [(4, "run_split_cf4"), (8, "run_split_cf8"), (11, "run_split_cf11"), (22, "run_split_nb25"), (33, "run_split_cf33")]
    labs = String[]
    for (cf, d) in runs
        f = joinpath(top, d, "julia_split_merged.pksc"); isfile(f) || continue
        l = "cf$(cf)"; cats[l] = cat_lagr(rd(joinpath(d, "julia_split_merged.pksc"))); push!(labs, l)
    end
    all = vcat(["JE", "F2"], labs)
    AM = Dict(l => am_masses(cats[l]) for l in all)
    dlin = linear_field() .* DZ; Dk = rfft(dlin)
    out = open(joinpath(@__DIR__, "..", "results", "matched_split_scan.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    X = Dict(l => (c = AM[l][1]; s = findall(>(5e12), AM[l][2]); xi(c.xE[s], c.yE[s], c.zE[s])) for l in all)
    MB = [5e12, 8e12, 1.3e13, 3.2e13, 7.9e13, 2.5e14]
    bE(l, i) = begin
        c, M = AM[l]; s = i == 0 ? findall(>(5e12), M) : findall(q -> MB[i] <= q < MB[i+1], M)
        h = deposit(c.xE[s], c.yE[s], c.zE[s]); (bx(rfft(h), Dk), bx_jk(h, dlin))
    end
    B = Dict((l, i) => bE(l, i) for l in all, i in 0:length(MB)-1)
    say("split path vs exact field (same δ): block = 1056/(4·cf); JS/JE and F2/JS (jackknife σ on b)")
    say(@sprintf("%-6s %6s %9s %9s | %-16s %-16s %-16s %-16s %-16s %-16s | %9s %-16s", "cf", "block", "merged", "ξ JS/JE",
                 "b>5e12 JS/JE", "5–8e12", "0.8–1.3e13", "1.3–3.2e13", "3.2–7.9e13", "0.8–2.5e14", "ξ F2/JS", "b>5e12 F2/JS"))
    r(a, b, i) = @sprintf("%.4f±%.4f", B[(a, i)][1] / B[(b, i)][1], jkerr(B[(a, i)][2] ./ B[(b, i)][2]))
    for l in labs
        cf = parse(Int, l[3:end])
        say(@sprintf("%-6s %6d %9d %9.4f | ", l, 1056 ÷ (4cf), length(cats[l].R), mean(X[l][WX]) / mean(X["JE"][WX])) *
            join([rpad(r(l, "JE", i), 17) for i in 0:length(MB)-1]) *
            @sprintf("| %9.4f %s", mean(X["F2"][WX]) / mean(X[l][WX]), r("F2", l, 0)))
    end
    say(@sprintf("reference: F2/JE ξ %.4f, b>5e12 %s;  full-sky Websky/ours v4fs ξ 0.948 ± 0.002", mean(X["F2"][WX]) / mean(X["JE"][WX]), r("F2", "JE", 0)))
    close(out)
end
scan_main()
