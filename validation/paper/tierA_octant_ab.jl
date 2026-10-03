#!/usr/bin/env julia
# Single-octant Tier-A A/B (MATCHED_FORTRAN_2026-10.md §9): the full-sky caps (tierA_fullsky.jl) that lie in one
# octant, measured on the FULL Websky catalogue and on two of our catalogues of that octant (same seed, same
# AM table): A = v4 (original splice), B = v5test (gaussian_split). Same-noise A/B, so B/A isolates the splice
# change; W/A and W/B are the Websky gap before and after.
#   env: OCT (default 000), CAT_A, CAT_B (AM'd pksc paths), LAB_A, LAB_B
include(joinpath(@__DIR__, "tierA_fullsky.jl"))       # caps, capof, helpers; fullsky_main() guarded

function octab_main()
    oc = get(ENV, "OCT", "000")
    cats = [(get(ENV, "LAB_A", "v4"), ENV["CAT_A"]), (get(ENV, "LAB_B", "v5test"), ENV["CAT_B"])]
    ks = findall(c -> octof(c) == oc, CENT); nc = length(ks)
    @info "caps in octant" oc nc
    capsW = Dict(k => Cap((0.0, 0.0, 0.0)) for k in ks)
    capsC = [Dict(k => Cap(obs_of(oc)) for k in ks) for _ in cats]
    nh = open(io -> Int(read(io, Int32)), WFULL); MAXH > 0 && (nh = min(nh, MAXH))
    open(WFULL) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        ch = 4_000_000; b = Vector{Float32}(undef, ch * 10); nd = 0
        while nd < nh
            m = min(ch, nh - nd); read!(io, view(b, 1:m*10))
            @inbounds for k in 1:m
                o = (k - 1) * 10
                X = HUB * Float64(b[o+1]); Y = HUB * Float64(b[o+2]); Z = HUB * Float64(b[o+3])
                r = sqrt(X^2 + Y^2 + Z^2); r > 0 || continue
                kc = capof(X / r, Y / r, Z / r); haskey(capsW, kc) || continue
                M = mass(HUB * Float64(b[o+7])); z = chi_to_z(CHI2Z, r)
                addhalo!(capsW[kc], X, Y, Z, Float64(b[o+4]), Float64(b[o+5]), Float64(b[o+6]), M, r, z)
            end
            nd += m
        end
    end
    ob = obs_of(oc)
    for (ic, (lab, path)) in enumerate(cats)
        C = capsC[ic]
        stream(path, (bf, b) -> begin
            R = Float64(bf[b+7]); R > 0 || return
            X = Float64(bf[b+1]); Y = Float64(bf[b+2]); Z = Float64(bf[b+3])
            dx = X - ob[1]; dy = Y - ob[2]; dz = Z - ob[3]; r = sqrt(dx^2 + dy^2 + dz^2)
            kc = capof(dx / r, dy / r, dz / r); haskey(C, kc) || return
            addhalo!(C[kc], X, Y, Z, Float64(bf[b+4]), Float64(bf[b+5]), Float64(bf[b+6]), mass(R), r, chi_to_z(CHI2Z, r))
        end)
        @info "streamed" lab
    end
    PW = Dict(k => pairstats(capsW[k], collect(CENT[k]), 10k + 1) for k in ks)
    PC = [Dict(k => pairstats(capsC[ic][k], collect(CENT[k]), 10k + 2) for k in ks) for ic in eachindex(cats)]
    out = joinpath(@__DIR__, "results", "tierA_octant_ab_oct$(oc)_$(cats[1][1])_$(cats[2][1]).txt")
    io = open(out, "w"); say(a...) = begin s = string(a...); println(s); println(io, s) end
    la, lb = cats[1][1], cats[2][1]
    # per-cap paired ratio for A vs B (same sky, same noise); W vs each as ratio of means
    sem(v) = begin ok = filter(isfinite, v); (mean(ok), std(ok) / sqrt(length(ok))) end
    function row(label, fW, fA, fB)
        w = [fW(k) for k in ks]; a = [fA(k) for k in ks]; b = [fB(k) for k in ks]
        (mw, sw), (ma, sa), (mb, sb) = sem(w), sem(a), sem(b)
        pr = filter(isfinite, b ./ a); mpr, spr = mean(pr), std(pr) / sqrt(length(pr))
        say(@sprintf("%-30s W %9.4g | %s %9.4g | %s %9.4g || W/%s %6.3f±%5.3f  W/%s %6.3f±%5.3f  %s/%s (paired) %6.4f±%6.4f",
                     label, mw, la, ma, lb, mb, la, mw / ma, sqrt((sw / ma)^2 + (mw * sa / ma^2)^2), lb, mw / mb,
                     sqrt((sw / mb)^2 + (mw * sb / mb^2)^2), lb, la, mpr, spr))
    end
    say("== single-octant Tier-A A/B, octant $oc, $nc caps; A = $la, B = $lb (same seed, same AM table)")
    row(@sprintf("ξ mean %.0f-%.0f Mpc/h", WIN_XI...), k -> winmean(PW[k].xi, XI_EDGES, WIN_XI), k -> winmean(PC[1][k].xi, XI_EDGES, WIN_XI), k -> winmean(PC[2][k].xi, XI_EDGES, WIN_XI))
    for (i, r) in enumerate(rc(XI_EDGES))
        row(@sprintf("ξ r=%.2f", r), k -> PW[k].xi[i], k -> PC[1][k].xi[i], k -> PC[2][k].xi[i])
    end
    for ib in 1:length(BM_EDGES)-1
        row(@sprintf("b∝√ξ M %.1e-%.1e", BM_EDGES[ib], BM_EDGES[ib+1]), k -> sqrt(max(PW[k].bxi[ib], 0)), k -> sqrt(max(PC[1][k].bxi[ib], 0)), k -> sqrt(max(PC[2][k].bxi[ib], 0)))
    end
    sig(v) = v[1] > 1 ? sqrt(v[3] / v[1] - (v[2] / v[1])^2) : NaN
    for iv in 1:length(VM_EDGES)-1
        row(@sprintf("σ_vr M %.1e-%.1e", VM_EDGES[iv], VM_EDGES[iv+1]), k -> sig(capsW[k].vs[:, iv]), k -> sig(capsC[1][k].vs[:, iv]), k -> sig(capsC[2][k].vs[:, iv]))
    end
    v12c = [(V12_EDGES[i] + V12_EDGES[i+1]) / 2 for i in 1:length(V12_EDGES)-1]
    wm(v) = mean(v[[i for i in eachindex(v12c) if WIN_V12[1] <= v12c[i] <= WIN_V12[2]]])
    row(@sprintf("v12 mean %.0f-%.0f Mpc/h", WIN_V12...), k -> wm(PW[k].v12), k -> wm(PC[1][k].v12), k -> wm(PC[2][k].v12))
    for (j, t) in enumerate(DNDZ_M), id in 1:length(DNDZ_EDGES)-1
        row(@sprintf("dN/dz M>%.0e z %.2f-%.2f", t, DNDZ_EDGES[id], DNDZ_EDGES[id+1]), k -> capsW[k].dndz[j, id], k -> capsC[1][k].dndz[j, id], k -> capsC[2][k].dndz[j, id])
    end
    close(io); @info "wrote" out
end
octab_main()
