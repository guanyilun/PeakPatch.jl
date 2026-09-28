#!/usr/bin/env julia
# Size of the Ω_m(a) 2LPT-coefficient fix (058207a): v3test (pre-fix) vs v3 (post-fix)
# oct000 RAW catalogs. Record order follows tile completion (not canonical), so halos are
# matched by fix-invariant attributes (RTHL, strain_11, strain_22, d2F — collapse does
# not use the ψ2 coefficient) and positions/velocities compared per matched pair.
using Printf, Statistics
const obs = -2618f0
cols = (1, 2, 3, 4, 5, 6, 7, 14, 15, 20)
function load(p)
    nh = open(io -> Int(read(io, Int32)), p)
    out = zeros(Float32, length(cols), nh)
    open(p) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        ch = 2_000_000; b = Vector{Float32}(undef, ch * 33); nd = 0
        while nd < nh
            m = min(ch, nh - nd); read!(io, view(b, 1:m*33))
            for k in 1:m, (j, c) in enumerate(cols); out[j, nd+k] = b[(k-1)*33+c]; end
            nd += m
        end
    end
    out
end
# fix-invariant key packed into a UInt128 (bit patterns of RTHL, strain_11, strain_22, d2F)
packkey(X) = [(UInt128(reinterpret(UInt32, X[7, i])) << 96) | (UInt128(reinterpret(UInt32, X[8, i])) << 64) |
              (UInt128(reinterpret(UInt32, X[9, i])) << 32) | UInt128(reinterpret(UInt32, X[10, i])) for i in axes(X, 2)]
function main()
    A = load("/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_v3test_oct000.pksc")
    B = load("/home/yguan/scratch/websky_6144/catalogs_v3/catalog_websky_6144_v3_oct000.pksc")
    ka = packkey(A); kb = packkey(B)
    pa = sortperm(ka); pb = sortperm(kb)
    ok = [k for k in eachindex(pa) if ka[pa[k]] == kb[pb[k]]]
    # keys must be unique for the pairing to be meaningful
    dup = count(k -> k > 1 && ka[pa[k]] == ka[pa[k-1]], eachindex(pa))
    @printf("halos %d / %d; matched invariant keys %d (%.6f); duplicate keys %d\n", size(A, 2), size(B, 2),
            length(ok), length(ok) / length(pa), dup)
    dp = [sqrt(sum(abs2, A[1:3, pa[k]] .- B[1:3, pb[k]])) for k in ok]
    dv = [sqrt(sum(abs2, A[4:6, pa[k]] .- B[4:6, pb[k]])) for k in ok]
    vr = [sqrt(sum(abs2, A[4:6, pa[k]])) for k in ok]
    r = [sqrt(sum(abs2, A[1:3, pa[k]] .- obs)) for k in ok]
    @printf("all: rms dpos %.4f Mpc/h (max %.3f); rms dvel %.3f km/s (max %.2f); rms |v| %.1f km/s\n",
            sqrt(mean(abs2, dp)), maximum(dp), sqrt(mean(abs2, dv)), maximum(dv), sqrt(mean(abs2, vr)))
    for (lo, hi) in ((0, 1300), (1300, 2300), (2300, 3500), (3500, 7000))
        m = (r .>= lo) .& (r .< hi)
        @printf("chi %4d-%4d n=%9d rms dpos %.4f Mpc/h  rms dvel %.3f km/s  rel dvel %.5f\n", lo, hi, count(m),
                sqrt(mean(abs2, dp[m])), sqrt(mean(abs2, dv[m])), sqrt(mean(abs2, dv[m])) / sqrt(mean(abs2, vr[m])))
    end
end
main()
