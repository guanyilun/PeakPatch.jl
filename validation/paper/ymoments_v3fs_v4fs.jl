#!/usr/bin/env julia
# oct000 halo counts and tSZ moment Σ(M/1e14)^{5/3} (∝ mean y) per mass and χ bin, v4fs vs v3fs
# ABUNDANCE-MATCHED catalogs (both _AMfs). Output: results/ymoments_oct000_v3fs_v4fs.txt (V4_RESULTS_2026-10-01.md).
using Printf
const RHO_M = 2.775e11 * 0.31
function mom(p; obs=-2618f0)
    nh = open(io -> Int(read(io, Int32)), p)
    zb = [0.0, 1000.0, 2000.0, 3000.0, 6000.0]       # comoving distance bins (Mpc/h)
    mb = [1e11, 1e12, 3e12, 1e13, 1e14, 1e16]
    S = zeros(length(mb) - 1, length(zb) - 1); N = zeros(Int, length(mb) - 1, length(zb) - 1)
    open(p) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        ch = 2_000_000; b = Vector{Float32}(undef, ch * 33); nd = 0
        while nd < nh
            m = min(ch, nh - nd); read!(io, view(b, 1:m*33))
            @inbounds for k in 1:m
                o = (k - 1) * 33; M = 4 / 3 * π * RHO_M * Float64(b[o+7])^3
                r = sqrt((b[o+1] - obs)^2 + (b[o+2] - obs)^2 + (b[o+3] - obs)^2)
                i = searchsortedlast(mb, M); j = searchsortedlast(zb, r)
                (1 <= i < length(mb) && 1 <= j < length(zb)) || continue
                S[i, j] += (M / 1e14)^(5 / 3); N[i, j] += 1
            end
            nd += m
        end
    end
    S, N, nh
end
S3, N3, n3 = mom("/home/yguan/scratch/websky_6144/catalogs_v3/catalog_websky_6144_v3_oct000_AMfs.pksc")
S4, N4, n4 = mom("/home/yguan/scratch/websky_6144/catalogs_v4/catalog_websky_6144_v4_oct000_AMfs.pksc")
@printf("halos v3fs %d v4fs %d\nΣ(M/1e14)^{5/3} (∝ mean-y) v4/v3, and N v4/v3; rows M 1e11-1e12,1e12-3e12,3e12-1e13,1e13-1e14,>1e14; cols χ 0-1,1-2,2-3,3-6 Gpc/h\n", n3, n4)
for i in 1:5
    @printf("M bin %d  S: %s   N: %s\n", i, join([@sprintf("%.3f", S4[i, j] / S3[i, j]) for j in 1:4], " "), join([@sprintf("%.3f", N4[i, j] / N3[i, j]) for j in 1:4], " "))
end
@printf("total Σ M^{5/3} v4/v3 = %.4f  (M>1e12: %.4f)\n", sum(S4) / sum(S3), sum(S4[2:end, :]) / sum(S3[2:end, :]))
