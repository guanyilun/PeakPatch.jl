#!/usr/bin/env julia
# Why does merge_catalog take ~4.6 h at production scale (333M pre-merge halos) but 51 s for
# 25.7M (bench)? Hypothesis: the fixed NC=256 spatial hash (src/Merger/Exclusion.jl) → at a
# 5236 Mpc/h box the hash cells are ~20.5 Mpc/h and hold ~20 halos each, so every halo scans
# ≥27 cells of linked-list entries scattered over GBs (cache-miss bound).
#
# Test on real data: halos of the production oct000 raw catalog inside the observer-corner
# sub-box (the densest lightcone region), exclusion timed with
#   (a) the src algorithm verbatim (linked-list hash) at the production cell size and finer;
#   (b) a counting-sorted cell list (CSR, halos stored cell-contiguously) at the same sizes.
# All variants must return the identical survivor set. src/ is not modified.
#
#   julia --project=validation -t 1 validation/performance/exclusion_scaling.jl [catalog] [subbox_Mpch]
using Mmap, Printf
using PeakPatch.Exclusion: SpatialHash, lagrangian_exclusion!

cat = length(ARGS) >= 1 ? ARGS[1] :
    "/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_prod_oct000.pksc"
Ls  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1309.0     # 1/4 of the box per axis
obs = (-2618.0, -2618.0, -2618.0)

# ---- load positions + RTHL of halos inside the sub-box ----
io = open(cat); nh = read(io, Int32); close(io)
nf = (filesize(cat) - 12) ÷ (4 * Int(nh))
A = Mmap.mmap(open(cat), Matrix{Float32}, (nf, Int(nh)), 12)
sel = Int[]
for i in 1:size(A, 2)
    (A[1, i] - obs[1] < Ls && A[2, i] - obs[2] < Ls && A[3, i] - obs[3] < Ls) && push!(sel, i)
end
x = Float64.(A[1, sel]); y = Float64.(A[2, sel]); z = Float64.(A[3, sel]); r = Float64.(A[7, sel])
n = length(x)
@printf("catalog: %d halos total; sub-box %.0f Mpc/h: %d halos (median RTHL %.2f, max %.1f Mpc/h)\n",
        nh, Ls, n, sort(r)[n ÷ 2], maximum(r))
order = sortperm(r; rev=true)
dmin = (minimum(x) - 1, minimum(y) - 1, minimum(z) - 1)
span = maximum((maximum(x) - dmin[1], maximum(y) - dmin[2], maximum(z) - dmin[3])) + 1

# ---- (a) src algorithm with a hash of arbitrary resolution (same struct, same loop) ----
function hash_ll(x, y, z, nc, dmin, span)
    cs = span / nc
    hoc = zeros(Int32, nc, nc, nc); ll = zeros(Int32, length(x))
    for m in eachindex(x)
        ix = clamp(floor(Int, (x[m] - dmin[1]) / cs) + 1, 1, nc)
        iy = clamp(floor(Int, (y[m] - dmin[2]) / cs) + 1, 1, nc)
        iz = clamp(floor(Int, (z[m] - dmin[3]) / cs) + 1, 1, nc)
        ll[m] = hoc[ix, iy, iz]; hoc[ix, iy, iz] = Int32(m)
    end
    SpatialHash(hoc, ll, nc, cs, dmin)
end

# ---- (b) CSR cell list: halos permuted cell-contiguously, coordinates packed per cell ----
function exclusion_csr(x, y, z, r, order, nc, dmin, span)
    cs = span / nc; n = length(x)
    cell = Vector{Int}(undef, n)
    @inbounds for m in 1:n
        ix = clamp(floor(Int, (x[m] - dmin[1]) / cs), 0, nc - 1)
        iy = clamp(floor(Int, (y[m] - dmin[2]) / cs), 0, nc - 1)
        iz = clamp(floor(Int, (z[m] - dmin[3]) / cs), 0, nc - 1)
        cell[m] = 1 + ix + nc * (iy + nc * iz)
    end
    start = zeros(Int, nc^3 + 1)
    @inbounds for m in 1:n; start[cell[m] + 1] += 1; end
    start[1] = 1; @inbounds for c in 1:nc^3; start[c + 1] += start[c]; end
    fill_ = copy(start); perm = Vector{Int}(undef, n)
    @inbounds for m in 1:n; perm[fill_[cell[m]]] = m; fill_[cell[m]] += 1; end
    px = x[perm]; py = y[perm]; pz = z[perm]; pr = r[perm]
    inv = similar(perm); inv[perm] = 1:n
    alive = fill(true, n)                                 # indexed in permuted order
    @inbounds for rank in 1:n
        i = inv[order[rank]]; alive[i] || continue
        ri = pr[i]; xi = px[i]; yi = py[i]; zi = pz[i]
        d = ceil(Int, ri / cs); c = cell[perm[i]] - 1
        ci = c % nc; cj = (c ÷ nc) % nc; ck = c ÷ (nc * nc)
        for kz in max(ck - d, 0):min(ck + d, nc - 1), jy in max(cj - d, 0):min(cj + d, nc - 1)
            base = 1 + nc * (jy + nc * kz)
            lo = start[base + max(ci - d, 0)]; hi = start[base + min(ci + d, nc - 1) + 1] - 1
            for j in lo:hi                                # one contiguous run per (y,z) row
                (j != i && alive[j] && ri >= pr[j]) || continue
                (px[j] - xi)^2 + (py[j] - yi)^2 + (pz[j] - zi)^2 < ri^2 && (alive[j] = false)
            end
        end
    end
    surv = fill(false, n); surv[perm[alive]] .= true
    return surv
end

prod_cs = (5236.0 + 4 * 30) / 256                          # production NC=256 cell ≈ 20.9 Mpc/h
println("variant               cell[Mpc/h]   nc     time[s]   survivors")
ref = nothing
for cs in (prod_cs, 10.0, 5.0, 2.5)
    nc = max(1, round(Int, span / cs))
    s = fill(true, n); sh = hash_ll(x, y, z, nc, dmin, span)
    t = @elapsed lagrangian_exclusion!(s, x, y, z, r, order, sh)
    global ref = ref === nothing ? s : ref
    @printf("src linked-list     %10.2f %6d %10.1f %11d %s\n", span / nc, nc, t, count(s), s == ref ? "" : "MISMATCH")
    nc^3 <= 2^31 || continue
    t = @elapsed (s2 = exclusion_csr(x, y, z, r, order, nc, dmin, span))
    @printf("CSR cell list       %10.2f %6d %10.1f %11d %s\n", span / nc, nc, t, count(s2), s2 == ref ? "" : "MISMATCH")
end
