#!/usr/bin/env julia
# Tier-A catalog statistics on the current catalogs (v3fs), with an error bar (plan A2/A3).
#
# The 2026-06/07 Tier-A numbers (BIAS_PAIRWISE_V12_2026-07-16.md etc.) compared ONE cap of the
# superseded oct000_finecell_AM catalog with Websky's public 10°x10° patch, without an error bar.
# Here the SAME cap (inscribed in the Websky patch, same shell and bins as those scripts) is
# placed at the centre of each of our 8 octants. The 8 caps give the cap-to-cap (sample-variance)
# scatter, and each statistic is reported as
#     Websky patch  vs  mean ± std over our 8 caps  ->  ratio, and z = (W − mean)/std
# i.e. "is Websky's patch consistent with a draw from our cap distribution".
# Statistics: N(>M|z) per sr, dN/dz, σ_vr(M), ξ(r) (Landy-Szalay), relative b(M), v12(r).
#
# Also (plan A4) the full-sky RAW and AM N(>M|z) against the Tinker08 expectation, from the
# same streaming pass over the 8 raw and 8 AM catalogs.
#
# v3fs catalogs: Eulerian positions (Mpc/h, cols 1-3), velocities km/s (4-6), RTHL (7).
# Websky patch: positions Mpc, velocities km/s, R Mpc; observer at the origin (×h here).
using Printf, Random, Statistics, LinearAlgebra
using PeakPatch
import PeakPatch.Cosmology: CosmologyParams, chi, build_chi_to_z, chi_to_z, growth_factor
import PeakPatch.MassFunction: precompute_sigma, tinker_dndlnM

const HUB = 0.68; const RHO_M = 2.775e11 * 0.31
const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
const CHI2Z = build_chi_to_z(COSMO; z_max=6.0)
const S3 = "/home/yguan/scratch/websky_6144/catalogs_v3"
const WSKY = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/halos_10x10.pksc"
const OCTS = split(get(ENV, "TIERA_OCTS", "000,001,010,011,100,101,110,111"), ",")   # smoke test: TIERA_OCTS=000
const MAXH = parse(Int, get(ENV, "TIERA_MAXH", "0"))                              # smoke test: read only the first MAXH halos per file
const OUT = MAXH > 0 ? "/tmp/tierA_smoke.txt" : joinpath(@__DIR__, "results", "tierA_v3fs.txt")
# octZYX bit = 1 -> observer at +2618 on that axis -> that axis looks toward −
obs_of(o) = (o[3] == '1' ? 2618.0 : -2618.0, o[2] == '1' ? 2618.0 : -2618.0, o[1] == '1' ? 2618.0 : -2618.0)
axis_of(o) = (a = [o[3] == '1' ? -1.0 : 1.0, o[2] == '1' ? -1.0 : 1.0, o[1] == '1' ? -1.0 : 1.0]; a ./ norm(a))
mass(R) = 4 / 3 * π * RHO_M * R^3

# ---- binning (same as compare_{bias,xi,pairwise_velocity,velocities}.jl) ----
const SH = (1600.0, 2000.0)                                   # clustering shell, Mpc/h (z≈0.6-0.8)
const MLOW_XI = 5e12; const MCUT_V12 = 1e13
const XI_EDGES = collect(10.0 .^ range(log10(1.0), log10(50.0); length=11))
const BM_EDGES = [10.0^l for l in (12.7, 13.1, 13.5, 13.9, 14.4)]
const V12_EDGES = collect(range(2.0, 60.0; length=16))
const WIN_XI = (3.0, 15.0); const WIN_B = (6.0, 18.0); const WIN_V12 = (5.0, 30.0)
const VM_EDGES = [10.0^l for l in 12.0:0.5:14.5]
const NM_THR = [1.2e12, 3e12, 1e13, 1e14]
const NZ_EDGES = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.5]
const DNDZ_EDGES = collect(0.0:0.25:4.5); const DNDZ_M = [3e12, 1e13]
const TK_Z = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.5]; const TK_M = [1e12, 3e12, 1e13, 1e14, 5e14, 1e15]

# ---- per-cap sample ----
mutable struct Cap
    x::Vector{Float64}; y::Vector{Float64}; z::Vector{Float64}          # shell halos M>MLOW_XI (Eulerian)
    vx::Vector{Float64}; vy::Vector{Float64}; vz::Vector{Float64}; M::Vector{Float64}
    obs::NTuple{3,Float64}
    nm::Matrix{Float64}          # N(>M_thr | z-bin), counts
    dndz::Matrix{Float64}        # counts per Δz bin for each DNDZ_M
    vs::Matrix{Float64}          # σ_vr accumulators per VM bin: n, Σvr, Σvr²
end
Cap(obs) = Cap(Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], obs,
               zeros(length(NM_THR), length(NZ_EDGES) - 1), zeros(length(DNDZ_M), length(DNDZ_EDGES) - 1),
               zeros(3, length(VM_EDGES) - 1))
binof(e, v) = (i = searchsortedlast(e, v); 1 <= i < length(e) ? i : 0)

function addhalo!(c::Cap, X, Y, Z, VX, VY, VZ, M, r, z)
    iz = binof(NZ_EDGES, z)
    iz > 0 && for (j, t) in enumerate(NM_THR); M > t && (c.nm[j, iz] += 1); end
    id = binof(DNDZ_EDGES, z)
    id > 0 && for (j, t) in enumerate(DNDZ_M); M > t && (c.dndz[j, id] += 1); end
    vr = ((X - c.obs[1]) * VX + (Y - c.obs[2]) * VY + (Z - c.obs[3]) * VZ) / r
    iv = binof(VM_EDGES, M)
    iv > 0 && (c.vs[1, iv] += 1; c.vs[2, iv] += vr; c.vs[3, iv] += vr^2)
    if SH[1] <= r <= SH[2] && M > MLOW_XI
        push!(c.x, X); push!(c.y, Y); push!(c.z, Z); push!(c.vx, VX); push!(c.vy, VY); push!(c.vz, VZ); push!(c.M, M)
    end
end

# ---- Websky patch and its inscribed cap (as in compare_bias.jl) ----
function load_wsky_all()
    io = open(WSKY); Nw = read(io, Int32); read(io, Int32); read(io, Int32)
    b = Vector{Float32}(undef, Int(Nw) * 10); read!(io, b); close(io); reshape(b, 10, :)
end
W = load_wsky_all()
wpos = HUB .* Float64.(W[1:3, :]); wr = vec(sqrt.(sum(abs2, wpos; dims=1)))
wax = vec(mean(wpos ./ wr'; dims=2)); wax ./= norm(wax)
cosmin = minimum((wpos' * wax) ./ wr)
const COSW = cos(acos(cosmin) / sqrt(2)); const OMEGA = 2π * (1 - COSW)
@printf("Websky patch: %d halos; inscribed cap half-angle %.3f°, Ω = %.5f sr\n", size(W, 2), acosd(COSW), OMEGA)
function stream(path, f)
    nh = open(io -> Int(read(io, Int32)), path); MAXH > 0 && (nh = min(nh, MAXH))
    open(path) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        ch = 2_000_000; bf = Vector{Float32}(undef, ch * 33); nd = 0
        while nd < nh
            m = min(ch, nh - nd); read!(io, view(bf, 1:m*33))
            @inbounds for k in 1:m; f(bf, (k - 1) * 33); end
            nd += m
        end
    end
end
# ---- pair statistics (chaining mesh; as compare_bias.jl / compare_pairwise_velocity.jl) ----
function pairhist(qx, qy, qz, tx, ty, tz, auto, edges; logbins=true, v=nothing)
    nb = length(edges) - 1; L = edges[end]
    xmin = min(minimum(qx), minimum(tx)); ymin = min(minimum(qy), minimum(ty)); zmin = min(minimum(qz), minimum(tz))
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    for j in eachindex(tx); push!(get!(grid, (floor(Int, (tx[j] - xmin) / L), floor(Int, (ty[j] - ymin) / L), floor(Int, (tz[j] - zmin) / L)), Int[]), j); end
    le = log10(edges[1]); lr = log10(edges[end])
    function worker(i0, i1)
        h = zeros(nb); s = zeros(nb)
        @inbounds for i in i0:i1
            cx = floor(Int, (qx[i] - xmin) / L); cy = floor(Int, (qy[i] - ymin) / L); cz = floor(Int, (qz[i] - zmin) / L)
            for ddx in -1:1, ddy in -1:1, ddz in -1:1
                cc = (cx + ddx, cy + ddy, cz + ddz); haskey(grid, cc) || continue
                for j in grid[cc]
                    auto && j <= i && continue
                    dx = qx[i] - tx[j]; dy = qy[i] - ty[j]; dz = qz[i] - tz[j]; d2 = dx^2 + dy^2 + dz^2
                    (0.0 < d2 < L^2) || continue; d = sqrt(d2)
                    bin = logbins ? floor(Int, (log10(d) - le) / ((lr - le) / nb)) + 1 : floor(Int, (d - edges[1]) / ((edges[end] - edges[1]) / nb)) + 1
                    1 <= bin <= nb || continue
                    h[bin] += 1
                    v === nothing || (s[bin] += ((v[1][i] - v[1][j]) * dx + (v[2][i] - v[2][j]) * dy + (v[3][i] - v[3][j]) * dz) / d)
                end
            end
        end
        (h, s)
    end
    nq = length(qx); nw = Threads.nthreads(); bn = [1 + div(nq * (w - 1), nw) for w in 1:nw+1]; bn[end] = nq + 1
    rs = fetch.([Threads.@spawn worker(bn[w], bn[w+1] - 1) for w in 1:nw])
    (reduce(+, first.(rs)), reduce(+, last.(rs)))
end
function randoms(c::Cap, ax, rng; mult=3)
    t = abs(ax[1]) < 0.9 ? [1.0, 0, 0] : [0, 1.0, 0]; e1 = normalize(cross(ax, t)); e2 = cross(ax, e1)
    rq = [sqrt((c.x[i] - c.obs[1])^2 + (c.y[i] - c.obs[2])^2 + (c.z[i] - c.obs[3])^2) for i in eachindex(c.x)]
    n = mult * length(c.x); x = zeros(n); y = zeros(n); z = zeros(n)
    for k in 1:n
        ct = COSW + rand(rng) * (1 - COSW); st = sqrt(1 - ct^2); ph = 2π * rand(rng); r = rq[rand(rng, eachindex(rq))]
        d = ct .* ax .+ st * cos(ph) .* e1 .+ st * sin(ph) .* e2
        x[k] = c.obs[1] + r * d[1]; y[k] = c.obs[2] + r * d[2]; z[k] = c.obs[3] + r * d[3]
    end
    (x, y, z)
end
function xi_ls(x, y, z, rx, ry, rz, edges)
    Nd = length(x); Nr = length(rx)
    DD, _ = pairhist(x, y, z, x, y, z, true, edges); RR, _ = pairhist(rx, ry, rz, rx, ry, rz, true, edges)
    DR, _ = pairhist(x, y, z, rx, ry, rz, false, edges)
    @. (DD / (Nd * (Nd - 1) / 2) - 2DR / (Nd * Nr) + RR / (Nr * (Nr - 1) / 2)) / max(RR / (Nr * (Nr - 1) / 2), 1e-30)
end
rc(e) = [sqrt(e[i] * e[i+1]) for i in 1:length(e)-1]
winmean(v, e, w) = (c = rc(e); mean(v[[i for i in eachindex(c) if w[1] <= c[i] <= w[2]]]))

function pairstats(c::Cap, ax, seed)
    nanres = (xi=fill(NaN, length(XI_EDGES) - 1), bxi=fill(NaN, length(BM_EDGES) - 1),
              v12=fill(NaN, length(V12_EDGES) - 1), n=length(c.x))
    length(c.x) < 50 && return nanres            # only in smoke tests (partial reads)
    rng = MersenneTwister(seed); rx, ry, rz = randoms(c, ax, rng)
    xi = xi_ls(c.x, c.y, c.z, rx, ry, rz, XI_EDGES)
    bxi = Float64[]
    for ib in 1:length(BM_EDGES)-1
        s = findall(m -> BM_EDGES[ib] <= m < BM_EDGES[ib+1], c.M)
        if length(s) < 50; push!(bxi, NaN); continue; end
        push!(bxi, winmean(xi_ls(c.x[s], c.y[s], c.z[s], rx, ry, rz, XI_EDGES), XI_EDGES, WIN_B))
    end
    s = findall(>(MCUT_V12), c.M)
    length(s) < 50 && return (xi=xi, bxi=bxi, v12=fill(NaN, length(V12_EDGES) - 1), n=length(c.x))
    n, sv = pairhist(c.x[s], c.y[s], c.z[s], c.x[s], c.y[s], c.z[s], true, V12_EDGES; logbins=false, v=(c.vx[s], c.vy[s], c.vz[s]))
    (xi=xi, bxi=bxi, v12=sv ./ max.(n, 1), n=length(c.x))
end
function main()
    tk_raw = zeros(length(TK_M), length(TK_Z) - 1); tk_am = zeros(length(TK_M), length(TK_Z) - 1)
wcap = Cap((0.0, 0.0, 0.0))
for i in axes(W, 2)
    r = wr[i]; (wpos[:, i]' * wax) / r >= COSW || continue
    addhalo!(wcap, wpos[1, i], wpos[2, i], wpos[3, i], Float64(W[4, i]), Float64(W[5, i]), Float64(W[6, i]),
             mass(HUB * Float64(W[7, i])), r, chi_to_z(CHI2Z, r))
end

# ---- stream our catalogs: caps (AM) + full-sky N(>M|z) (raw and AM) ----
caps = Cap[]
for o in OCTS
    ob = obs_of(o); ax = axis_of(o); c = Cap(ob)
    @info "octant $o" ob ax
    stream(joinpath(S3, "catalog_websky_6144_v3_oct$(o)_AMfs.pksc"), (bf, b) -> begin
        R = Float64(bf[b+7]); R > 0 || return
        X = Float64(bf[b+1]); Y = Float64(bf[b+2]); Z = Float64(bf[b+3])
        dx = X - ob[1]; dy = Y - ob[2]; dz = Z - ob[3]; r = sqrt(dx^2 + dy^2 + dz^2)
        M = mass(R); z = chi_to_z(CHI2Z, r)
        it = binof(TK_Z, z); it > 0 && for (j, t) in enumerate(TK_M); M > t && (tk_am[j, it] += 1); end
        (dx * ax[1] + dy * ax[2] + dz * ax[3]) / r >= COSW || return
        addhalo!(c, X, Y, Z, Float64(bf[b+4]), Float64(bf[b+5]), Float64(bf[b+6]), M, r, z)
    end)
    stream(joinpath(S3, "catalog_websky_6144_v3_oct$(o).pksc"), (bf, b) -> begin
        R = Float64(bf[b+7]); R > 0 || return
        r = sqrt((bf[b+1] - ob[1])^2 + (bf[b+2] - ob[2])^2 + (bf[b+3] - ob[3])^2)
        z = chi_to_z(CHI2Z, r); M = mass(R)
        it = binof(TK_Z, z); it > 0 && for (j, t) in enumerate(TK_M); M > t && (tk_raw[j, it] += 1); end
    end)
    push!(caps, c)
end

@info "pair statistics: Websky"; PW = pairstats(wcap, wax, 1)
PO = [(@info "pair statistics: oct $(OCTS[i])"; pairstats(caps[i], axis_of(OCTS[i]), 1 + i)) for i in eachindex(caps)]

# ---- report ----
io = open(OUT, "w")
function line(label, w, os)
    ok = filter(isfinite, os); m = mean(ok); sd = length(ok) > 1 ? std(ok) : NaN
    s = @sprintf("%-34s Websky %11.4g | ours %11.4g ± %-10.3g (8 caps) | W/ours %6.3f ± %5.3f | z %+5.1f",
                 label, w, m, sd, w / m, sd / m, (w - m) / sd)
    println(s); println(io, s)
end
hdr(t) = (println("\n== ", t); println(io, "\n== ", t))
hdr("header: Websky patch vs our 8 octant-centre caps (same inscribed cap, Ω=$(round(OMEGA; digits=5)) sr); v3fs catalogs")
hdr("N(>M) per sr in z-bins (after AM; largely enforced by AM)")
for (j, t) in enumerate(NM_THR), iz in 1:length(NZ_EDGES)-1
    line(@sprintf("N(>%.1e) z %.1f-%.1f", t, NZ_EDGES[iz], NZ_EDGES[iz+1]), wcap.nm[j, iz] / OMEGA, [c.nm[j, iz] / OMEGA for c in caps])
end
hdr("dN/dz per sr (Δz = 0.25)")
for (j, t) in enumerate(DNDZ_M), id in 1:length(DNDZ_EDGES)-1
    line(@sprintf("dN/dz M>%.0e z %.2f-%.2f", t, DNDZ_EDGES[id], DNDZ_EDGES[id+1]), wcap.dndz[j, id] / OMEGA, [c.dndz[j, id] / OMEGA for c in caps])
end
sig(v) = v[1] > 1 ? sqrt(v[3] / v[1] - (v[2] / v[1])^2) : NaN
hdr("σ_vr (km/s) per mass bin, whole cap (all z)")
for iv in 1:length(VM_EDGES)-1
    line(@sprintf("σ_vr M %.1e-%.1e", VM_EDGES[iv], VM_EDGES[iv+1]), sig(wcap.vs[:, iv]), [sig(c.vs[:, iv]) for c in caps])
end
hdr(@sprintf("ξ(r), shell %.0f-%.0f Mpc/h, M>%.0e (Landy-Szalay; N_wsky=%d, N_ours≈%d)", SH..., MLOW_XI, PW.n, round(Int, mean(p.n for p in PO))))
for (i, r) in enumerate(rc(XI_EDGES))
    line(@sprintf("ξ r=%.2f", r), PW.xi[i], [p.xi[i] for p in PO])
end
line(@sprintf("ξ mean over %.0f-%.0f Mpc/h", WIN_XI...), winmean(PW.xi, XI_EDGES, WIN_XI), [winmean(p.xi, XI_EDGES, WIN_XI) for p in PO])
hdr(@sprintf("relative bias: sqrt(<ξ_hh(M)> over %.0f-%.0f Mpc/h)", WIN_B...))
for ib in 1:length(BM_EDGES)-1
    line(@sprintf("b ∝ √ξ, M %.1e-%.1e", BM_EDGES[ib], BM_EDGES[ib+1]), sqrt(max(PW.bxi[ib], 0)), [sqrt(max(p.bxi[ib], 0)) for p in PO])
end
hdr(@sprintf("v12(r) km/s, M>%.0e, shell %.0f-%.0f Mpc/h", MCUT_V12, SH...))
v12c = [(V12_EDGES[i] + V12_EDGES[i+1]) / 2 for i in 1:length(V12_EDGES)-1]
for (i, r) in enumerate(v12c)
    line(@sprintf("v12 r=%.1f", r), PW.v12[i], [p.v12[i] for p in PO])
end
wm(v) = mean(v[[i for i in eachindex(v12c) if WIN_V12[1] <= v12c[i] <= WIN_V12[2]]])
line(@sprintf("v12 mean over %.0f-%.0f Mpc/h", WIN_V12...), wm(PW.v12), [wm(p.v12) for p in PO])

# ---- A4: full-sky raw / AM N(>M|z) vs Tinker08 ----
pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky_RAW_unnormalized.dat"))
Mg = 10 .^ range(11.5, 16; length=901); lnMg = log.(Mg); sg = precompute_sigma(Mg, pk, 0.31)
function ntinker(Mcut, za, zb)       # full sky
    ra, rb = chi(za, COSMO), chi(zb, COSMO); ns = 40; acc = 0.0
    for j in 1:ns
        r0 = ra + (j - 1) * (rb - ra) / ns; r1 = ra + j * (rb - ra) / ns
        z = chi_to_z(CHI2Z, (r0 + r1) / 2); Dg = growth_factor(z, COSMO); n = 0.0
        for i in 1:length(Mg)-1
            Mm = sqrt(Mg[i] * Mg[i+1]); Mm < Mcut && continue
            dls = (log(sg[i+1]) - log(sg[i])) / (lnMg[i+1] - lnMg[i])
            n += tinker_dndlnM(Mm, Dg * sqrt(sg[i] * sg[i+1]), dls, z, 0.31) * (lnMg[i+1] - lnMg[i])
        end
        acc += n * (4π / 3) * (r1^3 - r0^3)
    end
    acc
end
hdr("A4: full-sky N(>M | z) — RAW (pre-AM) and AM vs Tinker08 M200m (8 octants)")
s = @sprintf("%-26s %12s %12s %12s %8s %8s", "bin", "raw", "AM", "Tinker", "raw/Tk", "AM/Tk"); println(s); println(io, s)
for (j, t) in enumerate(TK_M), iz in 1:length(TK_Z)-1
    e = ntinker(t, TK_Z[iz], TK_Z[iz+1])
    s = @sprintf("M>%.0e z %.2f-%.2f   %12d %12d %12.1f %8.3f %8.3f", t, TK_Z[iz], TK_Z[iz+1], tk_raw[j, iz], tk_am[j, iz], e, tk_raw[j, iz] / e, tk_am[j, iz] / e)
    println(s); println(io, s)
end
close(io); @info "wrote" OUT
end
if abspath(PROGRAM_FILE) == (@__FILE__)   # include()-able for the helpers (tierA_diag.jl)
    main()
end
