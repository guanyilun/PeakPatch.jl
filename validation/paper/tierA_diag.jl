#!/usr/bin/env julia
# Why is v3fs ~7% more biased than Websky at fixed number density? (TIERA_V3FS_2026-09-28.md)
# Same caps/shell as tierA_v3fs.jl (Websky patch vs the identical cap at each octant centre).
#  D1 absolute bias b(M) = sqrt(ξ_hh / ξ_mm,lin) over 6–18 Mpc/h, both catalogs, vs Tinker10:
#     which side departs from theory?
#  D2 close-pair census: pairs per halo with u = d / (R_i + R_j) in bins (Eulerian d, AM R_TH):
#     does ours keep more overlapping neighbours (weaker exclusion)?
#  D3 ξ after one common, post-hoc Eulerian exclusion applied to BOTH catalogs (largest first,
#     drop a halo whose centre lies within f·R of a kept larger halo): does the ratio → 1?
#  D4 ours at fixed number density by RAW mass rank vs by AM mass: does AM change the selection?
# Raw and AM files are record-aligned (abundance_match keeps order), so both masses are known.
include(joinpath(@__DIR__, "tierA_v3fs.jl"))      # helpers only (main() is guarded)
import PeakPatch.Cosmology: Dlinear_ab, Dlinear_tables
using DelimitedFiles

# ---- linear theory at the shell (as compare_bias.jl) ----
pkd = readdlm(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky.dat"); comments=true, comment_char='#')
const KK = Float64.(pkd[:, 1]); const PK = Float64.(pkd[:, 2]) .* (2π)^3
Wth(x) = x < 1e-3 ? 1.0 - x^2 / 10 : 3 * (sin(x) - x * cos(x)) / x^3
spec(f) = sum(0.5 * (f(i) + f(i + 1)) * (KK[i+1] - KK[i]) for i in 1:length(KK)-1)
sigma2(R) = spec(i -> PK[i] * Wth(KK[i] * R)^2 * KK[i]^2) / (2π^2)
xi_mm0(r) = spec(i -> (kr = KK[i] * r; PK[i] * KK[i]^2 * (kr < 1e-3 ? 1.0 : sin(kr) / kr))) / (2π^2)
yT = log10(200.0); AT = 1 + 0.24yT * exp(-(4 / yT)^4); aT = 0.44yT - 0.88; BT = 0.183; bT = 1.5
CT = 0.019 + 0.107yT + 0.19exp(-(4 / yT)^4); cT = 2.4
tinker_b(M, Dz) = (ν = 1.686 / (Dz * sqrt(sigma2(cbrt(3M / (4π * RHO_M))))); 1 - AT * ν^aT / (ν^aT + 1.686^aT) + BT * ν^bT + CT * ν^cT)

const DB_EDGES = [5e12, 8e12, 1.3e13, 2e13, 3.2e13, 7.9e13, 2.5e14]
struct Samp; x::Vector{Float64}; y::Vector{Float64}; z::Vector{Float64}; M::Vector{Float64}; Mraw::Vector{Float64}; obs::NTuple{3,Float64}; end

function shell_ours(o)
    ob = obs_of(o); ax = axis_of(o)
    x = Float64[]; y = Float64[]; z = Float64[]; M = Float64[]; Mr = Float64[]
    fa = joinpath(S3, "catalog_websky_6144_v3_oct$(o)_AMfs.pksc"); fr = joinpath(S3, "catalog_websky_6144_v3_oct$(o).pksc")
    nh = open(io -> Int(read(io, Int32)), fa); MAXH > 0 && (nh = min(nh, MAXH))
    ia = open(fa); ir = open(fr); for f in (ia, ir); read(f, Int32); read(f, Float32); read(f, Float32); end
    ch = 2_000_000; ba = Vector{Float32}(undef, ch * 33); br = similar(ba); nd = 0
    while nd < nh
        m = min(ch, nh - nd); read!(ia, view(ba, 1:m*33)); read!(ir, view(br, 1:m*33))
        @inbounds for k in 1:m
            b = (k - 1) * 33
            X = Float64(ba[b+1]); Y = Float64(ba[b+2]); Z = Float64(ba[b+3])
            dx = X - ob[1]; dy = Y - ob[2]; dz = Z - ob[3]; r = sqrt(dx^2 + dy^2 + dz^2)
            SH[1] <= r <= SH[2] || continue
            (dx * ax[1] + dy * ax[2] + dz * ax[3]) / r >= COSW || continue
            Ma = mass(Float64(ba[b+7])); Mw = mass(Float64(br[b+7]))
            max(Ma, Mw) > 2e12 || continue
            push!(x, X); push!(y, Y); push!(z, Z); push!(M, Ma); push!(Mr, Mw)
        end
        nd += m
    end
    close(ia); close(ir); Samp(x, y, z, M, Mr, ob)
end
function shell_wsky()
    Wt = load_wsky_all(); x = Float64[]; y = Float64[]; z = Float64[]; M = Float64[]
    for i in axes(Wt, 2)
        X, Y, Z = HUB .* Float64.(Wt[1:3, i]); r = sqrt(X^2 + Y^2 + Z^2)
        SH[1] <= r <= SH[2] || continue
        (X * wax[1] + Y * wax[2] + Z * wax[3]) / r >= COSW || continue
        Mh = mass(HUB * Float64(Wt[7, i])); Mh > 2e12 || continue
        push!(x, X); push!(y, Y); push!(z, Z); push!(M, Mh)
    end
    Samp(x, y, z, M, fill(NaN, length(M)), (0.0, 0.0, 0.0))
end
sub(s::Samp, idx) = Samp(s.x[idx], s.y[idx], s.z[idx], s.M[idx], s.Mraw[idx], s.obs)
Rof(M) = cbrt(3M / (4π * RHO_M))

# ξ for a sample with randoms made from its own radii (tierA helpers expect a Cap-like object)
function xi_samp(s::Samp, ax, seed; edges=XI_EDGES)
    c = Cap(s.obs); append!(c.x, s.x); append!(c.y, s.y); append!(c.z, s.z)
    rx, ry, rz = randoms(c, ax, MersenneTwister(seed))
    xi_ls(s.x, s.y, s.z, rx, ry, rz, edges), (rx, ry, rz)
end
# D2: pairs per halo in u = d/(R_i+R_j) bins
const U_EDGES = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
function upairs(s::Samp)
    R = Rof.(s.M); L = 3 * 2 * maximum(R); h = zeros(length(U_EDGES) - 1)
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    for j in eachindex(s.x); push!(get!(grid, (floor(Int, s.x[j] / L), floor(Int, s.y[j] / L), floor(Int, s.z[j] / L)), Int[]), j); end
    for i in eachindex(s.x)
        cx = floor(Int, s.x[i] / L); cy = floor(Int, s.y[i] / L); cz = floor(Int, s.z[i] / L)
        for a in -1:1, b in -1:1, c in -1:1
            haskey(grid, (cx + a, cy + b, cz + c)) || continue
            for j in grid[(cx+a, cy+b, cz+c)]
                j <= i && continue
                d = sqrt((s.x[i] - s.x[j])^2 + (s.y[i] - s.y[j])^2 + (s.z[i] - s.z[j])^2); u = d / (R[i] + R[j])
                k = searchsortedlast(U_EDGES, u); 1 <= k < length(U_EDGES) && (h[k] += 1)
            end
        end
    end
    h ./ length(s.x)
end
# D3: common post-hoc Eulerian exclusion (largest first; drop centre within f·R_kept)
function exclude(s::Samp, f)
    o = sortperm(s.M; rev=true); R = Rof.(s.M); L = f * maximum(R) + 1e-9; keep = falses(length(o))
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    for i in o
        cx = floor(Int, s.x[i] / L); cy = floor(Int, s.y[i] / L); cz = floor(Int, s.z[i] / L); ok = true
        for a in -1:1, b in -1:1, c in -1:1
            haskey(grid, (cx + a, cy + b, cz + c)) || continue
            for j in grid[(cx+a, cy+b, cz+c)]
                ((s.x[i] - s.x[j])^2 + (s.y[i] - s.y[j])^2 + (s.z[i] - s.z[j])^2 < (f * R[j])^2) && (ok = false; break)
            end
            ok || break
        end
        ok && (keep[i] = true; push!(get!(grid, (cx, cy, cz), Int[]), i))
    end
    findall(keep)
end

function run()
    zsh = chi_to_z(CHI2Z, 0.5 * (SH[1] + SH[2])); gt = Dlinear_tables(COSMO); Dz, _, _ = Dlinear_ab(1 / (1 + zsh), gt)
    rcs = rc(XI_EDGES); wsel = [i for i in eachindex(rcs) if WIN_B[1] <= rcs[i] <= WIN_B[2]]
    ximm = mean(Dz^2 * xi_mm0(rcs[i]) for i in wsel)
    io = open(joinpath(@__DIR__, "results", MAXH > 0 ? "/tmp/tierA_diag_smoke.txt" : "tierA_diag_v3fs.txt"), "w")
    say(a...) = (s = string(a...); println(s); println(io, s))
    say(@sprintf("shell z=%.3f D=%.4f  <ξ_mm,lin>(6–18) = %.4f  σ8 check %.4f", zsh, Dz, ximm, sqrt(sigma2(8.0))))
    Ws = shell_wsky(); Os = [shell_ours(o) for o in OCTS]
    stat(v) = begin ok = filter(isfinite, v); (mean(ok), std(ok) * sqrt(1 + 1 / length(ok))) end
    fmt(w, os) = begin m, sd = stat(os); @sprintf("W %8.4g | ours %8.4g ± %-8.3g | W/ours %.3f  z %+5.1f", w, m, sd, w / m, (w - m) / sd) end
    # D1
    say("\n== D1 absolute bias b = sqrt(<ξ_hh>/<ξ_mm,lin>) over 6–18 Mpc/h, and Tinker10 b(M)")
    for i in 1:length(DB_EDGES)-1
        lo, hi = DB_EDGES[i], DB_EDGES[i+1]; Mc = sqrt(lo * hi)
        function bias_of(smp, ax, seed)
            ss = sub(smp, findall(m -> lo <= m < hi, smp.M))
            length(ss.x) < 50 && return NaN
            sqrt(max(winmean(first(xi_samp(ss, ax, seed)), XI_EDGES, WIN_B), 0) / ximm)
        end
        bw = bias_of(Ws, wax, 1)
        bo = [bias_of(Os[k], axis_of(OCTS[k]), 1 + k) for k in eachindex(Os)]
        bt = tinker_b(Mc, Dz)
        say(@sprintf("M %.1e-%.1e  n_W %6d  ", lo, hi, count(m -> lo <= m < hi, Ws.M)), fmt(bw, bo), @sprintf("  | Tinker10 %.3f  W/Tk %.3f  ours/Tk %.3f", bt, bw / bt, mean(filter(isfinite, bo)) / bt))
    end
    # D2
    say("\n== D2 close pairs per halo, M>5e12 (AM masses), in u = d/(R_i+R_j)")
    W5 = sub(Ws, findall(>(MLOW_XI), Ws.M)); O5 = [sub(s, findall(>(MLOW_XI), s.M)) for s in Os]
    hw = upairs(W5); ho = [upairs(s) for s in O5]
    for k in 1:length(U_EDGES)-1
        say(@sprintf("u %.2f-%.2f  ", U_EDGES[k], U_EDGES[k+1]), fmt(hw[k], [h[k] for h in ho]))
    end
    # D3
    say("\n== D3 ξ(3–15 Mpc/h) and ξ(1–3), M>5e12, after a common post-hoc Eulerian exclusion (centre within f·R of a larger kept halo)")
    for f in (0.0, 1.0, 2.0)
        sel(s) = f == 0 ? s : sub(s, exclude(s, f))
        ws = sel(W5); os = [sel(s) for s in O5]
        xw = first(xi_samp(ws, wax, 1)); xo = [first(xi_samp(os[k], axis_of(OCTS[k]), 1 + k)) for k in eachindex(os)]
        s13(x) = mean(x[[i for i in eachindex(rcs) if 1 <= rcs[i] <= 3]])
        say(@sprintf("f=%.0f  kept W %.4f ours %.4f  ", f, length(ws.x) / length(W5.x), mean(length(os[k].x) / length(O5[k].x) for k in eachindex(os))),
            "ξ(3–15): ", fmt(winmean(xw, XI_EDGES, WIN_XI), [winmean(x, XI_EDGES, WIN_XI) for x in xo]), "   ξ(1–3): ", fmt(s13(xw), [s13(x) for x in xo]))
    end
    # D4
    say("\n== D4 ours: same number of halos selected by RAW-mass rank vs by AM mass (M_AM > 5e12)")
    for k in eachindex(Os)
        s = Os[k]; n = count(>(MLOW_XI), s.M)
        ia = findall(>(MLOW_XI), s.M); ir = sortperm(s.Mraw; rev=true)[1:n]
        xa = first(xi_samp(sub(s, ia), axis_of(OCTS[k]), 100 + k)); xr = first(xi_samp(sub(s, ir), axis_of(OCTS[k]), 100 + k))
        say(@sprintf("oct %s  N %6d  overlap of selections %.4f  ξ(3–15) AM %.4f raw-rank %.4f  ratio %.4f", OCTS[k], n,
                     length(intersect(ia, ir)) / n, winmean(xa, XI_EDGES, WIN_XI), winmean(xr, XI_EDGES, WIN_XI),
                     winmean(xr, XI_EDGES, WIN_XI) / winmean(xa, XI_EDGES, WIN_XI)))
    end
    close(io)
end
run()
