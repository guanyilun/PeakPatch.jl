#!/usr/bin/env julia
# A/B test of the lightcone peak-candidate threshold (CLUSTERING_EXCESS_2026-09-28.md, open item):
#   A = constant fcrit = fsc_of_z(z_out = 0) everywhere (every catalog up to v3fs)
#   B = Fortran hpkvd: fcrit = fsc_of_z(z_tile), z_tile from the tile centre's distance
#       (hpkvd.f90:510-514 -> get_pks; [run] peak_threshold_per_tile = true)
# Both merged WITH the Fortran volume reduction (the production candidate). Lightcone box:
# v3 physics, N = 1536 (L = 1309 Mpc/h), observer at the box corner, z_max = 1.0. Compared:
# raw and merged counts per z bin, N(>M|z), and ξ(r) at fixed number density in the Tier-A
# shell 1600-2000 Mpc/h (z≈0.6-0.8), with randoms uniform in the same box ∩ shell region.
#   usage: julia --project=validation -t 16 validation/paper/threshold_ab.jl <config.toml>
using CUDA
using TOML, Random
include(joinpath(@__DIR__, "tierA_v3fs.jl"))          # helpers: pairhist, xi_ls, mass, rc, winmean

function region_randoms(n, obs, L, rng)
    x = Float64[]; y = Float64[]; z = Float64[]
    lo = obs .+ 0.0                                    # observer at the min corner of the box
    while length(x) < n
        p = (lo[1] + L * rand(rng), lo[2] + L * rand(rng), lo[3] + L * rand(rng))
        r = sqrt(sum(abs2, p .- obs)); SH[1] <= r <= SH[2] || continue
        push!(x, p[1]); push!(y, p[2]); push!(z, p[3])
    end
    x, y, z
end

function tmain()
    cfgd = TOML.parsefile(ARGS[1]); rcfg = cfgd["run"]; cosmo = COSMO
    zb = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    out = open(joinpath(@__DIR__, "results", "threshold_ab.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    res = Dict{String,Any}()
    for (lab, pt) in (("A_const", false), ("B_pertile", true))
        cfgd["run"]["peak_threshold_per_tile"] = pt
        cfg = PipelineConfig(cfgd); ntile = rcfg["ntile"]; nsub, N = grid_layout(cfg, ntile); L = N * cfg.boxsize / cfg.n
        obs = (Float64(cfg.cenx), Float64(cfg.ceny), Float64(cfg.cenz))
        raw = run_multitile_split(cfg; ntile=ntile, seed=rcfg["seed"], coarse_factor=rcfg["coarse_factor"], use_gpu=true,
                                  devices=collect(0:length(CUDA.devices())-1), verbose=true)
        m = merge_catalog(raw; verbose=true, volume_reduction=true)
        h = finalize_eulerian(m, cosmo, obs; ievol=1, z_out=0.0)
        zr(q) = chi_to_z(CHI2Z, sqrt((q.x - obs[1])^2 + (q.y - obs[2])^2 + (q.z - obs[3])^2))
        res[lab] = (raw=raw, h=h, zraw=[zr(q) for q in raw], zh=[zr(q) for q in h], obs=obs, L=L)
        say(@sprintf("[%s] raw %d  merged %d  (L=%.1f, obs=%s)", lab, length(raw), length(h), L, obs))
    end
    A = res["A_const"]; B = res["B_pertile"]
    say("\n== counts per z bin: raw peaks-with-halo / merged; B/A")
    for i in 1:length(zb)-1
        nz(v) = count(z -> zb[i] <= z < zb[i+1], v)
        say(@sprintf("z %.1f-%.1f  raw A %9d B %9d (B/A %.4f)   merged A %9d B %9d (B/A %.4f)", zb[i], zb[i+1],
                     nz(A.zraw), nz(B.zraw), nz(B.zraw) / nz(A.zraw), nz(A.zh), nz(B.zh), nz(B.zh) / nz(A.zh)))
    end
    say("\n== merged N(>M | z), B/A")
    for t in (1e12, 3e12, 1e13, 1e14), i in 1:length(zb)-1
        nm(H, Z) = count(k -> zb[i] <= Z[k] < zb[i+1] && mass(H[k].RTHL) > t, eachindex(H))
        a = nm(A.h, A.zh); b = nm(B.h, B.zh)
        say(@sprintf("M>%.0e z %.1f-%.1f   A %9d  B %9d  B/A %.4f", t, zb[i], zb[i+1], a, b, b / max(a, 1)))
    end
    # ξ at fixed number density in the Tier-A shell
    shell(R) = [k for k in eachindex(R.h) if SH[1] <= sqrt(sum(abs2, (R.h[k].x, R.h[k].y, R.h[k].z) .- R.obs)) <= SH[2]]
    sa = shell(A); sb = shell(B)
    nsel = count(k -> mass(A.h[k].RTHL) > MLOW_XI, sa)
    top(R, s) = s[sortperm([mass(R.h[k].RTHL) for k in s]; rev=true)[1:min(nsel, length(s))]]
    ta = top(A, sa); tb = top(B, sb)
    rx, ry, rz = region_randoms(3nsel, A.obs, A.L, MersenneTwister(5))
    xy(R, t) = ([Float64(R.h[k].x) for k in t], [Float64(R.h[k].y) for k in t], [Float64(R.h[k].z) for k in t])
    xa = xi_ls(xy(A, ta)..., rx, ry, rz, XI_EDGES); xb = xi_ls(xy(B, tb)..., rx, ry, rz, XI_EDGES)
    say(@sprintf("\n== shell %.0f-%.0f Mpc/h, top %d halos by mass in each (A's count above %.0e)", SH..., nsel, MLOW_XI))
    say(@sprintf("mass at the cut: A %.3e  B %.3e", mass(A.h[ta[end]].RTHL), mass(B.h[tb[end]].RTHL)))
    for (i, r) in enumerate(rc(XI_EDGES)); say(@sprintf("ξ r=%6.2f  A %.4f  B %.4f  B/A %.4f", r, xa[i], xb[i], xb[i] / xa[i])); end
    w13 = [i for (i, r) in enumerate(rc(XI_EDGES)) if 1 <= r <= 3]
    say(@sprintf("ξ(3–15) B/A %.4f   ξ(1–3) B/A %.4f   [residual to explain after reduction: Websky/ours ÷ reduction A/B = 0.877/0.916 = 0.957 ± 0.03]",
                 winmean(xb, XI_EDGES, WIN_XI) / winmean(xa, XI_EDGES, WIN_XI), mean(xb[w13]) / mean(xa[w13])))
    close(out)
end
tmain()
