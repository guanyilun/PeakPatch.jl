#!/usr/bin/env julia
# Halo velocity / displacement error of the multires split vs the exact global field, same δ
# (MATCHED_FORTRAN_2026-10.md §9). Halo velocities are v = a H f (ψ1 + 2ψ2) (merge_pkvd convention), so the
# 1LPT and 2LPT displacements are compared separately and as v ∝ ψ1 + 2ψ2, for raw peaks matched by
# Lagrangian cell + filter:
#   * per-halo: rms |Δ| / rms |ref|, correlation
#   * halo velocity field (x component, mass-unweighted mean per 4-cell, on matched halos): error power and
#     cross-correlation r(k) vs k — the large-scale coherent part that kSZ sees
#   usage: julia --project=validation -t 32 split_velocity.jl <rundir with fields/> <exact dir> <split dir>... (labels from dir names)
include(joinpath(@__DIR__, "matched_compare.jl"))

function vel_main()
    exd = ARGS[2]; spds = ARGS[3:end]
    rd(p) = read_pksc(p)[1]
    E = filter(q -> q.RTHL > 0 && inreg(q.x, q.y, q.z), rd(joinpath(exd, "julia_exact_raw.pksc")))
    out = open(joinpath(@__DIR__, "..", "results", "split_velocity.txt"), "w")
    say(a...) = begin s = string(a...); println(s); println(out, s) end
    say("halo displacement / velocity error, split vs exact (same δ), raw peaks matched by cell + filter, M > 5e12 (raw R)")
    Ec = cat_lagr(E)
    for spd in spds
        S = filter(q -> q.RTHL > 0 && inreg(q.x, q.y, q.z), rd(joinpath(spd, "julia_split_raw.pksc"))); Sc = cat_lagr(S)
        m = match(Ec, Sc; tol=0.5A)
        ok = findall(i -> m[i] > 0 && abs(Ec.Rf[i] - Sc.Rf[m[i]]) < 1e-3 && mass(Ec.R[i]) > 5e12, eachindex(m))
        say("\n== ", basename(spd), ": matched $(length(ok)) of $(count(i -> mass(Ec.R[i]) > 5e12, eachindex(Ec.R)))")
        for (lab, f) in (("ψ1", q -> (q.vx, q.vy, q.vz)), ("ψ2", q -> (q.vx2, q.vy2, q.vz2)),
                         ("v∝ψ1+2ψ2", q -> (q.vx + 2q.vx2, q.vy + 2q.vy2, q.vz + 2q.vz2)))
            r2 = 0.0; e2 = 0.0; x = 0.0; s2 = 0.0
            for i in ok
                a = f(E[i]); b = f(S[m[i]])
                for d in 1:3
                    r2 += a[d]^2; s2 += b[d]^2; e2 += (b[d] - a[d])^2; x += a[d] * b[d]
                end
            end
            say(@sprintf("  %-9s rms|Δ|/rms|ref| %.4f   corr %.5f   amplitude Σab/Σa² %.4f   rms|ref| %.3f Mpc/h", lab, sqrt(e2 / r2), x / sqrt(r2 * s2), x / r2, sqrt(r2 / length(ok))))
        end
        # large-scale velocity field (x), mean per 4-cell over matched halos, exact vs split on the same cells
        vr = zeros(NC, NC, NC); vs = zeros(NC, NC, NC); nn = zeros(NC, NC, NC); h = CB * A
        for i in ok
            q = E[i]; I = clamp(floor(Int, (q.x - XLO) / h) + 1, 1, NC); J = clamp(floor(Int, (q.y - XLO) / h) + 1, 1, NC)
            K = clamp(floor(Int, (q.z - XLO) / h) + 1, 1, NC)
            vr[I, J, K] += q.vx + 2q.vx2; p = S[m[i]]; vs[I, J, K] += p.vx + 2p.vx2; nn[I, J, K] += 1
        end
        R = rfft(vr); D = rfft(vs .- vr); kf = 2π / (NC * h)
        ed = [0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]; Pr = zeros(6); Pe = zeros(6)
        for k3 in 1:NC, k2 in 1:NC, k1 in 1:NC÷2+1
            q2 = k2 <= NC ÷ 2 + 1 ? k2 - 1 : k2 - 1 - NC; q3 = k3 <= NC ÷ 2 + 1 ? k3 - 1 : k3 - 1 - NC
            kk = kf * sqrt((k1 - 1)^2 + q2^2 + q3^2); b = searchsortedlast(ed, kk); (1 <= b <= 6 && kk > 0) || continue
            Pr[b] += abs2(R[k1, k2, k3]); Pe[b] += abs2(D[k1, k2, k3])
        end
        say("  halo velocity field (x, summed per cell over matched halos): error power / reference power by k [h/Mpc]")
        say("    " * join([@sprintf("%.2f–%.2f: %.2e   ", ed[b], ed[b+1], Pe[b] / Pr[b]) for b in 1:6]))
    end
    close(out)
end
vel_main()
