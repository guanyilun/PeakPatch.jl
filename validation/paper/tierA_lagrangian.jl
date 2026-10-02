#!/usr/bin/env julia
# Split the full-sky clustering offset (TIERA_FULLSKY_2026-10-02.md: Websky/ours ξ(3–15) = 0.948 ± 0.002)
# into SELECTION (which halos, at which Lagrangian positions) vs DISPLACEMENT (how far they move).
# Same 305 caps, shell (1600–2000 Mpc/h by Eulerian distance) and M > 5e12 sample as tierA_fullsky.jl.
# ξ(r) is measured with five position sets:
#   W_E  Websky Eulerian (cols 1–3)            W_Lt Websky stored Lagrangian xL (cols 8–10, merge_pkvd)
#   W_Lr Websky reconstructed  x − v/(a·100·E·f)  (= q − ψ2 for the merge_pkvd velocity convention)
#   O_E  ours Eulerian                          O_Lr ours reconstructed, same formula (our finalize uses
#                                                    v = a·100·E·f·(ψ1+2ψ2) at the Eulerian redshift)
# W_Lt vs W_Lr measures the reconstruction error; W_Lr vs O_Lr is the like-for-like Lagrangian comparison.
include(joinpath(@__DIR__, "tierA_fullsky.jl"))       # caps, capof, helpers; fullsky_main() guarded
import PeakPatch.Cosmology: Dlinear_ab, Dlinear_tables
const GT = Dlinear_tables(COSMO)
const OUTL = MAXH > 0 ? "/tmp/tierA_lagr_smoke.txt" : joinpath(@__DIR__, "results", "tierA_lagrangian_$(CAMP)fs.txt")
@inline function vfac(z)
    a = 1 / (1 + z); _, f, _ = Dlinear_ab(a, GT); a * 100 * sqrt(COSMO.Om * a^-3 + COSMO.OL) * f
end
mutable struct LS                                  # shell sample of one cap, several position sets
    P::Dict{Symbol,NTuple{3,Vector{Float64}}}; M::Vector{Float64}; obs::NTuple{3,Float64}
end
LS(obs, keys) = LS(Dict(k => (Float64[], Float64[], Float64[]) for k in keys), Float64[], obs)
push3!(t, x, y, z) = (push!(t[1], x); push!(t[2], y); push!(t[3], z))

function lagr_main()
    nc = length(CENT)
    LW = [LS((0.0, 0.0, 0.0), (:E, :Lt, :Lr)) for _ in 1:nc]
    LO = [LS(obs_of(octof(c)), (:E, :Lr)) for c in CENT]
    dW = zeros(3); nW = 0; dWr = zeros(2)              # |x − xL| and reconstruction-error statistics
    nh = open(io -> Int(read(io, Int32)), WFULL); MAXH > 0 && (nh = min(nh, MAXH))
    open(WFULL) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        ch = 4_000_000; b = Vector{Float32}(undef, ch * 10); nd = 0
        while nd < nh
            m = min(ch, nh - nd); read!(io, view(b, 1:m*10))
            @inbounds for k in 1:m
                o = (k - 1) * 10
                X = HUB * Float64(b[o+1]); Y = HUB * Float64(b[o+2]); Z = HUB * Float64(b[o+3])
                r = sqrt(X^2 + Y^2 + Z^2); (SH[1] <= r <= SH[2]) || continue
                M = mass(HUB * Float64(b[o+7])); M > MLOW_XI || continue
                kc = capof(X / r, Y / r, Z / r); kc > 0 || continue
                z = chi_to_z(CHI2Z, r); vf = vfac(z)
                L = LW[kc]; push!(L.M, M)
                push3!(L.P[:E], X, Y, Z)
                xl = HUB * Float64(b[o+8]); yl = HUB * Float64(b[o+9]); zl = HUB * Float64(b[o+10])
                push3!(L.P[:Lt], xl, yl, zl)
                xr = X - Float64(b[o+4]) / vf; yr = Y - Float64(b[o+5]) / vf; zr = Z - Float64(b[o+6]) / vf
                push3!(L.P[:Lr], xr, yr, zr)
                dW[1] += sqrt((X - xl)^2 + (Y - yl)^2 + (Z - zl)^2); dW[2] += sqrt((xr - xl)^2 + (yr - yl)^2 + (zr - zl)^2)
                dW[3] += sqrt((X - xr)^2 + (Y - yr)^2 + (Z - zr)^2); nW += 1
            end
            nd += m
        end
    end
    dO = 0.0; nO = 0
    for oc in OCTS
        ob = obs_of(oc)
        stream(joinpath(S3, "catalog_websky_6144_$(CAMP)_oct$(oc)_AMfs.pksc"), (bf, b) -> begin
            R = Float64(bf[b+7]); R > 0 || return
            M = mass(R); M > MLOW_XI || return
            X = Float64(bf[b+1]); Y = Float64(bf[b+2]); Z = Float64(bf[b+3])
            dx = X - ob[1]; dy = Y - ob[2]; dz = Z - ob[3]; r = sqrt(dx^2 + dy^2 + dz^2)
            (SH[1] <= r <= SH[2]) || return
            kc = capof(dx / r, dy / r, dz / r); (kc > 0 && octof(CENT[kc]) == oc) || return
            vf = vfac(chi_to_z(CHI2Z, r)); L = LO[kc]; push!(L.M, M)
            push3!(L.P[:E], X, Y, Z)
            xr = X - Float64(bf[b+4]) / vf; yr = Y - Float64(bf[b+5]) / vf; zr = Z - Float64(bf[b+6]) / vf
            push3!(L.P[:Lr], xr, yr, zr); dO += sqrt((X - xr)^2 + (Y - yr)^2 + (Z - zr)^2); nO += 1
        end)
    end
    # ξ per cap and position set (randoms from the Eulerian sample geometry, shared across sets)
    sets = [(:W, :E), (:W, :Lt), (:W, :Lr), (:O, :E), (:O, :Lr)]
    X3 = Dict(s => zeros(nc) for s in sets); XI = Dict(s => [zeros(length(XI_EDGES) - 1) for _ in 1:nc] for s in sets)
    for k in 1:nc
        ax = collect(CENT[k])
        for (side, L) in ((:W, LW[k]), (:O, LO[k]))
            length(L.M) < 50 && continue
            c = Cap(L.obs); E = L.P[:E]; append!(c.x, E[1]); append!(c.y, E[2]); append!(c.z, E[3])
            rx, ry, rz = randoms(c, ax, MersenneTwister(10k + (side == :W ? 1 : 2)))
            for (s2, key) in sets
                s2 == side || continue
                P = L.P[key]; xi = xi_ls(P[1], P[2], P[3], rx, ry, rz, XI_EDGES)
                XI[(s2, key)][k] = xi; X3[(s2, key)][k] = winmean(xi, XI_EDGES, WIN_XI)
            end
        end
    end
    io = open(OUTL, "w"); say(a...) = begin s = string(a...); println(s); println(io, s) end
    sem(v) = begin ok = filter(x -> isfinite(x) && x != 0, v); (mean(ok), std(ok) / sqrt(length(ok)), length(ok)) end
    function ratio(a, b, label)
        # per-cap paired ratio is not meaningful across skies; use ratio of means with SEMs
        ma, sa, na = sem(a); mb, sb, nb = sem(b)
        say(@sprintf("%-40s %9.4f ± %-8.4f / %9.4f ± %-8.4f = %6.3f ± %5.3f  (%d/%d caps)", label, ma, sa, mb, sb, ma / mb,
                     sqrt((sa / mb)^2 + (ma * sb / mb^2)^2), na, nb))
    end
    say("== ξ(3–15 Mpc/h), M>5e12, shell $(SH[1])–$(SH[2]) Mpc/h, $(nc) caps; ours = $(CAMP)fs")
    say(@sprintf("mean |x − xL| Websky (true) %.3f Mpc/h; |x − x_rec| Websky %.3f, ours %.3f; |x_rec − xL| Websky (reconstruction error) %.3f",
                 dW[1] / nW, dW[3] / nW, dO / nO, dW[2] / nW))
    ratio(X3[(:W, :E)], X3[(:O, :E)], "Eulerian        W_E / O_E")
    ratio(X3[(:W, :Lr)], X3[(:O, :Lr)], "Lagrangian rec  W_Lr / O_Lr")
    ratio(X3[(:W, :Lt)], X3[(:W, :Lr)], "control: W_Lt / W_Lr (reconstruction)")
    ratio(X3[(:W, :E)], X3[(:W, :Lr)], "Websky  E / L_rec")
    ratio(X3[(:O, :E)], X3[(:O, :Lr)], "ours    E / L_rec")
    say("\n== ξ(r) by radius: W_E/O_E | W_Lr/O_Lr | W_Lt/W_Lr")
    for (i, r) in enumerate(rc(XI_EDGES))
        f(s) = [XI[s][k][i] for k in 1:nc]
        a1 = mean(filter(isfinite, f((:W, :E)))) / mean(filter(isfinite, f((:O, :E))))
        a2 = mean(filter(isfinite, f((:W, :Lr)))) / mean(filter(isfinite, f((:O, :Lr))))
        a3 = mean(filter(isfinite, f((:W, :Lt)))) / mean(filter(isfinite, f((:W, :Lr))))
        say(@sprintf("r=%6.2f  %.4f  %.4f  %.4f", r, a1, a2, a3))
    end
    close(io); @info "wrote" OUTL
end
lagr_main()
