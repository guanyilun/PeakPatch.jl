#!/usr/bin/env julia
# Tier-A against the FULL Websky halo catalog (websky_ref/halos.pksc, 862.9M halos, observer at the
# origin, Mpc units), replacing the single 10°×10° patch (tierA_v3fs.jl). The same 5.00° cap and the
# same shell / bins / windows as tierA_v3fs.jl are laid out on a sky-wide grid (centres ≥ 10.5°
# apart and ≥ 5.5° from the octant planes, so each cap lies inside one octant of both skies), and every
# statistic is measured in every cap of BOTH catalogs at the same sky positions. Errors: standard error
# of the mean over caps on each side, combined in quadrature. Also full-sky N(>M|z) and dN/dz
# (no caps) for both catalogs, and Websky full-sky N(>M|z) vs Tinker08.
#   env: TIERA_CAMP (default v4) for ours; TIERA_MAXH for smoke tests
include(joinpath(@__DIR__, "tierA_v3fs.jl"))          # helpers (Cap, addhalo!, pairstats, binning); main() guarded
const WFULL = joinpath(get(ENV, "WS_REF", "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"), "halos.pksc")
const OUTFS = MAXH > 0 ? "/tmp/tierA_fullsky_smoke.txt" : joinpath(@__DIR__, "results", "tierA_fullsky_$(CAMP)fs.txt")

# ---- cap centres: greedy on a Fibonacci sphere ----
function cap_centres(; sep=deg2rad(10.5), plane=deg2rad(5.5), n=20000)
    out = NTuple{3,Float64}[]; ga = π * (3 - sqrt(5))
    for i in 0:n-1
        zc = 1 - 2(i + 0.5) / n; r = sqrt(1 - zc^2); ϕ = ga * i
        v = (r * cos(ϕ), r * sin(ϕ), zc)
        all(abs.(v) .>= sin(plane)) || continue
        all(u -> v[1] * u[1] + v[2] * u[2] + v[3] * u[3] < cos(sep), out) && push!(out, v)
    end
    out
end
const CENT = cap_centres()
# 1°×1° (θ, ϕ) lookup grid → candidate caps
const NT, NP = 180, 360
const GRID = [Int[] for _ in 1:NT, _ in 1:NP]
for it in 1:NT, ip in 1:NP
    θ = (it - 0.5) * π / NT; ϕ = (ip - 0.5) * 2π / NP; v = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))
    for (k, c) in enumerate(CENT)
        acos(clamp(v[1] * c[1] + v[2] * c[2] + v[3] * c[3], -1.0, 1.0)) < acos(COSW) + deg2rad(1.5) && push!(GRID[it, ip], k)
    end
end
@inline function capof(nx, ny, nz)
    it = clamp(floor(Int, acos(clamp(nz, -1.0, 1.0)) / π * NT) + 1, 1, NT)
    ip = clamp(floor(Int, mod(atan(ny, nx), 2π) / (2π) * NP) + 1, 1, NP)
    for k in GRID[it, ip]
        c = CENT[k]; nx * c[1] + ny * c[2] + nz * c[3] >= COSW && return k
    end
    0
end
octof(c) = string(c[3] < 0 ? '1' : '0', c[2] < 0 ? '1' : '0', c[1] < 0 ? '1' : '0')   # axis looks toward − ⇒ bit 1

function fullsky_main()
    nc = length(CENT)
    @info "caps" nc per_octant = [count(c -> octof(c) == o, CENT) for o in OCTS]
    capsW = [Cap((0.0, 0.0, 0.0)) for _ in 1:nc]
    capsO = [Cap(obs_of(octof(c))) for c in CENT]
    fsW = zeros(length(TK_M), length(TK_Z) - 1); fsO = zeros(length(TK_M), length(TK_Z) - 1)
    tally!(F, M, z) = begin it = binof(TK_Z, z); it > 0 && for (j, t) in enumerate(TK_M); M > t && (F[j, it] += 1); end end
    # Websky full sky: 10 floats per halo
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
                M = mass(HUB * Float64(b[o+7])); z = chi_to_z(CHI2Z, r)
                tally!(fsW, M, z)
                kc = capof(X / r, Y / r, Z / r)
                kc > 0 && addhalo!(capsW[kc], X, Y, Z, Float64(b[o+4]), Float64(b[o+5]), Float64(b[o+6]), M, r, z)
            end
            nd += m
        end
    end
    @info "Websky streamed" nh
    for oc in OCTS
        ob = obs_of(oc)
        stream(joinpath(S3, "catalog_websky_6144_$(CAMP)_oct$(oc)_AMfs.pksc"), (bf, b) -> begin
            R = Float64(bf[b+7]); R > 0 || return
            X = Float64(bf[b+1]); Y = Float64(bf[b+2]); Z = Float64(bf[b+3])
            dx = X - ob[1]; dy = Y - ob[2]; dz = Z - ob[3]; r = sqrt(dx^2 + dy^2 + dz^2)
            M = mass(R); z = chi_to_z(CHI2Z, r); tally!(fsO, M, z)
            kc = capof(dx / r, dy / r, dz / r)
            (kc > 0 && octof(CENT[kc]) == oc) && addhalo!(capsO[kc], X, Y, Z, Float64(bf[b+4]), Float64(bf[b+5]), Float64(bf[b+6]), M, r, z)
        end)
        @info "ours streamed" oc
    end
    @info "pair statistics" nc
    PW = Vector{Any}(undef, nc); PO = Vector{Any}(undef, nc)
    for k in 1:nc
        ax = collect(CENT[k]); PW[k] = pairstats(capsW[k], ax, 10k + 1); PO[k] = pairstats(capsO[k], ax, 10k + 2)
    end
    io = open(OUTFS, "w"); say(a...) = begin s = string(a...); println(s); println(io, s) end
    sem(v) = begin ok = filter(isfinite, v); (mean(ok), std(ok) / sqrt(length(ok)), length(ok)) end
    function line(label, w, o)
        mw, sw, nw = sem(w); mo, so, no = sem(o); s = sqrt(sw^2 + so^2)
        say(@sprintf("%-34s Websky %11.4g ± %-9.3g | ours %11.4g ± %-9.3g (%d/%d caps) | W/ours %6.3f ± %5.3f | z %+5.1f",
                     label, mw, sw, mo, so, nw, no, mw / mo, sqrt((sw / mo)^2 + (mw * so / mo^2)^2), (mw - mo) / s))
    end
    say("== Tier-A vs the FULL Websky catalog: same $(nc) caps (half-angle $(round(acosd(COSW); digits=3))°) in both skies; ours = $(CAMP)fs")
    say("   errors: standard error of the mean over caps per side, combined in quadrature")
    say("\n== N(>M) per cap, z-bins (after AM; largely enforced by AM)")
    for (j, t) in enumerate(NM_THR), iz in 1:length(NZ_EDGES)-1
        line(@sprintf("N(>%.1e) z %.1f-%.1f", t, NZ_EDGES[iz], NZ_EDGES[iz+1]), [c.nm[j, iz] for c in capsW], [c.nm[j, iz] for c in capsO])
    end
    say("\n== dN/dz per cap (Δz = 0.25)")
    for (j, t) in enumerate(DNDZ_M), id in 1:length(DNDZ_EDGES)-1
        line(@sprintf("dN/dz M>%.0e z %.2f-%.2f", t, DNDZ_EDGES[id], DNDZ_EDGES[id+1]), [c.dndz[j, id] for c in capsW], [c.dndz[j, id] for c in capsO])
    end
    sig(v) = v[1] > 1 ? sqrt(v[3] / v[1] - (v[2] / v[1])^2) : NaN
    say("\n== σ_vr (km/s) per mass bin")
    for iv in 1:length(VM_EDGES)-1
        line(@sprintf("σ_vr M %.1e-%.1e", VM_EDGES[iv], VM_EDGES[iv+1]), [sig(c.vs[:, iv]) for c in capsW], [sig(c.vs[:, iv]) for c in capsO])
    end
    say(@sprintf("\n== ξ(r), shell %.0f-%.0f Mpc/h, M>%.0e", SH..., MLOW_XI))
    for (i, r) in enumerate(rc(XI_EDGES)); line(@sprintf("ξ r=%.2f", r), [p.xi[i] for p in PW], [p.xi[i] for p in PO]); end
    line(@sprintf("ξ mean over %.0f-%.0f Mpc/h", WIN_XI...), [winmean(p.xi, XI_EDGES, WIN_XI) for p in PW], [winmean(p.xi, XI_EDGES, WIN_XI) for p in PO])
    say(@sprintf("\n== relative bias √<ξ_hh(M)> over %.0f-%.0f Mpc/h", WIN_B...))
    for ib in 1:length(BM_EDGES)-1
        line(@sprintf("b ∝ √ξ, M %.1e-%.1e", BM_EDGES[ib], BM_EDGES[ib+1]), [sqrt(max(p.bxi[ib], 0)) for p in PW], [sqrt(max(p.bxi[ib], 0)) for p in PO])
    end
    v12c = [(V12_EDGES[i] + V12_EDGES[i+1]) / 2 for i in 1:length(V12_EDGES)-1]
    wm(v) = mean(v[[i for i in eachindex(v12c) if WIN_V12[1] <= v12c[i] <= WIN_V12[2]]])
    say(@sprintf("\n== v12 (km/s), M>%.0e", MCUT_V12))
    for (i, r) in enumerate(v12c); line(@sprintf("v12 r=%.1f", r), [p.v12[i] for p in PW], [p.v12[i] for p in PO]); end
    line(@sprintf("v12 mean over %.0f-%.0f Mpc/h", WIN_V12...), [wm(p.v12) for p in PW], [wm(p.v12) for p in PO])
    # full sky, no caps; and Websky vs Tinker08
    pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky_RAW_unnormalized.dat"))
    Mg = 10 .^ range(11.5, 16; length=901); lnMg = log.(Mg); sgm = precompute_sigma(Mg, pk, 0.31)
    function ntk(Mcut, za, zb)
        ra, rb = chi(za, COSMO), chi(zb, COSMO); ns = 40; acc = 0.0
        for j in 1:ns
            r0 = ra + (j - 1) * (rb - ra) / ns; r1 = ra + j * (rb - ra) / ns
            z = chi_to_z(CHI2Z, (r0 + r1) / 2); Dg = growth_factor(z, COSMO); n = 0.0
            for i in 1:length(Mg)-1
                Mm = sqrt(Mg[i] * Mg[i+1]); Mm < Mcut && continue
                n += tinker_dndlnM(Mm, Dg * sqrt(sgm[i] * sgm[i+1]), (log(sgm[i+1]) - log(sgm[i])) / (lnMg[i+1] - lnMg[i]), z, 0.31) * (lnMg[i+1] - lnMg[i])
            end
            acc += n * (4π / 3) * (r1^3 - r0^3)
        end
        acc
    end
    say("\n== full sky (no caps): N(>M | z), Websky vs ours vs Tinker08")
    say(@sprintf("%-24s %12s %12s %12s %8s %8s %8s", "bin", "Websky", "ours", "Tinker", "W/ours", "W/Tk", "ours/Tk"))
    for (j, t) in enumerate(TK_M), iz in 1:length(TK_Z)-1
        e = ntk(t, TK_Z[iz], TK_Z[iz+1])
        say(@sprintf("M>%.0e z %.2f-%.2f %12d %12d %12.1f %8.4f %8.4f %8.4f", t, TK_Z[iz], TK_Z[iz+1], fsW[j, iz], fsO[j, iz], e,
                     fsW[j, iz] / fsO[j, iz], fsW[j, iz] / e, fsO[j, iz] / e))
    end
    close(io); @info "wrote" OUTFS
end
if abspath(PROGRAM_FILE) == (@__FILE__)   # include()-able (tierA_lagrangian.jl)
    fullsky_main()
end
