#!/usr/bin/env julia
# Production full-sky HEALPix halo painter for one frozen-campaign octant
# (docs/paper_comparison_plan_2026-09.md, item A1). Physics = the validated cap painters,
# ported verbatim in halo_profiles.jl; this file adds only HEALPix geometry and threading.
#
#   julia --project=validation -t N paint_octant.jl <OCT> <catalog.pksc> <outdir> [nside] [products] [tag]
#   julia --project=validation -t N paint_octant.jl --selftest
#
# products: comma list of kappa,tsz,ksz (default all). Catalog = finalized (Eulerian Mpc/h,
# km/s) AM catalog, observer from the octZYX bits (±2618 Mpc/h per axis). Outputs (Float32):
#   kappa_halo_comp  — construction B halo term (add the full-matter field κ map)
#   kappa_halo_plain — plain truncated NFW (component plots)
#   tsz_y            — Compton-y, all halos
#   ksz_halo_Wc      — μK, zero-net compensated (add T_CMB × full-matter field kSZ)
#   ksz_halo_W       — μK, uncompensated (component plots)
#
# Every halo is PER-HALO EXACT: pixel-centre sampling misses the cusp of halos near or
# below the pixel scale, so the residual vs the analytic solid-angle integral (plain) or
# vs zero (compensated) is deposited in the halo's centre pixel (TIER_B lesson #5).
# Threading: threads own disjoint HEALPix ring bands (no shared-pixel races); the
# per-halo painted sums for the deposits come from a read-only first pass.
using Healpix, Printf, Statistics, Random
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z, chi
include(joinpath(@__DIR__, "halo_profiles.jl"))

# ---------------------------------- ring geometry ----------------------------------------
struct Rings
    res::Healpix.Resolution
    nring::Int
    θ::Vector{Float64}; cθ::Vector{Float64}; sθ::Vector{Float64}
    first::Vector{Int}; np::Vector{Int}; shift::Vector{Float64}
end
function Rings(nside)
    res = Healpix.Resolution(nside); nr = 4nside - 1
    θ = zeros(nr); first = zeros(Int, nr); np = zeros(Int, nr); sh = zeros(nr)
    for i in 1:nr
        ri = Healpix.getringinfo(res, i)
        θ[i] = ri.colatitude_rad; first[i] = ri.firstPixIdx; np[i] = ri.numOfPixels
        sh[i] = ri.shifted ? 0.5 : 0.0
    end
    Rings(res, nr, θ, cos.(θ), sin.(θ), first, np, sh)
end
function ring_range(R::Rings, θ0, rad)
    ztop = cos(max(0.0, θ0 - rad)); zbot = cos(min(π, θ0 + rad))
    max(1, Healpix.ringAbove(R.res, ztop)), min(R.nring, Healpix.ringAbove(R.res, zbot) + 1)
end
ring_of_pix(R::Rings, p) = searchsortedlast(R.first, p)

# Walk every pixel centre within angular radius `rad` of unit vector v0 (θ0, ϕ0), on rings
# rlo:rhi only; f!(pix, θ) is called for each. Returns nothing (callers accumulate).
@inline function disc_walk(f!::F, R::Rings, θ0, ϕ0, v0, rad, rlo, rhi) where {F}
    c0 = cos(θ0); s0 = sin(θ0); crad = cos(rad)
    pole = (θ0 - rad <= 0) || (θ0 + rad >= π)
    d2max = 2 - 2crad
    @inbounds for ir in rlo:rhi
        np = R.np[ir]; dphi = 2π / np; sh = R.shift[ir]
        cr = R.cθ[ir]; sr = R.sθ[ir]
        full = pole
        jlo = 1; jhi = np
        if !full
            cdp = (crad - cr * c0) / (sr * s0)
            cdp > 1 && continue
            if cdp <= -1
                full = true
            else
                Δ = acos(cdp)
                jlo = floor(Int, (ϕ0 - Δ) / dphi - sh) - 1
                jhi = ceil(Int, (ϕ0 + Δ) / dphi - sh) + 2
                jhi - jlo + 1 >= np && (full = true; jlo = 1; jhi = np)
            end
        end
        for jj in jlo:jhi
            j = full ? jj : mod(jj - 1, np) + 1
            ϕ = (j - 1 + sh) * dphi
            px = sr * cos(ϕ); py = sr * sin(ϕ); pz = cr
            d2 = (px - v0[1])^2 + (py - v0[2])^2 + (pz - v0[3])^2
            d2 <= d2max || continue
            f!(R.first[ir] + j - 1, 2 * asin(min(sqrt(d2) / 2, 1.0)))
        end
    end
    nothing
end

# --------------------------------- halo parameters ---------------------------------------
struct Halo
    θ0::Float64; ϕ0::Float64; v::NTuple{3,Float64}
    z::Float64; χ::Float64
    M::Float64; W::Float64; Rk::Float64                 # κ: M200m Msun/h, kernel, radius[rad]
    mh::Float64; θv::Float64; yamp::Float64             # tSZ (all halos)
    ksz::Bool; kfac::Float64; tamp::Float64; ucomp::Float64; itk::Float64
end
function make_halo(x, y, z3, vx, vy, vz, R, obs, chi2z, TK)
    dx = x - obs[1]; dy = y - obs[2]; dz = z3 - obs[3]
    χ = sqrt(dx^2 + dy^2 + dz^2)
    v = (dx / χ, dy / χ, dz / χ)
    θ0 = acos(clamp(v[3], -1.0, 1.0)); ϕ0 = mod(atan(v[2], v[1]), 2π)
    z = chi_to_z(chi2z, χ)
    M = 4 / 3 * π * rho_mh * R^3                          # Msun/h
    W = kernel_kappa(z, χ)
    Rk = paint_radius_kappa(M) / χ
    mh = sqrt(deltacrit(z) / 200) * M / hub               # Msun, SIS M200c proxy
    θv = rvir_com_mpc(mh, z) * hub / χ
    yamp = y_amp(mh, z)
    ksz = mh > MMIN_KSZ && θv > THMIN_KSZ
    vr = (vx * v[1] + vy * v[2] + vz * v[3])
    kfac = -TCMB_UK * vr / C_KMS
    tamp = tau_amp(mh, z)
    itk = ksz ? itot_interp(TK, mh, z) : 0.0
    Halo(θ0, ϕ0, v, z, χ, M, W, Rk, mh, θv, yamp, ksz, kfac, tamp, itk * 3 / (256π), itk)
end
maxrad(h::Halo, P) = max(P.kappa ? h.Rk : 0.0, P.tsz ? 4h.θv : 0.0, (P.ksz && h.ksz) ? 4h.θv : 0.0)

# ------------------------------- per-product evaluation ----------------------------------
mutable struct Acc
    k0::Float64; kc::Float64; y::Float64; w::Float64; wc::Float64
end
Acc() = Acc(0.0, 0.0, 0.0, 0.0, 0.0)

# visit all products of one halo on rings rlo:rhi; if `maps === nothing` only accumulate
function visit!(acc::Acc, maps, R::Rings, h::Halo, P, TY, TK, rlo, rhi)
    if P.kappa
        a, b = ring_range(R, h.θ0, h.Rk); a = max(a, rlo); b = min(b, rhi)
        disc_walk(R, h.θ0, h.ϕ0, h.v, h.Rk, a, b) do pix, θ
            k0, kc = kappa_pair(θ * h.χ, h.M, h.W)
            acc.k0 += k0; acc.kc += kc
            if maps !== nothing
                maps.k0[pix] += k0; maps.kc[pix] += kc
            end
        end
    end
    if P.tsz
        rad = 4h.θv
        a, b = ring_range(R, h.θ0, rad); a = max(a, rlo); b = min(b, rhi)
        disc_walk(R, h.θ0, h.ϕ0, h.v, rad, a, b) do pix, θ
            yv = h.yamp * sigma_interp(TY, θ / h.θv, h.mh, h.z)
            acc.y += yv
            maps !== nothing && (maps.y[pix] += yv)
        end
    end
    if P.ksz && h.ksz
        rad = 4h.θv
        a, b = ring_range(R, h.θ0, rad); a = max(a, rlo); b = min(b, rhi)
        disc_walk(R, h.θ0, h.ϕ0, h.v, rad, a, b) do pix, θ
            xb = θ / h.θv
            Σ = sigma_interp(TK, xb, h.mh, h.z)
            w = h.kfac * h.tamp * Σ
            wc = h.kfac * h.tamp * (Σ - h.ucomp * 2 * sqrt(max(16.0 - xb^2, 0.0)))
            acc.w += w; acc.wc += wc
            if maps !== nothing
                maps.w[pix] += w; maps.wc[pix] += wc
            end
        end
    end
    acc
end

# analytic per-halo targets in pixel units (Σ pixel values after the deposit)
function targets(h::Halo, TY, Ωpix)
    k0 = kappa_integral(h.M, h.W, h.χ) / Ωpix
    y = h.yamp * itot_interp(TY, h.mh, h.z) * h.θv^2 / Ωpix
    w = h.ksz ? h.kfac * h.tamp * h.itk * h.θv^2 / Ωpix : 0.0
    (k0, 0.0, y, w, 0.0)
end

# ------------------------------------ painting -------------------------------------------
function paint_chunk!(maps, R::Rings, halos::Vector{Halo}, P, TY, TK, Ωpix, nband)
    n = length(halos)
    dep = zeros(5, n); cpix = zeros(Int, n); cring = zeros(Int, n)
    # pass 1 (read-only): painted sums → exactness deposits
    Threads.@threads :dynamic for i in 1:n
        h = halos[i]
        acc = visit!(Acc(), nothing, R, h, P, TY, TK, 1, R.nring)
        t = targets(h, TY, Ωpix)
        dep[1, i] = P.kappa ? t[1] - acc.k0 : 0.0
        dep[2, i] = P.kappa ? -acc.kc : 0.0
        dep[3, i] = P.tsz ? t[3] - acc.y : 0.0
        dep[4, i] = (P.ksz && h.ksz) ? t[4] - acc.w : 0.0
        dep[5, i] = (P.ksz && h.ksz) ? -acc.wc : 0.0
        cpix[i] = Healpix.ang2pixRing(R.res, h.θ0, h.ϕ0)
        cring[i] = ring_of_pix(R, cpix[i])
    end
    # band lists (serial; each halo touches few bands)
    rpb = cld(R.nring, nband)
    lists = [Int32[] for _ in 1:nband]
    for i in 1:n
        h = halos[i]
        a, b = ring_range(R, h.θ0, maxrad(h, P))
        a = min(a, cring[i]); b = max(b, cring[i])
        for bb in ((a - 1) ÷ rpb + 1):((b - 1) ÷ rpb + 1)
            push!(lists[bb], Int32(i))
        end
    end
    # pass 2: each band is painted by exactly one task
    Threads.@threads :dynamic for bb in 1:nband
        rlo = (bb - 1) * rpb + 1; rhi = min(bb * rpb, R.nring)
        for i32 in lists[bb]
            i = Int(i32); h = halos[i]
            visit!(Acc(), maps, R, h, P, TY, TK, rlo, rhi)
            if rlo <= cring[i] <= rhi
                p = cpix[i]
                P.kappa && (maps.k0[p] += dep[1, i]; maps.kc[p] += dep[2, i])
                P.tsz && (maps.y[p] += dep[3, i])
                (P.ksz && h.ksz) && (maps.w[p] += dep[4, i]; maps.wc[p] += dep[5, i])
            end
        end
    end
    dep
end

function newmaps(npix, P)
    z() = zeros(Float64, npix); e = Float64[]
    (k0=P.kappa ? z() : e, kc=P.kappa ? z() : e, y=P.tsz ? z() : e,
     w=P.ksz ? z() : e, wc=P.ksz ? z() : e)
end

# -------------------------------------- self-test ----------------------------------------
function selftest()
    rng = MersenneTwister(1)
    for nside in (64, 1024, 4096)
        R = Rings(nside); npix = 12nside^2
        # (1) our ring φ formula == pix2angRing
        bad = 0
        for _ in 1:20000
            p = rand(rng, 1:npix); ir = ring_of_pix(R, p)
            ϕ = (p - R.first[ir] + R.shift[ir]) * 2π / R.np[ir]
            θr, ϕr = Healpix.pix2angRing(R.res, p)
            (abs(θr - R.θ[ir]) < 1e-12 && abs(mod(ϕ - ϕr + π, 2π) - π) < 1e-9) || (bad += 1)
        end
        bad == 0 || error("ring geometry mismatch at nside $nside: $bad")
        # (2) disc pixel sets == Healpix.queryDiscRing (non-inclusive: centres within radius)
        nbad = 0
        radmax = nside == 64 ? 0.6 : nside == 1024 ? 0.05 : 0.01
        for _ in 1:(nside == 1024 ? 1000 : 300)
            θ0 = acos(2rand(rng) - 1); ϕ0 = 2π * rand(rng)
            rad = exp(log(1e-4) + rand(rng) * (log(radmax) - log(1e-4)))
            v0 = (sin(θ0) * cos(ϕ0), sin(θ0) * sin(ϕ0), cos(θ0))
            mine = Int[]
            a, b = ring_range(R, θ0, rad)
            disc_walk((p, θ) -> push!(mine, p), R, θ0, ϕ0, v0, rad, a, b)
            ref = Healpix.queryDiscRing(R.res, θ0, ϕ0, rad)
            # allow boundary ties (centre at |θ-rad| < 1e-9)
            d = symdiff(Set(mine), Set(ref))
            for p in d
                θp, ϕp = Healpix.pix2angRing(R.res, p)
                vp = (sin(θp) * cos(ϕp), sin(θp) * sin(ϕp), cos(θp))
                ang = acos(clamp(sum(vp .* v0), -1.0, 1.0))
                abs(ang - rad) < 1e-7 || (nbad += 1)
            end
            length(mine) == length(Set(mine)) || (nbad += 1)
        end
        nbad == 0 || error("disc pixel set mismatch vs queryDiscRing at nside $nside: $nbad")
        @printf("geometry self-test nside=%d passed\n", nside)
    end
    TY = BattTable(ptilde, 1e11, 2e16, 80); TK = BattTable(rhotilde, 1e13, 2e16, 60)
    profile_selftests(TY, TK)
    # (3) synthetic painting: exact totals; band-parallel == serial (nband=1)
    nside = 1024; R = Rings(nside); npix = 12nside^2; Ωpix = 4π / npix
    cosmo = CosmologyParams(Om, OmB, OL, hub, 0.965, 0.81)
    chi2z = build_chi_to_z(cosmo; z_max=6.0)
    obs = (-2618.0, -2618.0, -2618.0)
    halos = Halo[]
    for _ in 1:4000
        χ = 50 + 3000 * rand(rng)^0.5
        u = (rand(rng), rand(rng), rand(rng)); u = u ./ sqrt(sum(abs2, u))
        R_TH = exp(log(1.3) + rand(rng) * (log(9.0) - log(1.3)))   # 1e12..3e15 Msun/h
        push!(halos, make_halo(obs[1] + χ * u[1], obs[2] + χ * u[2], obs[3] + χ * u[3],
                               600randn(rng), 600randn(rng), 600randn(rng), R_TH, obs, chi2z, TK))
    end
    P = (kappa=true, tsz=true, ksz=true)
    mA = newmaps(npix, P); paint_chunk!(mA, R, halos, P, TY, TK, Ωpix, 64)
    mB = newmaps(npix, P); paint_chunk!(mB, R, halos, P, TY, TK, Ωpix, 1)
    tot = reduce(.+, (collect(targets(h, TY, Ωpix)) for h in halos))
    for (k, t) in zip((:k0, :kc, :y, :w, :wc), tot)
        s = sum(getfield(mA, k)); sB = sum(getfield(mB, k))
        scale = k in (:kc,) ? tot[1] : k == :wc ? sum(abs, [targets(h, TY, Ωpix)[4] for h in halos]) : abs(t)
        abs(s - t) <= 1e-8 * scale || error("exactness failed for $k: $s vs $t")
        maximum(abs.(getfield(mA, k) .- getfield(mB, k))) <= 1e-12 * maximum(abs.(getfield(mB, k))) ||
            error("band-parallel != serial for $k")
    end
    @printf("painting self-test passed (%d halos, %d threads, %d with kSZ)\n",
            length(halos), Threads.nthreads(), count(h -> h.ksz, halos))
end

# ---------------------------------------- main -------------------------------------------
if length(ARGS) >= 1 && ARGS[1] == "--selftest"
    selftest(); exit(0)
end

function main(ARGS)
    OCT = ARGS[1]; CAT = ARGS[2]; OUTD = ARGS[3]
    nside = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 4096
    prods = length(ARGS) >= 5 ? split(ARGS[5], ",") : ["kappa", "tsz", "ksz"]
    tag = length(ARGS) >= 6 ? ARGS[6] : "prod_oct$(OCT)_AMv2"          # output name tag
    P = (kappa="kappa" in prods, tsz="tsz" in prods, ksz="ksz" in prods)
    length(OCT) == 3 && all(in("01"), OCT) || error("OCT must be three 0/1 bits (ZYX)")
    obs = ntuple(d -> OCT[4-d] == '1' ? 2618.0 : -2618.0, 3)     # OCT = "ZYX"
    mkpath(OUTD)

    cosmo = CosmologyParams(Om, OmB, OL, hub, 0.965, 0.81)
    chi2z = build_chi_to_z(cosmo; z_max=6.0)
    chimax = chi(4.5, cosmo)
    @info "paint_octant" OCT obs CAT nside P chimax threads=Threads.nthreads()
    t0 = time()
    TY = BattTable(ptilde, 1e11, 2e16, 80); TK = BattTable(rhotilde, 1e13, 2e16, 60)
    profile_selftests(TY, TK)

    R = Rings(nside); npix = 12nside^2; Ωpix = 4π / npix
    maps = newmaps(npix, P)
    nband = 8 * Threads.nthreads()
    nh_tot = open(io -> Int(read(io, Int32)), CAT)
    tot = zeros(5); nsel = zeros(Int, 3); nused = 0
    open(CAT) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        nf = (filesize(CAT) - 12) ÷ (4nh_tot)
        nf in (11, 33) || error("unexpected record size $nf")
        chunk = 2_000_000; buf = Vector{Float32}(undef, chunk * nf); ndone = 0
        halos = Halo[]
        while ndone < nh_tot
            m = min(chunk, nh_tot - ndone); read!(io, view(buf, 1:m*nf))
            empty!(halos)
            for k in 1:m
                b = (k - 1) * nf
                dx = buf[b+1] - obs[1]; dy = buf[b+2] - obs[2]; dz = buf[b+3] - obs[3]
                r = sqrt(Float64(dx)^2 + Float64(dy)^2 + Float64(dz)^2)
                (30.0 <= r <= chimax) || continue
                push!(halos, make_halo(Float64(buf[b+1]), Float64(buf[b+2]), Float64(buf[b+3]),
                                       Float64(buf[b+4]), Float64(buf[b+5]), Float64(buf[b+6]),
                                       Float64(buf[b+7]), obs, chi2z, TK))
            end
            paint_chunk!(maps, R, halos, P, TY, TK, Ωpix, nband)
            for h in halos
                t = targets(h, TY, Ωpix)
                tot[1] += t[1]; tot[3] += t[3]; tot[4] += t[4]
            end
            nused += length(halos); nsel[3] += count(h -> h.ksz, halos)
            ndone += m
            ndone % 20_000_000 < chunk &&
                @info @sprintf("%.0f%% (%d/%d halos, %.1f min)", 100ndone / nh_tot, ndone, nh_tot, (time() - t0) / 60)
        end
    end

    @printf("\nhalos painted: %d of %d (30 ≤ χ ≤ %.1f Mpc/h); kSZ-selected: %d\n", nused, nh_tot, chimax, nsel[3])
    if P.kappa
        @printf("κ plain: Σpix = %.6e, analytic = %.6e (ratio %.9f); κ comp Σpix/plain = %.2e\n",
                sum(maps.k0), tot[1], sum(maps.k0) / tot[1], sum(maps.kc) / tot[1])
    end
    if P.tsz
        @printf("tSZ: Σpix/analytic = %.9f; mean y over octant (1/8 sky) = %.4e\n",
                sum(maps.y) / tot[3], sum(maps.y) / (npix / 8))
    end
    if P.ksz
        @printf("kSZ W: Σpix/analytic = %.9f; Wc net Σ/Σ|W| = %.2e\n",
                sum(maps.w) / tot[4], sum(maps.wc) / sum(abs, maps.w))
    end
    function save(name, v)
        isempty(v) && return
        m = HealpixMap{Float32,RingOrder}(nside); m.pixels .= Float32.(v)
        out = joinpath(OUTD, "$(name)_$(tag)_nside$(nside).fits")
        Healpix.saveToFITS(m, "!" * out, typechar="E")
        @info "wrote $out"
    end
    save("kappa_halo_comp", maps.kc); save("kappa_halo_plain", maps.k0)
    save("tsz_y", maps.y); save("ksz_halo_Wc", maps.wc); save("ksz_halo_W", maps.w)
    @printf("done in %.1f min\n", (time() - t0) / 60)
end

main(ARGS)
