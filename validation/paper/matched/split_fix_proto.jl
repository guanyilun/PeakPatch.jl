#!/usr/bin/env julia
# Prototype of the ψ fix for the multires split (docs/split_psi_fix_literature_2026-10.md). Same Threefry noise
# as the exact global field; tile fields compared in the core.
#   prod   : current production ψ1 = isolated(residual, √P·ψ-kernel) + interp(coarse ψ, D/T-compensated)
#   poisson: ψ1 = interp(ψ_L) + isolated Poisson[δ_tile − interp(δ_L)]          (MUSIC principle, Gaussian split)
#   treepm : ψ1 = interp(ψ_L) + isolated[(1 − G)·ψ-kernel ∗ δ_tile]             (TreePM-style short kernel)
# with δ_L, ψ_L = coarse fields × G(k) = exp(−k² r_s²) (same coarse noise, same compensation, same interpolator).
# δ_tile is the production tile δ (already accurate to ≲1% error power).
#   usage: julia --project=validation -t 32 split_fix_proto.jl <config.toml> <rs list, fine cells> <nbuff list>
using TOML, Printf, FFTW, Statistics
using PeakPatch
import PeakPatch.PowerSpectrum: load_pk
import PeakPatch.MultiResolution: _downsample_noise, _splice_compensation, _periodic_convolve!, _kernel_1lpt,
    _generate_extended_residual, _isolated_convolve, _interpolate_to_tile

# isolated (2×-padded) convolution of a real tile field with a k-space kernel (no √P)
function iso_kernel(f::Array{Float32,3}, Lb::Float64, kern)
    n = size(f, 1); n2 = 2n; dx = Lb / n
    p = zeros(Float32, n2, n2, n2); p[1:n, 1:n, 1:n] .= f
    pk = rfft(p); dk = 2π / (n2 * dx)
    kx = FFTW.rfftfreq(n2, n2 * dk); ky = FFTW.fftfreq(n2, n2 * dk)
    Threads.@threads for iz in 1:n2
        @inbounds for iy in 1:n2, ix in 1:size(pk, 1)
            k2 = kx[ix]^2 + ky[iy]^2 + ky[iz]^2
            pk[ix, iy, iz] = k2 == 0 ? 0 : pk[ix, iy, iz] * kern(kx[ix], ky[iy], ky[iz], k2)
        end
    end
    irfft(pk, n2)[1:n, 1:n, 1:n]
end

const SPEC = get(ENV, "SPEC", "0") == "1"
function errspec(r, s, a, kNc)
    e = s .- r; nsub = size(r, 1); R = rfft(r); E = rfft(e); kf = 2π / (nsub * a)
    ed = collect(10 .^ range(log10(kf), log10(π / a); length=15)); Pr = zeros(14); Pe = zeros(14)
    for k3 in 1:nsub, k2 in 1:nsub, k1 in 1:nsub÷2+1
        q2 = k2 <= nsub ÷ 2 + 1 ? k2 - 1 : k2 - 1 - nsub; q3 = k3 <= nsub ÷ 2 + 1 ? k3 - 1 : k3 - 1 - nsub
        b = searchsortedlast(ed, kf * sqrt((k1 - 1)^2 + q2^2 + q3^2)); 1 <= b <= 14 || continue
        Pr[b] += abs2(R[k1, k2, k3]); Pe[b] += abs2(E[k1, k2, k3])
    end
    ed, Pr, Pe
end
function errstats(r, s, a, kNc)
    e = s .- r; nsub = size(r, 1)
    R = rfft(r); E = rfft(e); kf = 2π / (nsub * a)
    band = (0.0, 0.0); low = (0.0, 0.0)
    for k3 in 1:nsub, k2 in 1:nsub, k1 in 1:nsub÷2+1
        q2 = k2 <= nsub ÷ 2 + 1 ? k2 - 1 : k2 - 1 - nsub; q3 = k3 <= nsub ÷ 2 + 1 ? k3 - 1 : k3 - 1 - nsub
        kk = kf * sqrt((k1 - 1)^2 + q2^2 + q3^2)
        if kNc <= kk < 2kNc
            band = (band[1] + abs2(R[k1, k2, k3]), band[2] + abs2(E[k1, k2, k3]))
        elseif 0 < kk < 0.5kNc
            low = (low[1] + abs2(R[k1, k2, k3]), low[2] + abs2(E[k1, k2, k3]))
        end
    end
    # rms error near the core edge (0–6 cells) vs interior (≥ 30)
    ed = (0.0, 0.0); it = (0.0, 0.0)
    for z in 1:nsub, y in 1:nsub, x in 1:nsub
        d = min(x, y, z, nsub + 1 - x, nsub + 1 - y, nsub + 1 - z) - 1
        d < 6 && (ed = (ed[1] + r[x, y, z]^2, ed[2] + e[x, y, z]^2))
        d >= 30 && (it = (it[1] + r[x, y, z]^2, it[2] + e[x, y, z]^2))
    end
    (sqrt(sum(abs2, e) / sum(abs2, r)), band[2] / band[1], low[2] / low[1], sqrt(ed[2] / ed[1]), sqrt(it[2] / it[1]))
end

function main()
    cfgd = TOML.parsefile(ARGS[1]); cfg = PipelineConfig(cfgd); rc = cfgd["run"]
    rss = parse.(Float64, split(ARGS[2], ",")); nbs = parse.(Int, split(ARGS[3], ","))
    ntile = rc["ntile"]; seed = rc["seed"]; nsub, N = grid_layout(cfg, ntile); @assert cfg.periodic_cores
    a = cfg.boxsize / cfg.n; L = N * a; M = ntile * rc["coarse_factor"]; block = N ÷ M; kNc = π / (block * a)
    pk = load_pk(cfg.pkfile)
    δg = generate_grf(N, pk, L, seed); ψg = displacements_1lpt(rfft(δg), N, L)[1]
    cn = _downsample_noise(N, M, seed); ck = rfft(cn)
    comp = cfg.coarse_compensation ? _splice_compensation(M, block) : nothing
    coarse(kern) = (q = copy(ck); _periodic_convolve!(q, pk, M, L; kernel_fn=kern, comp=comp); irfft(q, M))
    δc = coarse(nothing); ψc = coarse(_kernel_1lpt(1))
    out = open(joinpath(@__DIR__, "..", "results", "split_fix_proto_cf$(rc["coarse_factor"]).txt"), "w")
    say(x...) = begin s = string(x...); println(s); println(out, s) end
    say(@sprintf("N=%d nsub=%d block=%d (coarse cell %.2f Mpc/h), k_N,coarse=%.3f h/Mpc; ψ1 errors in tile cores (2 tiles)", N, nsub, block, block * a, kNc))
    say("variant           r_s[cells] nbuff | ψ rms err | Perr/Pref k∈[kN,2kN) | Perr/Pref k<kN/2 | rms err edge(0–6) | interior(≥30)")
    tiles = [(2, 2, 2), (3, 3, 3)]
    for nb in nbs
        nmesh = nsub + 2nb; Lb = nmesh * a; core = nb+1:nb+nsub
        # production δ_tile and ψ (baseline) for these tiles
        base = Dict{NTuple{3,Int},Tuple{Array{Float32,3},Array{Float32,3}}}()
        for t in tiles
            res = _generate_extended_residual(t..., nsub, nmesh, N, seed, cn, M, 0)
            δt = _isolated_convolve(res, pk, Lb, nmesh) .+ _interpolate_to_tile(δc, t..., nsub, nmesh, N, M)
            ψt = _isolated_convolve(res, pk, Lb, nmesh; kernel_fn=_kernel_1lpt(1)) .+ _interpolate_to_tile(ψc, t..., nsub, nmesh, N, M)
            base[t] = (δt, ψt)
        end
        function report(lab, rs, ψof)
            acc = zeros(5); ws = 0
            for t in tiles
                gi = (t[1]-1)*nsub+1:t[1]*nsub; gj = (t[2]-1)*nsub+1:t[2]*nsub; gk = (t[3]-1)*nsub+1:t[3]*nsub
                s = errstats(Float64.(ψg[gi, gj, gk]), Float64.(ψof(t)[core, core, core]), a, kNc)
                acc .+= collect(s); ws += 1
            end
            v = acc ./ ws
            if SPEC
                t = tiles[1]; gi = (t[1]-1)*nsub+1:t[1]*nsub; gj = (t[2]-1)*nsub+1:t[2]*nsub; gk = (t[3]-1)*nsub+1:t[3]*nsub
                ed, Pr, Pe = errspec(Float64.(ψg[gi, gj, gk]), Float64.(ψof(t)[core, core, core]), a, kNc)
                tp = sum(Pr)
                say("    k/k_N:  " * join([@sprintf("%6.2f ", sqrt(ed[b] * ed[b+1]) / kNc) for b in 1:14]))
                say("    Pe/Pr:  " * join([@sprintf("%6.4f ", Pe[b] / Pr[b]) for b in 1:14]))
                say("    Pe/ΣPr: " * join([@sprintf("%6.4f ", Pe[b] / tp) for b in 1:14]))
            end
            say(@sprintf("%-17s %8.1f %6d | %8.4f  | %12.2e         | %12.2e     | %8.4f          | %8.4f", lab, rs, nb, v[1], v[2], v[3], v[4], v[5]))
        end
        report("prod", NaN, t -> base[t][2])
        for rs in rss
            G(k2) = exp(-k2 * (rs * a)^2)
            δLc = coarse((kx, ky, kz, k2) -> G(k2)); ψLc = coarse((kx, ky, kz, k2) -> G(k2) * _kernel_1lpt(1)(kx, ky, kz, k2))
            ψpois = Dict{NTuple{3,Int},Array{Float32,3}}(); ψtpm = Dict{NTuple{3,Int},Array{Float32,3}}()
            for t in tiles
                δt = base[t][1]
                δL = _interpolate_to_tile(δLc, t..., nsub, nmesh, N, M); ψL = _interpolate_to_tile(ψLc, t..., nsub, nmesh, N, M)
                ψpois[t] = ψL .+ iso_kernel(δt .- δL, Lb, (kx, ky, kz, k2) -> im * kx / k2)
                ψtpm[t] = ψL .+ iso_kernel(δt, Lb, (kx, ky, kz, k2) -> (1 - G(k2)) * im * kx / k2)
            end
            report("poisson", rs, t -> ψpois[t]); report("treepm", rs, t -> ψtpm[t])
            GC.gc()
        end
    end
    close(out)
end
main()
