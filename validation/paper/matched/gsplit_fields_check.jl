#!/usr/bin/env julia
# Field-level check of the production gaussian_split helpers (MultiResolution._gsplit_setup/_gsplit_tile, CPU)
# for δ, ψ1 and the potential −δ/k² used by the field maps (κ, ISW), vs the exact global fields from the same
# noise, and the original splice for comparison (MATCHED_FORTRAN_2026-10.md §9).
#   usage: julia --project=validation -t 32 gsplit_fields_check.jl <config.toml (periodic)>
using TOML, Printf, FFTW, Statistics
using PeakPatch
import PeakPatch.PowerSpectrum: load_pk
import PeakPatch.MultiResolution: _downsample_noise, _splice_compensation, _periodic_convolve!, _kernel_1lpt, _kernel_pot,
    _generate_extended_residual, _isolated_convolve, _isolated_convolve_dispatch, _interpolate_to_tile,
    _gsplit_setup, _gsplit_tile, _gaussian_split_rs

function main()
    cfgd = TOML.parsefile(ARGS[1]); cfg = PipelineConfig(cfgd); rc = cfgd["run"]
    ntile = rc["ntile"]; seed = rc["seed"]; nsub, N = grid_layout(cfg, ntile); @assert cfg.periodic_cores
    nmesh = cfg.n; nbuff = cfg.nbuff; a = cfg.boxsize / nmesh; L = N * a; Lb = nmesh * a
    M = ntile * rc["coarse_factor"]; block = N ÷ M; kNc = π / (block * a)
    pk = load_pk(cfg.pkfile)
    δg = generate_grf(N, pk, L, seed); δk = rfft(δg)
    ψg = displacements_1lpt(δk, N, L)[1]
    # exact potential −δ_k/k²
    kf = 2π / L; kx = FFTW.rfftfreq(N, N * kf); ky = FFTW.fftfreq(N, N * kf)
    Threads.@threads for iz in 1:N
        @inbounds for iy in 1:N, ix in 1:N÷2+1
            k2 = kx[ix]^2 + ky[iy]^2 + ky[iz]^2; δk[ix, iy, iz] = k2 == 0 ? 0 : δk[ix, iy, iz] * (-1 / k2)
        end
    end
    φg = irfft(δk, N); δk = nothing
    cn = _downsample_noise(N, M, seed); ck = rfft(cn)
    comp = cfg.coarse_compensation ? _splice_compensation(M, block) : nothing
    coarse(kern) = (q = copy(ck); _periodic_convolve!(q, pk, M, L; kernel_fn=kern, comp=comp); irfft(q, M))
    δc = coarse(nothing); ψc = coarse(_kernel_1lpt(1)); φc = coarse(_kernel_pot())
    rs = _gaussian_split_rs(cfg.gaussian_split_rs, block, nbuff)
    gs = _gsplit_setup(cn, pk, M, L, comp, rs, a; need_pot=true)
    out = open(joinpath(@__DIR__, "..", "results", "gsplit_fields_check_cf$(rc["coarse_factor"]).txt"), "w")
    say(x...) = begin s = string(x...); println(s); println(out, s) end
    say(@sprintf("N=%d nsub=%d nbuff=%d block=%d r_s=%.1f; k_N,coarse=%.3f; 2³ central tiles, cores", N, nsub, nbuff, block, rs, kNc))
    acc = Dict((v, f) => [0.0, 0.0, 0.0, 0.0] for v in (:orig, :gsplit), f in (:δ, :ψ1, :φ))   # Σr², Σe², Σe² band, Σr² band
    t0 = (ntile - 2) ÷ 2
    for t in [(i, j, k) for k in t0+1:t0+2, j in t0+1:t0+2, i in t0+1:t0+2]
        res = _generate_extended_residual(t..., nsub, nmesh, N, seed, cn, M, 0)
        δself = _isolated_convolve(res, pk, Lb, nmesh)
        orig = (δ=δself .+ _interpolate_to_tile(δc, t..., nsub, nmesh, N, M),
                ψ1=_isolated_convolve(res, pk, Lb, nmesh; kernel_fn=_kernel_1lpt(1)) .+ _interpolate_to_tile(ψc, t..., nsub, nmesh, N, M),
                φ=_isolated_convolve_dispatch(false, res, pk, Lb, nmesh, 5, 0, 0) .+ _interpolate_to_tile(φc, t..., nsub, nmesh, N, M))
        δt, ψt, φt = _gsplit_tile(gs, false, δself, cn, t..., nsub, nmesh, N, M, Lb; need_pot=true)
        gsp = (δ=δt, ψ1=ψt[1], φ=φt)
        core = nbuff+1:nbuff+nsub
        gi = (t[1]-1)*nsub+1:t[1]*nsub; gj = (t[2]-1)*nsub+1:t[2]*nsub; gk = (t[3]-1)*nsub+1:t[3]*nsub
        for (f, G) in ((:δ, δg), (:ψ1, ψg), (:φ, φg))
            r = Float64.(G[gi, gj, gk])
            for (v, F) in ((:orig, orig), (:gsplit, gsp))
                e = Float64.(getfield(F, f)[core, core, core]) .- r
                R = rfft(r); E = rfft(e); kq = 2π / (nsub * a); bE = 0.0; bR = 0.0
                for k3 in 1:nsub, k2 in 1:nsub, k1 in 1:nsub÷2+1
                    q2 = k2 <= nsub ÷ 2 + 1 ? k2 - 1 : k2 - 1 - nsub; q3 = k3 <= nsub ÷ 2 + 1 ? k3 - 1 : k3 - 1 - nsub
                    kk = kq * sqrt((k1 - 1)^2 + q2^2 + q3^2)
                    kNc <= kk < 2kNc && (bE += abs2(E[k1, k2, k3]); bR += abs2(R[k1, k2, k3]))
                end
                acc[(v, f)] .+= [sum(abs2, r), sum(abs2, e), bE, bR]
            end
        end
        GC.gc()
    end
    say("field   splice   rms err/rms   err power k∈[kN,2kN)")
    for f in (:δ, :ψ1, :φ), v in (:orig, :gsplit)
        q = acc[(v, f)]
        say(@sprintf("%-6s  %-7s  %10.4f    %10.2e", f, v, sqrt(q[2] / q[1]), q[3] / q[4]))
    end
    close(out)
end
main()
