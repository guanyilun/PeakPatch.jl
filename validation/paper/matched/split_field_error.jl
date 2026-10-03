#!/usr/bin/env julia
# Field-level error of the multires split (MATCHED_FORTRAN_2026-10.md §8). This rebuilds, with the CPU branches
# of run_multitile_split, the tile fields that the production path uses for peak finding and shell analysis
# (δ, ψ1, ψ2) and compares them in the tile cores with the exact global fields from the same Threefry noise.
#   * correlation and rms error per field
#   * error power and cross-correlation r(k) vs k (core cube FFT)
#   * rms error vs distance from the core edge (missing residual from outside the tile)
#   * the same with an extended residual shell (nshell > 0) — is the missing outer residual the cause?
#   usage: julia --project=validation -t 32 split_field_error.jl <config.toml> <nshell list, e.g. 0,24> [ntiles_side]
using TOML, Printf, FFTW, Statistics
using PeakPatch
import PeakPatch.PowerSpectrum: load_pk
import PeakPatch.MultiResolution: _downsample_noise, _splice_compensation, _periodic_convolve!, _kernel_1lpt,
    _kernel_phi_ij, _kernel_2lpt, _generate_extended_residual, _isolated_convolve, _interpolate_to_tile,
    _apply_kernel_inplace!

function tile2lpt(δt, nmesh, Lb)
    dk = rfft(δt); src = δt .^ 2 .* 0.5f0
    for d in 1:3
        p = copy(dk); _apply_kernel_inplace!(p, nmesh, Lb, _kernel_phi_ij(d, d); zero_nyquist=false); src .-= irfft(p, nmesh) .^ 2 .* 0.5f0
    end
    for (a, b) in ((1, 2), (1, 3), (2, 3))
        p = copy(dk); _apply_kernel_inplace!(p, nmesh, Lb, _kernel_phi_ij(a, b); zero_nyquist=false); src .-= irfft(p, nmesh) .^ 2
    end
    sk = rfft(src); q = copy(sk); _apply_kernel_inplace!(q, nmesh, Lb, _kernel_2lpt(1)); irfft(q, nmesh)
end

function main()
    cfgd = TOML.parsefile(ARGS[1]); cfg = PipelineConfig(cfgd); rc = cfgd["run"]
    shells = parse.(Int, split(ARGS[2], ","))
    nside = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2
    ntile = rc["ntile"]; nmesh = cfg.n; nbuff = cfg.nbuff; seed = rc["seed"]
    nsub, N = grid_layout(cfg, ntile); @assert cfg.periodic_cores
    a = cfg.boxsize / nmesh; L = N * a; Lb = nmesh * a
    M = ntile * rc["coarse_factor"]; block = N ÷ M
    pk = load_pk(cfg.pkfile)
    @info "split field error" N nsub nmesh nbuff M block shells
    # exact global fields
    δg = generate_grf(N, pk, L, seed); δgk = rfft(δg)
    ψg = displacements_1lpt(δgk, N, L)[1]; ψ2g = displacements_2lpt(δgk, N, L)[1]; δgk = nothing
    # coarse fields as in run_multitile_split
    cn = _downsample_noise(N, M, seed); ck = rfft(cn)
    comp = cfg.coarse_compensation ? _splice_compensation(M, block) : nothing
    dck = copy(ck); _periodic_convolve!(dck, pk, M, L; comp=comp); δc = irfft(dck, M)
    pck = copy(ck); _periodic_convolve!(pck, pk, M, L; kernel_fn=_kernel_1lpt(1), comp=comp); ψc = irfft(pck, M)
    out = open(joinpath(@__DIR__, "..", "results", "split_field_error_cf$(rc["coarse_factor"]).txt"), "w")
    say(x...) = begin s = string(x...); println(s); println(out, s) end
    say(@sprintf("N=%d nsub=%d nmesh=%d nbuff=%d M=%d block=%d, cell %.4f Mpc/h; %d³ central tiles", N, nsub, nmesh, nbuff, M, block, a, nside))
    kny_c = π / (block * a)
    say(@sprintf("coarse Nyquist k = %.3f h/Mpc (λ = %.1f Mpc/h)", kny_c, 2π / kny_c))
    t0 = (ntile - nside) ÷ 2
    tiles = [(i, j, k) for k in t0+1:t0+nside, j in t0+1:t0+nside, i in t0+1:t0+nside]
    nb = 14; kmax = π / a; kedges = 10 .^ range(log10(2π / (nsub * a)), log10(kmax); length=nb + 1)
    dedges = [0, 3, 6, 10, 15, 20, 30, 50, 80, 1000]
    for ns in shells
        acc = Dict(f => (zeros(nb), zeros(nb), zeros(nb), zeros(Int, nb)) for f in (:δ, :ψ1, :ψ2))   # Pref, Perr, Pcross, n
        dacc = Dict(f => (zeros(length(dedges) - 1), zeros(length(dedges) - 1), zeros(Int, length(dedges) - 1)) for f in (:δ, :ψ1, :ψ2))
        tot = Dict(f => [0.0, 0.0, 0.0] for f in (:δ, :ψ1, :ψ2))   # Σref², Σerr², Σref·spl
        for (it, jt, kt) in tiles
            res = _generate_extended_residual(it, jt, kt, nsub, nmesh, N, seed, cn, M, ns)
            δt = _isolated_convolve(res, pk, Lb, nmesh; nshell=ns) .+ _interpolate_to_tile(δc, it, jt, kt, nsub, nmesh, N, M)
            ψt = _isolated_convolve(res, pk, Lb, nmesh; kernel_fn=_kernel_1lpt(1), nshell=ns) .+ _interpolate_to_tile(ψc, it, jt, kt, nsub, nmesh, N, M)
            ψ2t = tile2lpt(δt, nmesh, Lb)
            core = nbuff+1:nbuff+nsub
            gi = (it-1)*nsub+1:it*nsub; gj = (jt-1)*nsub+1:jt*nsub; gk = (kt-1)*nsub+1:kt*nsub
            for (f, T, G) in ((:δ, δt, δg), (:ψ1, ψt, ψg), (:ψ2, ψ2t, ψ2g))
                r = Float64.(G[gi, gj, gk]); s = Float64.(T[core, core, core]); e = s .- r
                tot[f][1] += sum(abs2, r); tot[f][2] += sum(abs2, e); tot[f][3] += sum(r .* s)
                R = rfft(r); E = rfft(e); Pr, Pe, Px, nn = acc[f]
                kf = 2π / (nsub * a)
                for k3 in 1:nsub, k2 in 1:nsub, k1 in 1:nsub÷2+1
                    q2 = k2 <= nsub ÷ 2 + 1 ? k2 - 1 : k2 - 1 - nsub; q3 = k3 <= nsub ÷ 2 + 1 ? k3 - 1 : k3 - 1 - nsub
                    kk = kf * sqrt((k1 - 1)^2 + q2^2 + q3^2); b = searchsortedlast(kedges, kk)
                    1 <= b <= nb || continue
                    Pr[b] += abs2(R[k1, k2, k3]); Pe[b] += abs2(E[k1, k2, k3]); Px[b] += real(R[k1, k2, k3] * conj(R[k1, k2, k3] + E[k1, k2, k3])); nn[b] += 1
                end
                De, Dr, Dn = dacc[f]
                for z in 1:nsub, y in 1:nsub, x in 1:nsub
                    d = min(x, y, z, nsub + 1 - x, nsub + 1 - y, nsub + 1 - z) - 1
                    b = searchsortedlast(dedges, d); De[b] += e[x, y, z]^2; Dr[b] += r[x, y, z]^2; Dn[b] += 1
                end
            end
            GC.gc()
        end
        say("\n== nshell = $ns (residual shell beyond the tile)")
        for f in (:δ, :ψ1, :ψ2)
            Σr, Σe, Σx = tot[f]
            say(@sprintf("  %-3s rms error / rms field %.4f   corr %.6f   amplitude Σrs/Σr² %.4f", f, sqrt(Σe / Σr),
                         Σx / sqrt(Σr * (Σr + Σe + 2 * (Σx - Σr))), Σx / Σr))
        end
        say("  k [h/Mpc]   " * join([@sprintf("%-22s", "$(f): Perr/Pref  r(k)") for f in (:δ, :ψ1, :ψ2)]))
        for b in 1:nb
            kc = sqrt(kedges[b] * kedges[b+1])
            say(@sprintf("  %8.4f   ", kc) * join([begin Pr, Pe, Px, nn = acc[f]; nn[b] == 0 ? @sprintf("%-22s", "-") :
                @sprintf("%-10.4f %-11.5f", Pe[b] / Pr[b], Px[b] / sqrt(Pr[b] * (Pr[b] + Pe[b] + 2 * (Px[b] - Pr[b])))) end for f in (:δ, :ψ1, :ψ2)]))
        end
        say("  distance from core edge [cells]: rms error / rms field")
        for b in 1:length(dedges)-1
            say(@sprintf("  %3d–%-4d  ", dedges[b], dedges[b+1]) * join([begin De, Dr, Dn = dacc[f]; Dn[b] == 0 ? "   -    " : @sprintf("%s %.4f   ", f, sqrt(De[b] / Dr[b])) end for f in (:δ, :ψ1, :ψ2)]))
        end
    end
    close(out)
end
main()
