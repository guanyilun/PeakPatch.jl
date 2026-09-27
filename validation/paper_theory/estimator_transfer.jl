#!/usr/bin/env julia
# Estimator transfer test for measure_fieldmap_cl.jl: Gaussian skies drawn from the THEORY
# spectra (kSZ: linear Doppler + OV; κ: linear Limber z<4.5), passed through the identical
# geometric-octant apodized mask, mean subtraction and banding. Output = band-mean(estimated)/
# band-mean(input). No pixel window (alm2map applies none), so no pixwin correction here.
#   julia --project=validation -t 8 validation/paper_theory/estimator_transfer.jl [nsim]
using Healpix, Printf, Statistics, Random, DelimitedFiles
const NS, LMAX = 2048, 3200
const R = joinpath(@__DIR__, "results")
nsim = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6
include_str = read(joinpath(@__DIR__, "measure_fieldmap_cl.jl"), String)
# reuse smooth_map / octant_weights / pseudo_cl verbatim from the measurement script
for fn in ("smooth_map", "octant_weights", "pseudo_cl")
    i0 = findfirst("function $fn", include_str)[1]
    i1 = findnext("\nend", include_str, i0)[end]
    eval(Meta.parse(include_str[i0:i1]))
end

function loadcol(file, col)
    rows = filter(l -> !startswith(l, "#") && !isempty(strip(l)), readlines(joinpath(R, file)))
    M = [parse.(Float64, split(r)) for r in rows]
    [m[1] for m in M], [m[col] for m in M]
end
interp_log(x, xs, ys) = exp(begin
    lx = log(x); i = clamp(searchsortedlast(log.(xs), lx), 1, length(xs) - 1)
    t = (lx - log(xs[i])) / (log(xs[i+1]) - log(xs[i]))
    (1 - t) * log(ys[i]) + t * log(ys[i+1])
end)

# theory C_ℓ (kSZ in μK², κ dimensionless)
ld, dfull = loadcol("ksz_doppler.txt", 4); _, dapprox = loadcol("ksz_doppler.txt", 5)
dop = [isnan(a) ? b : a for (a, b) in zip(dfull, dapprox)]
lo, dov = loadcol("ksz_ov.txt", 2)
clz = zeros(LMAX + 1); clk = zeros(LMAX + 1)
th = readdlm(joinpath(R, "fieldmaps_vs_theory_prod.txt"); comments=true)   # l_eff, Cl_linOURS in col 3
for l in 2:LMAX
    d = interp_log(l, ld, dop) + (l >= lo[1] ? interp_log(l, lo, dov) : 0.0)
    clz[l+1] = d * 2π / (l * (l + 1))
    clk[l+1] = interp_log(clamp(l, th[1, 1], th[end, 1]), th[:, 1], th[:, 3])
end

function sim_map(cl, rng)
    alm = Alm(LMAX, LMAX, zeros(ComplexF64, numberOfAlms(LMAX, LMAX)))
    for m in 0:LMAX, l in m:LMAX
        s = sqrt(cl[l+1])
        alm.alm[almIndex(alm, l, m)] = m == 0 ? s * randn(rng) : s * (randn(rng) + im * randn(rng)) / sqrt(2)
    end
    alm2map(alm, NS)
end

oct = "000"
dummy = HealpixMap{Float64,RingOrder}(NS); dummy.pixels .= 1.0
wts = octant_weights(dummy, oct)
edges = unique(round.(Int, exp.(range(log(20), log(3000), length=26))))
bands = [(edges[i], edges[i+1] - 1) for i in 1:length(edges)-1]
bm(cl, a, b; D=false) = mean(D ? cl[l+1] * l * (l + 1) / 2π : cl[l+1] for l in a:b)
tz = zeros(nsim, length(bands)); tk = zeros(nsim, length(bands))
for s in 1:nsim
    rng = MersenneTwister(1000 + s)
    for (cl, T) in ((clz, tz), (clk, tk))
        est = pseudo_cl(sim_map(cl, rng), wts)
        for (ib, (a, b)) in enumerate(bands)
            T[s, ib] = bm(est, a, b) / bm(cl, a, b)
        end
    end
    println("sim $s done"); flush(stdout)
end
open(joinpath(R, "estimator_transfer.txt"), "w") do io
    println(io, "# estimator transfer (estimated/input band means), $nsim Gaussian sims, geometric oct000 mask")
    println(io, "# l_lo l_hi l_eff  T_ksz  sd_ksz  T_kappa  sd_kappa")
    for (ib, (a, b)) in enumerate(bands)
        @printf(io, "%d %d %.1f %.4f %.4f %.4f %.4f\n", a, b, sqrt(a * b),
                mean(tz[:, ib]), std(tz[:, ib]), mean(tk[:, ib]), std(tk[:, ib]))
    end
end
println(read(joinpath(R, "estimator_transfer.txt"), String))
