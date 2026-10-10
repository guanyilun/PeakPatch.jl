#!/usr/bin/env julia
# Is the low-ℓ κ excess (field κ 4% above linear theory at ℓ = 50–150, v5) a property of our one realization?
# Our full sky is ONE periodic 5236 Mpc/h box seen from a corner by 8 octants (same seed), so its large-scale
# modes are few and shared across octants. This measures the realized-to-expected power R(k) of the production
# white noise (seed 12345, N = 6144) on the coarse grid, whose low-k modes are exactly the realization's,
# and predicts the realization-only κ ratio C_ℓ(with R) / C_ℓ(R = 1) with linear Limber over the z < 4.5 kernel
# (χ* = 14.2 Gpc, as the field maps). No GPU. env: M (coarse grid, default 512 = production)
using PeakPatch, FFTW, Printf, Statistics
import PeakPatch: CosmologyParams, chi, growth_factor, load_pk

const N = 6144; const L = 5236.0; const SEED = 12345
const M = parse(Int, get(ENV, "M", "512"))
const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
const CHISTAR = 14200.0 * 0.68
const OUT = joinpath(@__DIR__, "results", "kappa_realization_check.txt")

function main()
    t = time()
    noise = PeakPatch.MultiResolution._downsample_noise(N, M, SEED)            # unit-variance white noise per coarse cell
    @info "noise" M var = var(noise) secs = round(time() - t)
    nk = rfft(Float64.(noise))
    kf = 2π / L
    kx = FFTW.rfftfreq(M, M * kf); ky = FFTW.fftfreq(M, M * kf)
    edges = 10 .^ range(log10(kf * 0.999), log10(0.25), length=31)
    s = zeros(length(edges) - 1); c = zeros(Int, length(edges) - 1)
    for iz in 1:M, iy in 1:M, ix in 1:size(nk, 1)
        k = sqrt(kx[ix]^2 + ky[iy]^2 + ky[iz]^2); k == 0 && continue
        b = searchsortedlast(edges, k); 1 <= b < length(edges) || continue
        s[b] += abs2(nk[ix, iy, iz]) / M^3; c[b] += 1               # E = 1 per mode for white noise
    end
    kc = sqrt.(edges[1:end-1] .* edges[2:end]); R = s ./ max.(c, 1)
    io = open(OUT, "w"); say(a...) = (local line = string(a...); println(line); println(io, line))
    say("# realized/expected power of the production white noise (seed $SEED, N = $N, coarse M = $M, box $L Mpc/h)")
    say(@sprintf("%-10s %8s %10s %10s", "k[h/Mpc]", "modes", "R(k)", "1/sqrt(n)"))
    for b in eachindex(kc)
        c[b] > 0 && say(@sprintf("%-10.4f %8d %10.4f %10.4f", kc[b], c[b], R[b], 1 / sqrt(c[b])))
    end
    # linear Limber prediction of the realization factor for κ (z < 4.5)
    Rof(k) = (b = searchsortedlast(edges, k); 1 <= b < length(edges) && c[b] > 0 ? R[b] : 1.0)
    pk = load_pk(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky.dat"))
    zs = collect(range(0.005, 4.5, length=900)); χs = [chi(z, COSMO) for z in zs]; Ds = [growth_factor(z, COSMO) for z in zs]
    W(i) = (1 + zs[i]) * χs[i] * (1 - χs[i] / CHISTAR)
    say("\n# κ realization factor C_ℓ(R)/C_ℓ(R=1), linear Limber over z < 4.5")
    for (lo, hi) in ((50, 150), (150, 300), (300, 600), (600, 1000))
        num = 0.0; den = 0.0
        for l in lo:hi, i in 2:length(zs)
            dχ = χs[i] - χs[i-1]; k = (l + 0.5) / χs[i]; w = (2l + 1) * dχ * W(i)^2 / χs[i]^2 * Ds[i]^2 * pk(k)
            num += w * Rof(k); den += w
        end
        say(@sprintf("ell %4d-%-4d  predicted realization factor %.4f", lo, hi, num / den))
    end
    close(io)
end
main()
