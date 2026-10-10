#!/usr/bin/env julia
# Low-ℓ κ: is the field-κ excess over linear theory at ℓ ≈ 60–115 (v5: +3–11% per bin) the sky variance of our
# one periodic realization? kappa_realization_check.jl tested only the 3D mode amplitudes |δ_k|² (factor 1.002);
# the full-sky C_ℓ of one realization also depends on the phases (cross terms between 3D modes), which only a map
# captures. Here a linear Born κ map is ray-traced from the production coarse noise (seed 12345, the exact long
# modes of the run, observer at the box corner as for all 8 octants) and from NENS independent white-noise seeds
# on the same grid, the same box and the same code:
#   R_ℓ = C_ℓ(seed 12345) / ⟨C_ℓ(other seeds)⟩   = our realization's linear-sky factor (grid smoothing cancels)
#   r_ℓ = corr(toy 12345, production field κ)     ≈ 1 checks that the toy is the same sky
#   field/lin                                      from the production field map (Nside 4096 → 256) and Limber
# If R_ℓ ≈ field/lin bin by bin, the low-ℓ excess is our realization, not the pipeline.
#   env: WS_ROOT, NENS (default 8), NSIDE (default 256), NOISE_CACHE (raw Float32 coarse noise)
using PeakPatch, FFTW, Healpix, Printf, Statistics, Random, DelimitedFiles
import PeakPatch: CosmologyParams, chi, load_pk
import PeakPatch.Cosmology: Dlinear_tables, Dlinear_ab

const N = 6144; const M = 512; const L = 5236.0; const SEED = 12345
const DC = L / M                                                     # coarse cell [Mpc/h]
const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
const CHISTAR = chi(1089.0, COSMO)                                   # FieldMap default (no radiation)
const ZMAX = 4.5; const DCHI = 5.0
const NENS = parse(Int, get(ENV, "NENS", "8"))
const NSIDE = parse(Int, get(ENV, "NSIDE", "256")); const LMAX = 2NSIDE
const WS = get(ENV, "WS_ROOT", "/mnt/ceph/users/yguan/projects/uoft/peakpatch/websky_6144")
const CACHE = get(ENV, "NOISE_CACHE", joinpath(WS, "scratch", "coarse_noise_seed$(SEED)_M$(M).f32"))
const OUT = joinpath(@__DIR__, "results", "kappa_linear_toy.txt")
include(joinpath(@__DIR__, "spectra.jl"))

function coarse_noise()
    if isfile(CACHE)
        a = Array{Float32,3}(undef, M, M, M); read!(CACHE, a); return a
    end
    a = PeakPatch.MultiResolution._downsample_noise(N, M, SEED)
    mkpath(dirname(CACHE)); write(CACHE, a); a
end

# linear δ(z=0) on the coarse grid, the run's own convolution (no Gaussian split, no compensation)
function delta0(noise, pk)
    q = rfft(Float64.(noise)); PeakPatch.MultiResolution._periodic_convolve!(q, pk, M, L); irfft(q, M)
end

# Born κ = Σ W(χ)Δχ D(z) δ0(χ n̂), trilinear on the periodic grid; observer at the box corner (cell centers at
# (I − ½)·DC from it). The 8 production octants all see this one point of the periodic box.
function kappamap(d)
    zs = range(0.0, ZMAX, length=20001); cs = [chi(z, COSMO) for z in zs]
    cmax = cs[end]; nstep = floor(Int, cmax / DCHI)
    gt = Dlinear_tables(COSMO)                                       # the growth FieldMap uses
    w = zeros(nstep); χs = [(i - 0.5) * DCHI for i in 1:nstep]
    for i in 1:nstep
        j = clamp(searchsortedlast(cs, χs[i]), 1, length(cs) - 1)
        z = zs[j] + (zs[j+1] - zs[j]) * (χs[i] - cs[j]) / (cs[j+1] - cs[j])
        w[i] = 1.5 * COSMO.Om / 2997.92458^2 * (1 + z) * χs[i] * (1 - χs[i] / CHISTAR) * Dlinear_ab(1 / (1 + z), gt)[1] * DCHI
    end
    res = Healpix.Resolution(NSIDE); m = HealpixMap{Float64,RingOrder}(NSIDE)
    Threads.@threads for p in 1:12NSIDE^2
        θ, ϕ = Healpix.pix2angRing(res, p)
        nx, ny, nz = sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)
        acc = 0.0
        @inbounds for i in 1:nstep
            gx = χs[i] * nx / DC - 0.5; gy = χs[i] * ny / DC - 0.5; gz = χs[i] * nz / DC - 0.5
            fx, fy, fz = floor(gx), floor(gy), floor(gz); tx, ty, tz = gx - fx, gy - fy, gz - fz
            i0 = mod(Int(fx), M) + 1; j0 = mod(Int(fy), M) + 1; k0 = mod(Int(fz), M) + 1
            i1 = i0 == M ? 1 : i0 + 1; j1 = j0 == M ? 1 : j0 + 1; k1 = k0 == M ? 1 : k0 + 1
            v = (1 - tz) * ((1 - ty) * ((1 - tx) * d[i0, j0, k0] + tx * d[i1, j0, k0]) + ty * ((1 - tx) * d[i0, j1, k0] + tx * d[i1, j1, k0])) +
                tz * ((1 - ty) * ((1 - tx) * d[i0, j0, k1] + tx * d[i1, j0, k1]) + ty * ((1 - tx) * d[i0, j1, k1] + tx * d[i1, j1, k1]))
            acc += w[i] * v
        end
        m.pixels[p] = acc
    end
    m
end

function main()
    io = open(OUT, "w"); say(a...) = (local line = string(a...); println(line); println(io, line))
    pk = load_pk(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky.dat"))
    edges = lbins(LMAX; lmin=20)
    cl(a, b=a) = binned(clof(a, b), edges)[1]
    t = time(); k0 = kappamap(delta0(coarse_noise(), pk)); a0 = almof(k0, LMAX)
    @info "toy seed $SEED" secs = round(time() - t)
    c0 = cl(a0)
    ens = Vector{Vector{Float64}}()
    for s in 1:NENS
        n = randn(Xoshiro(1000 + s), Float32, M, M, M)
        push!(ens, cl(almof(kappamap(delta0(n, pk)), LMAX))); @info "ensemble seed" s
    end
    E = reduce(hcat, ens); μ = vec(mean(E; dims=2)); σμ = vec(std(E; dims=2)) ./ sqrt(NENS)
    f = loadmap(joinpath(WS, "fullsky_v5fs", "kappa_field_v5fs_fullsky_nside4096.fits"); nside=NSIDE)
    af = almof(f, LMAX); cf = cl(af); cx = cl(a0, af)
    pw = [mean(pixwin2_gauss(l, NSIDE) for l in edges[b]:edges[b+1]-1) for b in 1:length(edges)-1]
    kl = readdlm(joinpath(@__DIR__, "..", "paper_theory", "results", "kappa_limber.txt"); comments=true)
    lin(l) = (i = clamp(searchsortedlast(kl[:, 1], l), 1, size(kl, 1) - 1); u = (l - kl[i, 1]) / (kl[i+1, 1] - kl[i, 1]);
              kl[i, 11] * (1 - u) + kl[i+1, 11] * u)                     # linOURS_0-4.5 (χ* no radiation)
    tl = [sum((2l + 1) * lin(l) for l in edges[b]:edges[b+1]-1) / sum(2l + 1 for l in edges[b]:edges[b+1]-1) for b in 1:length(edges)-1]
    nm = [sum(2l + 1 for l in edges[b]:edges[b+1]-1) for b in 1:length(edges)-1]
    say("# linear Born κ toy (z < $ZMAX) from the production coarse noise, seed $SEED, M = $M, vs $NENS other seeds; Nside $NSIDE")
    say("# R = C(toy $SEED)/<C(ensemble)> (±err of the mean); ens/lin = <C ensemble>/Limber (grid smoothing);")
    say("# field/lin = production field κ (Nside 4096→$NSIDE, Gaussian pixwin removed) / Limber linear; r = corr(toy, field); σG = sqrt(2/Nmodes)")
    say(@sprintf("%-7s %7s %7s %7s %9s %7s %6s", "ell", "R", "±", "ens/lin", "field/lin", "r", "σG"))
    for b in eachindex(c0)
        ℓ = (edges[b] + edges[b+1] - 1) / 2; ℓ > 0.6LMAX && break
        say(@sprintf("%-7.1f %7.3f %7.3f %7.3f %9.3f %7.3f %6.3f", ℓ, c0[b] / μ[b], c0[b] / μ[b] * σμ[b] / μ[b],
                     μ[b] / tl[b], cf[b] / pw[b] / tl[b], cx[b] / sqrt(c0[b] * cf[b]), sqrt(2 / nm[b])))
    end
    for (lo, hi) in ((50, 150), (60, 115), (150, 300))
        bs = [b for b in eachindex(c0) if edges[b] >= lo && edges[b+1] - 1 <= hi]
        wsum(x) = sum(x[b] * nm[b] for b in bs) / sum(nm[b] for b in bs)
        say(@sprintf("band ℓ %d-%d: R = %.3f   field/lin = %.3f   r = %.3f", lo, hi, wsum(c0 ./ μ), wsum(cf ./ pw ./ tl),
                     wsum(cx ./ sqrt.(c0 .* cf))))
    end
    close(io)
end
main()
