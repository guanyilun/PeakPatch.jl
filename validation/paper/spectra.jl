# Full-sky angular power spectra for the paper comparisons (paper_comparison_plan A2/A3).
# include() this file. Conventions (A3, adopted 2026-09-26):
#   * full-sky maps, no mask; map2alm with niter=0 (both maps same Nside -> quadrature and
#     pixel-window effects cancel in ours/ref ratios; compare at the COARSER Nside when
#     they differ, e.g. the Nside-2048 released tSZ map).
#   * bins: linear Δℓ=2 at ℓ<10, then logarithmic Δℓ/ℓ ≈ 0.1, unit weight per mode.
#   * errors: Gaussian/Knox per realization; octant jackknife (apodized octant masks,
#     f_sky-corrected pseudo-C_ℓ) for non-Gaussian products; the larger is used.
using Healpix, Printf, Statistics

function loadmap(path; scale=1.0, nside=0)
    m = Healpix.readMapFromFITS(path, 1, Float64)
    scale != 1 && (m.pixels .*= scale)
    if nside > 0 && m.resolution.nside != nside
        m = Healpix.udgrade(m, nside)
    end
    m
end

almof(m::HealpixMap, lmax) = Healpix.map2alm(m; lmax=lmax, mmax=lmax, niter=0)
clof(a, b=a) = Healpix.alm2cl(a, b)                       # index ℓ+1

function lbins(lmax; lmin=2, dlog=0.1)
    edges = lmin < 10 ? collect(lmin:2:10) : [lmin]
    while edges[end] < lmax
        push!(edges, max(edges[end] + 2, round(Int, edges[end] * (1 + dlog))))
    end
    edges[end] = lmax + 1
    edges                                                   # bin b = [edges[b], edges[b+1])
end

# mode-weighted bin average of C_ℓ (ℓ from index), effective ℓ, number of modes
function binned(cl, edges)
    nb = length(edges) - 1
    cb = zeros(nb); le = zeros(nb); nm = zeros(nb)
    for b in 1:nb, l in edges[b]:edges[b+1]-1
        l + 1 > length(cl) && continue
        w = 2l + 1
        cb[b] += w * cl[l+1]; le[b] += w * l; nm[b] += w
    end
    cb ./ max.(nm, 1), le ./ max.(nm, 1), nm
end

# Gaussian error of a binned auto spectrum (per realization): σ = C·√(2/Nmodes)
gauss_auto(cb, nm) = cb .* sqrt.(2 ./ max.(nm, 1))
# Knox error of a binned cross spectrum
knox_cross(cab, caa, cbb, nm) = sqrt.((cab .^ 2 .+ caa .* cbb) ./ max.(nm, 1))

# ------------------------------- octant jackknife ------------------------------------
# apodized mask of the SKY octant seen from observer bits OCT ("ZYX"; bit 1 = observer at
# +2618 → that axis looks toward −), cosine taper of width
# θap [rad] measured from the three bounding great circles (angular distance asin|x_i|)
function octant_mask(nside, oct::AbstractString; θap=deg2rad(2.0))
    res = Healpix.Resolution(nside); npix = 12nside^2
    sgn = ntuple(d -> oct[4-d] == '1' ? -1.0 : 1.0, 3)
    w = zeros(npix)
    Threads.@threads for p in 1:npix
        v = Healpix.pix2vecRing(res, p)
        inside = all(d -> sgn[d] * v[d] > 0, 1:3)
        inside || continue
        dmin = minimum(asin(min(abs(v[d]), 1.0)) for d in 1:3)
        w[p] = dmin >= θap ? 1.0 : 0.5 - 0.5cos(π * dmin / θap)
    end
    m = HealpixMap{Float64,RingOrder}(nside); m.pixels .= w
    m
end

const OCTS = ["000", "001", "010", "011", "100", "101", "110", "111"]

# per-octant binned pseudo-C_ℓ / ⟨w²⟩ for a set of (label => map) and pairs; returns
# Dict(pair => nb×8 matrix)
function octant_spectra(maps::Dict, pairs, lmax, edges; θap=deg2rad(2.0))
    nside = first(values(maps)).resolution.nside
    out = Dict(p => zeros(length(edges) - 1, 8) for p in pairs)
    for (io, oct) in enumerate(OCTS)
        w = octant_mask(nside, oct; θap=θap)
        w2 = mean(abs2, w.pixels)
        alms = Dict{String,Any}()
        for k in unique(vcat([collect(p) for p in pairs]...))
            mm = HealpixMap{Float64,RingOrder}(nside); mm.pixels .= maps[k].pixels .* w.pixels
            alms[k] = almof(mm, lmax)
        end
        for p in pairs
            out[p][:, io] = binned(clof(alms[p[1]], alms[p[2]]) ./ w2, edges)[1]
        end
    end
    out
end
jk_sigma(M) = vec(std(M; dims=2)) ./ sqrt(size(M, 2))       # error of the full-sky mean

# Gaussian pixel-window approximation, for THEORY comparisons only (ratios at equal Nside
# need none): w_ℓ² ≈ exp(−ℓ(ℓ+1)σ²), σ = √(4π/Npix)/√12
pixwin2_gauss(l, nside) = exp(-l * (l + 1) * (4π / (12nside^2)) / 12)

function write_table(path, header, cols...)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", header)
        for i in eachindex(cols[1])
            println(io, join((@sprintf("%.6e", c[i]) for c in cols), " "))
        end
    end
    @info "wrote $path"
end
