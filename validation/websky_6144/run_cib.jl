#!/usr/bin/env julia
# CIB painting via XGPaint (user's fork, ~/work/XGPaint.jl) — no reinvention:
# CIB_Planck2013 IS the Websky CIB model (Shang HOD + greybody; Td0=20.7, α=0.2
# z-dependent dust temperature, η=2.4, β=1.6 = the 08-OCT-2019 Websky update).
# We feed our octant AM catalog (Eulerian positions -> Mpc, masses -> Msun) and
# compare mean intensity + whole-octant pseudo-C_ell against cib_nu*.fits.
#
# Run with the XGPaint environment:
#   julia --project=/home/yguan/work/XGPaint.jl -t N run_cib.jl [freqGHz] [nside] [nmax]
# nmax (testing): cap on halos streamed (0 = all).
using XGPaint, Healpix, Printf, Statistics
import Cosmology: comoving_radial_dist
import Unitful

const CAT = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_finecell_AM.pksc"
const WREF = "/home/yguan/scratch/websky_6144/websky_ref"
const OUTD = "/home/yguan/scratch/websky_6144/fieldmaps"
const hub = 0.68
const rho_mh = 2.775e11 * 0.31
const OBS = -2618.0
const MMIN_MSUN = 1.0e12          # CIB model min_mass

freq_ghz = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 545.0
nside = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4096
nmax = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0

cosmo = get_cosmology(Float32; h=Float32(hub), OmegaM=0.31f0)

# chi[Mpc/h] -> a interpolation (for the old-convention Eulerian reconstruction)
zg = collect(0.0:0.002:6.0)
chig = [Unitful.ustrip(comoving_radial_dist(cosmo, z)) * hub for z in zg]   # Mpc/h
function chi2a(chi)
    i = clamp(searchsortedlast(chig, chi), 1, length(zg) - 1)
    t = (chi - chig[i]) / (chig[i+1] - chig[i])
    z = zg[i] * (1 - t) + zg[i+1] * t
    1.0 / (1.0 + z)
end

# stream catalog -> (3,N) positions [Mpc, observer at origin] + masses [Msun]
function load_halos(path, nmax)
    nh_tot = open(path) do io Int(read(io, Int32)) end
    nh = nmax > 0 ? min(nmax, nh_tot) : nh_tot
    px = Float32[]; py = Float32[]; pz = Float32[]; mm = Float32[]
    sizehint!(px, nh ÷ 2); sizehint!(py, nh ÷ 2); sizehint!(pz, nh ÷ 2); sizehint!(mm, nh ÷ 2)
    open(path) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        nf = 33; chunk = 1_000_000; buf = Vector{Float32}(undef, chunk * nf); ndone = 0
        while ndone < nh
            m = min(chunk, nh - ndone); read!(io, view(buf, 1:m*nf))
            @inbounds for k in 1:m
                b = (k - 1) * nf
                R = Float64(buf[b+7])
                M = 4 / 3 * π * rho_mh * R^3 / hub          # Msun (M200m proxy)
                M > MMIN_MSUN || continue
                q1 = Float64(buf[b+1]); q2 = Float64(buf[b+2]); q3 = Float64(buf[b+3])
                rq = sqrt((q1 - OBS)^2 + (q2 - OBS)^2 + (q3 - OBS)^2)
                a = chi2a(rq)
                x = (q1 + Float64(buf[b+4]) * a - Float64(buf[b+8]) * a^2 - OBS) / hub
                y = (q2 + Float64(buf[b+5]) * a - Float64(buf[b+9]) * a^2 - OBS) / hub
                z = (q3 + Float64(buf[b+6]) * a - Float64(buf[b+10]) * a^2 - OBS) / hub
                r = sqrt(x^2 + y^2 + z^2)
                (40.0 <= r <= 7600.0) || continue           # Mpc; z<~4.6
                push!(px, Float32(x)); push!(py, Float32(y)); push!(pz, Float32(z))
                push!(mm, Float32(M))
            end
            ndone += m
        end
    end
    pos = Matrix{Float32}(undef, 3, length(mm))
    pos[1, :] = px; pos[2, :] = py; pos[3, :] = pz
    pos, mm
end

@info "streaming catalog..." CAT nmax
t0 = time()
halo_pos, halo_mass = load_halos(CAT, nmax)
@info "halos loaded" n = length(halo_mass) minM = minimum(halo_mass) maxM = maximum(halo_mass) min = round(time() - t0; digits=1)

model = CIB_Planck2013{Float32}(nside=nside)
@info "generating sources (centrals + satellites)..."
t0 = time()
sources = generate_sources(model, cosmo, halo_pos, halo_mass)
@info "sources" N_cen = sources.N_cen N_sat = sources.N_sat min = round((time() - t0) / 60; digits=1)

m = HealpixMap{Float32,RingOrder}(nside)
@info "painting at $(freq_ghz) GHz..."
t0 = time()
paint!(m, Float32(freq_ghz * 1e9), model, sources)
@info "painted" min = round((time() - t0) / 60; digits=1)

# XGPaint healpix paint adds source FLUXES (Jy) per pixel -> intensity MJy/sr
omega_pix = 4π / nside2npix(nside)
m.pixels ./= Float32(omega_pix * 1e6)

out = joinpath(OUTD, "cib_nu$(lpad(Int(freq_ghz),4,'0'))_websky_6144_oct000_AM_nside$(nside).fits")
Healpix.saveToFITS(m, "!" * out, typechar="E")
@info "wrote $out"

# ---- comparison vs reference (octant mean; full-sky ref mean) ----
refpath = joinpath(WREF, "cib_nu$(lpad(Int(freq_ghz),4,'0')).fits")
if isfile(refpath) && nmax == 0
    ref = Healpix.readMapFromFITS(refpath, 1, Float32)
    inoct = m.pixels .> 0
    @printf("octant mean intensity: ours %.4e MJy/sr;  ref (same pixels) %.4e;  ref full-sky %.4e\n",
            mean(m.pixels[inoct]), mean(ref.pixels[inoct]), mean(ref.pixels))
    @printf("ratio ours/ref (same pixels): %.3f\n", mean(m.pixels[inoct]) / mean(ref.pixels[inoct]))
else
    @info "reference not present or test mode — skipping comparison" refpath
end
