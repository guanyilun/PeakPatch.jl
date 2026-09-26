#!/usr/bin/env julia
# Production CIB maps for one frozen-campaign octant (docs/paper_comparison_plan_2026-09.md
# A1), generalizing ../run_cib.jl: XGPaint CIB_Planck2013 (= the Websky CIB model), sources
# generated ONCE, painted at every requested frequency.
#
#   julia --project=/home/yguan/work/XGPaint.jl -t N paint_cib.jl <OCT> <catalog_AMv2.pksc> <outdir> [nside] [freqs] [wcut|nocut] [tag]
#
# Catalog = finalized (Eulerian Mpc/h, km/s) AM catalog; observer from the octZYX bits.
# `wcut` (default) applies Websky's z-dependent completeness (halo_mass_completion.txt) —
# the Websky-comparison product; `nocut` keeps the Tinker-complete catalog.
# Pixel accumulation is done here SERIALLY: XGPaint's HEALPix paint! adds fluxes into
# shared pixels from Threads.@threads (a rare lost-update race), so we only use its
# fill_fluxes!. Maps are MJy/sr (fluxes are MJy; divide by Ω_pix once).
using XGPaint, Healpix, Printf, Statistics, DelimitedFiles
import Cosmology: comoving_radial_dist
import Unitful

const hub = 0.68
const rho_mh = 2.775e11 * 0.31
const MMIN_MSUN = 1.0e12                     # CIB model min_mass (as ../run_cib.jl)
const WREF = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref"

OCT = ARGS[1]; CAT = ARGS[2]; OUTD = ARGS[3]
nside = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 4096
freqs = length(ARGS) >= 5 ? parse.(Float64, split(ARGS[5], ",")) : [100.0, 143, 217, 353, 545, 857]
WCUT = !(length(ARGS) >= 6 && ARGS[6] == "nocut")
TAG = length(ARGS) >= 7 ? ARGS[7] : "prod_oct$(OCT)_AMv2"               # output name tag
length(OCT) == 3 && all(in("01"), OCT) || error("OCT must be three 0/1 bits (ZYX)")
obs = ntuple(d -> OCT[4-d] == '1' ? 2618.0 : -2618.0, 3)     # OCT = "ZYX", Mpc/h
mkpath(OUTD)

_mc = readdlm(joinpath(WREF, "halo_mass_completion.txt"); comments=true, comment_char='#')
const MC_Z = Float64.(_mc[:, 1]); const MC_M = Float64.(_mc[:, 2])
function mmin_of_z(z)
    i = clamp(searchsortedlast(MC_Z, z), 1, length(MC_Z) - 1)
    t = clamp((z - MC_Z[i]) / (MC_Z[i+1] - MC_Z[i]), 0.0, 1.0)
    MC_M[i] * (1 - t) + MC_M[i+1] * t
end

cosmo = get_cosmology(Float32; h=Float32(hub), OmegaM=0.31f0)
zg = collect(0.0:0.002:6.0)
chig = [Unitful.ustrip(comoving_radial_dist(cosmo, z)) * hub for z in zg]   # Mpc/h
function chi2z(chi)
    i = clamp(searchsortedlast(chig, chi), 1, length(zg) - 1)
    t = (chi - chig[i]) / (chig[i+1] - chig[i])
    zg[i] * (1 - t) + zg[i+1] * t
end
const CHIMAX = chig[searchsortedfirst(zg, 4.5)]              # z<4.5, as the other maps

function load_halos(path)
    nh = open(io -> Int(read(io, Int32)), path)
    nf = (filesize(path) - 12) ÷ (4nh)
    px = Float32[]; py = Float32[]; pz = Float32[]; mm = Float32[]
    open(path) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        chunk = 1_000_000; buf = Vector{Float32}(undef, chunk * nf); ndone = 0
        while ndone < nh
            m = min(chunk, nh - ndone); read!(io, view(buf, 1:m*nf))
            @inbounds for k in 1:m
                b = (k - 1) * nf
                M = 4 / 3 * π * rho_mh * Float64(buf[b+7])^3 / hub      # Msun (M200m proxy)
                M > MMIN_MSUN || continue
                x = Float64(buf[b+1]) - obs[1]; y = Float64(buf[b+2]) - obs[2]; z = Float64(buf[b+3]) - obs[3]
                r = sqrt(x^2 + y^2 + z^2)                               # Mpc/h
                (30.0 <= r <= CHIMAX) || continue
                WCUT && (M * hub > mmin_of_z(chi2z(r)) || continue)
                push!(px, Float32(x / hub)); push!(py, Float32(y / hub)); push!(pz, Float32(z / hub))
                push!(mm, Float32(M))
            end
            ndone += m
        end
    end
    pos = Matrix{Float32}(undef, 3, length(mm))
    pos[1, :] = px; pos[2, :] = py; pos[3, :] = pz
    pos, mm
end

t0 = time()
@info "paint_cib" OCT obs CAT nside freqs WCUT threads=Threads.nthreads()
halo_pos, halo_mass = load_halos(CAT)
@info "halos loaded" n=length(halo_mass) min=round((time() - t0) / 60; digits=1)
model = CIB_Planck2013{Float32}(nside=nside)
sources = generate_sources(model, cosmo, halo_pos, halo_mass)
halo_pos = nothing; halo_mass = nothing; GC.gc()
@info "sources" N_cen=sources.N_cen N_sat=sources.N_sat min=round((time() - t0) / 60; digits=1)

fc = Vector{Float32}(undef, sources.N_cen); fs = Vector{Float32}(undef, sources.N_sat)
acc = zeros(Float64, 12nside^2)
Ωpix = Healpix.nside2pixarea(nside)
tag = WCUT ? "wcut" : "nocut"
for ν in freqs
    XGPaint.fill_fluxes!(Float32(ν * 1e9), model, sources, fc, fs)
    fill!(acc, 0.0)
    @inbounds for i in 1:sources.N_cen; acc[sources.hp_ind_cen[i]] += fc[i]; end
    @inbounds for i in 1:sources.N_sat; acc[sources.hp_ind_sat[i]] += fs[i]; end
    tot = sum(Float64, fc) + sum(Float64, fs)
    abs(sum(acc) / tot - 1) < 1e-9 || error("flux bookkeeping failed at $ν GHz")
    m = HealpixMap{Float32,RingOrder}(nside); m.pixels .= Float32.(acc ./ Ωpix)
    out = joinpath(OUTD, "cib_nu$(lpad(Int(ν), 4, '0'))_$(tag)_$(TAG)_nside$(nside).fits")
    Healpix.saveToFITS(m, "!" * out, typechar="E")
    @printf("%4.0f GHz: total flux %.4e MJy, octant mean %.4e MJy/sr -> %s (%.1f min)\n",
            ν, tot, tot / (4π / 8), out, (time() - t0) / 60)
end
