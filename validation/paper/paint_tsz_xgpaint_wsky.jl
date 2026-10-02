#!/usr/bin/env julia
# Independent-painter test: paint Websky's OWN halos (halos_10x10.pksc) with upstream XGPaint's
# Battaglia16ThermalSZProfile (Battaglia+12 AGN fit, the same parameters as pks2map bbps model 1)
# at Nside 2048, for comparison with the released tsz_2048.fits and with our port of the Fortran
# physics (check_tsz_painter_wsky.jl / compare_tsz_three_painters.jl).
# Two mass conventions:
#   asis : the Websky (M200m top-hat) mass passed as "M200c" — the convention our painter and the
#          Fortran maptable use (y amplitude ∝ M, radius from 200 ρ_crit)
#   nfw7 : M200c from M200m assuming NFW with c = 7 w.r.t. r200m (Stein+2020 §3.1.1)
#   usage: julia --project=/home/yguan/work/XGPaint.jl -t 32 paint_tsz_xgpaint_wsky.jl <outdir>
using XGPaint, Healpix, Printf
const H = 0.68; const OM = 0.31; const OB = 0.049
const WSKY = "/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/halos_10x10.pksc"
outdir = ARGS[1]; mkpath(outdir)

# halos: positions Mpc (observer at origin), R_TH Mpc; M = 4π/3 ρ_m R³ in Msun
io = open(WSKY); N = read(io, Int32); read(io, Int32); read(io, Int32)
b = Vector{Float32}(undef, Int(N) * 10); read!(io, b); close(io); b = reshape(b, 10, :)
rho_m = 2.775e11 * OM * H^2                                     # Msun / Mpc³
x = Float64.(b[1, :]); y = Float64.(b[2, :]); z = Float64.(b[3, :])
r = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
M200m = 4π / 3 * rho_m .* Float64.(b[7, :]) .^ 3
ra = atan.(y, x); dec = asin.(z ./ r)

# comoving distance (Mpc) -> z, flat ΛCDM
zg = collect(range(0, 6; length=60001)); E(zz) = sqrt(OM * (1 + zz)^3 + 1 - OM)
chig = zeros(length(zg)); for i in 2:length(zg); chig[i] = chig[i-1] + 299792.458 / (100H) * 0.5 * (1 / E(zg[i]) + 1 / E(zg[i-1])) * (zg[i] - zg[i-1]); end
zof(c) = begin i = searchsortedlast(chig, c); i >= length(chig) ? zg[end] : zg[i] + (c - chig[i]) / (chig[i+1] - chig[i]) * (zg[i+1] - zg[i]) end
red = zof.(r)

# M200m -> M200c for NFW with c200m = 7 (r_s = r200m / 7): find r with mean density 200 ρ_crit(z)
μ(q) = log(1 + q) - q / (1 + q)
function m200c_of(Mm, zz)
    rhoc = 2.775e11 * H^2 * E(zz)^2                              # physical Msun/Mpc³
    rhom = rhoc * OM * (1 + zz)^3 / E(zz)^2
    r200m = cbrt(3Mm / (4π * 200rhom)); rs = r200m / 7
    f(rr) = 3 * Mm * μ(rr / rs) / μ(7.0) / (4π * rr^3) - 200rhoc   # decreasing in rr
    lo, hi = 0.05r200m, r200m
    for _ in 1:80; mid = (lo + hi) / 2; f(mid) > 0 ? (lo = mid) : (hi = mid); end
    rr = (lo + hi) / 2; Mm * μ(rr / rs) / μ(7.0)
end
M200c = m200c_of.(M200m, red)
@printf("halos %d; z %.3f–%.3f; M200c/M200m median %.3f\n", N, minimum(red), maximum(red), sort(M200c ./ M200m)[end÷2])

model = Battaglia16ThermalSZProfile(Omega_c=OM - OB, Omega_b=OB, h=H)
interp = build_interpolator(model; cache_file=joinpath(outdir, "b16_websky_cosmo.jld2"), overwrite=false)
res = Healpix.Resolution(2048)
for (lab, M) in (("asis", M200m), ("nfw7", M200c))
    o = sortperm(dec); m = HealpixMap{Float64,RingOrder}(2048); w = HealpixRingProfileWorkspace{Float64}(res)
    @time paint!(m, w, interp, M[o], red[o], ra[o], dec[o])
    out = joinpath(outdir, "tsz_y_xgpaint_$(lab)_wskypatch_nside2048.fits")
    Healpix.saveToFITS(m, "!" * out, typechar="D"); @info "wrote" out
end
