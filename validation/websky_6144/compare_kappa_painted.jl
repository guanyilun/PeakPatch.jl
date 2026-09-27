#!/usr/bin/env julia
# Paint κ with XGPaint's NFWKappaProfile (fork branch lensing-kappa) from BOTH halo catalogs
# and compare against Websky's released kap.fits — the first actual paper-figure-level test.
#
# Run with:  julia --project=$HOME/work/XGPaint.jl -t 8 compare_kappa_painted.jl
#
# Three maps on the matched footprint (square of half-width θc/√2 inscribed in the shared
# 5° cap), painted/sampled identically at 0.25′ in a cap-centered rotated frame (CAR at the
# equator → uniform pixels):
#   1. ours      : oct000 finecell_AM halos, NFW c=7 profiles × Born CMB kernel
#   2. wsky_halo : Websky 10×10 patch halos, SAME painting
#   3. kap.fits  : Websky released map sampled on the SAME footprint in the WEBSKY frame —
#                  i.e. the SAME SKY as (2), so (3)−(2) is exactly their field component
#                  (+ unresolved halos + z>4.5 tail).
# Readouts: wsky_halo/kap (halo fraction of the full map, catalog-independent physics);
# ours/wsky_halo (catalog test; realization scatter applies — different sky, see
# KAPPA_2026-07-16.md). kap.fits carries the Nside=4096 pixel window (~0.86′), our painted
# maps a 0.25′ one — restrict conclusions to ℓ ≲ 2000.
using XGPaint, Pixell, Healpix, FFTW, Printf, Statistics, LinearAlgebra

const hub = 0.68
const rho_mh = 2.775e11 * 0.31            # Msun/h per (Mpc/h)³ (peak-patch convention)
const Mcut_h = 2e12                        # Msun/h, both catalogs complete above this
const RES_ARCMIN = 0.25
const Dr = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const WREF = "/home/yguan/scratch/websky_6144/websky_ref"

model = NFWKappaProfile(Omega_c=0.31 - 0.049, Omega_b=0.049, h=hub)
r2z = XGPaint.build_r2z_interpolator(0.0, 6.0, model.cosmo)   # physical Mpc → z
ortho_basis(a) = (t = abs(a[1]) < 0.9 ? [1.0,0,0] : [0.0,1,0]; e1 = normalize(cross(a,t)); e2 = cross(a,e1); (e1,e2))

# ---------- loaders → (α, δ, z, M_Msun) in the cap-centered frame ----------
function load_wsky(ax, e1, e2, s)
    io = open(joinpath(WREF, "halos_10x10.pksc"))
    Nw = read(io, Int32); read(io, Int32); read(io, Int32); nf = 10
    buf = Vector{Float32}(undef, Int(Nw)*nf); read!(io, buf); close(io)
    α = Float64[]; δ = Float64[]; zz = Float64[]; M = Float64[]
    for i in 1:Nw
        b = (i-1)*nf
        R = Float64(buf[b+7]) * hub                      # Mpc → Mpc/h
        Mh = 4/3 * π * rho_mh * R^3; Mh > Mcut_h || continue
        X = Float64(buf[b+1]); Y = Float64(buf[b+2]); Z = Float64(buf[b+3])   # Mpc (physical)
        r = sqrt(X^2 + Y^2 + Z^2); (30.0 <= r <= 7600.0) || continue
        na = (X*ax[1] + Y*ax[2] + Z*ax[3]) / r; na > 0 || continue
        a1 = atan((X*e1[1] + Y*e1[2] + Z*e1[3]) / r, na)
        d1 = asin(clamp((X*e2[1] + Y*e2[2] + Z*e2[3]) / r, -1, 1))
        (abs(a1) <= s && abs(d1) <= s) || continue
        push!(α, a1); push!(δ, d1); push!(zz, r2z(r)); push!(M, Mh / hub)     # Msun
    end
    (α, δ, zz, M)
end

function load_ours(ax, e1, e2, s; obs=-2618.0)
    α = Float64[]; δ = Float64[]; zz = Float64[]; M = Float64[]
    io = open(joinpath(Dr, "catalog_websky_6144_oct000_finecell_AM.pksc"))
    N = read(io, Int32); read(io, Float32); read(io, Float32)
    nf = 33; chunk = 4_000_000; buf = Vector{Float32}(undef, chunk*nf); ndone = 0
    while ndone < N
        m = min(chunk, Int(N) - ndone); read!(io, view(buf, 1:m*nf))
        @inbounds for k in 1:m
            b = (k-1)*nf
            R = Float64(buf[b+7]); R > 0 || continue
            Mh = 4/3 * π * rho_mh * R^3; Mh > Mcut_h || continue
            q1 = Float64(buf[b+1]); q2 = Float64(buf[b+2]); q3 = Float64(buf[b+3])   # Mpc/h
            rq = sqrt((q1-obs)^2 + (q2-obs)^2 + (q3-obs)^2)
            (30.0 <= rq <= 5170.0) || continue
            zq = r2z(rq / hub)
            a = 1.0 / (1.0 + zq)
            # old-convention catalog: Eulerian = q + d1·a − d2·a² (stored 2LPT carries +3/7)
            ex = q1 + Float64(buf[b+4])*a - Float64(buf[b+8])*a^2  - obs
            ey = q2 + Float64(buf[b+5])*a - Float64(buf[b+9])*a^2  - obs
            ez = q3 + Float64(buf[b+6])*a - Float64(buf[b+10])*a^2 - obs
            r = sqrt(ex^2 + ey^2 + ez^2); (30.0 <= r <= 5170.0) || continue    # Mpc/h
            na = (ex*ax[1] + ey*ax[2] + ez*ax[3]) / r; na > 0 || continue
            a1 = atan((ex*e1[1] + ey*e1[2] + ez*e1[3]) / r, na)
            d1 = asin(clamp((ex*e2[1] + ey*e2[2] + ez*e2[3]) / r, -1, 1))
            (abs(a1) <= s && abs(d1) <= s) || continue
            push!(α, a1); push!(δ, d1); push!(zz, r2z(r / hub)); push!(M, Mh / hub)
        end
        ndone += m
    end
    close(io); (α, δ, zz, M)
end

# ---------- footprint from the Websky patch (inscribed cap, as all previous scripts) ----------
@info "scanning websky patch for footprint..."
io = open(joinpath(WREF, "halos_10x10.pksc"))
Nw = read(io, Int32); read(io, Int32); read(io, Int32); nf = 10
buf = Vector{Float32}(undef, Int(Nw)*nf); read!(io, buf); close(io)
sx = 0.0; sy = 0.0; sz = 0.0
for i in 1:Nw; b = (i-1)*nf; global sx += buf[b+1]; global sy += buf[b+2]; global sz += buf[b+3]; end
wn = sqrt(sx^2 + sy^2 + sz^2); wax = [sx, sy, sz] ./ wn
cosmin = minimum((Float64(buf[(i-1)*nf+1])*wax[1] + Float64(buf[(i-1)*nf+2])*wax[2] +
                  Float64(buf[(i-1)*nf+3])*wax[3]) /
                 sqrt(Float64(buf[(i-1)*nf+1])^2 + Float64(buf[(i-1)*nf+2])^2 + Float64(buf[(i-1)*nf+3])^2)
                 for i in 1:Nw)
θc = acos(cosmin) / sqrt(2)
s = θc / sqrt(2)                                # square half-width, radians
@info "footprint" cap_deg=round(rad2deg(θc); digits=2) square_deg=round(2*rad2deg(s); digits=2)
we1, we2 = ortho_basis(wax)
oax = [1.0, 1.0, 1.0] ./ sqrt(3.0); oe1, oe2 = ortho_basis(oax)

# ---------- paint both catalogs ----------
sdeg = rad2deg(s)
box = [sdeg -sdeg; -sdeg sdeg] * Pixell.degree
shape, wcs = geometry(Pixell.CarClenshawCurtis{Float64}, box, RES_ARCMIN * Pixell.arcminute)
ws = profileworkspace(shape, wcs)
@info "map" shape=shape

@info "loading + painting websky halos..."
wα, wδ, wz, wM = load_wsky(wax, we1, we2, s)
mw = Enmap(zeros(shape), wcs); paint!(mw, ws, model, wM, wz, wα, wδ)
@info "websky painted" n=length(wM)

@info "loading + painting our catalog (26 GB stream)..."
oα, oδ, oz, oM = load_ours(oax, oe1, oe2, s)
mo = Enmap(zeros(shape), wcs); paint!(mo, ws, model, oM, oz, oα, oδ)
@info "ours painted" n=length(oM)

# ---------- sample kap.fits on the same footprint (WEBSKY frame → same sky as mw) ----------
@info "reading kap.fits..."
kap = Healpix.readMapFromFITS(joinpath(WREF, "kap.fits"), 1, Float32)
res = kap.resolution
αmap, δmap = posmap(shape, wcs)
mk = zeros(shape)
Threads.@threads for j in 1:size(mk, 2)
    for i in 1:size(mk, 1)
        ca = cos(αmap[i,j]); sa = sin(αmap[i,j]); cd = cos(δmap[i,j]); sd = sin(δmap[i,j])
        nx = cd*ca*wax[1] + cd*sa*we1[1] + sd*we2[1]
        ny = cd*ca*wax[2] + cd*sa*we1[2] + sd*we2[2]
        nz = cd*ca*wax[3] + cd*sa*we1[3] + sd*we2[3]
        mk[i,j] = kap.pixels[Healpix.vec2pixRing(res, nx, ny, nz)]
    end
end

# ---------- flat-sky C_ℓ ----------
function cl_map(A, L, ledges)
    m = A .- mean(A)
    F = fft(m); n1, n2 = size(A)
    f1 = fftfreq(n1, n1/L) .* 2π; f2 = fftfreq(n2, n2/L) .* 2π
    nb = length(ledges) - 1; S = zeros(nb); Nm = zeros(nb)
    for j in 1:n2, i in 1:n1
        l = sqrt(f1[i]^2 + f2[j]^2); (ledges[1] <= l < ledges[end]) || continue
        b = searchsortedlast(ledges, l); S[b] += abs2(F[i,j]); Nm[b] += 1
    end
    [Nm[b] > 0 ? S[b]/Nm[b] * L^2 / (n1*n2)^2 : 0.0 for b in 1:nb]
end

L = 2s
ledges = [10.0^l for l in range(log10(150.0), log10(8000.0); length=14)]
lc = [sqrt(ledges[i]*ledges[i+1]) for i in 1:length(ledges)-1]
co = cl_map(Array(parent(mo)), L, ledges)
cw = cl_map(Array(parent(mw)), L, ledges)
ck = cl_map(mk, L, ledges)

@printf("\nmap rms κ: ours=%.4f  wsky_halo=%.4f  kap.fits=%.4f\n",
        std(parent(mo)), std(parent(mw)), std(mk))
@printf("%-7s %-11s %-11s %-11s %-11s %-11s\n","ell","Cl_ours","Cl_wskyhalo","Cl_kapfits","ours/whalo","whalo/kap")
for i in eachindex(lc)
    (co[i] > 0 && cw[i] > 0 && ck[i] > 0) || continue
    @printf("%-7.0f %-11.3e %-11.3e %-11.3e %-11.3f %-11.3f\n",
            lc[i], co[i], cw[i], ck[i], co[i]/cw[i], cw[i]/ck[i])
end
sel = findall(l -> 300 <= l <= 2000, lc)
@printf("\nmeans over 300<=ell<=2000: ours/wsky_halo=%.3f   wsky_halo/kap.fits=%.3f\n",
        mean(co[sel]./cw[sel]), mean(cw[sel]./ck[sel]))
@printf("wsky_halo and kap.fits are the SAME SKY: their ratio is the halo fraction of the map\n")
@printf("(deficit = field component + halos below the cut/completeness + z>4.5 tail).\n")
@printf("ours vs wsky_halo is the catalog test, subject to different-realization scatter at high ell.\n")
