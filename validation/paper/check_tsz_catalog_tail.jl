#!/usr/bin/env julia
# tSZ deficit, catalog side (V2_RESULTS_2026-09-27.md open item 1): the painter reproduces
# Websky's map from their halos (check_tsz_painter_wsky.jl), so compare the high-mass,
# low-z tail — N(>M) per steradian in z bins — between Websky's 10°×10° patch and our v2
# oct000 AM catalog (1/8 sky), plus the tSZ-weighted moment Σ mh^{5/3}, Σ mh^{10/3}
# (the mean-y and 1-halo proxies) per steradian.
using Printf, Statistics
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z
const hub = 0.68; const rho_mh = 2.775e11 * 0.31
cosmo = CosmologyParams(0.31, 0.049, 0.69, hub, 0.965, 0.81)
chi2z = build_chi_to_z(cosmo; z_max=6.0)
zb = [0.0, 0.1, 0.25, 0.5, 1.0]; Ms = [1e14, 2e14, 5e14, 1e15]          # Msun/h (M200m)

function tally!(C, S5, S10, M, z)
    b = searchsortedlast(zb, z); (1 <= b < length(zb)) || return
    for (j, m) in enumerate(Ms); M > m && (C[b, j] += 1); end
    M > 1e13 && (S5[b] += (M / 1e14)^(5 / 3); S10[b] += (M / 1e14)^(10 / 3))
end

# Websky patch: 10 floats, positions Mpc (observer at origin), R Mpc
io = open("/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref/halos_10x10.pksc")
Nw = read(io, Int32); read(io, Int32); read(io, Int32)
buf = Vector{Float32}(undef, Int(Nw) * 10); read!(io, buf); close(io)
Cw = zeros(4, 4); S5w = zeros(4); S10w = zeros(4)
cmin = 1.0
for i in 1:Nw
    b = (i - 1) * 10
    x, y, z = Float64(buf[b+1]), Float64(buf[b+2]), Float64(buf[b+3])
    r = sqrt(x^2 + y^2 + z^2) * hub
    global cmin = min(cmin, x * hub / r)
    M = 4 / 3 * π * rho_mh * (Float64(buf[b+7]) * hub)^3
    tally!(Cw, S5w, S10w, M, chi_to_z(chi2z, r))
end
# patch solid angle: |y|,|z| ≤ 5° square in gnomonic-ish coords ≈ (10°)² (flat approx)
Ωw = deg2rad(10.0)^2
@printf("Websky patch: %d halos, min cos to +x = %.4f (half-angle %.2f°), Ω ≈ %.4f sr\n", Nw, cmin, acosd(cmin), Ωw)

# ours: v2 oct000 AM (33 floats, Eulerian Mpc/h, observer (-2618)^3)
path = "/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_v2_oct000_AM.pksc"
Co = zeros(4, 4); S5o = zeros(4); S10o = zeros(4)
nh = open(io -> Int(read(io, Int32)), path)
open(path) do io
    read(io, Int32); read(io, Float32); read(io, Float32)
    chunk = 2_000_000; bf = Vector{Float32}(undef, chunk * 33); nd = 0
    while nd < nh
        m = min(chunk, nh - nd); read!(io, view(bf, 1:m*33))
        for k in 1:m
            b = (k - 1) * 33
            x = bf[b+1] + 2618.0; y = bf[b+2] + 2618.0; z = bf[b+3] + 2618.0
            r = sqrt(x^2 + y^2 + z^2); r < 3000 || continue
            tally!(Co, S5o, S10o, 4 / 3 * π * rho_mh * Float64(bf[b+7])^3, chi_to_z(chi2z, r))
        end
        nd += m
    end
end
Ωo = 4π / 8
@printf("\nN(>M) per sr: ours (oct000, v2 AM) | Websky patch | ratio ours/websky ± Poisson(websky)\n")
for b in 1:length(zb)-1, (j, m) in enumerate(Ms)
    no = Co[b, j] / Ωo; nw = Cw[b, j] / Ωw
    @printf("z %.2f-%.2f  M>%.0e: %9.1f | %9.1f  (n_w=%4d)  ratio %.3f ± %.3f\n", zb[b], zb[b+1], m,
            no, nw, Cw[b, j], no / nw, Cw[b, j] > 0 ? no / nw / sqrt(Cw[b, j]) : NaN)
end
@printf("\ntSZ moments per sr (M>1e13):  Σ(M/1e14)^{5/3}  ours|websky|ratio    Σ(M/1e14)^{10/3} ours|websky|ratio\n")
for b in 1:length(zb)-1
    @printf("z %.2f-%.2f   %9.1f | %9.1f | %.3f      %9.1f | %9.1f | %.3f\n", zb[b], zb[b+1],
            S5o[b] / Ωo, S5w[b] / Ωw, (S5o[b] / Ωo) / (S5w[b] / Ωw), S10o[b] / Ωo, S10w[b] / Ωw, (S10o[b] / Ωo) / (S10w[b] / Ωw))
end
