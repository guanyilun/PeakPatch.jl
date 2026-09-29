#!/usr/bin/env julia
# High-mass tail after abundance matching (tSZ C_ℓ deficit, V2_RESULTS_2026-09-27.md item 1).
# Our AM maps rank k in each Δz=0.1 bin to M(N_Tinker=k) on a 10^4-bin grid → the top
# halo of each bin is pinned at M(N=1). Compare v2 oct000 RAW vs AM vs the Tinker
# expectation of N(>M) for the most massive halos, z<1, one octant.
using PeakPatch, Printf
import PeakPatch.Cosmology: CosmologyParams, chi, build_chi_to_z, chi_to_z, growth_factor
import PeakPatch.MassFunction: precompute_sigma, tinker_dndlnM
const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
const Om = 0.31; const rho_m = 2.775e11 * Om
pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "..", "websky_6144", "data", "pk_websky_RAW_unnormalized.dat"))
chi2z = build_chi_to_z(COSMO; z_max=3.0)
const D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const OCT = get(ENV, "OCT", "000")
# octZYX bit = 1 -> observer at +2618 on that axis
const obsv = (OCT[3] == '1' ? 2618.0 : -2618.0, OCT[2] == '1' ? 2618.0 : -2618.0, OCT[1] == '1' ? 2618.0 : -2618.0)

function tail(path; Mcut=3e14, zmax=1.0)
    out = Tuple{Float64,Float64}[]
    nh = open(io -> Int(read(io, Int32)), path)
    open(path) do io
        read(io, Int32); read(io, Float32); read(io, Float32)
        chunk = 2_000_000; bf = Vector{Float32}(undef, chunk * 33); nd = 0
        while nd < nh
            m = min(chunk, nh - nd); read!(io, view(bf, 1:m*33))
            for k in 1:m
                b = (k - 1) * 33
                M = 4 / 3 * π * rho_m * Float64(bf[b+7])^3
                M > Mcut || continue
                r = sqrt((bf[b+1] - obsv[1])^2 + (bf[b+2] - obsv[2])^2 + (bf[b+3] - obsv[3])^2)
                z = chi_to_z(chi2z, r); z < zmax && push!(out, (M, z))
            end
            nd += m
        end
    end
    out
end
const TAG = get(ENV, "TAG", "v2")   # v2 (nbuff 16) | v3test (nbuff 25)
const S3 = "/home/yguan/scratch/websky_6144/catalogs_v3"      # v3 raw and v3fs AM catalogs live on scratch
const RAWD = TAG in ("v3", "v3fs") ? S3 : D
raw = tail(joinpath(RAWD, "catalog_websky_6144_$(TAG == "v3fs" ? "v3" : TAG)_oct$(OCT).pksc"))
am = tail(TAG == "v3fs" ? joinpath(S3, "catalog_websky_6144_v3_oct$(OCT)_AMfs.pksc") :
          joinpath(D, "catalog_websky_6144_$(TAG)_oct$(OCT)_AM.pksc"))

Mg = 10 .^ range(12, 16; length=801); lnMg = log.(Mg); sg = precompute_sigma(Mg, pk, Om)
function Ntinker(Mcut, za, zb)
    ra, rb = chi(za, COSMO), chi(zb, COSMO); ns = 60; acc = 0.0
    for j in 1:ns
        r0 = ra + (j - 1) * (rb - ra) / ns; r1 = ra + j * (rb - ra) / ns
        z = chi_to_z(chi2z, (r0 + r1) / 2); Dg = growth_factor(z, COSMO)
        dV = (1 / 8) * (4π / 3) * (r1^3 - r0^3); n = 0.0
        for i in 1:length(Mg)-1
            Mm = sqrt(Mg[i] * Mg[i+1]); Mm < Mcut && continue
            dls = (log(sg[i+1]) - log(sg[i])) / (lnMg[i+1] - lnMg[i])
            n += tinker_dndlnM(Mm, Dg * sqrt(sg[i] * sg[i+1]), dls, z, Om) * (lnMg[i+1] - lnMg[i])
        end
        acc += n * dV
    end
    acc
end
cnt(v, M, za, zb) = count(t -> t[1] > M && za <= t[2] < zb, v)
@printf("[%s] oct%s (1/8 sky), N(>M) in z bins: RAW | AM | Tinker expectation\n", TAG, OCT)
for (za, zb) in ((0.0, 0.25), (0.25, 0.5), (0.5, 1.0), (0.0, 1.0))
    for M in (3e14, 5e14, 1e15, 2e15)
        @printf("z %.2f-%.2f  M>%.0e: %6d | %6d | %8.2f\n", za, zb, M, cnt(raw, M, za, zb), cnt(am, M, za, zb), Ntinker(M, za, zb))
    end
end
s(v, za, zb) = sum((t[1] / 1e14)^(10 / 3) for t in v if za <= t[2] < zb; init=0.0)
@printf("\nΣ(M/1e14)^{10/3} over M>3e14: z<0.5 RAW %.1f AM %.1f (AM/RAW %.3f); z<1 RAW %.1f AM %.1f\n",
        s(raw, 0, 0.5), s(am, 0, 0.5), s(am, 0, 0.5) / s(raw, 0, 0.5), s(raw, 0, 1), s(am, 0, 1))
@printf("max M: RAW %.3e  AM %.3e\n", maximum(first, raw), maximum(first, am))
# pile-up check at the R_TH clamp: raw-mass histogram near the top (R ∝ M^{1/3})
rth(M) = cbrt(3M / (4π * rho_m))
@printf("\nraw R_TH distribution of M>3e14, z<1 (pile-up at the cap shows as an excess in the top bin):\n")
Rs = sort(rth.(first.(raw)))
for Rlo in 14.0:1.0:22.0
    @printf("  R_TH in [%.0f,%.0f) Mpc/h: %6d\n", Rlo, Rlo + 1, count(r -> Rlo <= r < Rlo + 1, Rs))
end
@printf("max raw R_TH %.2f Mpc/h\n", maximum(Rs))
