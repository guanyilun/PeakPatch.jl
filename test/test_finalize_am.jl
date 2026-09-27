using PeakPatch
using Test
using Random
import PeakPatch.Cosmology: CosmologyParams, Dlinear_tables, Dlinear_ab, chi, build_chi_to_z, chi_to_z
import PeakPatch.MassFunction: precompute_sigma, tinker_dndlnM

# finalize_eulerian and abundance matching are the two post-processing steps where
# the Websky validation found real bugs (D/a vs D, 2LPT sign, Om_total, missing
# Eulerian conversion). These tests check them against independent physics, not
# against their own formulas.

const COSMO = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)   # 1st arg = Om TOTAL

@testset "finalize_eulerian" begin
    # Linear theory at z=0: v = H0·f·ψ = 100·f km/s per Mpc/h of displacement, with the
    # growth rate f ≈ Ω_m^0.55 (Linder 2005; ≲1% for ΛCDM).
    f0 = 0.31^0.55
    q = (10.0, -20.0, 30.0); d1 = (1.0, -2.0, 0.5); d2 = (0.1, 0.05, -0.2)
    h = HaloRecord(Float32.(q)..., Float32.(d1)..., 3f0, Float32.(d2)..., 0.7f0)
    out = finalize_eulerian([h], COSMO, (0.0, 0.0, 0.0); ievol=0, z_out=0.0)[1]

    @testset "Eulerian position = q + ψ₁ + ψ₂" begin
        @test out.x ≈ q[1] + d1[1] + d2[1] rtol=1e-6
        @test out.y ≈ q[2] + d1[2] + d2[2] rtol=1e-6
        @test out.z ≈ q[3] + d1[3] + d2[3] rtol=1e-6
    end
    @testset "z=0 velocity = 100·f·(ψ₁ + 2ψ₂) km/s, f ≈ Ωm^0.55" begin
        @test out.vx ≈ 100 * f0 * (d1[1] + 2 * d2[1]) rtol=0.01
        @test out.vy ≈ 100 * f0 * (d1[2] + 2 * d2[2]) rtol=0.01
        @test out.vz ≈ 100 * f0 * (d1[3] + 2 * d2[3]) rtol=0.02
    end
    @testset "bookkeeping" begin
        @test (out.vx2, out.vy2, out.vz2) == (0f0, 0f0, 0f0)   # 2LPT consumed
        @test out.RTHL == h.RTHL && out.overdensity == h.overdensity
        z = finalize_eulerian([HaloRecord(1f0, 2f0, 3f0, 0f0, 0f0, 0f0, 1f0, 0f0, 0f0, 0f0, 0f0)],
                              COSMO, (0.0, 0.0, 0.0); ievol=0)[1]
        @test (z.x, z.y, z.z, z.vx, z.vy, z.vz) == (1f0, 2f0, 3f0, 0f0, 0f0, 0f0)
        @test isempty(finalize_eulerian(HaloRecord[], COSMO, (0.0, 0.0, 0.0)))
    end
    @testset "lightcone: redshift from Eulerian distance to the observer" begin
        # Halo at comoving distance χ(z=1) along x from an off-origin observer: the
        # velocity factor must be a·H(a)·f(a) at z=1, i.e. 100·E(1)·f(1)/2.
        obs = (-500.0, 100.0, 50.0)
        r1 = chi(1.0, COSMO)
        hz = HaloRecord(Float32(obs[1] + r1), Float32(obs[2]), Float32(obs[3]),
                        1f0, 0f0, 0f0, 3f0, 0f0, 0f0, 0f0, 0f0)
        o = finalize_eulerian([hz], COSMO, obs; ievol=1)[1]
        a = 0.5; E = sqrt(0.31 / a^3 + 0.69); Om_a = 0.31 / a^3 / E^2
        @test o.vx ≈ a * 100 * E * Om_a^0.55 * (1 + 0) rtol=0.01   # ψ = 1 Mpc/h (+ Eulerian shift is 1 Mpc/h: negligible Δz)
        # and it must differ from the snapshot (z_out=0) answer by the growth of a·H·f
        o0 = finalize_eulerian([hz], COSMO, obs; ievol=0, z_out=0.0)[1]
        @test o.vx / o0.vx ≈ a * E * Om_a^0.55 / f0 rtol=0.01
    end
    @testset "ExtHaloRecord keeps its extra fields" begin
        e = ExtHaloRecord(Float32.(q)..., Float32.(d1)..., 3f0, Float32.(d2)..., 0.7f0,
                          Float32.(1:22)...)
        eo = finalize_eulerian([e], COSMO, (0.0, 0.0, 0.0); ievol=0)[1]
        @test eo.x ≈ out.x && eo.vx ≈ out.vx
        @test (eo.e_v, eo.Rf, eo.gradrf_z) == (e.e_v, e.Rf, e.gradrf_z)
    end
end

@testset "abundance_match" begin
    # Synthetic octant (fsky = 1/8) catalog in two redshift bins. After AM the k-th most
    # massive halo of a bin must sit where the INDEPENDENTLY integrated Tinker target
    # gives N(>M) = k, i.e. AM reproduces the target cumulative counts. Halos are put on
    # the outer side of each bin centre so the (linear-in-z) table lookup is exactly
    # that bin's mapping.
    pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "..", "validation", "websky_6144",
                                                  "data", "pk_websky_RAW_unnormalized.dat"))
    zlo, zhi, nz = 0.2, 0.4, 2
    zedge = range(zlo, zhi; length=nz + 1); zc = [(zedge[i] + zedge[i+1]) / 2 for i in 1:nz]
    rng = MersenneTwister(1)
    Om = COSMO.Om; rho_m = 2.775e11 * Om
    halos = HaloRecord[]; binof = Int[]
    for ib in 1:nz
        za, zb = ib == 1 ? (zedge[1], zc[1]) : (zc[2], zedge[3])    # outer half of the bin
        ra, rb = chi(za, COSMO), chi(zb, COSMO)
        for _ in 1:3000
            r = cbrt(ra^3 + rand(rng) * (rb^3 - ra^3))                # uniform in volume
            u = abs.(randn(rng, 3)); u ./= sqrt(sum(abs2, u))          # +++ octant direction
            M = 10^(12.5 + 2.5 * rand(rng)); R = cbrt(3M / (4π * rho_m))
            push!(halos, HaloRecord(Float32.(r .* u)..., 0f0, 0f0, 0f0, Float32(R), 0f0, 0f0, 0f0, 0f0))
            push!(binof, ib)
        end
    end
    table = build_abundance_table(halos, COSMO, pk; hmf=:tinker, z_min=zlo, z_max=zhi,
                                  nzbins=nz, nMbins=4000, Mmin=1e12, Mmax=1e16,
                                  obs=(0.0, 0.0, 0.0), fsky=1/8)
    am = abundance_match(halos, table, COSMO; obs=(0.0, 0.0, 0.0))
    mass(h) = (4π / 3) * rho_m * Float64(h.RTHL)^3

    # Independent target: fsky·∫dV ∫_M dn/dlnM dlnM with its own σ(M) and shell sums.
    Mg = 10 .^ range(12, 16; length=801); lnMg = log.(Mg); sg = precompute_sigma(Mg, pk, Om)
    chi2z = build_chi_to_z(COSMO; z_max=2.0)
    function Ntarget(Mcut, ib)
        ra, rb = chi(zedge[ib], COSMO), chi(zedge[ib+1], COSMO); ns = 40; acc = 0.0
        for j in 1:ns
            r0 = ra + (j - 1) * (rb - ra) / ns; r1 = ra + j * (rb - ra) / ns
            z = chi_to_z(chi2z, (r0 + r1) / 2); D = PeakPatch.Cosmology.growth_factor(z, COSMO)
            dV = (1 / 8) * (4π / 3) * (r1^3 - r0^3)
            n = 0.0
            for i in 1:length(Mg)-1
                Mm = sqrt(Mg[i] * Mg[i+1]); Mm < Mcut && continue
                dls = (log(sg[i+1]) - log(sg[i])) / (lnMg[i+1] - lnMg[i])
                n += tinker_dndlnM(Mm, D * sqrt(sg[i] * sg[i+1]), dls, z, Om) * (lnMg[i+1] - lnMg[i])
            end
            acc += n * dV
        end
        return acc
    end

    for ib in 1:nz
        idx = findall(==(ib), binof)
        Mraw = mass.(halos[idx]); Mam = mass.(am[idx])
        @testset "bin $ib (z $(zedge[ib])–$(zedge[ib+1]))" begin
            # rank-preserving: AM mass is non-decreasing in raw mass (ties from Float32
            # RTHL rounding allowed). The top halo of each z-bin used to map BELOW the
            # next ones (identity default above N_PP=0; fixed in build_abundance_table).
            @test issorted(Mam[sortperm(Mraw)])
            Ms = sort(Mam; rev=true)
            for k in (30, 300, 3000)
                # N_target(>M_k) = k up to table/bin discretization (few %)
                @test Ntarget(Ms[k], ib) ≈ k rtol=0.06
            end
        end
    end

    @testset "fsky scales the target volume" begin
        # 8× the target volume at fixed catalog counts → each rank maps to a larger mass
        t1 = build_abundance_table(halos, COSMO, pk; hmf=:tinker, z_min=zlo, z_max=zhi, nzbins=nz,
                                   nMbins=4000, Mmin=1e12, Mmax=1e16, obs=(0.0, 0.0, 0.0), fsky=1.0)
        am1 = abundance_match(halos, t1, COSMO; obs=(0.0, 0.0, 0.0))
        @test all(mass(a1) > mass(a8) for (a1, a8) in zip(am1, am))
    end
end
