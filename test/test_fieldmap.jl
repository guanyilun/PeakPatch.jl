using PeakPatch
using Test
using Random
import Healpix
import PeakPatch.Cosmology: CosmologyParams, build_chi_to_z, chi_to_z, chi

# CPU checks of run_multitile_fieldmap on a small lightcone box (from
# validation/websky_6144/test_fieldmap_smoke.jl, which also covers the GPU paths):
# exact bookkeeping identities plus analytic anchors for the mean κ and τ.

@testset "FieldMap" begin
    datadir = joinpath(@__DIR__, "..", "validation", "websky_6144", "data")
    alatt = 8.0; nmesh = 64; nbuff = 8; ntile = 2
    nsub = nmesh - 2nbuff; N = nsub * ntile + 2nbuff; Nc = nsub * ntile      # N = 112
    obs = (-384.0, -384.0, -384.0)                                            # grid corner
    config = Dict{String,Any}(
        "cosmology" => Dict{String,Any}("Om" => 0.31, "OB" => 0.049, "OL" => 0.69, "h" => 0.68),
        "grid" => Dict{String,Any}("n" => nmesh, "boxsize" => nmesh * alatt, "nbuff" => nbuff,
                                   "cenx" => obs[1], "ceny" => obs[2], "cenz" => obs[3]),
        "run" => Dict{String,Any}("ievol" => 1, "z_max" => 0.1, "z_out" => 0.0, "ilpt" => 2,
                                  "ioutshear" => 0),
        "files" => Dict{String,Any}("pk" => joinpath(datadir, "pk_websky.dat"),
                                    "filterbank" => joinpath(datadir, "filters_websky.dat"),
                                    "homeltab" => joinpath(datadir, "HomelTab_websky.dat"),
                                    "output" => joinpath(mktempdir(), "fieldmap_test.pksc")),
    )
    cfg = PipelineConfig(config)
    nside = 64; res = Healpix.Resolution(nside); npix = Healpix.nside2npix(nside)
    v2p = (x, y, z) -> Healpix.vec2pixRing(res, x, y, z)
    cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.808)
    chimax = chi(0.1, cosmo); rho_m = 2.775e11 * 0.31; rmin_eff = 2 * alatt
    cellpos(i) = (i - (Nc + 1) / 2) * alatt
    rdist(x, y, z) = sqrt((x - obs[1])^2 + (y - obs[2])^2 + (z - obs[3])^2)

    @testset "ang2pix_ring (GPU painter) == Healpix.jl" begin
        rng = MersenneTwister(7); nbad = 0
        for _ in 1:200_000
            x, y, z = randn(rng), randn(rng), randn(rng)
            PeakPatch.MultiResolution.ang2pix_ring(nside, x, y, z) == Healpix.vec2pixRing(res, x, y, z) || (nbad += 1)
        end
        @test nbad == 0
    end

    K = [:kappa, :mass, :tau, :ksz, :isw]
    run(; kw...) = run_multitile_fieldmap(cfg; ntile=ntile, seed=12345, npix=npix, vec2pix=v2p,
                                          use_gpu=false, verbose=false, kw...)
    maps  = run(; kernels=K, subdiv_max=1)
    maps3 = run(; kernels=K, subdiv_max=3)

    @testset "mass conservation (exact cell bookkeeping)" begin
        ncell = count(rmin_eff <= rdist(cellpos(i), cellpos(j), cellpos(k)) <= chimax
                      for k in 1:Nc, j in 1:Nc, i in 1:Nc)
        @test sum(maps[:mass]) ≈ rho_m * alatt^3 * ncell rtol=1e-6
        # sub-cell splitting truncates boundary sub-cells: ~1% low, never more
        @test 0.98 < sum(maps3[:mass]) / sum(maps[:mass]) <= 1.0
    end

    @testset "mean κ and τ vs analytic line-of-sight integrals" begin
        chi2z = build_chi_to_z(cosmo; z_max=1.0); chistar = chi(1089.0, cosmo)
        nint = 2000; dr = (chimax - rmin_eff) / nint; kacc = 0.0; tacc = 0.0
        sigT_ne0 = 6.65246e-29 * 11.2299 * 0.9 * 0.049 * 0.68^2 * 3.0857e22 / 0.68
        for ii in 0:nint-1
            r = rmin_eff + (ii + 0.5) * dr; z = chi_to_z(chi2z, r)
            kacc += 1.5 * 0.31 * (1 / 2997.92458)^2 * (1 + z) * (1 - r / chistar) * r * dr
            x_e = z < 3 ? (1 - 0.245 / 2) : (1 - 3 * 0.245 / 4)
            tacc += sigT_ne0 * x_e * (1 + z)^2 * dr
        end
        # octant of sky; the painted map is grid-quantized at the χ boundaries
        @test sum(maps[:kappa]) / npix ≈ kacc / 8 rtol=0.005
        @test sum(maps[:tau]) / npix ≈ tacc / 8 rtol=0.005
    end

    @testset "kSZ / ISW sanity" begin
        cov = findall(!iszero, maps[:tau])
        kszm = sum(maps[:ksz][cov]) / length(cov); kszr = sqrt(sum(abs2, maps[:ksz][cov]) / length(cov))
        @test abs(kszm) / kszr < 0.5          # signed map: mean ≪ rms (0.24 at this tiny volume)
        iswr = sqrt(sum(abs2, maps[:isw][cov]) / length(cov))
        @test isfinite(iswr) && iswr > 0
    end

    @testset "octant containment" begin
        on = findall(>(0), maps[:mass])
        bad = count(p -> any(<(-0.15), Healpix.pix2vecRing(res, p)), on)
        @test bad <= 0.001 * length(on)       # boundary pixels only
        @test 0.10 < length(on) / npix < 0.15 # ≈ 1/8 of the sky
    end

    @testset "multi-worker dispatch == single worker" begin
        # Painted totals are ψ-independent → sums must agree to summation roundoff. Guards
        # the aliased-accumulator bug that painted exactly 2^(n-1) = 8× (job 4280770).
        maps4 = run(; kernels=K, subdiv_max=3, cpu_workers=4)
        for k in (:mass, :kappa, :isw)
            @test sum(maps4[k]) ≈ sum(maps3[k]) rtol=1e-9
        end
    end

    @testset "exclude_halos removes exactly the enclosed cells" begin
        halos = (x=[-250.0, -150.0], y=[-250.0, -320.0], z=[-250.0, -300.0], R=[30.0, 20.0])
        mex = run(; kernels=[:mass], subdiv_max=1, exclude_halos=halos)
        nex = 0
        for n in 1:2, k in 1:Nc, j in 1:Nc, i in 1:Nc
            x, y, z = cellpos(i), cellpos(j), cellpos(k)
            (x - halos.x[n])^2 + (y - halos.y[n])^2 + (z - halos.z[n])^2 <= halos.R[n]^2 || continue
            rmin_eff <= rdist(x, y, z) <= chimax && (nex += 1)
        end
        @test nex > 0
        @test sum(maps[:mass]) - sum(mex[:mass]) ≈ rho_m * alatt^3 * nex rtol=1e-6
    end

    @testset "cross-tile exclusion mask == brute force" begin
        hb = (x=[0.0], y=[-200.0], z=[-200.0], R=[30.0])     # straddles the x=0 tile boundary
        dcore = nsub * alatt; x0 = -(ntile / 2) * dcore
        bt = ntuple(d -> clamp(floor(Int, ((hb.x[1], hb.y[1], hb.z[1])[d] - x0) / dcore) + 1, 1, ntile), 3)
        bins = Dict{NTuple{3,Int},Vector{Int}}(bt => [1])
        nmask = 0
        for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
            xb, yb, zb = PeakPatch.MultiResolution.tile_center(it, jt, kt, ntile, dcore)
            nmask += count(PeakPatch.MultiResolution._build_exclusion_mask(
                hb, bins, it, jt, kt, ntile, dcore, nmesh, nbuff, alatt, xb, yb, zb))
        end
        nbrute = count((cellpos(i) - hb.x[1])^2 + (cellpos(j) - hb.y[1])^2 + (cellpos(k) - hb.z[1])^2 <= hb.R[1]^2
                       for k in 1:Nc, j in 1:Nc, i in 1:Nc)
        @test nmask == nbrute
    end
end
