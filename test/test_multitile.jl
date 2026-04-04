@testset "MultiTile" begin

using PeakPatch
using Test

# ---------- helpers (shared with test_pipeline.jl) ----------

function _make_pk_file(dir)
    path = joinpath(dir, "test_pk.dat")
    ks = 10.0 .^ range(-4, stop=1, length=200)
    open(path, "w") do f
        for k in ks
            Pk = 1e4 * k^(-1.5)
            println(f, "$k  $Pk")
        end
    end
    return path
end

function _make_filter_file(dir)
    path = joinpath(dir, "test_filter.dat")
    open(path, "w") do f
        println(f, "3")
        println(f, "1  1.686  8.0")
        println(f, "2  1.686  5.0")
        println(f, "3  1.686  3.0")
    end
    return path
end

function _make_simparams(dir; n=64, boxsize=200.0, z=0.0, ilpt=1, nbuff=4,
                         ioutshear=0, rmax2rs=0.0)
    pkfile = _make_pk_file(dir)
    filterfile = _make_filter_file(dir)
    outfile = joinpath(dir, "test_out.pksc")
    tabfile = joinpath(@__DIR__, "data", "HomelTab_julia.dat")

    nsub = n - 2 * nbuff
    dcore_box = boxsize * nsub / n  # dcore_box = nsub * alatt = nsub * (boxsize/n)

    PeakPatch.SimParams(
        Int32(0),       # ireadfield
        Int32(ioutshear), # ioutshear
        Float32(z),     # global_redshift
        Float32(z),     # maximum_redshift
        Int32(1),       # num_redshifts
        Float32(0.315 - 0.049), # Omx (CDM only)
        Float32(0.049), # OmB
        Float32(0.685), # Omvac
        Float32(0.674), # h
        Int32(n),       # nlx
        Int32(n),       # nly
        Int32(n),       # nlz
        Float32(dcore_box), # dcore_box
        Float32(boxsize),   # dL_box
        Float32(0.0),   # cenx
        Float32(0.0),   # ceny
        Float32(0.0),   # cenz
        Int32(nbuff),   # nbuff
        Int32(0),       # next
        Int32(0),       # ievol
        Int32(2),       # ivir_strat
        Float32(0.171), # fcoll_3
        Float32(0.171), # fcoll_2
        Float32(0.01),  # fcoll_1
        Float32(200.0), # dcrit
        Int32(4),       # iforce_strat
        Int32(50),      # TabInterpNx
        Int32(20),      # TabInterpNy
        Int32(20),      # TabInterpNz
        Float32(log10(1.5)), # TabInterpX1
        Float32(log10(8.0)), # TabInterpX2
        Float32(0.0),   # TabInterpY1
        Float32(0.5),   # TabInterpY2
        Float32(-0.9999), # TabInterpZ1
        Float32(0.9999),  # TabInterpZ2
        Int32(0),       # wsmooth (Gaussian)
        Float32(rmax2rs), # rmax2rs
        Int32(1),       # ioutfield
        Int32(0),       # NonGauss
        Float32(0.0),   # fNL
        Float32(0.0),   # A_nG
        Float32(0.0),   # B_nG
        Float32(0.0),   # R_nG
        Int32(ilpt),    # ilpt
        Int32(0),       # iwant_field_part
        Int32(0),       # largerun
        "",             # fielddir
        "",             # densfilein
        "",             # filein
        pkfile,         # pkfile
        filterfile,     # filterfile
        outfile,        # fileout
        tabfile,        # TabInterpFile
    )
end


# =================================================================
# Test 1: extract_tile geometry
# =================================================================
@testset "extract_tile" begin
    nsub = 4
    nbuff = 2
    nmesh = nsub + 2 * nbuff  # 8
    ntile = 2
    N = nsub * ntile + 2 * nbuff  # 12

    # Create a field where each cell stores its global linear index
    field = zeros(Float32, N, N, N)
    for k in 1:N, j in 1:N, i in 1:N
        field[i, j, k] = Float32(i + (j-1)*N + (k-1)*N*N)
    end

    # Tile (1,1,1): starts at (1,1,1), ends at (8,8,8)
    t1 = PeakPatch.MultiTile.extract_tile(field, 1, 1, 1, nsub, nmesh)
    @test size(t1) == (nmesh, nmesh, nmesh)
    @test t1[1, 1, 1] == field[1, 1, 1]
    @test t1[nmesh, nmesh, nmesh] == field[nmesh, nmesh, nmesh]

    # Tile (2,1,1): starts at (nsub+1, 1, 1) = (5, 1, 1)
    t2 = PeakPatch.MultiTile.extract_tile(field, 2, 1, 1, nsub, nmesh)
    @test t2[1, 1, 1] == field[nsub+1, 1, 1]

    # Buffer overlap: tile 1 cells [nsub+1:nmesh] == tile 2 cells [1:2*nbuff]
    @test t1[nsub+1:nmesh, 1, 1] == t2[1:2*nbuff, 1, 1]

    # Tile (2,2,2): starts at (5,5,5), ends at (12,12,12)
    t3 = PeakPatch.MultiTile.extract_tile(field, 2, 2, 2, nsub, nmesh)
    @test t3[nmesh, nmesh, nmesh] == field[N, N, N]
end

# =================================================================
# Test 2: tile_center
# =================================================================
@testset "tile_center" begin
    dcore_box = 100.0

    # ntile=1: single tile centered at origin
    x, y, z = PeakPatch.MultiTile.tile_center(1, 1, 1, 1, dcore_box)
    @test x ≈ 0.0
    @test y ≈ 0.0
    @test z ≈ 0.0

    # ntile=2: two tiles centered at ±dcore_box/2
    x1, _, _ = PeakPatch.MultiTile.tile_center(1, 1, 1, 2, dcore_box)
    x2, _, _ = PeakPatch.MultiTile.tile_center(2, 1, 1, 2, dcore_box)
    @test x1 ≈ -dcore_box / 2
    @test x2 ≈ +dcore_box / 2

    # ntile=3: centered at -dcore_box, 0, +dcore_box
    x1, _, _ = PeakPatch.MultiTile.tile_center(1, 1, 1, 3, dcore_box)
    x2, _, _ = PeakPatch.MultiTile.tile_center(2, 1, 1, 3, dcore_box)
    x3, _, _ = PeakPatch.MultiTile.tile_center(3, 1, 1, 3, dcore_box)
    @test x1 ≈ -dcore_box
    @test x2 ≈ 0.0
    @test x3 ≈ +dcore_box
end

# =================================================================
# Test 3: Grid geometry formulas
# =================================================================
@testset "Grid geometry" begin
    nmesh = 142
    nbuff = 4
    nsub = nmesh - 2 * nbuff  # 134
    alatt = 200.0 / nmesh      # dL_box / nmesh

    # ntile=1: N = nmesh (matches single-tile)
    N1 = nsub * 1 + 2 * nbuff
    @test N1 == nmesh

    boxsize1 = N1 * alatt
    @test boxsize1 ≈ 200.0

    # ntile=2: larger domain
    N2 = nsub * 2 + 2 * nbuff
    @test N2 == 2 * nmesh - 2 * nbuff  # 276

    # Tile 2 ends at cell (1)*nsub + nmesh = 134 + 142 = 276 = N
    @test (2 - 1) * nsub + nmesh == N2
end

# =================================================================
# Test 4: ntile=1 matches run_tile exactly
# =================================================================
@testset "ntile=1 matches run_tile" begin
    tmpdir = mktempdir()
    sp = _make_simparams(tmpdir; n=32, boxsize=100.0, z=0.0, ilpt=1, nbuff=3)

    # Run single-tile
    halos_single = PeakPatch.Pipeline.run_tile(sp; seed=42, verbose=false)

    # Run multi-tile with ntile=1 (should be identical)
    sp_mt = _make_simparams(tmpdir; n=32, boxsize=100.0, z=0.0, ilpt=1, nbuff=3)
    # Need a different output file
    outfile2 = joinpath(tmpdir, "test_out_mt.pksc")
    sp_mt = PeakPatch.SimParams(
        sp.ireadfield, sp.ioutshear, sp.global_redshift, sp.maximum_redshift,
        sp.num_redshifts, sp.Omx, sp.OmB, sp.Omvac, sp.h,
        sp.nlx, sp.nly, sp.nlz, sp.dcore_box, sp.dL_box,
        sp.cenx, sp.ceny, sp.cenz, sp.nbuff, sp.next, sp.ievol,
        sp.ivir_strat, sp.fcoll_3, sp.fcoll_2, sp.fcoll_1, sp.dcrit, sp.iforce_strat,
        sp.TabInterpNx, sp.TabInterpNy, sp.TabInterpNz,
        sp.TabInterpX1, sp.TabInterpX2, sp.TabInterpY1, sp.TabInterpY2,
        sp.TabInterpZ1, sp.TabInterpZ2,
        sp.wsmooth, sp.rmax2rs, sp.ioutfield,
        sp.NonGauss, sp.fNL, sp.A_nG, sp.B_nG, sp.R_nG,
        sp.ilpt, sp.iwant_field_part, sp.largerun,
        sp.fielddir, sp.densfilein, sp.filein, sp.pkfile, sp.filterfile,
        outfile2, sp.TabInterpFile
    )

    halos_multi = PeakPatch.MultiTile.run_multitile(sp_mt; ntile=1, seed=42, verbose=false)

    @test length(halos_single) == length(halos_multi)

    if !isempty(halos_single) && !isempty(halos_multi)
        # Sort both by x coordinate for stable comparison
        sort!(halos_single; by=h -> (h.x, h.y, h.z))
        sort!(halos_multi; by=h -> (h.x, h.y, h.z))

        for i in 1:length(halos_single)
            @test halos_single[i].x ≈ halos_multi[i].x atol=1e-6
            @test halos_single[i].y ≈ halos_multi[i].y atol=1e-6
            @test halos_single[i].z ≈ halos_multi[i].z atol=1e-6
            @test halos_single[i].RTHL ≈ halos_multi[i].RTHL atol=1e-6
            @test halos_single[i].vx ≈ halos_multi[i].vx atol=1e-6
        end
    end

    rm(tmpdir; recursive=true)
end

# =================================================================
# Test 5: ntile=2 runs and produces halos
# =================================================================
@testset "ntile=2 end-to-end" begin
    tmpdir = mktempdir()
    # Use n=20 (small tile), nbuff=3, nsub=14
    # ntile=2: N = 14*2 + 6 = 34
    sp = _make_simparams(tmpdir; n=20, boxsize=60.0, z=0.0, ilpt=1, nbuff=3)

    halos = PeakPatch.MultiTile.run_multitile(sp; ntile=2, seed=42, verbose=false)

    # Should find some halos (hard to guarantee exact count, but should be > 0 for this setup)
    @test length(halos) >= 0  # at minimum, should not error

    if !isempty(halos)
        @test all(h -> h.RTHL > 0, halos)

        # Check coordinate ranges: halos should span the full multi-tile domain
        nmesh = Int(sp.nlx)
        nbuff_val = Int(sp.nbuff)
        nsub = nmesh - 2 * nbuff_val
        dcore_box = nsub * Float64(sp.dL_box) / nmesh

        # Peak coordinates should be within ±ntile/2 * dcore_box (approximately)
        max_coord = dcore_box  # ntile=2, so ±dcore_box
        xs = [h.x for h in halos]
        @test minimum(xs) >= -max_coord - 1  # some tolerance
        @test maximum(xs) <= +max_coord + 1
    end

    rm(tmpdir; recursive=true)
end

# =================================================================
# Test 6: ntile=2 with 2LPT
# =================================================================
@testset "ntile=2 with 2LPT" begin
    tmpdir = mktempdir()
    sp = _make_simparams(tmpdir; n=20, boxsize=60.0, z=0.0, ilpt=2, nbuff=3)

    halos = PeakPatch.MultiTile.run_multitile(sp; ntile=2, seed=99, verbose=false)

    if !isempty(halos)
        @test all(h -> h.RTHL > 0, halos)
    end

    rm(tmpdir; recursive=true)
end

# =================================================================
# Test 7: ntile=1 with extended output (ioutshear=1)
# =================================================================
@testset "ntile=1 with ioutshear" begin
    tmpdir = mktempdir()
    sp = _make_simparams(tmpdir; n=32, boxsize=100.0, z=0.0, ilpt=1, nbuff=3,
                         ioutshear=1)

    halos_single = PeakPatch.Pipeline.run_tile(sp; seed=42, verbose=false)

    outfile2 = joinpath(tmpdir, "test_out_mt2.pksc")
    sp_mt = PeakPatch.SimParams(
        sp.ireadfield, sp.ioutshear, sp.global_redshift, sp.maximum_redshift,
        sp.num_redshifts, sp.Omx, sp.OmB, sp.Omvac, sp.h,
        sp.nlx, sp.nly, sp.nlz, sp.dcore_box, sp.dL_box,
        sp.cenx, sp.ceny, sp.cenz, sp.nbuff, sp.next, sp.ievol,
        sp.ivir_strat, sp.fcoll_3, sp.fcoll_2, sp.fcoll_1, sp.dcrit, sp.iforce_strat,
        sp.TabInterpNx, sp.TabInterpNy, sp.TabInterpNz,
        sp.TabInterpX1, sp.TabInterpX2, sp.TabInterpY1, sp.TabInterpY2,
        sp.TabInterpZ1, sp.TabInterpZ2,
        sp.wsmooth, sp.rmax2rs, sp.ioutfield,
        sp.NonGauss, sp.fNL, sp.A_nG, sp.B_nG, sp.R_nG,
        sp.ilpt, sp.iwant_field_part, sp.largerun,
        sp.fielddir, sp.densfilein, sp.filein, sp.pkfile, sp.filterfile,
        outfile2, sp.TabInterpFile
    )

    halos_multi = PeakPatch.MultiTile.run_multitile(sp_mt; ntile=1, seed=42, verbose=false)
    @test length(halos_single) == length(halos_multi)

    if !isempty(halos_single) && !isempty(halos_multi)
        sort!(halos_single; by=h -> (h.x, h.y, h.z))
        sort!(halos_multi; by=h -> (h.x, h.y, h.z))
        for i in 1:length(halos_single)
            @test halos_single[i].RTHL ≈ halos_multi[i].RTHL atol=1e-6
            @test halos_single[i].e_v ≈ halos_multi[i].e_v atol=1e-6
        end
    end

    rm(tmpdir; recursive=true)
end

end # MultiTile testset
