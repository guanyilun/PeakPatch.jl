using PeakPatch
using PeakPatch: Exclusion, Merger, Catalog
using Test

@testset "Merger" begin

    @testset "sphere_overlap" begin
        # No overlap: spheres too far apart
        v1, v2 = sphere_overlap(10.0, 2.0, 3.0)
        @test v1 == 0.0
        @test v2 == 0.0

        # Touching: d = r1 + r2
        v1, v2 = sphere_overlap(5.0, 2.0, 3.0)
        @test v1 == 0.0
        @test v2 == 0.0

        # One inside the other: d < |r1 - r2| → returns (0,0), handled by exclusion
        v1, v2 = sphere_overlap(0.5, 3.0, 1.0)
        @test v1 == 0.0
        @test v2 == 0.0

        # Equal spheres at d=0 → one inside other
        v1, v2 = sphere_overlap(0.0, 2.0, 2.0)
        @test v1 == 0.0
        @test v2 == 0.0

        # Partial overlap: two equal spheres at d = r
        # Each cap height h = r/2, cap volume = (π/3)(r/2)²(3r - r/2) = (π/3)(r²/4)(5r/2)
        r = 2.0
        d = r  # d = r, so overlap region exists
        v1, v2 = sphere_overlap(d, r, r)
        @test v1 ≈ v2  # symmetric
        h = r / 2.0
        expected = (π / 3.0) * h^2 * (3.0 * r - h)
        @test v1 ≈ expected

        # Total overlap volume should be v1 + v2
        # For equal spheres at d=r: known formula
        total_overlap = v1 + v2
        @test total_overlap > 0.0
        @test total_overlap < (4π/3) * r^3  # less than full sphere volume

        # Asymmetric case: r1=3, r2=1, d=3 (r2 center on r1 surface)
        v1, v2 = sphere_overlap(3.0, 3.0, 1.0)
        @test v1 > 0.0
        @test v2 > 0.0
        @test v1 < v2  # smaller sphere loses more fraction
    end

    @testset "SpatialHash build and query" begin
        # 4 halos in a small box
        x = [1.0, 2.0, 1.1, 50.0]
        y = [1.0, 2.0, 1.1, 50.0]
        z = [1.0, 2.0, 1.1, 50.0]

        sh = Exclusion.build_hash(x, y, z, 4;
            domain_min=(0.0, 0.0, 0.0), domain_max=(100.0, 100.0, 100.0))
        @test sh.nc == 256
        @test sh.cell_size ≈ 100.0 / 256

        # Halos 1 and 3 should be in the same cell (positions ~1.0)
        ix1, iy1, iz1 = Exclusion._cell_idx(sh, 1.0, 1.0, 1.0)
        ix3, iy3, iz3 = Exclusion._cell_idx(sh, 1.1, 1.1, 1.1)
        @test ix1 == ix3
        @test iy1 == iy3
        @test iz1 == iz3

        # Halo 4 should be in a different cell
        ix4, iy4, iz4 = Exclusion._cell_idx(sh, 50.0, 50.0, 50.0)
        @test ix4 != ix1
    end

    @testset "lagrangian_exclusion! — basic" begin
        # 3 halos: A (large, at origin), B (small, inside A), C (far away)
        x = [0.0, 0.3, 10.0]
        y = [0.0, 0.3, 10.0]
        z = [0.0, 0.3, 10.0]
        r = [2.0, 0.5, 1.0]

        order = sortperm(r; rev=true)  # [1, 3, 2]
        sh = Exclusion.build_hash(x, y, z, 3;
            domain_min=(-5.0, -5.0, -5.0), domain_max=(15.0, 15.0, 15.0))
        survived = [true, true, true]

        Exclusion.lagrangian_exclusion!(survived, x, y, z, r, order, sh)

        @test survived[1] == true   # A survives (largest)
        @test survived[2] == false  # B killed (center at 0.52 < r_A=2.0)
        @test survived[3] == true   # C survives (far away)
    end

    @testset "lagrangian_exclusion! — equal radius" begin
        # Two equal halos overlapping: larger-index processed second
        x = [0.0, 0.5]
        y = [0.0, 0.0]
        z = [0.0, 0.0]
        r = [2.0, 2.0]

        order = sortperm(r; rev=true)  # [1, 2] (stable sort, 1 first)
        sh = Exclusion.build_hash(x, y, z, 2;
            domain_min=(-5.0, -5.0, -5.0), domain_max=(5.0, 5.0, 5.0))
        survived = [true, true]

        Exclusion.lagrangian_exclusion!(survived, x, y, z, r, order, sh)

        # Halo 2 center at 0.5 < r_1=2.0, and r_1 >= r_2, so halo 2 is killed
        @test survived[1] == true
        @test survived[2] == false
    end

    @testset "volume_reduction! — overlapping pair" begin
        # Two equal halos, partially overlapping (not center-inside)
        r_val = 2.0
        d = 3.0  # > r but < 2r, so overlap but no exclusion
        x = [0.0, d]
        y = [0.0, 0.0]
        z = [0.0, 0.0]
        r = [r_val, r_val]

        order = sortperm(r; rev=true)
        sh = Exclusion.build_hash(x, y, z, 2;
            domain_min=(-5.0, -5.0, -5.0), domain_max=(10.0, 10.0, 10.0))
        survived = [true, true]

        new_r = Exclusion.volume_reduction!(survived, x, y, z, r, order, sh)

        @test survived[1] == true
        @test survived[2] == true
        @test new_r[1] < r_val  # radius reduced
        @test new_r[2] < r_val
        @test new_r[1] ≈ new_r[2]  # symmetric → equal reduction
    end

    @testset "volume_reduction! — no overlap" begin
        x = [0.0, 100.0]
        y = [0.0, 0.0]
        z = [0.0, 0.0]
        r = [2.0, 3.0]

        order = sortperm(r; rev=true)
        sh = Exclusion.build_hash(x, y, z, 2;
            domain_min=(-5.0, -5.0, -5.0), domain_max=(105.0, 105.0, 105.0))
        survived = [true, true]

        new_r = Exclusion.volume_reduction!(survived, x, y, z, r, order, sh)

        @test new_r[1] ≈ 2.0
        @test new_r[2] ≈ 3.0
    end

    @testset "merge_catalog — synthetic HaloRecords" begin
        # Create synthetic halos:
        # H1: large halo at origin
        # H2: small halo inside H1 → should be excluded
        # H3: halo overlapping H1 but not center-inside → survives with reduced radius
        # H4: isolated halo → survives unchanged
        halos = [
            HaloRecord(0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 5.0f0, 0f0, 0f0, 0f0, 1.0f0),   # H1: r=5
            HaloRecord(1f0, 1f0, 1f0, 0f0, 0f0, 0f0, 1.0f0, 0f0, 0f0, 0f0, 0.5f0),   # H2: r=1, dist=1.73 < 5
            HaloRecord(7f0, 0f0, 0f0, 0f0, 0f0, 0f0, 3.0f0, 0f0, 0f0, 0f0, 0.8f0),   # H3: r=3, dist=7, overlap (7 < 5+3)
            HaloRecord(50f0, 50f0, 50f0, 0f0, 0f0, 0f0, 2.0f0, 0f0, 0f0, 0f0, 0.3f0), # H4: isolated
        ]

        merged = merge_catalog(halos)

        @test length(merged) == 3  # H2 excluded
        # Check H2 is gone (the one with overdensity 0.5 at position 1,1,1)
        positions = [(h.x, h.y, h.z) for h in merged]
        @test !((1f0, 1f0, 1f0) in positions)
        # H4 should be unchanged
        h4 = merged[findfirst(h -> h.x == 50f0, merged)]
        @test h4.RTHL == 2.0f0
        # H1 should have reduced radius (overlap with H3)
        h1 = merged[findfirst(h -> h.x == 0f0, merged)]
        @test h1.RTHL < 5.0f0
        @test h1.RTHL > 4.0f0  # not drastically reduced
    end

    @testset "merge_catalog — empty input" begin
        @test merge_catalog(HaloRecord[]) == HaloRecord[]
        @test merge_catalog(ExtHaloRecord[]) == ExtHaloRecord[]
    end

    @testset "merge_catalog — single halo" begin
        h = HaloRecord(0f0, 0f0, 0f0, 1f0, 2f0, 3f0, 5.0f0, 0f0, 0f0, 0f0, 1.0f0)
        result = merge_catalog([h])
        @test length(result) == 1
        @test result[1].RTHL == 5.0f0
        @test result[1].vx == 1f0  # fields preserved
    end

    @testset "merge_catalog — no overlaps" begin
        halos = [
            HaloRecord(0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 1.0f0, 0f0, 0f0, 0f0, 1.0f0),
            HaloRecord(100f0, 0f0, 0f0, 0f0, 0f0, 0f0, 1.0f0, 0f0, 0f0, 0f0, 1.0f0),
            HaloRecord(0f0, 100f0, 0f0, 0f0, 0f0, 0f0, 1.0f0, 0f0, 0f0, 0f0, 1.0f0),
        ]
        result = merge_catalog(halos)
        @test length(result) == 3
        for i in 1:3
            @test result[i].RTHL == 1.0f0
        end
    end

    @testset "merge_catalog — chain exclusion" begin
        # A > B > C, all nested: A at origin, B inside A, C inside B
        halos = [
            HaloRecord(0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 10.0f0, 0f0, 0f0, 0f0, 1.0f0),
            HaloRecord(1f0, 0f0, 0f0, 0f0, 0f0, 0f0, 5.0f0, 0f0, 0f0, 0f0, 1.0f0),
            HaloRecord(1.5f0, 0f0, 0f0, 0f0, 0f0, 0f0, 2.0f0, 0f0, 0f0, 0f0, 1.0f0),
        ]
        result = merge_catalog(halos)
        @test length(result) == 1
        @test result[1].RTHL ≈ 10.0f0  # only the largest survives, no overlap partner
    end

    @testset "merge_catalog — ExtHaloRecord" begin
        h = ExtHaloRecord(
            0f0, 0f0, 0f0, 1f0, 2f0, 3f0, 5.0f0, 0f0, 0f0, 0f0, 1.0f0,
            0.1f0, 0.2f0,
            0.01f0, 0.02f0, 0.03f0, 0.04f0, 0.05f0, 0.06f0,
            0.5f0, 2.0f0,
            0.1f0, 0.2f0, 0.3f0,
            0.4f0, 0.5f0, 0.6f0,
            1.0f0, 1.5f0, 0.7f0,
            0.01f0, 0.02f0, 0.03f0
        )
        result = merge_catalog([h])
        @test length(result) == 1
        @test result[1].e_v == 0.1f0
        @test result[1].strain_11 == 0.01f0
        @test result[1].Rf == 1.0f0
    end

    @testset "merge_catalog — preserves non-overlapping fields" begin
        # Two halos with distinct velocities, no overlap
        h1 = HaloRecord(0f0, 0f0, 0f0, 10f0, 20f0, 30f0, 1.0f0, 4f0, 5f0, 6f0, 0.5f0)
        h2 = HaloRecord(50f0, 0f0, 0f0, -1f0, -2f0, -3f0, 2.0f0, -4f0, -5f0, -6f0, 0.8f0)
        result = merge_catalog([h1, h2])
        @test length(result) == 2
        # Find h2 in output
        r2 = result[findfirst(h -> h.x == 50f0, result)]
        @test r2.vx == -1f0
        @test r2.vy == -2f0
        @test r2.overdensity == 0.8f0
    end

end
