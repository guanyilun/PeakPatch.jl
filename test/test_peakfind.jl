@testset "PeakFind" begin
    using PeakPatch: PeakCandidate, find_peaks

    @testset "single peak" begin
        n = 11
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        # Place a peak at center
        delta[6, 6, 6] = 2.0f0
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test length(peaks) == 1
        @test peaks[1].i == 6 && peaks[1].j == 6 && peaks[1].k == 6
        @test peaks[1].delta ≈ 2.0f0
        @test peaks[1].Rsmooth == 3.0
    end

    @testset "no peaks below threshold" begin
        n = 11
        delta = fill(0.5f0, n, n, n)
        mask = zeros(Int8, n, n, n)
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test isempty(peaks)
    end

    @testset "no local max in gradient field" begin
        # Monotonically increasing field: every interior cell has a neighbor
        # with a strictly greater value, so no strict local max exists
        n = 11
        delta = zeros(Float32, n, n, n)
        for i in 1:n, j in 1:n, k in 1:n
            delta[i, j, k] = Float32(i + j + k)
        end
        mask = zeros(Int8, n, n, n)
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test isempty(peaks)
    end

    @testset "mask skip" begin
        n = 11
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        delta[6, 6, 6] = 2.0f0
        mask[6, 6, 6] = 1
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test isempty(peaks)
    end

    @testset "buffer zone exclusion" begin
        n = 11
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        # Peak at boundary (i=1, within nbuff=1)
        delta[1, 6, 6] = 2.0f0
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test isempty(peaks)
    end

    @testset "multiple peaks" begin
        n = 15
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        delta[4, 4, 4] = 2.0f0
        delta[12, 12, 12] = 3.0f0
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test length(peaks) == 2
        # Order: k-outer, j-middle, i-inner -> (4,4,4) first, then (12,12,12)
        @test peaks[1].i == 4 && peaks[1].j == 4 && peaks[1].k == 4
        @test peaks[2].i == 12 && peaks[2].j == 12 && peaks[2].k == 12
    end

    @testset "flat index" begin
        n = 11
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        delta[3, 5, 7] = 2.0f0
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test length(peaks) == 1
        expected_ipp = 3 + (5 - 1) * 11 + (7 - 1) * 11 * 11
        @test peaks[1].ipp == expected_ipp
    end

    @testset "physical coordinates" begin
        n = 21
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        delta[11, 11, 11] = 2.0f0
        alatt = 2.0
        xbx, ybx, zbx = 10.0, 20.0, 30.0
        peaks = find_peaks(delta, mask, xbx, ybx, zbx, alatt, 1, 1.0f0, 3.0)
        @test length(peaks) == 1
        # cen = 0.5*(21+1) = 11, so coord = xbx + alatt*(11-11) = xbx
        @test peaks[1].x ≈ xbx
        @test peaks[1].y ≈ ybx
        @test peaks[1].z ≈ zbx
    end

    @testset "max peaks limit" begin
        n = 15
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        # Place 4 well-separated peaks
        for (i, j, k) in [(4, 4, 4), (4, 4, 12), (12, 12, 4), (12, 12, 12)]
            delta[i, j, k] = 2.0f0
        end
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0;
                           max_peaks=2)
        @test length(peaks) == 2
    end

    @testset "mask mutation" begin
        n = 11
        delta = zeros(Float32, n, n, n)
        mask = zeros(Int8, n, n, n)
        delta[6, 6, 6] = 2.0f0
        peaks = find_peaks(delta, mask, 0.0, 0.0, 0.0, 1.0, 1, 1.0f0, 3.0)
        @test length(peaks) == 1
        @test mask[6, 6, 6] == 1
    end
end
