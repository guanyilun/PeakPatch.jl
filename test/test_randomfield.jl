@testset "RandomField" begin
    using FFTW
    using Statistics: mean, var, std

    @testset "GRF basic statistics" begin
        pk_uniform(k) = 1.0
        n = 64
        boxsize = 100.0
        seed = 42

        field = generate_grf(n, pk_uniform, boxsize, seed)

        # Output should be n×n×n
        @test size(field) == (n, n, n)
        # Element type should be Float32
        @test eltype(field) == Float32

        # Mean should be approximately zero (DC mode removed)
        @test abs(mean(field)) < 0.1

        # Variance should be positive and finite
        @test var(field) > 0.0
        @test isfinite(var(field))
    end

    @testset "GRF reproducibility (same seed)" begin
        pk_uniform(k) = 1.0
        f1 = generate_grf(32, pk_uniform, 100.0, 12345)
        f2 = generate_grf(32, pk_uniform, 100.0, 12345)
        @test f1 ≈ f2
    end

    @testset "GRF different seeds → different fields" begin
        pk_uniform(k) = 1.0
        f1 = generate_grf(32, pk_uniform, 100.0, 1)
        f2 = generate_grf(32, pk_uniform, 100.0, 2)
        @test !(f1 ≈ f2)
    end

    @testset "GRF power spectrum matches input" begin
        # Generate two fields with different P(k) amplitudes and verify
        # that the field variance scales linearly with the amplitude.
        n = 64
        boxsize = 200.0

        pk_A(k) = 100.0
        pk_B(k) = 400.0  # 4× amplitude

        field_A = generate_grf(n, pk_A, boxsize, 42)
        field_B = generate_grf(n, pk_B, boxsize, 42)

        var_A = var(field_A)
        var_B = var(field_B)

        # Variance should scale linearly with P(k) amplitude
        @test var_B / var_A ≈ 4.0 rtol=0.1
    end

    @testset "GRF absolute variance ≈ A × n³ × dk³" begin
        # For a white-noise field P(k) = A, the expected real-space variance is:
        #   〈δ²〉 = A × n³ × dk³
        # Derivation (Parseval with the n³-amplitude convention):
        #   δ(x) = (1/n³) × IFFT[noise_k × sqrt(A × dk³ × n³)]
        #   〈|δ|^2〉 = (1/n^6) × Σ_k A × dk³ × n³ × E[|noise_k|^2]
        # Since E[Σ_k |noise_k|^2] = n³ (Parseval for unit-variance white noise):
        #   〈δ²〉 = (1/n^6) × A × dk³ × n³ × n³ = A × dk³ × n³
        # This is confirmed numerically: ratio measured/(A×n³×dk³) ≈ 1.0.
        n = 32
        boxsize = 100.0
        A = 50.0
        pk_const(k) = A

        # Average over several seeds to beat down sampling variance
        vars = [Float64(var(generate_grf(n, pk_const, boxsize, UInt64(s)))) for s in 1:8]
        measured_var = mean(vars)

        dk = 2π / boxsize
        expected_var = A * Float64(n)^3 * dk^3

        # Expect within 15% (sampling variance)
        @test measured_var ≈ expected_var rtol=0.15
    end

    @testset "GRF Gaussianity" begin
        pk_uniform(k) = 1.0
        n = 64
        field = generate_grf(n, pk_uniform, 200.0, 99)

        # A Gaussian field should have near-zero skewness and kurtosis ≈ 3
        m = mean(field)
        s = std(field)
        skewness = mean(((field .- m) ./ s) .^ 3)
        kurt = mean(((field .- m) ./ s) .^ 4)

        @test abs(skewness) < 0.1     # Gaussian → skew ≈ 0
        @test abs(kurt - 3.0) < 0.3   # Gaussian → kurt ≈ 3
    end

    # ================================================================
    # Threefry counter-based RNG noise generation
    # ================================================================

    @testset "Threefry determinism" begin
        n = 32
        a = zeros(Float64, n, n, n)
        b = zeros(Float64, n, n, n)
        PeakPatch.fill_noise_threefry!(a, n, 42)
        PeakPatch.fill_noise_threefry!(b, n, 42)
        @test a == b
    end

    @testset "Threefry different seeds → different fields" begin
        n = 32
        a = zeros(Float64, n, n, n)
        b = zeros(Float64, n, n, n)
        PeakPatch.fill_noise_threefry!(a, n, 42)
        PeakPatch.fill_noise_threefry!(b, n, 99)
        @test a != b
    end

    @testset "Threefry Gaussian statistics" begin
        n = 128
        noise = zeros(Float64, n, n, n)
        PeakPatch.fill_noise_threefry!(noise, n, 42)
        μ = mean(noise)
        σ = std(noise)
        # n=128 → 2M samples, SE of mean ≈ 0.0007
        @test abs(μ) < 0.01
        @test abs(σ - 1.0) < 0.01
    end

    @testset "Threefry Gaussianity (skewness/kurtosis)" begin
        n = 128
        noise = zeros(Float64, n, n, n)
        PeakPatch.fill_noise_threefry!(noise, n, 77)
        m = mean(noise)
        s = std(noise)
        skewness = mean(((noise .- m) ./ s) .^ 3)
        kurt = mean(((noise .- m) ./ s) .^ 4)
        @test abs(skewness) < 0.05
        @test abs(kurt - 3.0) < 0.1
    end

    @testset "Threefry region matches full array" begin
        n = 64
        full = zeros(Float64, n, n, n)
        PeakPatch.fill_noise_threefry!(full, n, 42)

        # Test several sub-regions including edge cases
        regions = [
            (1:32, 1:32, 1:32),
            (33:64, 1:64, 1:64),
            (10:50, 20:60, 5:45),
            (1:1, 1:1, 1:1),      # single cell
            (1:64, 1:64, 1:64),    # full array
        ]
        for (r1, r2, r3) in regions
            region = zeros(Float64, length(r1), length(r2), length(r3))
            PeakPatch.fill_noise_threefry_region!(region, (r1, r2, r3), n, 42)
            @test region ≈ full[r1, r2, r3]
        end
    end

    @testset "Threefry decomposition-independent" begin
        # Splitting the grid into non-overlapping pencils should produce
        # the same noise as a single fill — this is the key MPI property
        n = 32
        full = zeros(Float64, n, n, n)
        PeakPatch.fill_noise_threefry!(full, n, 42)

        reassembled = zeros(Float64, n, n, n)
        # Simulate 4-way pencil decomposition along k
        nk_per = n ÷ 4
        for p in 0:3
            k_start = p * nk_per + 1
            k_end = (p + 1) * nk_per
            piece = zeros(Float64, n, n, nk_per)
            PeakPatch.fill_noise_threefry_region!(piece, (1:n, 1:n, k_start:k_end), n, 42)
            reassembled[:, :, k_start:k_end] .= piece
        end
        @test reassembled == full
    end

    @testset "Threefry Float32 output" begin
        n = 16
        noise = zeros(Float32, n, n, n)
        PeakPatch.fill_noise_threefry!(noise, n, 42)
        @test eltype(noise) == Float32
        @test all(isfinite, noise)
    end

    @testset "generate_grf basic" begin
        pk_uniform(k) = 1.0
        n = 64
        boxsize = 100.0

        field = generate_grf(n, pk_uniform, boxsize, 42)
        @test size(field) == (n, n, n)
        @test eltype(field) == Float32
        @test abs(mean(field)) < 0.1  # DC zeroed
        @test var(field) > 0.0
        @test isfinite(var(field))
    end

    @testset "generate_grf reproducibility" begin
        pk_uniform(k) = 1.0
        f1 = generate_grf(32, pk_uniform, 100.0, 42)
        f2 = generate_grf(32, pk_uniform, 100.0, 42)
        @test f1 ≈ f2
    end

    @testset "generate_grf power spectrum scaling" begin
        n = 64
        boxsize = 200.0
        pk_A(k) = 100.0
        pk_B(k) = 400.0

        field_A = generate_grf(n, pk_A, boxsize, 42)
        field_B = generate_grf(n, pk_B, boxsize, 42)

        @test var(field_B) / var(field_A) ≈ 4.0 rtol=0.1
    end

    @testset "FFT round-trip" begin
        n = 32
        arr = randn(Float32, n, n, n)
        plan = plan_rfft(arr)
        arr_k = plan * arr
        recovered = irfft(arr_k, n)
        @test maximum(abs.(arr .- recovered)) < 1e-5
    end
end
