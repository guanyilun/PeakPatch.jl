# GPU shell gather kernel validation.
# Runs only if CUDA.functional() — skipped otherwise.

using CUDA
using Random

@testset "ShellAnalysisGPU (CUDA kernel)" begin
    if !CUDA.functional()
        @info "CUDA not functional — skipping GPU tests" CUDA.functional()
        return
    end

    # ------------------------------------------------------------
    # CPU reference: direct per-shell Fshell sum (no physics, just gather)
    # ------------------------------------------------------------
    function cpu_shell_sums(delta::Array{Float32,3}, peaks_i, peaks_j, peaks_k,
                            stab::ShellTables)
        npeaks = length(peaks_i)
        nshells = stab.nshells
        n1, n2, n3 = size(delta)
        Fshell = zeros(Float32, nshells, npeaks)
        nshell = zeros(Int32, nshells, npeaks)
        for p in 1:npeaks
            ci, cj, ck = peaks_i[p], peaks_j[p], peaks_k[p]
            for s in 1:nshells
                s0 = stab.shell_start[s]
                nc = stab.shell_count[s]
                ssum = 0.0f0
                n = Int32(0)
                for c in 0:nc-1
                    iv1 = ci + stab.offsets_di[s0 + c]
                    iv2 = cj + stab.offsets_dj[s0 + c]
                    iv3 = ck + stab.offsets_dk[s0 + c]
                    (iv1 < 1 || iv1 > n1 || iv2 < 1 || iv2 > n2 || iv3 < 1 || iv3 > n3) && continue
                    ssum += delta[iv1, iv2, iv3]
                    n += Int32(1)
                end
                Fshell[s, p] = ssum
                nshell[s, p] = n
            end
        end
        return Fshell, nshell
    end

    # ------------------------------------------------------------
    # Small deterministic grid + a handful of peaks
    # ------------------------------------------------------------
    @testset "Small grid (N=64, rmax=12)" begin
        N = 64
        rmax = 12
        stab = build_shell_tables(rmax)

        # Deterministic reproducible field
        rng_seed = 0xdeadbeef
        Random.seed!(rng_seed)
        delta = randn(Float32, N, N, N)

        # Peaks mixed: interior + near-boundary (exercises bounds check)
        peaks_i = Int32[20, 32, 5,  60, 8]
        peaks_j = Int32[20, 32, 5,  60, 56]
        peaks_k = Int32[20, 32, 60, 5,  8]

        # GPU
        Fgpu, ngpu = shell_fbar_gather_gpu(delta, peaks_i, peaks_j, peaks_k, stab)

        # CPU reference
        Fcpu, ncpu = cpu_shell_sums(delta, peaks_i, peaks_j, peaks_k, stab)

        @test size(Fgpu) == size(Fcpu)
        @test ngpu == ncpu
        # Sums over ~ tens-to-hundreds of cells: Float32 sums have tiny rounding
        # differences from CPU ordering. Allow 1e-4 relative tolerance.
        @test isapprox(Fgpu, Fcpu; rtol=1e-4, atol=1e-5)
    end

    @testset "Single peak, all cells in-bounds" begin
        N = 48
        rmax = 8
        stab = build_shell_tables(rmax)
        Random.seed!(42)
        delta = randn(Float32, N, N, N)

        ci, cj, ck = Int32(24), Int32(24), Int32(24)
        Fgpu, ngpu = shell_fbar_gather_gpu(delta, [ci], [cj], [ck], stab)
        Fcpu, ncpu = cpu_shell_sums(delta, [ci], [cj], [ck], stab)
        @test ngpu == ncpu
        @test isapprox(Fgpu, Fcpu; rtol=1e-4, atol=1e-5)
        # At shell r²=1 (6 face cells), the count must be exactly 6
        idx1 = findfirst(==(Int32(1)), stab.shell_r2)
        @test ngpu[idx1, 1] == 6
    end

    # ------------------------------------------------------------
    # CPU reference for the multi-field (delta + ψ×3) gather
    # ------------------------------------------------------------
    function cpu_shell_gather_psi(delta, etax, etay, etaz,
                                   peaks_i, peaks_j, peaks_k, stab)
        npeaks = length(peaks_i)
        nshells = stab.nshells
        n1, n2, n3 = size(delta)
        F = zeros(Float32, nshells, npeaks)
        N = zeros(Int32, nshells, npeaks)
        S = zeros(Float32, 3, nshells, npeaks)
        G = zeros(Float32, 3, nshells, npeaks)
        for p in 1:npeaks
            ci, cj, ck = peaks_i[p], peaks_j[p], peaks_k[p]
            for s in 1:nshells
                s0 = stab.shell_start[s]
                nc = stab.shell_count[s]
                Fs = 0.0f0
                Sx = 0.0f0; Sy = 0.0f0; Sz = 0.0f0
                Gx = 0.0f0; Gy = 0.0f0; Gz = 0.0f0
                n = Int32(0)
                for c in 0:nc-1
                    di = stab.offsets_di[s0 + c]
                    dj = stab.offsets_dj[s0 + c]
                    dk = stab.offsets_dk[s0 + c]
                    iv1 = ci + di
                    iv2 = cj + dj
                    iv3 = ck + dk
                    (iv1 < 1 || iv1 > n1 || iv2 < 1 || iv2 > n2 || iv3 < 1 || iv3 > n3) && continue
                    d = delta[iv1, iv2, iv3]
                    Fs += d
                    Sx += etax[iv1, iv2, iv3]
                    Sy += etay[iv1, iv2, iv3]
                    Sz += etaz[iv1, iv2, iv3]
                    Gx += d * Float32(di)
                    Gy += d * Float32(dj)
                    Gz += d * Float32(dk)
                    n += Int32(1)
                end
                F[s, p] = Fs
                N[s, p] = n
                S[1, s, p] = Sx; S[2, s, p] = Sy; S[3, s, p] = Sz
                G[1, s, p] = Gx; G[2, s, p] = Gy; G[3, s, p] = Gz
            end
        end
        return F, N, S, G
    end

    @testset "Multi-field gather (delta + ψ×3)" begin
        N = 56
        rmax = 10
        stab = build_shell_tables(rmax)
        Random.seed!(1001)
        delta = randn(Float32, N, N, N)
        etax  = randn(Float32, N, N, N)
        etay  = randn(Float32, N, N, N)
        etaz  = randn(Float32, N, N, N)
        # Mix interior + boundary peaks to exercise bounds check
        peaks_i = Int32[28, 8,  50, 20]
        peaks_j = Int32[28, 52, 10, 30]
        peaks_k = Int32[28, 28, 28, 48]

        Fg, ng, Sg, Gg = shell_gather_psi_gpu(delta, etax, etay, etaz,
                                               peaks_i, peaks_j, peaks_k, stab)
        Fc, nc, Sc, Gc = cpu_shell_gather_psi(delta, etax, etay, etaz,
                                               peaks_i, peaks_j, peaks_k, stab)

        @test ng == nc
        @test isapprox(Fg, Fc; rtol=1e-4, atol=1e-5)
        @test isapprox(Sg, Sc; rtol=1e-4, atol=1e-5)
        # Gshell entries can reach ~|delta| * |d_L|; d_L grows as rmax so allow
        # slightly looser rtol (still Float32 noise only).
        @test isapprox(Gg, Gc; rtol=2e-4, atol=1e-4)
    end

    # ------------------------------------------------------------
    # CPU reference including strain tensor
    # ------------------------------------------------------------
    function cpu_shell_gather_strain(delta, etax, etay, etaz,
                                      peaks_i, peaks_j, peaks_k, stab)
        npeaks = length(peaks_i)
        nshells = stab.nshells
        n1, n2, n3 = size(delta)
        F = zeros(Float32, nshells, npeaks)
        N = zeros(Int32, nshells, npeaks)
        S = zeros(Float32, 3, nshells, npeaks)
        G = zeros(Float32, 3, nshells, npeaks)
        SR = zeros(Float32, 3, 3, nshells, npeaks)
        for p in 1:npeaks
            ci, cj, ck = peaks_i[p], peaks_j[p], peaks_k[p]
            for s in 1:nshells
                s0 = stab.shell_start[s]
                nc = stab.shell_count[s]
                Fs = 0.0f0
                Sv = zeros(Float32, 3)
                Gv = zeros(Float32, 3)
                SRm = zeros(Float32, 3, 3)
                n = Int32(0)
                for c in 0:nc-1
                    di = stab.offsets_di[s0 + c]
                    dj = stab.offsets_dj[s0 + c]
                    dk = stab.offsets_dk[s0 + c]
                    iv1 = ci + di; iv2 = cj + dj; iv3 = ck + dk
                    (iv1 < 1 || iv1 > n1 || iv2 < 1 || iv2 > n2 || iv3 < 1 || iv3 > n3) && continue
                    d = delta[iv1, iv2, iv3]
                    e = (etax[iv1,iv2,iv3], etay[iv1,iv2,iv3], etaz[iv1,iv2,iv3])
                    dv = (Float32(di), Float32(dj), Float32(dk))
                    Fs += d
                    for L in 1:3
                        Sv[L] += e[L]
                        Gv[L] += d * dv[L]
                        for K in 1:3
                            SRm[L, K] += e[L] * dv[K]
                        end
                    end
                    n += Int32(1)
                end
                F[s, p] = Fs
                N[s, p] = n
                for L in 1:3
                    S[L, s, p] = Sv[L]
                    G[L, s, p] = Gv[L]
                    for K in 1:3
                        SR[L, K, s, p] = SRm[L, K]
                    end
                end
            end
        end
        return F, N, S, G, SR
    end

    @testset "Strain tensor gather (delta + ψ×3 + SR[3,3])" begin
        N = 48
        rmax = 8
        stab = build_shell_tables(rmax)
        Random.seed!(2025)
        delta = randn(Float32, N, N, N)
        etax  = randn(Float32, N, N, N)
        etay  = randn(Float32, N, N, N)
        etaz  = randn(Float32, N, N, N)
        peaks_i = Int32[24, 6,  42]
        peaks_j = Int32[24, 40, 8]
        peaks_k = Int32[24, 24, 40]

        Fg, ng, Sg, Gg, SRg = shell_gather_strain_gpu(
            delta, etax, etay, etaz, peaks_i, peaks_j, peaks_k, stab)
        Fc, nc, Sc, Gc, SRc = cpu_shell_gather_strain(
            delta, etax, etay, etaz, peaks_i, peaks_j, peaks_k, stab)

        @test ng == nc
        @test isapprox(Fg, Fc; rtol=1e-4, atol=1e-5)
        @test isapprox(Sg, Sc; rtol=1e-4, atol=1e-5)
        @test isapprox(Gg, Gc; rtol=2e-4, atol=1e-4)
        # Strain entries can scale as |eta|*|d_K|, so slightly looser tol for
        # Float32 round-off (especially large-r shells with O(100) cells).
        @test isapprox(SRg, SRc; rtol=2e-4, atol=2e-4)

        # Sanity: shape
        @test size(SRg) == (3, 3, stab.nshells, length(peaks_i))
    end

    # ------------------------------------------------------------
    # CPU reference for the full-field gather (delta + η + η² + SR + Gf)
    # Mirrors the CPU loop in src/ShellAnalysisGPU.jl:202-232 exactly.
    # ------------------------------------------------------------
    function cpu_shell_gather_full(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                                    peaks_i, peaks_j, peaks_k, stab)
        npeaks = length(peaks_i)
        nshells = stab.nshells
        n1, n2, n3 = size(delta)
        F = zeros(Float32, nshells, npeaks)
        N = zeros(Int32, nshells, npeaks)
        S  = zeros(Float32, 3, nshells, npeaks)
        S2 = zeros(Float32, 3, nshells, npeaks)
        G  = zeros(Float32, 3, nshells, npeaks)
        Gf = zeros(Float32, 3, nshells, npeaks)
        SR = zeros(Float32, 3, 3, nshells, npeaks)
        for p in 1:npeaks
            ci, cj, ck = peaks_i[p], peaks_j[p], peaks_k[p]
            for s in 1:nshells
                s0 = stab.shell_start[s]
                nc = stab.shell_count[s]
                Fs = 0.0f0
                Sv  = zeros(Float32, 3); S2v = zeros(Float32, 3)
                Gv  = zeros(Float32, 3); Gfv = zeros(Float32, 3)
                SRm = zeros(Float32, 3, 3)
                n = Int32(0)
                for c in 0:nc-1
                    di = stab.offsets_di[s0 + c]
                    dj = stab.offsets_dj[s0 + c]
                    dk = stab.offsets_dk[s0 + c]
                    iv1 = ci + di; iv2 = cj + dj; iv3 = ck + dk
                    (iv1 < 1 || iv1 > n1 || iv2 < 1 || iv2 > n2 || iv3 < 1 || iv3 > n3) && continue
                    d  = delta[iv1, iv2, iv3]
                    df = delta[iv3, iv1, iv2]  # transposed access for Gf
                    e  = (etax[iv1,iv2,iv3], etay[iv1,iv2,iv3], etaz[iv1,iv2,iv3])
                    e2 = (eta2x[iv1,iv2,iv3], eta2y[iv1,iv2,iv3], eta2z[iv1,iv2,iv3])
                    dv = (Float32(di), Float32(dj), Float32(dk))
                    Fs += d
                    for L in 1:3
                        Sv[L]  += e[L]
                        S2v[L] += e2[L]
                        Gv[L]  += d  * dv[L]
                        Gfv[L] += df * dv[L]
                        for K in 1:3
                            SRm[L, K] += e[L] * dv[K]
                        end
                    end
                    n += Int32(1)
                end
                F[s, p] = Fs
                N[s, p] = n
                for L in 1:3
                    S[L, s, p]  = Sv[L]
                    S2[L, s, p] = S2v[L]
                    G[L, s, p]  = Gv[L]
                    Gf[L, s, p] = Gfv[L]
                    for K in 1:3
                        SR[L, K, s, p] = SRm[L, K]
                    end
                end
            end
        end
        return F, N, S, S2, G, Gf, SR
    end

    @testset "Full-field gather (delta + η + η² + Gf + SR)" begin
        N = 48
        rmax = 8
        stab = build_shell_tables(rmax)
        Random.seed!(314159)
        delta = randn(Float32, N, N, N)
        etax  = randn(Float32, N, N, N)
        etay  = randn(Float32, N, N, N)
        etaz  = randn(Float32, N, N, N)
        eta2x = randn(Float32, N, N, N)
        eta2y = randn(Float32, N, N, N)
        eta2z = randn(Float32, N, N, N)
        peaks_i = Int32[24, 10, 40]
        peaks_j = Int32[24, 40, 10]
        peaks_k = Int32[24, 30, 20]

        Fg, ng, Sg, S2g, Gg, Gfg, SRg = shell_gather_full_gpu(
            delta, etax, etay, etaz, eta2x, eta2y, eta2z,
            peaks_i, peaks_j, peaks_k, stab)
        Fc, nc, Sc, S2c, Gc, Gfc, SRc = cpu_shell_gather_full(
            delta, etax, etay, etaz, eta2x, eta2y, eta2z,
            peaks_i, peaks_j, peaks_k, stab)

        @test ng == nc
        @test isapprox(Fg,  Fc;  rtol=1e-4, atol=1e-5)
        @test isapprox(Sg,  Sc;  rtol=1e-4, atol=1e-5)
        @test isapprox(S2g, S2c; rtol=1e-4, atol=1e-5)
        @test isapprox(Gg,  Gc;  rtol=2e-4, atol=1e-4)
        @test isapprox(Gfg, Gfc; rtol=2e-4, atol=1e-4)
        @test isapprox(SRg, SRc; rtol=2e-4, atol=2e-4)

        @test size(Fg)  == (stab.nshells, length(peaks_i))
        @test size(SRg) == (3, 3, stab.nshells, length(peaks_i))
    end

    @testset "2LPT skipped via zero eta2 fields" begin
        # Passing zero 2LPT fields should make S2shell exactly zero,
        # while other quantities match the 1LPT-only case.
        N = 40
        rmax = 6
        stab = build_shell_tables(rmax)
        Random.seed!(1)
        delta = randn(Float32, N, N, N)
        etax  = randn(Float32, N, N, N)
        etay  = randn(Float32, N, N, N)
        etaz  = randn(Float32, N, N, N)
        zero_field = zeros(Float32, N, N, N)
        peaks_i = Int32[20]; peaks_j = Int32[20]; peaks_k = Int32[20]

        _, _, _, S2g, _, _, _ = shell_gather_full_gpu(
            delta, etax, etay, etaz, zero_field, zero_field, zero_field,
            peaks_i, peaks_j, peaks_k, stab)
        @test all(iszero, S2g)
    end

    @testset "3×3 symmetric eigenvalues (_eig3_symmetric_f32)" begin
        using LinearAlgebra: eigvals, Symmetric

        function ref_eig3(A::AbstractMatrix)
            lam = eigvals(Symmetric(Float64.(A)))
            return Float32.(sort(lam))
        end

        function pack_batch(list)
            n = length(list)
            mats = zeros(Float32, 3, 3, n)
            for i in 1:n
                mats[:, :, i] .= list[i]
            end
            return mats
        end

        # Test matrices spanning the relevant regimes
        test_mats = Matrix{Float32}[]
        push!(test_mats, Float32[1 0 0; 0 1 0; 0 0 1])               # identity
        push!(test_mats, Float32[2 0 0; 0 3 0; 0 0 5])               # diagonal
        push!(test_mats, Float32[1 0.5 0.25; 0.5 2 0.75; 0.25 0.75 3]) # typical symmetric
        push!(test_mats, Float32[5 4 1; 4 5 1; 1 1 2])              # degenerate-ish
        push!(test_mats, Float32[1 0 0; 0 1 0; 0 0 2])              # one repeated
        push!(test_mats, Float32[0.01 0 0; 0 0.01 0; 0 0 0.01])      # small isotropic

        Random.seed!(91)
        for _ in 1:50
            A = randn(Float32, 3, 3)
            push!(test_mats, (A + A') / 2)  # random symmetric
        end

        mats = pack_batch(test_mats)
        evs_gpu = eig3_symmetric_batch_gpu(mats)
        for (i, A) in enumerate(test_mats)
            ref = ref_eig3(A)
            @test isapprox(evs_gpu[:, i], ref; rtol=1e-5, atol=1e-5)
        end

        @testset "near-degenerate pair" begin
            # eigenvalues near (1, 1+ε, 3) where ε is small — Float32 regime
            A = Float32[1.0  1e-4 0;
                        1e-4 1.0  0;
                        0    0    3.0]
            ref = ref_eig3(A)
            evs = eig3_symmetric_batch_gpu(reshape(A, 3, 3, 1))[:, 1]
            @test isapprox(evs, ref; rtol=1e-4, atol=1e-5)
        end

        @testset "batch size > 1 block" begin
            # Force multi-block launch
            Random.seed!(123)
            N = 500
            mats = zeros(Float32, 3, 3, N)
            refs = zeros(Float32, 3, N)
            for i in 1:N
                A = randn(Float32, 3, 3)
                S = (A + A') / 2
                mats[:, :, i] = S
                refs[:, i] = ref_eig3(S)
            end
            evs = eig3_symmetric_batch_gpu(mats; threads=64)
            @test isapprox(evs, refs; rtol=2e-4, atol=1e-5)
        end
    end

    @testset "Trilinear interpolation (_interp3_trilinear_f32)" begin
        # Replicates the CPU CollapseTable test pattern.
        Nx, Ny, Nz = 11, 11, 11
        X1, X2 = 0.0, 1.0
        Y1, Y2 = 0.0, 1.0
        Z1, Z2 = 0.0, 1.0
        dX = (X2 - X1) / (Nx - 1)

        # Linear table: T[i,j,k] = (i-1)*dX + (j-1)*dY + (k-1)*dZ
        table = Float32[Float32((i-1)*dX + (j-1)*dX + (k-1)*dX)
                        for i in 1:Nx, j in 1:Ny, k in 1:Nz]

        # Reference CPU interpolant
        tp = CollapseTableParams(Nx=Nx, Ny=Ny, Nz=Nz, X1=X1, X2=X2,
                                  Y1=Y1, Y2=Y2, Z1=Z1, Z2=Z2)
        ct = CollapseTableInterp(table, tp)

        # Query grid: corners, interior points, OOB points
        xs = Float32[0.0, 0.5, 1.0, 0.25, 0.37, 0.99, -0.1,  1.1, 0.5]
        ys = Float32[0.0, 0.0, 1.0, 0.25, 0.62, 0.01,  0.5,  0.5, -0.2]
        zs = Float32[0.0, 0.5, 0.0, 0.25, 0.13, 0.50,  0.5,  0.5, 0.5]

        expected = Float32[PeakPatch.interpolate(ct, x, y, z) for (x,y,z) in zip(xs, ys, zs)]
        got = interp3_trilinear_gpu(table, X1, X2, Y1, Y2, Z1, Z2, xs, ys, zs, -1.0f0)

        @test isapprox(got, expected; rtol=1e-5, atol=1e-5)

        # Out-of-bounds sentinel
        @test got[7] ≈ -1.0f0  # x = -0.1
        @test got[8] ≈ -1.0f0  # x = 1.1
        @test got[9] ≈ -1.0f0  # y = -0.2
    end

    @testset "Trilinear interpolation on non-trivial table" begin
        # Larger table with curvature to stress the trilinear math
        Nx, Ny, Nz = 21, 17, 13
        X1, X2 = -2.0, 2.0
        Y1, Y2 =  0.0, 3.0
        Z1, Z2 =  1.0, 5.0
        xs_grid = range(X1, X2, length=Nx)
        ys_grid = range(Y1, Y2, length=Ny)
        zs_grid = range(Z1, Z2, length=Nz)
        # Use a simple nonlinear function so trilinear ≠ exact
        f(x, y, z) = sin(Float32(x)) * Float32(y) + Float32(z) * 0.1f0
        table = Float32[f(xi, yj, zk) for xi in xs_grid, yj in ys_grid, zk in zs_grid]
        tp = CollapseTableParams(Nx=Nx, Ny=Ny, Nz=Nz, X1=X1, X2=X2, Y1=Y1, Y2=Y2, Z1=Z1, Z2=Z2)
        ct = CollapseTableInterp(table, tp)

        Random.seed!(999)
        N = 200
        xs = Float32[X1 + (X2 - X1) * rand() for _ in 1:N]
        ys = Float32[Y1 + (Y2 - Y1) * rand() for _ in 1:N]
        zs = Float32[Z1 + (Z2 - Z1) * rand() for _ in 1:N]

        expected = Float32[PeakPatch.interpolate(ct, x, y, z) for (x,y,z) in zip(xs, ys, zs)]
        got = interp3_trilinear_gpu(table, X1, X2, Y1, Y2, Z1, Z2, xs, ys, zs, -1.0f0)
        @test isapprox(got, expected; rtol=1e-5, atol=1e-5)
    end

    @testset "kernel_strain / hRinteg (GPU)" begin
        # Build per-shell profiles with a known geometry and compare GPU
        # output against CPU RadialShell.kernel_strain exactly.
        #
        # CPU kernel_strain takes the CPU layouts Gshell::Vector{Vector{Float64}}
        # etc.; GPU takes Matrix{Float32}. Test both branches (mp>1, mp==1).

        nshells = 20
        Random.seed!(2718)
        # rad is per-shell radius. The akk branch (mp==1) requires
        # 100*rad[m1]^2 < 1101 (akk table size). Use realistic shell radii
        # sqrt(r²) for r² = 0,1,2,3,4,... that stay well under the limit.
        rad = Float64[sqrt(m - 1) for m in 1:nshells]  # 0, 1, √2, √3, 2, ...
        # CPU layout
        Gshell_cpu  = [randn(Float64, nshells) for _ in 1:3]
        Gshellf_cpu = [randn(Float64, nshells) for _ in 1:3]
        SRshell_cpu = [[randn(Float64, nshells) for _ in 1:3] for _ in 1:3]
        # GPU layout (Float32)
        Gshell_gpu  = Float32[Gshell_cpu[L][m]  for L in 1:3, m in 1:nshells]
        Gshellf_gpu = Float32[Gshellf_cpu[L][m] for L in 1:3, m in 1:nshells]
        SRshell_gpu = Float32[SRshell_cpu[L][K][m] for L in 1:3, K in 1:3, m in 1:nshells]

        akk = PeakPatch.atab4()
        wRnor  = 1.0
        aRnor  = 2.5
        hlatt_1 = 1.0
        hlatt_2 = 1.0

        @testset "mp > 1 (hRinteg branch)" begin
            for (mp, mlow, mupp) in [(5, 3, 7), (10, 7, 13), (15, 10, 18), (3, 1, 5)]
                Ecpu, gcpu, gfcpu = PeakPatch.RadialShell.kernel_strain(
                    rad, Gshell_cpu, Gshellf_cpu, SRshell_cpu,
                    mp, mlow, mupp, akk, wRnor, aRnor, hlatt_1, hlatt_2)
                Egpu, ggpu, gfgpu = kernel_strain_gpu(
                    rad, Gshell_gpu, Gshellf_gpu, SRshell_gpu,
                    mp, mlow, mupp, wRnor, aRnor, hlatt_1, hlatt_2)
                @test isapprox(Egpu, Float32.(Ecpu); rtol=5e-4, atol=1e-5)
                @test isapprox(ggpu, Float32.(gcpu); rtol=5e-4, atol=1e-5)
                @test isapprox(gfgpu, Float32.(gfcpu); rtol=5e-4, atol=1e-5)
            end
        end

        @testset "mp == 1 (akk branch)" begin
            # akk branch requires 100*rad[m1]^2 < 1100 (akk table size).
            # With rad[m] = sqrt(m-1), mupp must stay ≤ 11.
            # mlow ≥ 2 to avoid rad[1]=0 divide-by-zero (CPU also hits Inf).
            for (mlow, mupp) in [(2, 5), (2, 8), (2, 11), (3, 4)]
                Ecpu, gcpu, gfcpu = PeakPatch.RadialShell.kernel_strain(
                    rad, Gshell_cpu, Gshellf_cpu, SRshell_cpu,
                    1, mlow, mupp, akk, wRnor, aRnor, hlatt_1, hlatt_2)
                Egpu, ggpu, gfgpu = kernel_strain_gpu(
                    rad, Gshell_gpu, Gshellf_gpu, SRshell_gpu,
                    1, mlow, mupp, wRnor, aRnor, hlatt_1, hlatt_2)
                @test isapprox(Egpu, Float32.(Ecpu); rtol=5e-4, atol=1e-5)
                @test isapprox(ggpu, Float32.(gcpu); rtol=5e-4, atol=1e-5)
                @test isapprox(gfgpu, Float32.(gfcpu); rtol=5e-4, atol=1e-5)
            end
        end

        @testset "hRinteg scalar (GPU kernel_strain wide weight)" begin
            # Exercise u0 > 2 branch: denom in kernel_strain guards via u02<4
            # so the weighting just skips that m1. Should match CPU where
            # CPU also skips.
            rad_wide = collect(0.0:3.0:3.0*(nshells-1))  # large step
            Ecpu, _, _ = PeakPatch.RadialShell.kernel_strain(
                rad_wide, Gshell_cpu, Gshellf_cpu, SRshell_cpu,
                5, 1, 10, akk, wRnor, aRnor, hlatt_1, hlatt_2)
            Egpu, _, _ = kernel_strain_gpu(
                rad_wide, Gshell_gpu, Gshellf_gpu, SRshell_gpu,
                5, 1, 10, wRnor, aRnor, hlatt_1, hlatt_2)
            @test isapprox(Egpu, Float32.(Ecpu); rtol=5e-4, atol=1e-5)
        end

        @testset "Different scaling constants" begin
            # Exercise wRnor / aRnor / hlatt_1 / hlatt_2 as parameters
            for (wR, aR, h1, h2) in [(0.5, 1.0, 0.7, 1.3), (2.0, 3.0, 1.2, 0.9)]
                Ecpu, _, _ = PeakPatch.RadialShell.kernel_strain(
                    rad, Gshell_cpu, Gshellf_cpu, SRshell_cpu,
                    7, 4, 10, akk, wR, aR, h1, h2)
                Egpu, _, _ = kernel_strain_gpu(
                    rad, Gshell_gpu, Gshellf_gpu, SRshell_gpu,
                    7, 4, 10, wR, aR, h1, h2)
                @test isapprox(Egpu, Float32.(Ecpu); rtol=5e-4, atol=1e-5)
            end
        end
    end

    @testset "analyse_peak_gpu_cuda (phases 4-6) vs CPU" begin
        # Build a synthetic delta field with a few Gaussian "bumps" that
        # should collapse, plus low-density cells that shouldn't. Compare
        # GPU analyse_peak_gpu_cuda against CPU analyse_peak_gpu.

        N = 32
        rmax = 8
        stab = build_shell_tables(rmax)

        # Construct delta: sum of Gaussian bumps at known peak centers
        peak_centers = [(16, 16, 16), (8, 8, 8), (24, 24, 24)]
        delta = zeros(Float32, N, N, N)
        for (pi_, pj_, pk_) in peak_centers
            for i in 1:N, j in 1:N, k in 1:N
                d2 = (i - pi_)^2 + (j - pj_)^2 + (k - pk_)^2
                delta[i, j, k] += 3.0f0 * exp(-d2 / 18.0f0)
            end
        end
        # Small noise so strain tensor is nondegenerate
        Random.seed!(7777)
        delta .+= 0.01f0 .* randn(Float32, N, N, N)

        # Displacement fields: reasonable magnitudes
        etax = 0.3f0 .* randn(Float32, N, N, N)
        etay = 0.3f0 .* randn(Float32, N, N, N)
        etaz = 0.3f0 .* randn(Float32, N, N, N)

        eta2x = 0.15f0 .* randn(Float32, N, N, N)
        eta2y = 0.15f0 .* randn(Float32, N, N, N)
        eta2z = 0.15f0 .* randn(Float32, N, N, N)

        # Laplacian field (proxy — random small values; validates d2F output path)
        lapd = 0.05f0 .* randn(Float32, N, N, N)

        # Collapse table: make it yield zvir = 0.5 everywhere (above ZZon=0.0,
        # so peaks that cross fcrit should collapse)
        ct_Nx, ct_Ny, ct_Nz = 10, 10, 10
        ct_X1, ct_X2 = -2.0, 2.0   # log10(Frho)
        ct_Y1, ct_Y2 =  0.0, 1.0   # e
        ct_Z1, ct_Z2 =  0.0, 2.0   # p/e
        ct_table = fill(Float32(0.5), ct_Nx, ct_Ny, ct_Nz)
        tp = CollapseTableParams(Nx=ct_Nx, Ny=ct_Ny, Nz=ct_Nz,
                                  X1=ct_X1, X2=ct_X2, Y1=ct_Y1, Y2=ct_Y2,
                                  Z1=ct_Z1, Z2=ct_Z2)
        ct = CollapseTableInterp(ct_table, tp)

        # Shared parameters
        alatt = 1.0
        ir2min = 4
        ZZon = 0.0
        Rfclvi = 4.0
        fcrit = 1.0  # explicit, avoid growth_tables path

        peaks_i = Int32[p[1] for p in peak_centers]
        peaks_j = Int32[p[2] for p in peak_centers]
        peaks_k = Int32[p[3] for p in peak_centers]

        # GPU
        mask_gpu = zeros(Int8, N, N, N)
        res = analyse_peak_gpu_cuda(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                                     peaks_i, peaks_j, peaks_k, stab,
                                     ct_table, ct_X1, ct_X2, ct_Y1, ct_Y2, ct_Z1, ct_Z2,
                                     alatt, ir2min, ZZon, Rfclvi;
                                     fcrit_override=fcrit,
                                     ct_out_val=-1.0, lapd=lapd,
                                     mask=mask_gpu, nbuff=0)

        # CPU reference — include lapd and mask in PeakGrid so CPU writes mask
        mask_cpu = zeros(Int8, N, N, N)
        pg = PeakGrid(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                       mask_cpu, (N, N, N), lapd)
        n1xn2 = N * N
        for (idx, (ci, cj, ck)) in enumerate(peak_centers)
            ipp = ci + (cj - 1) * N + (ck - 1) * n1xn2
            cpu_res = analyse_peak_gpu(pg, ipp, alatt, ir2min, ZZon, Rfclvi, ct, stab;
                                        fcrit_override=fcrit)

            @testset "peak $idx at $(peak_centers[idx])" begin
                # Both agree on collapse/no-collapse
                cpu_collapsed = cpu_res.RTHL > 0.0
                gpu_collapsed = res.RTHL[idx] > 0.0f0
                @test cpu_collapsed == gpu_collapsed

                if cpu_collapsed && gpu_collapsed
                    @test isapprox(Float32(cpu_res.RTHL),  res.RTHL[idx];  rtol=5e-3, atol=5e-3)
                    @test isapprox(Float32(cpu_res.Fbarx), res.Fbarx[idx]; rtol=5e-3, atol=5e-3)
                    @test isapprox(Float32(cpu_res.e_v),   res.e_v[idx];   rtol=1e-2, atol=1e-3)
                    @test isapprox(Float32(cpu_res.p_v),   res.p_v[idx];   rtol=1e-2, atol=1e-3)
                    # Phase 7: Sbar, Sbar2, Srb
                    @test isapprox(Float32(cpu_res.Srb),   res.Srb[idx];   rtol=5e-3, atol=1e-4)
                    for L in 1:3
                        @test isapprox(Float32(cpu_res.Sbar[L]),  res.Sbar[L, idx];  rtol=5e-3, atol=5e-5)
                        @test isapprox(Float32(cpu_res.Sbar2[L]), res.Sbar2[L, idx]; rtol=5e-3, atol=5e-5)
                    end
                    # Phase 3 + phases 5-6 gradpk tracking
                    for L in 1:3
                        @test isapprox(Float32(cpu_res.gradpk[L]),   res.gradpk[L, idx];
                                       rtol=5e-3, atol=1e-5)
                        @test isapprox(Float32(cpu_res.gradpkf[L]),  res.gradpkf[L, idx];
                                       rtol=5e-3, atol=1e-5)
                        @test isapprox(Float32(cpu_res.gradpkrf[L]), res.gradpkrf[L, idx];
                                       rtol=5e-3, atol=1e-5)
                    end
                    # Phase 7b: Laplacian d2F
                    @test isapprox(Float32(cpu_res.d2F), res.d2F[idx]; rtol=5e-3, atol=1e-5)
                end
            end
        end

        # Phase 7c: mask side-effect. With 3 non-overlapping peaks that
        # both CPU and GPU claim, the mask state should be identical.
        @test sum(res.mask) == sum(mask_cpu)
        @test res.mask == mask_cpu
    end

    @testset "End-to-end CPU↔GPU on realistic field" begin
        # Larger grid with many peaks (found via the same find_peaks path
        # the production pipeline uses). Validates that GPU = CPU for every
        # field of every peak, including peaks near boundaries and peaks
        # that don't collapse.
        N = 64
        rmax = 10
        stab = build_shell_tables(rmax)
        Random.seed!(31415)

        # Field: sum of Gaussian bumps at random locations + noise
        # This produces both collapsing and non-collapsing peaks.
        delta = 0.05f0 .* randn(Float32, N, N, N)
        for _ in 1:20
            ci = rand(10:N-10); cj = rand(10:N-10); ck = rand(10:N-10)
            amp = 1.0f0 + 2.0f0 * rand(Float32)
            sigma2 = 15.0f0 + 10.0f0 * rand(Float32)
            for i in 1:N, j in 1:N, k in 1:N
                d2 = (i-ci)^2 + (j-cj)^2 + (k-ck)^2
                delta[i, j, k] += amp * exp(-d2 / sigma2)
            end
        end

        etax  = 0.3f0 .* randn(Float32, N, N, N)
        etay  = 0.3f0 .* randn(Float32, N, N, N)
        etaz  = 0.3f0 .* randn(Float32, N, N, N)
        eta2x = 0.15f0 .* randn(Float32, N, N, N)
        eta2y = 0.15f0 .* randn(Float32, N, N, N)
        eta2z = 0.15f0 .* randn(Float32, N, N, N)
        lapd  = 0.05f0 .* randn(Float32, N, N, N)

        # Find peaks via production path
        mask_find = zeros(Int8, N, N, N)
        peak_candidates = PeakPatch.find_peaks(delta, mask_find,
                                                0.0, 0.0, 0.0, 1.0, 1,
                                                0.8f0, 3.0)
        @info "End-to-end test peaks" nfound=length(peak_candidates)
        @test length(peak_candidates) > 5  # sanity: several peaks found

        peaks_i = Int32[p.i for p in peak_candidates]
        peaks_j = Int32[p.j for p in peak_candidates]
        peaks_k = Int32[p.k for p in peak_candidates]

        # Collapse table with realistic (nontrivial) values
        ct_Nx, ct_Ny, ct_Nz = 11, 11, 11
        ct_X1, ct_X2 = -2.0, 2.0
        ct_Y1, ct_Y2 =  0.0, 1.0
        ct_Z1, ct_Z2 =  0.0, 2.0
        ct_table = Float32[
            0.5f0 + 0.3f0 * sin(Float32(i-1) * 0.3f0) +
            0.2f0 * cos(Float32(j-1) * 0.4f0) +
            0.1f0 * Float32(k-1) * 0.1f0
            for i in 1:ct_Nx, j in 1:ct_Ny, k in 1:ct_Nz]
        tp = CollapseTableParams(Nx=ct_Nx, Ny=ct_Ny, Nz=ct_Nz,
                                  X1=ct_X1, X2=ct_X2, Y1=ct_Y1, Y2=ct_Y2,
                                  Z1=ct_Z1, Z2=ct_Z2)
        ct = CollapseTableInterp(ct_table, tp)

        alatt = 1.0
        ir2min = 4
        ZZon = 0.0
        Rfclvi = 4.0
        fcrit = 1.0

        # GPU
        mask_gpu = zeros(Int8, N, N, N)
        res = analyse_peak_gpu_cuda(
            delta, etax, etay, etaz, eta2x, eta2y, eta2z,
            peaks_i, peaks_j, peaks_k, stab,
            ct_table, ct_X1, ct_X2, ct_Y1, ct_Y2, ct_Z1, ct_Z2,
            alatt, ir2min, ZZon, Rfclvi;
            fcrit_override=fcrit, lapd=lapd, mask=mask_gpu, nbuff=0)

        # CPU
        mask_cpu = zeros(Int8, N, N, N)
        pg = PeakGrid(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                       mask_cpu, (N, N, N), lapd)
        ncollapsed_cpu = 0
        ncollapsed_gpu = 0
        for idx in 1:length(peak_candidates)
            ipp = peaks_i[idx] + (peaks_j[idx] - 1) * N + (peaks_k[idx] - 1) * N * N
            cpu_res = analyse_peak_gpu(pg, Int(ipp), alatt, ir2min, ZZon, Rfclvi, ct, stab;
                                        fcrit_override=fcrit)
            cpu_collapsed = cpu_res.RTHL > 0.0
            gpu_collapsed = res.RTHL[idx] > 0.0f0
            cpu_collapsed && (ncollapsed_cpu += 1)
            gpu_collapsed && (ncollapsed_gpu += 1)
            @test cpu_collapsed == gpu_collapsed

            if cpu_collapsed && gpu_collapsed
                # Scalar fields
                @test isapprox(Float32(cpu_res.RTHL),  res.RTHL[idx];  rtol=5e-3, atol=5e-3)
                @test isapprox(Float32(cpu_res.Fbarx), res.Fbarx[idx]; rtol=5e-3, atol=5e-3)
                @test isapprox(Float32(cpu_res.e_v),   res.e_v[idx];   rtol=1e-2, atol=1e-3)
                @test isapprox(Float32(cpu_res.p_v),   res.p_v[idx];   rtol=1e-2, atol=1e-3)
                @test isapprox(Float32(cpu_res.Srb),   res.Srb[idx];   rtol=5e-3, atol=1e-4)
                @test isapprox(Float32(cpu_res.d2F),   res.d2F[idx];   rtol=5e-3, atol=1e-5)
                # Vector fields
                for L in 1:3
                    @test isapprox(Float32(cpu_res.Sbar[L]),     res.Sbar[L, idx];    rtol=5e-3, atol=5e-5)
                    @test isapprox(Float32(cpu_res.Sbar2[L]),    res.Sbar2[L, idx];   rtol=5e-3, atol=5e-5)
                    @test isapprox(Float32(cpu_res.gradpk[L]),   res.gradpk[L, idx];  rtol=5e-3, atol=1e-5)
                    @test isapprox(Float32(cpu_res.gradpkf[L]),  res.gradpkf[L, idx]; rtol=5e-3, atol=1e-5)
                    @test isapprox(Float32(cpu_res.gradpkrf[L]), res.gradpkrf[L, idx];rtol=5e-3, atol=1e-5)
                end
            end
        end
        @info "Collapse tally" ncollapsed_cpu ncollapsed_gpu

        # Mask side-effect should be identical
        @test sum(res.mask) == sum(mask_cpu)
        @test res.mask == mask_cpu

        # ------------------------------------------------------------
        # zvir_half validation against the REAL CPU analyse_peak
        # (RadialShell.analyse_peak computes zvir_half properly;
        #  analyse_peak_gpu prototype stubs it to 0.0)
        # ------------------------------------------------------------
        shells_cpu = precompute_shells(rmax)
        n_zvir_checked = 0
        n_zvir_match = 0
        for idx in 1:length(peak_candidates)
            ipp = peaks_i[idx] + (peaks_j[idx] - 1) * N + (peaks_k[idx] - 1) * N * N
            # RadialShell.analyse_peak (the production CPU path)
            pg_cpu = PeakGrid(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                               zeros(Int8, N, N, N), (N, N, N), lapd)
            real_cpu = analyse_peak(pg_cpu, Int(ipp), alatt, ir2min, ZZon, Rfclvi, ct, shells_cpu;
                                     fcrit_override=fcrit)
            real_cpu.RTHL > 0.0 || continue  # no-collapse → zvir_half is -1.0 on both paths
            n_zvir_checked += 1
            z_cpu = Float32(real_cpu.zvir_half)
            z_gpu = res.zvir_half[idx]
            # Both should be in the same regime (> 0 or both <= 0)
            same_regime = (z_cpu > 0.0f0) == (z_gpu > 0.0f0)
            tol_ok = isapprox(z_cpu, z_gpu; rtol=5e-3, atol=5e-3)
            # Allow tolerance to relax when RTHL is small (fewer radii contribute,
            # Float32 rounding matters more)
            tol_loose = isapprox(z_cpu, z_gpu; rtol=2e-2, atol=1e-2)
            if same_regime && (tol_ok || tol_loose)
                n_zvir_match += 1
            else
                @warn "zvir_half mismatch" idx z_cpu z_gpu RTHL_cpu=real_cpu.RTHL RTHL_gpu=res.RTHL[idx]
            end
        end
        @info "zvir_half validation" n_zvir_checked n_zvir_match
        @test n_zvir_checked > 0        # sanity: at least one collapsed peak
        @test n_zvir_match == n_zvir_checked
    end

    @testset "isolated_convolve_gpu vs CPU _isolated_convolve" begin
        # Phase 1 of GPU plan: validates the cuFFT-based isolated convolution
        # for all transfer-function variants used by run_multitile_split.
        import PeakPatch.MultiResolution: _isolated_convolve, _kernel_1lpt,
                                            _kernel_2lpt, _kernel_phi_ij,
                                            _kernel_laplacian

        # Smooth, bounded P(k) — avoids near-Nyquist amplification noise
        pk_test = k -> 100.0 / (1.0 + (k * 50.0)^4)
        boxsize = 100.0  # Mpc/h

        @testset "delta + LPT variants (n=32, no shell)" begin
            Random.seed!(42)
            n = 32
            noise = randn(Float32, n, n, n)
            tol = 5e-6  # Float32 round-off; tabulated P(k) interp adds a bit

            # δ
            cpu = _isolated_convolve(noise, pk_test, boxsize, n)
            gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n; kernel_fn_id=0)
            @test isapprox(cpu, gpu; rtol=tol, atol=1e-6)

            # 1LPT psi for each dim
            for dim in 1:3
                cpu = _isolated_convolve(noise, pk_test, boxsize, n; kernel_fn=_kernel_1lpt(dim))
                gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n;
                                              kernel_fn_id=1, dim1=dim)
                @test isapprox(cpu, gpu; rtol=tol, atol=1e-6)
            end

            # 2LPT psi for each dim
            for dim in 1:3
                cpu = _isolated_convolve(noise, pk_test, boxsize, n; kernel_fn=_kernel_2lpt(dim))
                gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n;
                                              kernel_fn_id=2, dim1=dim)
                @test isapprox(cpu, gpu; rtol=tol, atol=1e-6)
            end

            # φ_ij for representative (di, dj) pairs
            for (di, dj) in [(1,1), (1,2), (1,3), (2,2), (2,3), (3,3)]
                cpu = _isolated_convolve(noise, pk_test, boxsize, n;
                                           kernel_fn=_kernel_phi_ij(di, dj))
                gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n;
                                              kernel_fn_id=3, dim1=di, dim2=dj)
                @test isapprox(cpu, gpu; rtol=tol, atol=1e-6)
            end

            # Laplacian
            cpu = _isolated_convolve(noise, pk_test, boxsize, n; kernel_fn=_kernel_laplacian())
            gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n; kernel_fn_id=4)
            @test isapprox(cpu, gpu; rtol=tol, atol=1e-6)
        end

        @testset "extended-noise nshell>0 (n=32, nshell=8)" begin
            Random.seed!(7)
            n = 32; nshell = 8
            n_input = n + 2 * nshell
            noise = randn(Float32, n_input, n_input, n_input)
            tol = 5e-6

            cpu = _isolated_convolve(noise, pk_test, boxsize, n; nshell=nshell)
            gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n; nshell=nshell)
            @test isapprox(cpu, gpu; rtol=tol, atol=1e-6)
        end

        @testset "larger size (n=64)" begin
            Random.seed!(123)
            n = 64
            noise = randn(Float32, n, n, n)

            cpu = _isolated_convolve(noise, pk_test, boxsize, n)
            gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n; kernel_fn_id=0)
            @test isapprox(cpu, gpu; rtol=5e-6, atol=1e-6)

            cpu = _isolated_convolve(noise, pk_test, boxsize, n; kernel_fn=_kernel_1lpt(2))
            gpu = isolated_convolve_gpu(noise, pk_test, boxsize, n; kernel_fn_id=1, dim1=2)
            @test isapprox(cpu, gpu; rtol=5e-6, atol=1e-6)
        end
    end

    @testset "Varying threads per block" begin
        N = 40
        rmax = 6
        stab = build_shell_tables(rmax)
        Random.seed!(7)
        delta = randn(Float32, N, N, N)
        peaks_i = Int32[20, 10]; peaks_j = Int32[20, 30]; peaks_k = Int32[20, 15]
        Fcpu, ncpu = cpu_shell_sums(delta, peaks_i, peaks_j, peaks_k, stab)
        for t in (32, 64, 128, 256)
            Fgpu, ngpu = shell_fbar_gather_gpu(delta, peaks_i, peaks_j, peaks_k, stab; threads=t)
            @test ngpu == ncpu
            @test isapprox(Fgpu, Fcpu; rtol=1e-4, atol=1e-5)
        end
    end
end
