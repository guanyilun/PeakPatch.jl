# Post-hoc pooled permutation test of lag correlations (C_orig + C_confirm samples).
# Run: julia --project=validation validation/rng/posthoc_lag_pooled.jl  (reads rng_test/ on scratch)
using DelimitedFiles, Statistics, Random, Printf
D = "/home/yguan/scratch/websky_6144/rng_test"
ld(f) = (m = readdlm(joinpath(D, f), ','); (String.(m[:, 1]), Float64.(m[:, 2:end])))
n1, T1 = ld("samples__orig_threefry.csv"); _, X1 = ld("samples__orig_xoshiro.csv")
n2, T2 = ld("samples__confirm_threefry.csv"); _, X2 = ld("samples__confirm_xoshiro.csv")
@assert n1 == n2
T = hcat(T1, T2); X = hcat(X1, X2); Z = hcat(T, X); nt = size(T, 2)
chi(A, B, rows) = sum(((mean(A[rows, :]; dims=2) .- mean(B[rows, :]; dims=2)) ./
                       sqrt.(var(A[rows, :]; dims=2) ./ size(A, 2) .+ var(B[rows, :]; dims=2) ./ size(B, 2))) .^ 2)
rng = Xoshiro(7)
for ax in ("x", "y", "z")
    rows = findall(s -> startswith(s, "xi_$(ax)("), n1)
    c0 = chi(T, X, rows); ge = 0; np = 20000
    for _ in 1:np
        p = randperm(rng, size(Z, 2)); ge += chi(Z[:, p[1:nt]], Z[:, p[nt+1:end]], rows) >= c0
    end
    @printf("xi_%s lags (%d): pooled %d Threefry vs %d Xoshiro  chi2 = %.1f  permutation p = %.3f\n",
            ax, length(rows), nt, size(X, 2), c0, (ge + 1) / (np + 1))
end
rows = findall(s -> startswith(s, "xi_x("), n1)
@printf("mean xi_x over lags: Threefry %+.3f  Xoshiro %+.3f (white-noise σ units)\n", mean(T[rows, :]), mean(X[rows, :]))
