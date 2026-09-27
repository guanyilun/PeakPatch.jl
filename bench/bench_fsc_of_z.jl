using PeakPatch

# Load collapse table
const BASEDIR = dirname(dirname(@__FILE__))
ct_array, ct_params = PeakPatch.read_homeltab(joinpath(BASEDIR, "validation/websky_6144/data/HomelTab_websky.dat"))
ct = PeakPatch.CollapseTableInterp(ct_array, ct_params)

println("Table params: X1=$(ct_params.X1), X2=$(ct_params.X2), Nx=$(ct_params.Nx)")
println("Walk range: $(ct.x2) down to $(ct_params.X1), step 1e-4 => max $(round(Int, (ct.x2 - ct_params.X1)/1e-4)) steps")

# Also get growth tables for comparison
cosmo = PeakPatch.CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.808)
growth_tables = PeakPatch.Dlinear_tables(cosmo)

# Benchmark new fsc_of_z (table walk)
N = 100_000
zvals = rand(N) .* 4.0  # z from 0 to 4

# Warmup
for z in zvals[1:100]
    PeakPatch.fsc_of_z(z, ct)
end

t0 = time()
for z in zvals
    PeakPatch.fsc_of_z(z, ct)
end
t_new = time() - t0
println("\nNew fsc_of_z (table walk): $(round(t_new*1e6/N; digits=1)) us/call, $(round(t_new; digits=3))s for $N calls")

# Benchmark old formula for comparison
for z in zvals[1:100]
    a = 1.0 / (1.0 + z)
    dlin, = PeakPatch.Cosmology.Dlinear_ab(a, growth_tables)
    _ = 1.686 * dlin
end

t0 = time()
for z in zvals
    a = 1.0 / (1.0 + z)
    dlin, = PeakPatch.Cosmology.Dlinear_ab(a, growth_tables)
    _ = 1.686 * dlin
end
t_old = time() - t0
println("Old fsc_of_z (D(z) mult): $(round(t_old*1e6/N; digits=1)) us/call, $(round(t_old; digits=3))s for $N calls")
println("Slowdown: $(round(t_new/t_old; digits=1))x")

# Estimate impact
npk_total = 968_000
println("\nEstimated overhead for $npk_total peaks: $(round(t_new/N * npk_total; digits=1))s")
println("Estimated overhead for 10M peaks: $(round(t_new/N * 10_000_000; digits=1))s")

# Bisection approach
function fsc_of_z_bisect(z, ct, x1)
    target = 1.0 + z
    lo = x1
    hi = ct.x2
    for _ in 1:50
        mid = (lo + hi) / 2
        zvir = PeakPatch.CollapseTable.interpolate(ct, mid, 0.0, 0.0)
        if zvir > target
            hi = mid
        else
            lo = mid
        end
    end
    return 10.0^((lo + hi) / 2)
end

x1 = ct_params.X1

# Warmup
for z in zvals[1:100]
    fsc_of_z_bisect(z, ct, x1)
end

t0 = time()
for z in zvals
    fsc_of_z_bisect(z, ct, x1)
end
t_bisect = time() - t0
println("\nBisection fsc_of_z: $(round(t_bisect*1e6/N; digits=1)) us/call, $(round(t_bisect; digits=3))s for $N calls")
println("Bisection vs linear walk: $(round(t_new/t_bisect; digits=1))x faster")

# Verify bisection accuracy
let max_err = 0.0
    for z in zvals[1:1000]
        v1 = PeakPatch.fsc_of_z(z, ct)
        v2 = fsc_of_z_bisect(z, ct, x1)
        max_err = max(max_err, abs(v1 - v2) / v1)
    end
    println("Max relative error (bisection vs walk): $(max_err)")
end

# How many steps does the walk actually take on average?
function count_walk_steps(z, ct)
    target = 1.0 + z
    fv = ct.x2
    dfv = 1e-4
    steps = 0
    zvir = 100.0
    while zvir > target
        zvir = PeakPatch.CollapseTable.interpolate(ct, fv, 0.0, 0.0)
        fv -= dfv
        steps += 1
    end
    return steps
end

avg_steps = sum(count_walk_steps(z, ct) for z in zvals[1:1000]) / 1000
println("\nAverage walk steps (z uniform 0-4): $(round(avg_steps; digits=0))")
println("Walk steps at z=0: $(count_walk_steps(0.0, ct))")
println("Walk steps at z=2: $(count_walk_steps(2.0, ct))")
println("Walk steps at z=4: $(count_walk_steps(4.0, ct))")
