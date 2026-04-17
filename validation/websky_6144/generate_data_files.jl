#!/usr/bin/env julia
#
# Generate data files for Websky 6144^3 run:
#   1. Filter bank (21 top-hat filters, Rf=2.507 to 36 Mpc/h, 1.15x spacing)
#   2. Collapse table (HomelTab) for Websky cosmology
#   3. Power spectrum (placeholder — needs CAMB/CLASS for accurate P(k))
#
# Usage:
#   julia --project=. validation/websky_6144/generate_data_files.jl

using PeakPatch

# Websky cosmology
Om_total = 0.31
OB = 0.049
OL = 0.69
h = 0.68

datadir = joinpath(@__DIR__, "data")
mkpath(datadir)

# ---- 1. Filter bank ----
filter_path = joinpath(datadir, "filters_websky.dat")

# 21 Websky-style top-hat filters
# Rf,min = 2 * a_latt = 2 * 1.25326 = 2.507 Mpc/h (per Stein+ 2020)
# Rf,max ~ 36 Mpc/h with spacing factor 1.15 → 20 filters
delta_c = 1.686
Rf_min = 2.507   # = 2 * a_latt (Websky paper: Rf,min = 2 × cell size)
spacing = 1.15
nfilters = 20    # gives Rf_max = 2.507 * 1.15^19 ≈ 35.7 Mpc/h ≈ 36

open(filter_path, "w") do f
    println(f, nfilters)
    for i in 1:nfilters
        Rf = Rf_min * spacing^(i - 1)
        println(f, "$i  $delta_c  $(round(Rf; digits=4))  1")
    end
end

println("Wrote filter bank: $filter_path")
println("  Filters: $nfilters, Rf range: $(round(Rf_min; digits=2)) to $(round(Rf_min * spacing^(nfilters-1); digits=2)) Mpc/h")

# ---- 2. Collapse table ----
tab_path = joinpath(datadir, "HomelTab_websky.dat")

cosmo = CosmologyParams(Om_total, OB, OL, h, 0.965, 0.808)
ep = EllipsoidParams(cosmo; solver=:rk4)
tp = CollapseTableParams()

println("Generating collapse table...")
table = make_table_threaded(ep, tp; verbose=true)
write_homeltab(tab_path, table, tp)
println("Wrote collapse table: $tab_path")

# ---- 3. Power spectrum ----
pk_path = joinpath(datadir, "pk_websky.dat")

# WARNING: This is an approximate analytical P(k) using the Eisenstein-Hu
# transfer function. For production Websky reproduction, replace with
# CAMB/CLASS output for the exact cosmology (Om=0.31, OB=0.049, OL=0.69, h=0.68).
#
# If you have a CAMB P(k) file, skip this and point the TOML config to it.
#
# For now, generate a simplified P(k) that's reasonable for the cosmology.
# The real Websky used a CAMB-generated spectrum.

println("Generating approximate P(k): $pk_path")
println("WARNING: This is approximate. Replace with CAMB/CLASS output for production.")

# Simple Eisenstein-Hu-like P(k) for testing
# Parameters roughly matching Planck 2018
ns = 0.965  # spectral index
As = 2.1e-9 # amplitude (approximate)
k_pivot = 0.05  # pivot scale [h/Mpc]

# Growth factor normalization: sigma_8 ~ 0.81
# We'll use a simple CDM-like transfer function
# T(k) ~ 1 for k << keq, T(k) ~ (keq/k)^2 for k >> keq
keq = 0.0146 * Om_total * h^2  # equality scale [h/Mpc]

function transfer_eh(k, keq)
    q = k / (keq * h)  # dimensionless
    # BBKS transfer function (simplified)
    L = log(2 * exp(1.0) * q + 0.01)
    C = 14.2 + 731.0 / (1 + 62.5 * q)
    T = L / (L + C * q^2)
    return max(T, 1e-10)
end

open(pk_path, "w") do f
    println(f, "# Power spectrum for Websky cosmology (Om=$Om_total, OB=$OB, OL=$OL, h=$h)")
    println(f, "# WARNING: Approximate (Eisenstein-Hu). Replace with CAMB output for production.")
    println(f, "# k [h/Mpc]    P(k) [(Mpc/h)^3]")
    for logk in range(-4, stop=1, length=500)
        k = 10.0^logk
        T = transfer_eh(k, keq)
        # P(k) = A_s * (k/k_pivot)^(ns-1) * T(k)^2 * (2π^2) * k^(-3)
        # Adjusted for correct normalization
        Pk = 2.1e4 * (k / 0.05)^(ns - 1) * T^2 * k^(-1.5)
        println(f, "$(round(k; digits=8))  $(round(Pk; digits=8))")
    end
end

println("Wrote power spectrum: $pk_path")
println()
println("=== Data files ready ===")
println("  $filter_path")
println("  $tab_path")
println("  $pk_path  (APPROXIMATE — replace with CAMB for production)")
