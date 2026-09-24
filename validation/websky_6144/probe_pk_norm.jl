#!/usr/bin/env julia
# Find the exact P(k) normalization factor. Build the global field from pk_websky
# divided by candidate factors and compare field sigma(R) to the THEORY sigma(R)
# (standard integral of the raw physical P(k)). The divisor giving ratio ~1.0 wins.
#
# Usage: julia --project=validation -t 8 validation/websky_6144/probe_pk_norm.jl

using PeakPatch, FFTW, Statistics, Printf, QuadGK

N = 256; alatt = 1.25326; box = N * alatt; seed = 13579
pkpath = joinpath(@__DIR__, "data", "pk_websky.dat")
pk_raw = PeakPatch.PowerSpectrum.load_pk(pkpath)   # physical (Mpc/h)^3, standard convention

# theory sigma(R): standard integral  sigma^2 = ∫ k^2 P(k) W(kR)^2 dk / (2π^2)
W(x) = x < 1e-4 ? 1.0 : 3*(sin(x) - x*cos(x))/x^3
function theory_sigma(R)
    f(k) = k^2 * pk_raw(k) * W(k*R)^2 / (2π^2)
    s2, _ = quadgk(f, 1e-4, 50.0; rtol=1e-4)
    sqrt(s2)
end

function field_sigma(divisor)
    pk = k -> pk_raw(k) / divisor
    d = PeakPatch.RandomField.generate_grf(N, pk, box, seed)
    dk = rfft(d)
    [std(PeakPatch.Filters.smooth_field(dk, N, box, R, 1)) for R in (2.0,4.0,8.0,16.0,32.0)]
end

Rs = (2.0,4.0,8.0,16.0,32.0)
th = [theory_sigma(R) for R in Rs]
h = 0.68
candidates = [("raw (÷1)", 1.0), ("÷(2π)³", (2π)^3), ("÷(2πh)³", (2π*h)^3)]

@printf("\n%-10s", "R[Mpc/h]"); for R in Rs; @printf(" | %8.1f", R); end; println()
@printf("%-10s", "theory σ"); for t in th; @printf(" | %8.4f", t); end; println()
println("-"^70)
for (name, div) in candidates
    fs = field_sigma(div)
    @printf("%-10s", name); for s in fs; @printf(" | %8.4f", s); end; println()
    @printf("%-10s", "  ratio");  for (s,t) in zip(fs,th); @printf(" | %8.2f", s/t); end; println()
end
println("\nThe divisor whose 'ratio' row is ~1.0 across all R is the correct factor.")
@printf("(2π)³ = %.3f   (2πh)³ = %.3f   (2π)^1.5 = %.3f\n", (2π)^3, (2π*h)^3, (2π)^1.5)
