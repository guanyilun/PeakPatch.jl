module MassFunction

using QuadGK
import ..Cosmology: CosmologyParams, growth_factor

# ---- Mass-radius conversions ----

"""Mean matter density ρ̄_m in M_sun/h / (Mpc/h)^3."""
rho_mean(Om::Real) = 2.775e11 * Om

"""Lagrangian radius R [Mpc/h] enclosing mass M [M_sun/h]."""
R_of_M(M::Real, Om::Real) = (3.0 * M / (4π * rho_mean(Om)))^(1.0/3.0)

"""Mass [M_sun/h] enclosed in Lagrangian radius R [Mpc/h]."""
M_of_R(R::Real, Om::Real) = (4π/3.0) * rho_mean(Om) * R^3

# ---- Variance σ(R) ----

"""Top-hat window function in Fourier space: W(kR) = 3(sin(kR) - kR cos(kR)) / (kR)³."""
function W_tophat(kR::Real)
    kR < 1e-6 && return 1.0
    return 3.0 * (sin(kR) - kR * cos(kR)) / kR^3
end

"""
    sigma_R(R, pk; kmin=1e-5, kmax=100.0, nk=10000)

Compute σ(R) from a tabulated P(k) interpolator using log-spaced numerical integration.

σ²(R) = 1/(2π²) ∫ dk k² P(k) W²(kR)
"""
function sigma_R(R::Real, pk; kmin::Real=1e-5, kmax::Real=100.0, nk::Int=10000)
    logk = range(log(kmin), log(kmax); length=nk)
    dk = logk[2] - logk[1]
    s2 = 0.0
    for i in eachindex(logk)
        k = exp(logk[i])
        W = W_tophat(k * R)
        s2 += k^3 * pk(k) * W^2 * dk
    end
    return sqrt(s2 / (2π^2))
end

"""
    sigma_M(M, pk, Om; kwargs...)

Convenience: σ(M) = σ(R(M)) at z=0.
"""
sigma_M(M::Real, pk, Om::Real; kwargs...) = sigma_R(R_of_M(M, Om), pk; kwargs...)

"""
    dlnsigma_dlnM(M, pk, Om; dlogM=0.01, kwargs...)

Numerical derivative d ln σ / d ln M via central finite difference.
"""
function dlnsigma_dlnM(M::Real, pk, Om::Real; dlogM::Real=0.01, kwargs...)
    dM = M * dlogM
    sig_plus = sigma_M(M + dM, pk, Om; kwargs...)
    sig_minus = sigma_M(M - dM, pk, Om; kwargs...)
    return (log(sig_plus) - log(sig_minus)) / (log(M + dM) - log(M - dM))
end

# ---- Mass functions ----

"""
    tinker_dndlnM(M, sigma, dlnsig_dlnM, z, Om)

Tinker et al. (2008) mass function dn/dlnM [(Mpc/h)⁻³] for Δ=200.

Uses Eqs. 3, 5-8 with parameters from Table 2 (Δ=200 row).
Redshift evolution from Eqs. 5-8.
"""
function tinker_dndlnM(M::Real, sigma::Real, dlnsig_dlnM::Real, z::Real, Om::Real)
    # Table 2, Δ=200 parameters (z=0 values)
    A0 = 0.186
    a0 = 1.47
    b0 = 2.57
    c0 = 1.19

    # Redshift evolution (Eqs. 5-8)
    A = A0 * (1.0 + z)^(-0.14)
    a = a0 * (1.0 + z)^(-0.06)
    alpha = 10.0^(-1.0 * (0.75 / log10(200.0 / 75.0))^1.2)
    b = b0 * (1.0 + z)^(-alpha)
    c = c0  # no z-evolution for c

    # f(σ) from Eq. 3
    f_sigma = A * ((sigma / b)^(-a) + 1.0) * exp(-c / sigma^2)

    ρ_m = rho_mean(Om)
    return (ρ_m / M) * f_sigma * abs(dlnsig_dlnM)
end

"""
    sheth_tormen_dndlnM(M, sigma, dlnsig_dlnM, z, Om)

Sheth-Tormen (1999) mass function dn/dlnM [(Mpc/h)⁻³].

Parameters: A=0.3222, a=0.707, p=0.3, δ_c=1.686.
"""
function sheth_tormen_dndlnM(M::Real, sigma::Real, dlnsig_dlnM::Real, z::Real, Om::Real)
    A = 0.3222
    a = 0.707
    p = 0.3
    δ_c = 1.686

    ν = δ_c / sigma
    ν2 = a * ν^2
    f_nu = A * sqrt(2.0 * a / π) * (1.0 + ν2^(-p)) * exp(-ν2 / 2.0)

    ρ_m = rho_mean(Om)
    return (ρ_m / M) * f_nu * abs(dlnsig_dlnM)
end

"""
    cumulative_ngtm(Mcent, dndlnM, dlnM)

Compute cumulative N(>M) from dn/dlnM by summing from high M to low M.
Returns a vector of the same length as `Mcent`.
"""
function cumulative_ngtm(dndlnM::AbstractVector, dlnM::Real)
    ngtm = similar(dndlnM)
    ngtm[end] = dndlnM[end] * dlnM
    for i in (length(dndlnM)-1):-1:1
        ngtm[i] = ngtm[i+1] + dndlnM[i] * dlnM
    end
    return ngtm
end

"""
    precompute_sigma(Medge, pk, Om; kwargs...)

Precompute σ(M) at z=0 for mass bin edges. Returns vector of σ values.
Scale by D(z) to get σ(M,z) = D(z) × σ(M,z=0).
"""
function precompute_sigma(Medge::AbstractVector, pk, Om::Real; kwargs...)
    sigma_edge = Vector{Float64}(undef, length(Medge))
    Threads.@threads for i in eachindex(Medge)
        sigma_edge[i] = sigma_M(Medge[i], pk, Om; kwargs...)
    end
    return sigma_edge
end

end # module MassFunction
