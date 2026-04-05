module Cosmology

using QuadGK

struct CosmologyParams
    Om::Float64  # total matter Ω_m (CDM + baryons)
    OB::Float64  # baryon density Ω_b
    OL::Float64  # dark energy Ω_Λ
    h::Float64   # H0 / (100 km/s/Mpc)
    ns::Float64  # scalar spectral index
    s8::Float64  # σ_8
end

# E²(z) = H²(z)/H₀² for flat ΛCDM
E2(z, c::CosmologyParams) = c.Om * (1 + z)^3 + c.OL

# Hubble parameter H(z) in km/s/Mpc
H(z, c::CosmologyParams) = 100.0 * c.h * sqrt(E2(z, c))

# Comoving distance χ(z) in Mpc/h
function chi(z, c::CosmologyParams)
    integral, _ = quadgk(zp -> 1.0 / sqrt(E2(zp, c)), 0.0, z)
    return (2.998e5 / 100.0) * c.h * integral  # Mpc/h
end

# Linear growth factor D(z), normalized to D(0) = 1
# Carroll, Press & Turner (1992) approximation
function growth_factor(z, c::CosmologyParams)
    function _g(zp)
        Omz = c.Om * (1 + zp)^3 / E2(zp, c)
        OLz = c.OL / E2(zp, c)
        return 2.5 * Omz / (Omz^(4.0/7.0) - OLz + (1.0 + Omz / 2.0) * (1.0 + OLz / 70.0))
    end
    return _g(z) / (1.0 + z) / _g(0.0)
end

# Linear growth rate f(z) = d ln D / d ln a
function growth_rate(z, c::CosmologyParams)
    Omz = c.Om * (1 + z)^3 / E2(z, c)
    return Omz^0.55  # Luther approximation
end

# Critical overdensity δ_c(z) for spherical collapse.
# Nakamura & Suto (1997) fitting formula for flat ΛCDM:
#   δ_c = (3/20)(12π)^{2/3} × (1 + 0.0123 log₁₀ Ω_m(z))
function delta_c(z, c::CosmologyParams)
    Omz = c.Om * (1 + z)^3 / E2(z, c)
    return (3.0 / 20.0) * (12π)^(2.0 / 3.0) * (1.0 + 0.0123 * log10(Omz))
end

# -------------------------------------------------------------------------
# DlinearTables: exact linear growth factor via RK4 integration
# Port of Dlinear_setup from psubs_Dlinear.f90
# -------------------------------------------------------------------------

struct DlinearTables
    atab::Vector{Float64}          # scale factor grid (ascending)
    D_atab::Vector{Float64}        # D(a)/a at each grid point (normalized)
    dlnD_dlnatab::Vector{Float64}  # d ln D / d ln a at each grid point
    amintab::Float64               # min scale factor in table
    amaxtab::Float64               # max scale factor in table
    D_ainfty::Float64              # D(a)/a at smallest a
    D_a0::Float64                  # D(a)/a at largest a
    HD_Hainfty::Float64            # dlnD/dlna at smallest a
    HD_Ha0::Float64                # dlnD/dlna at largest a
end

# 4-point or 2-point Lagrange interpolation (matches lagrint_one)
function _lagrint_one(xa, val, x)
    n = length(xa)
    if n == 4
        a1 = (x - xa[2]) / (xa[1] - xa[2]) * (x - xa[3]) / (xa[1] - xa[3]) * (x - xa[4]) / (xa[1] - xa[4])
        a2 = (x - xa[1]) / (xa[2] - xa[1]) * (x - xa[3]) / (xa[2] - xa[3]) * (x - xa[4]) / (xa[2] - xa[4])
        a3 = (x - xa[1]) / (xa[3] - xa[1]) * (x - xa[2]) / (xa[3] - xa[2]) * (x - xa[4]) / (xa[3] - xa[4])
        a4 = (x - xa[1]) / (xa[4] - xa[1]) * (x - xa[2]) / (xa[4] - xa[2]) * (x - xa[3]) / (xa[4] - xa[3])
        return a1 * val[1] + a2 * val[2] + a3 * val[3] + a4 * val[4]
    else
        a1 = (x - xa[2]) / (xa[1] - xa[2])
        a2 = 1.0 - a1
        return a1 * val[1] + a2 * val[2]
    end
end

# Lagrange interpolation for 2 output columns (matches lagrint_Dlin for D_a and HD_Ha)
function _lagrint_2col(xa, val1, val2, x)
    n = length(xa)
    if n == 4
        a1 = (x - xa[2]) / (xa[1] - xa[2]) * (x - xa[3]) / (xa[1] - xa[3]) * (x - xa[4]) / (xa[1] - xa[4])
        a2 = (x - xa[1]) / (xa[2] - xa[1]) * (x - xa[3]) / (xa[2] - xa[3]) * (x - xa[4]) / (xa[2] - xa[4])
        a3 = (x - xa[1]) / (xa[3] - xa[1]) * (x - xa[2]) / (xa[3] - xa[2]) * (x - xa[4]) / (xa[3] - xa[4])
        a4 = (x - xa[1]) / (xa[4] - xa[1]) * (x - xa[2]) / (xa[4] - xa[2]) * (x - xa[3]) / (xa[4] - xa[3])
        p1 = a1 * val1[1] + a2 * val1[2] + a3 * val1[3] + a4 * val1[4]
        p2 = a1 * val2[1] + a2 * val2[2] + a3 * val2[3] + a4 * val2[4]
        return (p1, p2)
    else
        a1 = (x - xa[2]) / (xa[1] - xa[2])
        a2 = 1.0 - a1
        p1 = a1 * val1[1] + a2 * val1[2]
        p2 = a1 * val2[1] + a2 * val2[2]
        return (p1, p2)
    end
end

# Binary search (matches HUNT from psubs_Dlinear.f90)
function _hunt(xx, x)
    n = length(xx)
    ascnd = xx[n] > xx[1]
    jlo = 0
    jhi = n + 1
    while jhi - jlo > 1
        jm = (jhi + jlo) ÷ 2
        if (x > xx[jm]) == ascnd
            jlo = jm
        else
            jhi = jm
        end
    end
    return jlo
end

# Derivatives for the growth ODE. Integration variable is ln(a).
# State: y = [time, ln(D), dlnD/dlna, Omnrd*a, tau]
function _deriv_Dlin!(y, d, a, h_val, Omnr, Omvac, Omcurv, Omer, OmX, OmB, Omhdm, fhdmclus, Omnrd, taud, ldecay, Omnrdi)
    Omerd = ldecay == 1 ? Omnrdi * exp(-y[1] / taud) * a : 0.0

    hub = Omer / a + Omnr + Omvac * a^3 + Omcurv * a + Omerd / a + Omnrd
    d[1] = 1.0 / (h_val * sqrt(hub / a^3))

    fb = y[3]
    d[2] = fb

    q = Omer / a + Omnr + Omvac * a^3 + Omcurv * a + Omerd / a + Omnrd
    qs = 0.5 / q
    q_val = qs * (2.0 * Omer / a + Omnr - 2.0 * Omvac * a^3 + 2.0 * Omerd / a + Omnrd)
    fnr = (OmX + OmB + fhdmclus * Omhdm) / Omnr
    if OmB != 0.0 && a <= 1.0e-3
        fnr = (OmX + fhdmclus * Omhdm) / Omnr
    end
    fer = 0.0  # flat ΛCDM
    qf = qs * (2.0 * Omer / a * fer + Omnr * fnr + Omnrd * fnr)

    d[3] = 3.0 * qf - fb * fb - (1.0 - q_val) * fb
    d[4] = Omnrd * a * d[1] / taud
    d[5] = d[1] / a
    return nothing
end

function Dlinear_tables(cosmo::CosmologyParams)
    Omnr = cosmo.Om
    Omvac = cosmo.OL
    Omcurv = 1.0 - cosmo.Om - cosmo.OL
    Omer = 0.0
    OmX = cosmo.Om - cosmo.OB
    Omhdm = 0.0
    fhdmclus = 0.0
    ldecay = 0
    Omnrdi = 0.0
    taud = 1.0
    Omnrd = 0.0

    amaxtab = 1.0 / 0.5  # = 2.0
    amintab_target = 1.0 / 1001.0
    ainit = 1.0e-7
    nmax = 2000

    a = ainit
    alogn = log(ainit)
    amlog = log(amaxtab)
    hh = (amlog - alogn) / nmax

    t0 = 1.0 / cosmo.h / sqrt(Omnr)

    # Initial conditions
    aeq = 0.0
    yaeq = sqrt(a + aeq)
    yaeq0 = sqrt(aeq)
    y = Vector{Float64}(undef, 5)
    y[1] = t0 * (2.0 / 3.0 * (yaeq^3 - yaeq0^3) - 2.0 * aeq * (yaeq - yaeq0))
    y[2] = log(a)
    fbclus = 0.0
    fb = 0.25 * (sqrt(24.0 * (OmX + fbclus * cosmo.OB + fhdmclus * Omhdm) / Omnr + 1.0) - 1.0)
    y[3] = fb
    y[4] = 0.0
    y[5] = t0 * 2.0 * (yaeq - yaeq0)

    atab_all = Vector{Float64}(undef, nmax)
    D_atab_all = Vector{Float64}(undef, nmax)
    dlnD_dlnatab_all = Vector{Float64}(undef, nmax)

    ip = 0
    lnorm = 1
    anorm = 0.001
    bnorm = 1.0
    ip0 = 0

    deriv_args = (cosmo.h, Omnr, Omvac, Omcurv, Omer, OmX, cosmo.OB, Omhdm, fhdmclus, Omnrd, taud, ldecay, Omnrdi)

    for i in 1:nmax
        w1 = similar(y)
        w2 = similar(y)
        w3 = similar(y)
        w4 = similar(y)

        _deriv_Dlin!(y, w1, a, deriv_args...)

        alogn_half = alogn + hh * 0.5
        a_half = exp(alogn_half)

        yt = y .+ hh .* w1 .* 0.5
        _deriv_Dlin!(yt, w2, a_half, deriv_args...)

        yt = y .+ hh .* w2 .* 0.5
        _deriv_Dlin!(yt, w3, a_half, deriv_args...)

        alogn_new = alogn + hh
        a_new = exp(alogn_new)
        yt = y .+ hh .* w3
        _deriv_Dlin!(yt, w4, a_new, deriv_args...)

        y .= y .+ hh ./ 6.0 .* (w1 .+ w4 .+ 2.0 .* (w2 .+ w3))

        alogn = alogn_new
        a = a_new

        if lnorm == 1 && a >= anorm
            lnorm = 0
            b = exp(y[2])
            bnorm = b / a
        end

        if a >= amintab_target && a <= amaxtab
            ip += 1
            b = exp(y[2])
            fb = y[3]
            atab_all[ip] = a
            D_atab_all[ip] = b / a
            dlnD_dlnatab_all[ip] = fb
            if a >= 1.0 && ip0 == 0
                ip0 = ip
            end
        end
    end

    ntab = ip

    # Normalize so D(a=1) = 1
    if ip0 > 0 && ip0 > 1
        ainterp = (atab_all[ip0] - 1.0) / (atab_all[ip0] - atab_all[ip0 - 1])
        D_agrowth = (D_atab_all[ip0] * (1.0 - ainterp) + D_atab_all[ip0 - 1] * ainterp) / bnorm
        D_agrowth_1fac = 1.0 / (D_agrowth * bnorm)
    else
        D_agrowth_1fac = 1.0 / bnorm
    end

    for i in 1:ntab
        D_atab_all[i] *= D_agrowth_1fac
    end

    DlinearTables(
        atab_all[1:ntab],
        D_atab_all[1:ntab],
        dlnD_dlnatab_all[1:ntab],
        atab_all[1],
        atab_all[ntab],
        D_atab_all[1],
        D_atab_all[ntab],
        dlnD_dlnatab_all[1],
        dlnD_dlnatab_all[ntab]
    )
end

# Lookup: given scale factor a_b, return (Dlinear, HD_Ha, D_a)
function Dlinear_ab(ai, tables::DlinearTables)
    ntab = length(tables.atab)
    if ai > tables.amintab && ai < tables.amaxtab
        kk = _hunt(tables.atab, ai)
        if kk < 2 || kk >= ntab
            nlint = 2
            kkl = max(kk - 1, 1)
        else
            nlint = 4
            kkl = kk - 2
        end
        idx = kkl:kkl + nlint - 1
        xa = tables.atab[idx]
        D_a, HD_Ha = _lagrint_2col(xa, tables.D_atab[idx], tables.dlnD_dlnatab[idx], ai)
        return (D_a * ai, HD_Ha, D_a)
    elseif ai >= tables.amaxtab
        return (tables.D_a0 * ai, tables.HD_Ha0, tables.D_a0)
    else
        return (tables.D_ainfty * ai, tables.HD_Hainfty, tables.D_ainfty)
    end
end

# Dfnofa: D(a) normalized so D(a=1)=1
function Dfnofa(ai, tables::DlinearTables)
    ntab = length(tables.atab)
    if ai > tables.amintab && ai < tables.amaxtab
        kk = _hunt(tables.atab, ai)
        if kk < 2 || kk >= ntab
            nlint = 2
            kkl = max(kk - 1, 1)
        else
            nlint = 4
            kkl = kk - 2
        end
        idx = kkl:kkl + nlint - 1
        xa = tables.atab[idx]
        pint = _lagrint_one(xa, tables.D_atab[idx], ai)
        return ai * pint
    elseif ai >= tables.amaxtab
        return ai * tables.D_a0
    else
        return ai * tables.D_ainfty
    end
end

# -------------------------------------------------------------------------
# Chi-to-z inversion: comoving distance → redshift lookup table
# -------------------------------------------------------------------------

struct ChiToZTable
    zs::Vector{Float64}
    chis::Vector{Float64}
    z_max::Float64
end

"""
    build_chi_to_z(cosmo; z_max=10.0, npts=2000)

Build a monotonic lookup table for inverting χ(z) → z via bisection.
"""
function build_chi_to_z(cosmo::CosmologyParams; z_max::Float64=10.0, npts::Int=2000)
    zs = collect(range(0.0, z_max; length=npts))
    chis = [chi(z, cosmo) for z in zs]
    return ChiToZTable(zs, chis, z_max)
end

"""
    chi_to_z(table, chi_val)

Convert comoving distance [Mpc/h] to redshift using a precomputed lookup table.
"""
function chi_to_z(table::ChiToZTable, chi_val::Float64)
    chi_val <= table.chis[1] && return 0.0
    chi_val >= table.chis[end] && return table.z_max
    lo, hi = 1, length(table.zs)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        if table.chis[mid] <= chi_val
            lo = mid
        else
            hi = mid
        end
    end
    frac = (chi_val - table.chis[lo]) / (table.chis[hi] - table.chis[lo])
    return table.zs[lo] + frac * (table.zs[hi] - table.zs[lo])
end

"""
    peak_redshift(obs_x, obs_y, obs_z, pk_x, pk_y, pk_z, chi2z)

Compute redshift of a peak at (pk_x, pk_y, pk_z) as seen by an observer at
(obs_x, obs_y, obs_z), using the precomputed chi-to-z table.
"""
function peak_redshift(obs_x, obs_y, obs_z, pk_x, pk_y, pk_z, chi2z::ChiToZTable)
    dx = Float64(pk_x) - Float64(obs_x)
    dy = Float64(pk_y) - Float64(obs_y)
    dz = Float64(pk_z) - Float64(obs_z)
    r = sqrt(dx^2 + dy^2 + dz^2)
    return chi_to_z(chi2z, r)
end

end # module Cosmology
