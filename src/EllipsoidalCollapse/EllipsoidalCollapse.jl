module EllipsoidalCollapse

import ..Cosmology: CosmologyParams, DlinearTables, Dlinear_tables, Dlinear_ab, Dfnofa
using StaticArrays
using OrdinaryDiffEq: ODEProblem, Tsit5, DiscreteCallback, solve, terminate!

const ONE_THIRD = 1.0 / 3.0

# Configuration for ellipsoidal collapse
struct EllipsoidParams
    cosmo::CosmologyParams
    iforce_strat::Int    # 4 (stbg+Lstrain)
    ivir_strat::Int      # 2 (per-axis fcoll)
    fcoll_1::Float64     # 0.01
    fcoll_2::Float64     # 0.171
    fcoll_3::Float64     # 0.171
    tfac::Float64        # 0.01
    nstepmax::Int        # 10000
    zinit_fac::Float64   # 20.0
    iwant_rd::Int        # 1 (use Carlson RD)
    dlin_tables::DlinearTables
    solver::Symbol       # :diffeq (default, adaptive Tsit5) or :rk4 (Fortran-compatible)
end

function EllipsoidParams(cosmo::CosmologyParams;
    iforce_strat::Int = 4,
    ivir_strat::Int = 2,
    fcoll_1::Float64 = 0.01,
    fcoll_2::Float64 = 0.171,
    fcoll_3::Float64 = 0.171,
    tfac::Float64 = 0.01,
    nstepmax::Int = 10000,
    zinit_fac::Float64 = 20.0,
    iwant_rd::Int = 1,
    solver::Symbol = :diffeq)
    tables = Dlinear_tables(cosmo)
    EllipsoidParams(cosmo, iforce_strat, ivir_strat,
        fcoll_1, fcoll_2, fcoll_3, tfac, nstepmax, zinit_fac, iwant_rd, tables, solver)
end

# Pre-computed cosmology ratios
struct CosmoCache
    Rvac_nr::Float64   # Omvac/Omnr
    Rcurv_nr::Float64  # Omcurv/Omnr
end

function CosmoCache(ep::EllipsoidParams)
    Omnr = ep.cosmo.Om
    Omvac = ep.cosmo.OL
    Omcurv = 1.0 - ep.cosmo.Om - ep.cosmo.OL
    CosmoCache(Omvac / Omnr, Omcurv / Omnr)
end

# Ha = sqrt((1 + Rvac*a^3 + Rcurv*a) / a)
function Ha_b_nr(a_b, cc::CosmoCache)
    sqrt((1.0 + cc.Rvac_nr * a_b^3 + cc.Rcurv_nr * a_b) / a_b)
end

# -------------------------------------------------------------------------
# Carlson's elliptic integral RD(x,y,z)
# Port from Solvers.f90 lines 36-76
# -------------------------------------------------------------------------
function _elliptic_rd(x, y, z)
    ERRTOL = 0.05
    TINY = 1.0e-25
    BIG = 4.5e21
    C1 = 3.0 / 14.0
    C2 = 1.0 / 6.0
    C3 = 9.0 / 22.0
    C4 = 3.0 / 26.0
    C5 = 0.25 * C3
    C6 = 1.5 * C4

    xt = x
    yt = y
    zt = z
    sum = 0.0
    fac = 1.0
    delx = 0.0
    dely = 0.0
    delz = 0.0
    ave = 0.0

    while true
        sqrtx = sqrt(xt)
        sqrty = sqrt(yt)
        sqrtz = sqrt(zt)
        alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz
        sum += fac / (sqrtz * (zt + alamb))
        fac *= 0.25
        xt = 0.25 * (xt + alamb)
        yt = 0.25 * (yt + alamb)
        zt = 0.25 * (zt + alamb)
        ave = 0.2 * (xt + yt + 3.0 * zt)
        delx = (ave - xt) / ave
        dely = (ave - yt) / ave
        delz = (ave - zt) / ave
        max(abs(delx), abs(dely), abs(delz)) <= ERRTOL && break
    end

    ea = delx * dely
    eb = delz * delz
    ec = ea - eb
    ed = ea - 6.0 * eb
    ee = ed + ec + ec

    return 3.0 * sum + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * delz * ee) +
        delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea))) / (ave * sqrt(ave))
end

# -------------------------------------------------------------------------
# Shape integrals b_i
# Port from HomogeneousEllipsoid.f90 lines 373-430
# -------------------------------------------------------------------------
function get_b_2(a1, a2, a3, iwant_rd::Int)
    r_2 = a2 / a1
    r_3 = a3 / a1

    if r_3 >= r_2
        # Oblate or spherical
        if r_3 >= 0.999
            return SVector(ONE_THIRD, ONE_THIRD, ONE_THIRD)
        elseif r_3 <= 0.001
            ecc2 = 1.0 - r_3^2
            one_ecc2 = r_3^2
            one_ecc = 0.5 * r_3^2
            ecc = sqrt(ecc2)
            b3 = 0.5 * (one_ecc2 / ecc2 * (1.0 / one_ecc2 - log((1.0 + ecc) / one_ecc) * 0.5 / ecc))
            return SVector(1.0 - 2.0 * b3, b3, b3)
        else
            ecc2 = 1.0 - r_3^2
            ecc = sqrt(ecc2)
            b3 = 0.5 * ((1.0 - ecc2) / ecc2 * (1.0 / (1.0 - ecc2) - log((1.0 + ecc) / (1.0 - ecc)) * 0.5 / ecc))
            return SVector(1.0 - 2.0 * b3, b3, b3)
        end
    end

    # Triaxial case: r_3 < r_2
    if iwant_rd == 1
        r_22 = r_2 * r_2
        r_32 = r_3 * r_3
        qrat = ONE_THIRD * r_2 * r_3
        b3 = qrat * _elliptic_rd(r_22, 1.0, r_32)
        b2 = qrat * _elliptic_rd(r_32, 1.0, r_22)
        return SVector(1.0 - (b2 + b3), b2, b3)
    else
        error("Legendre form (iwant_rd != 1) not implemented in Julia port")
    end
end

# -------------------------------------------------------------------------
# RK4 integrator (matches Solvers.f90 RK4)
# -------------------------------------------------------------------------
function rk4_homel!(y, dy, t, h, yout, derivs_fn, derivs_state, yt, dyt, dym)
    hh = h * 0.5
    h6 = h / 6.0
    xh = t + hh
    n = length(y)

    @inbounds for i in 1:n
        yt[i] = y[i] + hh * dy[i]
    end
    derivs_fn(xh, yt, dyt, derivs_state)

    @inbounds for i in 1:n
        yt[i] = y[i] + hh * dyt[i]
    end
    derivs_fn(xh, yt, dym, derivs_state)

    @inbounds for i in 1:n
        yt[i] = y[i] + h * dym[i]
        dym[i] = dyt[i] + dym[i]
    end
    derivs_fn(t + h, yt, dyt, derivs_state)

    @inbounds for i in 1:n
        yout[i] = y[i] + h6 * (dy[i] + dyt[i] + 2.0 * dym[i])
    end
end

# -------------------------------------------------------------------------
# Mutable state for the derivative function
# -------------------------------------------------------------------------
mutable struct DerivState
    lvirv::Vector{Int}      # virialization flags (3 elements)
    a_3eq::Float64          # frozen axis 3 value
    a_2eq::Float64          # frozen axis 2 value
    a_1eq::Float64          # frozen axis 1 value
    a_3eq2::Float64         # threshold for axis 2 vir
    a_3eq1::Float64         # threshold for axis 1 vir
    Ha_b_nr_val::Float64    # last computed Ha_b_nr
    # Strain eigenvalues (constant for a given run)
    aLam_1::Float64
    aLam_2::Float64
    aLam_3::Float64
    Frho_3::Float64         # Frho/3
end

function DerivState(aLam_1, aLam_2, aLam_3, Frho)
    DerivState(
        [0, 0, 0],           # lvirv
        0.0, 0.0, 0.0,       # a_3eq, a_2eq, a_1eq
        0.0, 0.0,            # a_3eq2, a_3eq1
        0.0,                 # Ha_b_nr_val
        aLam_1, aLam_2, aLam_3,
        Frho / 3.0
    )
end

# -------------------------------------------------------------------------
# Derivatives get_derivs!
# Port from HomogeneousEllipsoid.f90 lines 231-370
# -------------------------------------------------------------------------
function get_derivs!(t, y, dy, state::Tuple{EllipsoidParams, CosmoCache, DerivState})
    ep, cc, ds = state
    dy .= 0.0

    a_b = t
    a_b3 = a_b^3
    Ha = Ha_b_nr(a_b, cc)
    ds.Ha_b_nr_val = Ha
    Ha_inv = 1.0 / Ha

    avec = @view y[1:3]

    # Virialization checks (axis 3 first)
    if ds.lvirv[3] == 0 && avec[3] <= ep.fcoll_3 * a_b
        ds.lvirv[3] = 1
        ds.a_3eq = ep.fcoll_3 * a_b
        if ep.ivir_strat == 1
            ds.a_3eq2 = ds.a_3eq * 1.001
            ds.a_3eq1 = ds.a_3eq * 1.0003
        end
    end
    if ep.ivir_strat == 2 && ds.lvirv[2] == 0
        ds.a_3eq2 = ep.fcoll_2 * a_b
    end
    if ep.ivir_strat == 2 && ds.lvirv[1] == 0
        ds.a_3eq1 = ep.fcoll_1 * a_b
    end

    if ds.lvirv[3] == 1 && avec[2] <= ds.a_3eq2 && ds.lvirv[2] == 0
        ds.lvirv[2] = 1
        ds.a_2eq = ds.a_3eq2
    end
    if ds.lvirv[3] == 1 && avec[1] <= ds.a_3eq1 && ds.lvirv[1] == 0
        ds.lvirv[1] = 1
        ds.a_1eq = ds.a_3eq1
    end

    # Freeze virialized axes
    if ds.lvirv[3] == 1
        y[3] = ds.a_3eq
        y[6] = 0.0
        dy[3] = 0.0
        dy[6] = 0.0
    end
    if ds.lvirv[2] == 1
        y[2] = ds.a_2eq
        y[5] = 0.0
        dy[2] = 0.0
        dy[5] = 0.0
    end
    if ds.lvirv[1] == 1
        y[1] = ds.a_1eq
        y[4] = 0.0
        dy[1] = 0.0
        dy[4] = 0.0
    end

    # Recompute avec after possible freezing
    a1 = y[1]; a2 = y[2]; a3 = y[3]
    delta1_e = a_b3 / (a1 * a2 * a3)

    # Force computation
    if ep.iforce_strat == 0
        d0 = -cc.Rvac_nr * Ha_inv
        d1_int = 1.5 * delta1_e / a_b3 * Ha_inv
    else
        d0 = (0.5 / a_b3 - cc.Rvac_nr) * Ha_inv
        d1_int = 1.5 * (delta1_e - 1.0) / a_b3 * Ha_inv
        if ep.iforce_strat == 4 || ep.iforce_strat == 5 || ep.iforce_strat == 6
            dlin, HD_Ha, _ = Dlinear_ab(a_b, ep.dlin_tables)
            d1_ext = 1.5 / a_b3 * Ha_inv * dlin
        elseif ep.iforce_strat == 3
            d1_ext = 1.5 * 2.5 / a_b3 * Ha_inv
        end
    end

    # Shape integrals (for iforce_strat != 5,6)
    if ep.iforce_strat != 5 && ep.iforce_strat != 6
        bvec_2 = get_b_2(a1, a2, a3, ep.iwant_rd)
    end

    # Equations of motion for each axis
    # Axis 3
    if ds.lvirv[3] == 0
        dy[3] = y[6] * Ha_inv
        if ep.iforce_strat == 0
            dy[6] = -y[3] * d1_int * bvec_2[3]
        elseif ep.iforce_strat == 1
            dy[6] = -y[3] * (d1_int * bvec_2[3] + d0)
        elseif ep.iforce_strat == 3
            dy[6] = -y[3] * (d1_int * bvec_2[3] + d0 + d1_ext * (bvec_2[3] - ONE_THIRD))
        elseif ep.iforce_strat == 4 || ep.iforce_strat == 6
            dy[6] = -y[3] * (d1_int * bvec_2[3] + d0 + d1_ext * (ds.aLam_3 - ds.Frho_3))
        elseif ep.iforce_strat == 5
            dy[6] = -y[3] * d0 - d1_ext * a_b * ds.aLam_3
        end
    end

    # Axis 2
    if ds.lvirv[2] == 0
        dy[2] = y[5] * Ha_inv
        if ep.iforce_strat == 0
            dy[5] = -y[2] * d1_int * bvec_2[2]
        elseif ep.iforce_strat == 1
            dy[5] = -y[2] * (d1_int * bvec_2[2] + d0)
        elseif ep.iforce_strat == 3
            dy[5] = -y[2] * (d1_int * bvec_2[2] + d0 + d1_ext * (bvec_2[2] - ONE_THIRD))
        elseif ep.iforce_strat == 4 || ep.iforce_strat == 6
            dy[5] = -y[2] * (d1_int * bvec_2[2] + d0 + d1_ext * (ds.aLam_2 - ds.Frho_3))
        elseif ep.iforce_strat == 5
            dy[5] = -y[2] * d0 - d1_ext * a_b * ds.aLam_2
        end
    end

    # Axis 1
    if ds.lvirv[1] == 0
        dy[1] = y[4] * Ha_inv
        if ep.iforce_strat == 0
            dy[4] = -y[1] * d1_int * bvec_2[1]
        elseif ep.iforce_strat == 1
            dy[4] = -y[1] * (d1_int * bvec_2[1] + d0)
        elseif ep.iforce_strat == 3
            dy[4] = -y[1] * (d1_int * bvec_2[1] + d0 + d1_ext * (bvec_2[1] - ONE_THIRD))
        elseif ep.iforce_strat == 4 || ep.iforce_strat == 6
            dy[4] = -y[1] * (d1_int * bvec_2[1] + d0 + d1_ext * (ds.aLam_1 - ds.Frho_3))
        elseif ep.iforce_strat == 5
            dy[4] = -y[1] * d0 - d1_ext * a_b * ds.aLam_1
        end
    end

    return nothing
end

# -------------------------------------------------------------------------
# Adaptive ODE solver via OrdinaryDiffEq (Tsit5)
# -------------------------------------------------------------------------

# Mutable state for tracking virialization events during DiffEq integration
mutable struct VirState
    lvirv::Vector{Int}       # virialization flags per axis (0 or 1)
    a_eq::Vector{Float64}    # frozen axis values at virialization
    zdynv::Vector{Float64}   # redshift at virialization for each axis (-1 = not yet)
end

VirState() = VirState([0, 0, 0], [0.0, 0.0, 0.0], [-1.0, -1.0, -1.0])

"""
ODE right-hand side for ellipsoidal collapse (DiffEq interface).
Includes virialization logic inline (same physics as the RK4's get_derivs!),
freezing axes when they reach their collapse threshold.

State: u = [a1, a2, a3, da1, da2, da3]
Time:  t = a_b (scale factor)
Parameters: p = (ep, cc, aLam, vir)
"""
function _diffeq_rhs!(du, u, p, t)
    ep, cc, aLam, vir = p

    a_b = t
    a_b3 = a_b^3
    Ha = Ha_b_nr(a_b, cc)
    Ha_inv = 1.0 / Ha

    fcoll = (ep.fcoll_1, ep.fcoll_2, ep.fcoll_3)

    # Virialization checks (matching get_derivs! logic exactly)
    if vir.lvirv[3] == 0 && u[3] <= fcoll[3] * a_b
        vir.lvirv[3] = 1
        vir.a_eq[3] = fcoll[3] * a_b
        vir.zdynv[3] = 1.0 / a_b
    end

    a_3eq2 = if ep.ivir_strat == 2 && vir.lvirv[2] == 0
        fcoll[2] * a_b
    elseif ep.ivir_strat == 1 && vir.lvirv[3] == 1
        vir.a_eq[3] * 1.001
    else
        0.0
    end

    a_3eq1 = if ep.ivir_strat == 2 && vir.lvirv[1] == 0
        fcoll[1] * a_b
    elseif ep.ivir_strat == 1 && vir.lvirv[3] == 1
        vir.a_eq[3] * 1.0003
    else
        0.0
    end

    if vir.lvirv[3] == 1 && vir.lvirv[2] == 0 && u[2] <= a_3eq2
        vir.lvirv[2] = 1
        vir.a_eq[2] = a_3eq2
        vir.zdynv[2] = 1.0 / a_b
    end

    if vir.lvirv[3] == 1 && vir.lvirv[1] == 0 && u[1] <= a_3eq1
        vir.lvirv[1] = 1
        vir.a_eq[1] = a_3eq1
        vir.zdynv[1] = 1.0 / a_b
    end

    # Freeze virialized axes
    if vir.lvirv[3] == 1; u[3] = vir.a_eq[3]; u[6] = 0.0; end
    if vir.lvirv[2] == 1; u[2] = vir.a_eq[2]; u[5] = 0.0; end
    if vir.lvirv[1] == 1; u[1] = vir.a_eq[1]; u[4] = 0.0; end

    a1 = u[1]; a2 = u[2]; a3 = u[3]
    delta1_e = a_b3 / (a1 * a2 * a3)

    # Force computation
    Frho_3 = aLam[4]
    if ep.iforce_strat == 0
        d0 = -cc.Rvac_nr * Ha_inv
        d1_int = 1.5 * delta1_e / a_b3 * Ha_inv
        d1_ext = 0.0
    else
        d0 = (0.5 / a_b3 - cc.Rvac_nr) * Ha_inv
        d1_int = 1.5 * (delta1_e - 1.0) / a_b3 * Ha_inv
        if ep.iforce_strat == 4 || ep.iforce_strat == 5 || ep.iforce_strat == 6
            dlin, _, _ = Dlinear_ab(a_b, ep.dlin_tables)
            d1_ext = 1.5 / a_b3 * Ha_inv * dlin
        elseif ep.iforce_strat == 3
            d1_ext = 1.5 * 2.5 / a_b3 * Ha_inv
        else
            d1_ext = 0.0
        end
    end

    bvec_2 = if ep.iforce_strat != 5 && ep.iforce_strat != 6
        get_b_2(a1, a2, a3, ep.iwant_rd)
    else
        nothing
    end

    aLam_vals = (aLam[1], aLam[2], aLam[3])
    for i in 1:3
        if vir.lvirv[i] == 1
            du[i] = 0.0; du[i+3] = 0.0
        else
            du[i] = u[i+3] * Ha_inv
            if ep.iforce_strat == 0
                du[i+3] = -u[i] * d1_int * bvec_2[i]
            elseif ep.iforce_strat == 1
                du[i+3] = -u[i] * (d1_int * bvec_2[i] + d0)
            elseif ep.iforce_strat == 3
                du[i+3] = -u[i] * (d1_int * bvec_2[i] + d0 + d1_ext * (bvec_2[i] - ONE_THIRD))
            elseif ep.iforce_strat == 4 || ep.iforce_strat == 6
                du[i+3] = -u[i] * (d1_int * bvec_2[i] + d0 + d1_ext * (aLam_vals[i] - Frho_3))
            elseif ep.iforce_strat == 5
                du[i+3] = -u[i] * d0 - d1_ext * a_b * aLam_vals[i]
            end
        end
    end
    return nothing
end

function _evolve_diffeq(Frho, e_v, p_v, ep::EllipsoidParams; idynax::Int = 1)
    cc = CosmoCache(ep)

    aLam_3 = Frho / 3.0 * (1.0 + 3.0 * e_v + p_v)
    aLam_2 = Frho / 3.0 * (1.0 - 2.0 * p_v)
    aLam_1 = Frho / 3.0 * (1.0 - 3.0 * e_v + p_v)
    aLam = (aLam_1, aLam_2, aLam_3, Frho / 3.0)

    vir = VirState()

    zinit = ep.zinit_fac * Frho
    a_b = 1.0 / zinit
    dlin, HD_Ha, _ = Dlinear_ab(a_b, ep.dlin_tables)
    Ha = Ha_b_nr(a_b, cc)

    u0 = Vector{Float64}(undef, 6)
    u0[3] = a_b * (1.0 - dlin * aLam_3)
    u0[2] = a_b * (1.0 - dlin * aLam_2)
    u0[1] = a_b * (1.0 - dlin * aLam_1)
    u0[6] = Ha * (1.0 - (1.0 + HD_Ha) * dlin * aLam_3)
    u0[5] = Ha * (1.0 - (1.0 + HD_Ha) * dlin * aLam_2)
    u0[4] = Ha * (1.0 - (1.0 + HD_Ha) * dlin * aLam_1)

    prob = ODEProblem(_diffeq_rhs!, u0, (a_b, 2.0), (ep, cc, aLam, vir))

    terminate_cb = DiscreteCallback(
        (u, t, integrator) -> integrator.p[4].lvirv[1] == 1,
        integrator -> terminate!(integrator)
    )

    sol = solve(prob, Tsit5();
                callback = terminate_cb,
                abstol = 1e-8, reltol = 1e-6,
                maxiters = 100_000,
                save_everystep = false, save_start = false, save_end = true)

    zdynax = vir.zdynv[idynax]
    if zdynax == -1.0
        return (-1.0, -1.0, -1.0)
    else
        Ddynax = Dfnofa(1.0 / zdynax, ep.dlin_tables)
        return (zdynax, Ddynax, Frho * Ddynax)
    end
end

# -------------------------------------------------------------------------
# RK4 solver (Fortran-compatible, for validation)
# Port from HomogeneousEllipsoid.f90 lines 10-228
# -------------------------------------------------------------------------
function evolve_ellipse_full(Frho, e_v, p_v, ep::EllipsoidParams; idynax::Int = 1)
    if ep.solver == :diffeq
        return _evolve_diffeq(Frho, e_v, p_v, ep; idynax=idynax)
    end

    cc = CosmoCache(ep)

    # Strain eigenvalues
    aLam_3 = Frho / 3.0 * (1.0 + 3.0 * e_v + p_v)
    aLam_2 = Frho / 3.0 * (1.0 - 2.0 * p_v)
    aLam_1 = Frho / 3.0 * (1.0 - 3.0 * e_v + p_v)

    # Derivative state
    ds = DerivState(aLam_1, aLam_2, aLam_3, Frho)
    derivs_state = (ep, cc, ds)

    # Initial conditions
    zinit = ep.zinit_fac * Frho
    a_b = 1.0 / zinit
    t = a_b

    dlin, HD_Ha, D_a = Dlinear_ab(a_b, ep.dlin_tables)

    Ha = Ha_b_nr(a_b, cc)
    ds.Ha_b_nr_val = Ha

    y = Vector{Float64}(undef, 6)
    y[3] = a_b * (1.0 - dlin * aLam_3)
    y[2] = a_b * (1.0 - dlin * aLam_2)
    y[1] = a_b * (1.0 - dlin * aLam_1)
    y[6] = Ha * (1.0 - (1.0 + HD_Ha) * dlin * aLam_3)
    y[5] = Ha * (1.0 - (1.0 + HD_Ha) * dlin * aLam_2)
    y[4] = Ha * (1.0 - (1.0 + HD_Ha) * dlin * aLam_1)

    # Event tracking arrays (7 events: 3 vir, 3 turn, 1 density)
    nzdyn = 7
    zdynv = fill(-1.0, nzdyn)
    ldynv = fill(0, nzdyn)
    delta1_ev = fill(-1.0, nzdyn)
    aturnv = fill(0.0, 3)

    dy = Vector{Float64}(undef, 6)
    yout = Vector{Float64}(undef, 6)

    # Pre-allocated RK4 workspace (avoids allocation per step)
    yt = Vector{Float64}(undef, 6)
    dyt = Vector{Float64}(undef, 6)
    dym = Vector{Float64}(undef, 6)

    istep = 0

    while true
        istep += 1
        a_3 = y[3]
        dtstep = ep.tfac * a_3 * sqrt(a_3) * ds.Ha_b_nr_val

        # NOTE: "funny RK4" requires first call to derivs
        get_derivs!(t, y, dy, derivs_state)
        rk4_homel!(y, dy, t, dtstep, yout, get_derivs!, derivs_state, yt, dyt, dym)

        # Copy yout to y (in-place for the frozen axes that may have been set in derivs)
        copyto!(y, yout)
        t += dtstep

        delta1_e = t^3 / (y[1] * y[2] * y[3])

        # Track virialization events (indices 1-3)
        for ii in 1:3
            if ds.lvirv[ii] != ldynv[ii]
                zdynv[ii] = 1.0 / t
                ldynv[ii] = 1
                delta1_ev[ii] = delta1_e
            end
        end

        # Track turnaround events (indices 4-6)
        for ii in 4:6
            if ldynv[ii] != 1
                iax = ii - 3
                if y[iax] < aturnv[iax]
                    zdynv[ii] = 1.0 / t
                    ldynv[ii] = 1
                    delta1_ev[ii] = delta1_e
                else
                    aturnv[iax] = y[iax]
                end
            end
        end

        # Track density event (index 7)
        if ldynv[7] != 1
            # dcrit not used in standard case, skip
        end

        # Termination conditions
        if ds.lvirv[1] == 1 || istep >= ep.nstepmax || t >= 2.0
            zdynax = zdynv[idynax]
            if zdynax == -1.0
                return (-1.0, -1.0, -1.0)
            else
                Ddynax = Dfnofa(1.0 / zdynax, ep.dlin_tables)
                fcdynax = Frho * Ddynax
                return (zdynax, Ddynax, fcdynax)
            end
        end
    end
end

end # module EllipsoidalCollapse
