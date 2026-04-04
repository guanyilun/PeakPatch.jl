module DiffEqExt

using OrdinaryDiffEq
using PeakPatch

import PeakPatch: _evolve_diffeq
import PeakPatch.Cosmology: CosmologyParams, DlinearTables, Dlinear_ab, Dfnofa
import PeakPatch.EllipsoidalCollapse: EllipsoidParams, CosmoCache, Ha_b_nr,
    get_b_2

const ONE_THIRD = 1.0 / 3.0

# Mutable state for tracking virialization events during integration
mutable struct VirState
    lvirv::Vector{Int}       # virialization flags per axis (0 or 1)
    a_eq::Vector{Float64}    # frozen axis values at virialization
    zdynv::Vector{Float64}   # redshift at virialization for each axis (-1 = not yet)
end

VirState() = VirState([0, 0, 0], [0.0, 0.0, 0.0], [-1.0, -1.0, -1.0])

"""
ODE right-hand side for ellipsoidal collapse.
Includes virialization logic inline (like the RK4 version's get_derivs!),
freezing axes when they reach their collapse threshold.

State: u = [a1, a2, a3, da1, da2, da3]
Time:  t = a_b (scale factor)
Parameters: p = (ep, cc, aLam, vir)
"""
function ellipsoid_rhs!(du, u, p, t)
    ep, cc, aLam, vir = p

    a_b = t
    a_b3 = a_b^3
    Ha = Ha_b_nr(a_b, cc)
    Ha_inv = 1.0 / Ha

    fcoll = (ep.fcoll_1, ep.fcoll_2, ep.fcoll_3)

    # Virialization checks (matching get_derivs! logic exactly)
    # Axis 3 first
    if vir.lvirv[3] == 0 && u[3] <= fcoll[3] * a_b
        vir.lvirv[3] = 1
        vir.a_eq[3] = fcoll[3] * a_b
        vir.zdynv[3] = 1.0 / a_b
    end

    # Axis 2 threshold (ivir_strat == 2: uses fcoll_2 * a_b directly)
    a_3eq2 = if ep.ivir_strat == 2 && vir.lvirv[2] == 0
        fcoll[2] * a_b
    elseif ep.ivir_strat == 1 && vir.lvirv[3] == 1
        vir.a_eq[3] * 1.001
    else
        0.0
    end

    # Axis 1 threshold
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

    # Freeze virialized axes (modify u in place, same as Fortran/RK4)
    if vir.lvirv[3] == 1
        u[3] = vir.a_eq[3]; u[6] = 0.0
    end
    if vir.lvirv[2] == 1
        u[2] = vir.a_eq[2]; u[5] = 0.0
    end
    if vir.lvirv[1] == 1
        u[1] = vir.a_eq[1]; u[4] = 0.0
    end

    # Read axes after freezing
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

    # Shape integrals
    bvec_2 = if ep.iforce_strat != 5 && ep.iforce_strat != 6
        get_b_2(a1, a2, a3, ep.iwant_rd)
    else
        nothing
    end

    # Equations of motion for each axis
    aLam_vals = (aLam[1], aLam[2], aLam[3])
    for i in 1:3
        if vir.lvirv[i] == 1
            du[i] = 0.0
            du[i+3] = 0.0
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

function PeakPatch._evolve_diffeq(Frho, e_v, p_v, ep::EllipsoidParams; idynax::Int = 1)
    cc = CosmoCache(ep)

    # Strain eigenvalues
    aLam_3 = Frho / 3.0 * (1.0 + 3.0 * e_v + p_v)
    aLam_2 = Frho / 3.0 * (1.0 - 2.0 * p_v)
    aLam_1 = Frho / 3.0 * (1.0 - 3.0 * e_v + p_v)
    aLam = (aLam_1, aLam_2, aLam_3, Frho / 3.0)

    vir = VirState()

    # Initial conditions (same as RK4 version)
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

    tspan = (a_b, 2.0)
    p = (ep, cc, aLam, vir)

    prob = ODEProblem(ellipsoid_rhs!, u0, tspan, p)

    # Termination callback: stop when axis 1 virializes
    terminate_cb = DiscreteCallback(
        (u, t, integrator) -> integrator.p[4].lvirv[1] == 1,
        integrator -> terminate!(integrator)
    )

    sol = solve(prob, Tsit5();
                callback = terminate_cb,
                abstol = 1e-8, reltol = 1e-6,
                maxiters = 100_000,
                save_everystep = false,
                save_start = false,
                save_end = true)

    # Extract result
    zdynax = vir.zdynv[idynax]
    if zdynax == -1.0
        return (-1.0, -1.0, -1.0)
    else
        Ddynax = Dfnofa(1.0 / zdynax, ep.dlin_tables)
        fcdynax = Frho * Ddynax
        return (zdynax, Ddynax, fcdynax)
    end
end

end # module DiffEqExt
