module RadialShell

using LinearAlgebra
using StaticArrays
import ..CollapseTable: CollapseTableInterp, interpolate
import ..Cosmology: Dlinear_ab

# ---------- constants ----------
const PI = Float64(pi)
const FOURPI = 4.0 * PI
const ONE_THIRD = 1.0 / 3.0
const TWO_THIRDS = 2.0 / 3.0

# ---------- dump-reason diagnostics (thread-safe) ----------
const _dump_premask    = Threads.Atomic{Int}(0)
const _dump_no_fcrit   = Threads.Atomic{Int}(0)
const _dump_m0_one     = Threads.Atomic{Int}(0)
const _dump_no_collapse = Threads.Atomic{Int}(0)
const _dump_rthl_neg   = Threads.Atomic{Int}(0)

function reset_dump_counters!()
    Threads.atomic_xchg!(_dump_premask, 0)
    Threads.atomic_xchg!(_dump_no_fcrit, 0)
    Threads.atomic_xchg!(_dump_m0_one, 0)
    Threads.atomic_xchg!(_dump_no_collapse, 0)
    Threads.atomic_xchg!(_dump_rthl_neg, 0)
    return nothing
end

function get_dump_counts()
    return (premask=_dump_premask[], no_fcrit=_dump_no_fcrit[],
            m0_one=_dump_m0_one[], no_collapse=_dump_no_collapse[],
            rthl_neg=_dump_rthl_neg[])
end

# ---------- data structures ----------

struct ShellCell
    di::Int32
    dj::Int32
    dk::Int32
    r2::Int32
end

struct PeakGrid{T<:AbstractFloat}
    delta::Array{T,3}
    etax::Array{T,3}
    etay::Array{T,3}
    etaz::Array{T,3}
    eta2x::Union{Nothing,Array{T,3}}
    eta2y::Union{Nothing,Array{T,3}}
    eta2z::Union{Nothing,Array{T,3}}
    mask::Union{Nothing,Array{Int8,3}}
    nn::NTuple{3,Int}
    lapd::Union{Nothing,Array{T,3}}  # Laplacian of unsmoothed delta (for ioutshear≥1)
end

struct PeakResult
    RTHL::Float64
    Srb::Float64
    Sbar::SVector{3,Float64}
    Sbar2::SVector{3,Float64}
    strain_mat::SMatrix{3,3,Float64,9}
    gradpk::SVector{3,Float64}
    gradpkrf::SVector{3,Float64}
    Fbarx::Float64
    e_v::Float64
    p_v::Float64
    zvir_half::Float64      # formation redshift (collapse redshift at RTHL/2)
    gradpkf::SVector{3,Float64}  # gradient of filtered overdensity at RTHL
    d2F::Float64               # Laplacian of overdensity at RTHL
end

# Convenience constructor accepting AbstractVector/AbstractMatrix (backward compat)
function PeakResult(RTHL, Srb, Sbar::AbstractVector, Sbar2::AbstractVector,
                    strain_mat::AbstractMatrix, gradpk::AbstractVector,
                    gradpkrf::AbstractVector, Fbarx, e_v, p_v, zvir_half,
                    gradpkf::AbstractVector, d2F)
    PeakResult(Float64(RTHL), Float64(Srb),
               SVector{3,Float64}(Sbar), SVector{3,Float64}(Sbar2),
               SMatrix{3,3,Float64,9}(strain_mat),
               SVector{3,Float64}(gradpk), SVector{3,Float64}(gradpkrf),
               Float64(Fbarx), Float64(e_v), Float64(p_v), Float64(zvir_half),
               SVector{3,Float64}(gradpkf), Float64(d2F))
end

const _zero_svec3 = @SVector zeros(3)
const _zero_smat3 = @SMatrix zeros(3,3)

function no_collapse()
    PeakResult(-1.0, 0.0, _zero_svec3, _zero_svec3, _zero_smat3,
               _zero_svec3, _zero_svec3, 0.0, 0.0, 0.0, -1.0,
               _zero_svec3, 0.0)
end

# ---------- SPH kernel integral (hRinteg) ----------
# Integral of SPH kernel 4pi*W(x)*x*(x^2 - u12) dx from |u0| to 2
# Matches Fortran hRinteg in peakvoidsubs.f90 lines 1007-1040

function hRinteg(u0::Float64, u12::Float64)::Float64
    x0 = abs(u0)
    x0 >= 2.0 && return 0.0

    one7 = 1.0 / 7.0
    result = 1.6 * u12 - 6.4 * one7  # value at upper boundary x=2

    # Note: do NOT return early for x0==0. The x0<1 branch with the
    # -0.2*(u12-1/7) correction gives 1.4*u12-31/35, matching Fortran.
    # Returning 1.6*u12-6.4/7 here would be wrong.

    x02 = x0 * x0
    x04 = x02 * x02

    if x0 >= 1.0
        h1 = x04 * (2.0 - x0 * (2.4 - x0 * (1.0 - x0 * one7)))
        h2 = x02 * (4.0 - x0 * (4.0 - x0 * (1.5 - x0 * 0.2)))
        result = result - h2 * u12 + h1
    else
        result = result - 0.2 * (u12 - one7)
        if x0 > 0.0
            h1 = x04 * (1.0 - x02 * (1.0 - 3.0 * x0 * one7))
            h2 = x02 * (2.0 - x02 * (1.5 - 0.6 * x0))
            result = result - h2 * u12 + h1
        end
    end
    return result
end

# ---------- akk kernel table (atab4) ----------
# Precompute cubic spline (Schoenberg M4) SPH kernel
# akk[i] for i=1..1100, u = (i-1)*0.01, ru = sqrt(u)

function _build_atab4()
    akk = Vector{Float64}(undef, 1101)
    ca = 3.0 / (2.0 * PI)
    cb = 1.0 / (4.0 * PI)
    for i in 1:1100
        u = (i - 1) * 0.01
        ru = sqrt(u)
        if ru < 1.0
            akk[i] = ca * (TWO_THIRDS - u + 0.5 * u * ru)
        elseif ru <= 2.0
            akk[i] = cb * (2.0 - ru)^3
        else
            akk[i] = 0.0
        end
    end
    akk[1101] = 0.0
    return akk
end

const _AKK_TAB = _build_atab4()

"""Return the cached SPH kernel table (1101-element Vector{Float64})."""
atab4() = _AKK_TAB

# ---------- precompute_shells (icloud) ----------
# Generate all lattice cells within radius rmax, sorted by r^2

function precompute_shells(rmax::Int)
    cells = ShellCell[]
    r2max = rmax * rmax
    for dk in -rmax:rmax
        for dj in -rmax:rmax
            for di in -rmax:rmax
                r2 = di*di + dj*dj + dk*dk
                if r2 <= r2max
                    push!(cells, ShellCell(Int32(di), Int32(dj), Int32(dk), Int32(r2)))
                end
            end
        end
    end
    sort!(cells; by = c -> c.r2)
    return cells
end

# ---------- get_ijk ----------
# Convert 1-d index ipp to 3-d indices (i,j,k) in grid nn

function get_ijk(n1xn2::Int, n1::Int, ipp::Int)
    k = (ipp - 1) ÷ n1xn2 + 1
    j = ((ipp - 1) % n1xn2) ÷ n1 + 1
    i = (ipp - 1) % n1 + 1
    return (i, j, k)
end

# ---------- fsc_of_z ----------
# Critical overdensity at redshift z (linear growth factor normalization)
# Uses delta_c = 1.686 as the spherical collapse threshold

function fsc_of_z(z::Float64, tables)
    # D(1/(1+z)) / D(1) = growth factor at redshift z, normalized to 1 at z=0
    a = 1.0 / (1.0 + z)
    dlin, = Dlinear_ab(a, tables)
    return Float64(1.686 * dlin)
end

# ---------- normalize_strain ----------
# Rescale strain so trace = F

function normalize_strain!(strain::Matrix{Float64}, F::Float64)
    tr = strain[1,1] + strain[2,2] + strain[3,3]
    if tr != 0.0
        strain .*= F / tr
    end
    return strain
end

function normalize_strain(strain::SMatrix{3,3,Float64,9}, F::Float64)
    tr = strain[1,1] + strain[2,2] + strain[3,3]
    tr != 0.0 ? strain * (F / tr) : strain
end

# ---------- eigenvalue helper ----------
# Returns (Lam, iflag) where Lam sorted ascending

function get_evals(mat::AbstractMatrix{Float64})
    F = Symmetric(mat)
    lam = eigvals(F)       # sorted ascending
    vecs = eigvecs(F)      # columns are eigenvectors
    return (lam, vecs, 0)  # iflag always 0 for Julia eigendecomposition
end

# ---------- kernel-weighted strain tensor ----------
# Sum over shells mlow..mupp with SPH kernel centered at rad(mp)
# Returns (Ebar, grad, gradf) — all 3-vectors / 3x3 matrices

function kernel_strain(rad, Gshell, Gshellf, SRshell, mp, mlow, mupp, akk_tab,
                       wRnor::Float64, aRnor::Float64, hlatt_1::Float64, hlatt_2::Float64)
    Ebar = @MMatrix zeros(3, 3)
    grad = @MVector zeros(3)
    gradf = @MVector zeros(3)

    if mp > 1
        for m1 in mlow:mupp
            u0 = hlatt_1 * (rad[mp] - rad[m1])
            u02 = u0 * u0
            if u02 < 4.0
                u12 = hlatt_2 * (rad[mp]^2 + rad[m1]^2)
                wt = hRinteg(u0, u12) / (u12 - u02)^2
                for L in 1:3
                    grad[L] += wt * Gshell[L][m1]
                    gradf[L] += wt * Gshellf[L][m1]
                    for K in 1:3
                        Ebar[L,K] -= 0.5 * wt * (SRshell[L][K][m1] + SRshell[K][L][m1])
                    end
                end
            end
        end
        scale = wRnor / rad[mp]
        Ebar .*= scale
        grad .*= scale
        gradf .*= scale
    else
        # mp == 1: use akk kernel table, normalize by aRnor
        for m1 in mlow:mupp
            u0 = hlatt_1 * rad[m1]
            con = 100.0 * u0 * u0
            icon = floor(Int, con)
            diff = con - icon
            i2c = icon + 1
            aww = akk_tab[i2c] + diff * (akk_tab[i2c+1] - akk_tab[i2c])
            wt = 0.5 * aww / rad[m1]
            for L in 1:3
                grad[L] += wt * Gshell[L][m1]
                gradf[L] += wt * Gshellf[L][m1]
                for K in 1:3
                    Ebar[L,K] -= 0.5 * wt * (SRshell[L][K][m1] + SRshell[K][L][m1])
                end
            end
        end
        Ebar .*= aRnor
        grad .*= aRnor
        gradf .*= aRnor
    end
    return (SMatrix{3,3}(Ebar), SVector{3}(grad), SVector{3}(gradf))
end

# ---------- analyse_peak ----------
# Top-level orchestrator. Returns PeakResult.

function analyse_peak(pg::PeakGrid, ipp::Int, alatt::Float64, ir2min::Int,
                      ZZon::Float64, Rfclvi::Float64, ct::CollapseTableInterp,
                      shells::Vector{ShellCell};
                      nbuff::Int = 0, growth_tables = nothing,
                      rmax2rs::Float64 = 0.0, fcrit_override = nothing,
                      fortran_compat::Bool = false)
    # Limit npart per peak when rmax2rs > 0 (matches Fortran hpkvd.f90:650-656)
    npartmax = length(shells)
    if rmax2rs > 0.0
        f3p = fortran_compat ? Float32(4.0f0/3.0f0 * Float32(pi)) : 4.0/3.0 * pi
        npart = min(floor(Int, f3p * (rmax2rs * Rfclvi / alatt)^3), npartmax)
    else
        npart = npartmax
    end

    # normalization constants
    hlatt = 1.0
    anor  = 1.0 / hlatt^3
    aRnor = 3.0 * anor / alatt
    wnor  = anor / FOURPI
    wRnor = 3.0 * wnor / alatt

    n1, n2, n3 = pg.nn
    n1xn2 = n1 * n2

    # peak position on lattice
    ipk = get_ijk(n1xn2, n1, ipp)

    # rmax2rs pre-mask check: skip if already masked by larger halo
    if rmax2rs > 0.0 && pg.mask !== nothing
        if pg.mask[ipk[1], ipk[2], ipk[3]] > floor(Int, Rfclvi / alatt * rmax2rs)
            Threads.atomic_add!(_dump_premask, 1)
            return no_collapse()
        end
    end

    # fcrit
    fcrit = fcrit_override !== nothing ? Float64(fcrit_override) :
            (growth_tables !== nothing ? fsc_of_z(ZZon - 1.0, growth_tables) : 1.686)

    # allocate profile arrays
    maxm = npart + 2
    rad    = Vector{Float64}(undef, maxm)
    Fbar   = Vector{Float64}(undef, maxm)
    nshell = Vector{Int}(undef, maxm)
    Sshell  = [Vector{Float64}(undef, maxm) for _ in 1:3]
    S2shell = [Vector{Float64}(undef, maxm) for _ in 1:3]
    Gshell  = [Vector{Float64}(undef, maxm) for _ in 1:3]
    Gshellf = [Vector{Float64}(undef, maxm) for _ in 1:3]
    SRshell = [[Vector{Float64}(undef, maxm) for _ in 1:3] for _ in 1:3]

    # ---- Phase 1+2: Build radial profiles ----
    m = 1
    m0 = 1
    rad[1] = 0.0
    Fbar[1] = pg.delta[ipk[1], ipk[2], ipk[3]]
    nshell[1] = 1

    ilpt2 = pg.eta2x !== nothing
    etavec = @MVector zeros(3)
    eta2vec = @MVector zeros(3)
    etavec[1] = pg.etax[ipk[1], ipk[2], ipk[3]]
    etavec[2] = pg.etay[ipk[1], ipk[2], ipk[3]]
    etavec[3] = pg.etaz[ipk[1], ipk[2], ipk[3]]
    if ilpt2
        eta2vec[1] = pg.eta2x[ipk[1], ipk[2], ipk[3]]
        eta2vec[2] = pg.eta2y[ipk[1], ipk[2], ipk[3]]
        eta2vec[3] = pg.eta2z[ipk[1], ipk[2], ipk[3]]
    end

    for L in 1:3
        Sshell[L][1] = etavec[L]
        S2shell[L][1] = eta2vec[L]
        Sshell[L][2] = 0.0
        S2shell[L][2] = 0.0
        Gshell[L][1] = 0.0
        Gshell[L][2] = 0.0
        Gshellf[L][1] = 0.0
        Gshellf[L][2] = 0.0
        for K in 1:3
            SRshell[L][K][1] = 0.0
            SRshell[L][K][2] = 0.0
        end
    end

    Fshell = 0.0
    nsh = 0
    ir2p = shells[2].r2
    dFbarp = Fbar[1]

    # initial rupp for ir2min
    rupp_init = sqrt(Float64(ir2min)) + 2.0 * hlatt
    ir2upp = floor(Int, rupp_init * rupp_init) + 1

    ifcrit = 1
    mupp = 1

    for jp in 2:npart
        ir2 = shells[jp].r2

        if ir2 != ir2p
            m += 1
            rad[m] = sqrt(Float64(ir2p))
            nshell[m] = nsh
            dFbar = nsh > 0 ? Fshell / nsh : 0.0
            rad3p = rad[m-1]^3
            rad3  = rad[m]^3
            Fbar[m] = (rad3p * Fbar[m-1] + 0.5 * (dFbarp + dFbar) * (rad3 - rad3p)) / rad3
            dFbarp = dFbar

            if ifcrit == 1
                m0 = m
                if ir2p > ir2min
                    rupp_init = rad[m0] + 2.0 * hlatt
                    ir2upp = floor(Int, rupp_init * rupp_init) + 1
                    ifcrit = 0
                    mupp = m
                end
            else
                if rad[m] < rupp_init
                    mupp = m
                end
            end

            ir2p = ir2
            Fshell = 0.0
            nsh = 0
            for L in 1:3
                Sshell[L][m+1] = 0.0
                S2shell[L][m+1] = 0.0
                Gshell[L][m] /= rad[m]
                Gshellf[L][m] /= rad[m]
                Gshell[L][m+1] = 0.0
                Gshellf[L][m+1] = 0.0
                for K in 1:3
                    SRshell[L][K][m] /= rad[m]
                    SRshell[L][K][m+1] = 0.0
                end
            end
        end

        # bounds check
        iv = (ipk[1] + shells[jp].di, ipk[2] + shells[jp].dj, ipk[3] + shells[jp].dk)
        if iv[1] < 1 || iv[1] > n1 || iv[2] < 1 || iv[2] > n2 || iv[3] < 1 || iv[3] > n3
            continue
        end

        nsh += 1
        delta_val = pg.delta[iv[1], iv[2], iv[3]]
        Fshell += delta_val
        etavec[1] = pg.etax[iv[1], iv[2], iv[3]]
        etavec[2] = pg.etay[iv[1], iv[2], iv[3]]
        etavec[3] = pg.etaz[iv[1], iv[2], iv[3]]
        if ilpt2
            eta2vec[1] = pg.eta2x[iv[1], iv[2], iv[3]]
            eta2vec[2] = pg.eta2y[iv[1], iv[2], iv[3]]
            eta2vec[3] = pg.eta2z[iv[1], iv[2], iv[3]]
        else
            eta2vec .= 0.0
        end

        ixsvec = (shells[jp].di, shells[jp].dj, shells[jp].dk)
        for L in 1:3
            Sshell[L][m+1] += etavec[L]
            S2shell[L][m+1] += eta2vec[L]
            Gshell[L][m+1] += delta_val * ixsvec[L]
            Gshellf[L][m+1] += pg.delta[iv[3], iv[1], iv[2]] * ixsvec[L]
            for K in 1:3
                SRshell[L][K][m+1] += etavec[L] * ixsvec[K]
            end
        end
    end
    # m is now the last shell index

    # ---- Phase 3: Gradient at Rf ----
    rmrf = 10.0
    mrf = 1
    ir2p_rf = shells[2].r2
    mrfi = 1
    for jp in 2:npart
        ir2 = shells[jp].r2
        if abs(Float64(ir2) - (Rfclvi / alatt)^2) < rmrf
            mrf = mrfi
            rmrf = min(rmrf, abs(Float64(ir2) - (Rfclvi / alatt)^2))
        end
        if ir2 != ir2p_rf
            mrfi += 1
            ir2p_rf = ir2
        end
    end

    # mlow/mupp for Rf gradient
    rlow_rf = max(rad[mrf] - 2.0, 0.0)
    mlow_rf = mrf
    if mrf > 1
        for m1 in mrf-1:-1:1
            rad[m1] > rlow_rf && (mlow_rf = m1)
        end
    end
    rupp_rf = rad[mrf] + 2.0
    mupp_rf = mrf
    for m1 in mrf+1:m
        rad[m1] < rupp_rf && (mupp_rf = m1)
    end

    akk_tab = _AKK_TAB
    _, gradpkrf, _ = kernel_strain(rad, Gshell, Gshellf, SRshell,
                                    mrf, mlow_rf, mupp_rf, akk_tab,
                                    wRnor, aRnor, 1.0, 1.0)

    # ---- Phase 4: Find fcrit crossing ----
    # All shells are built (no early break), so we search among all computed Fbar values.
    if Fbar[m0] >= fcrit
        # Search outward from m0 until Fbar < fcrit
        found = false
        for mm in m0:m
            if Fbar[mm] < fcrit
                m0 = mm
                found = true
                break
            end
        end
        if !found
            # Fbar >= fcrit everywhere: peak may have collapsed at outermost shell
            # Fortran sets mupp=m and continues (peak IS collapsed)
            mupp = m
        end
        if found
            rupp_m0 = rad[m0] + 2.0
            mupp = m0
            for mm in m0+1:m
                rad[mm] < rupp_m0 && (mupp = mm)
            end
        end
    else
        # step inward to find crossing
        mstart = max(m0 - 1, 1)
        found = false
        for mp in mstart:-1:1
            if Fbar[mp] >= fcrit
                found = true
                break
            end
            m0 = mp
        end
        if !found
            Threads.atomic_add!(_dump_no_fcrit, 1)
            return no_collapse()
        end
        rupp_m0 = rad[m0] + 2.0
        mupp_new = m0
        for mm in m0+1:m
            rad[mm] < rupp_m0 && (mupp_new = mm)
        end
        mupp = mupp_new
    end
    # Now Fbar[m0-1] >= fcrit > Fbar[m0]

    # ---- Phase 5: Step inward checking virialization ----
    # mlow for m0
    rlow = max(rad[m0] - 2.0, 0.0)
    mlow = m0
    if m0 > 1
        for m1 in m0-1:-1:1
            rad[m1] > rlow && (mlow = m1)
        end
    end

    # Pre-declare variables needed for the collapse interpolation
    # SMatrix/SVector have value semantics — no copy() needed
    zvir1 = -1.0
    zvir1p = -1.0
    Ebar_last = _zero_smat3
    gradpk_last = _zero_svec3
    Frhoh = 0.0
    Frho = 0.0
    e_v = 0.0; p_v = 0.0
    strain_mat = _zero_smat3
    gradpk_prev = _zero_svec3
    gradpkf_last = _zero_svec3
    gradpkf_prev = _zero_svec3
    Frhpk = 0.0; Fnupk = 0.0; Fevpk = 0.0; Fpvpk = 0.0
    collapsed = false

    # Compute strain at m0
    Ebar_m0, gradpk_m0, gradpkf_m0 = kernel_strain(rad, Gshell, Gshellf, SRshell,
                                           m0, mlow, mupp, akk_tab,
                                           wRnor, aRnor, 1.0, 1.0)

    Lam_m0, _, iflag_m0 = get_evals(Ebar_m0)
    Frho_m0 = Lam_m0[1] + Lam_m0[2] + Lam_m0[3]
    e_v_m0 = 0.0; p_v_m0 = 0.0
    if Frho_m0 > 0.0
        e_v_m0 = 0.5 * (Lam_m0[3] - Lam_m0[1]) / Frho_m0
        p_v_m0 = 0.5 * (Lam_m0[3] + Lam_m0[1] - 2.0 * Lam_m0[2]) / Frho_m0
    end

    Frhoh_m0 = Frho_m0
    Frho_m0 = Fbar[m0]
    # Fortran uses uninitialized Frhoc at this point, effectively returning
    # zvir1p = -1. Replicate this to match Fortran output.
    zvir1p_m0 = -1.0

    if zvir1p_m0 >= ZZon
        # Already virialized at m0
        zvir1 = -1.0
        zvir1p = zvir1p_m0
        Ebar_last = Ebar_m0
        gradpk_last = gradpk_m0
        gradpkf_last = gradpkf_m0
        strain_mat = Ebar_m0
        gradpk_prev = gradpk_m0
        gradpkf_prev = gradpkf_m0
        Frhoh = Frhoh_m0
        Frho = Frho_m0
        e_v = e_v_m0; p_v = p_v_m0
        Frhpk = Frhoh_m0
        Fnupk = Frho_m0
        Fevpk = e_v_m0
        Fpvpk = p_v_m0
        collapsed = true
    else
        zvir1 = zvir1p_m0
        Frhpk = Frhoh_m0
        Fnupk = Frho_m0
        Fevpk = e_v_m0
        Fpvpk = p_v_m0
        gradpk_prev = gradpk_m0
        gradpkf_prev = gradpkf_m0
        strain_mat = Ebar_m0
        Frhoh = Frhoh_m0
        Frho = Frho_m0
        e_v = e_v_m0; p_v = p_v_m0
        Ebar_last = Ebar_m0
        gradpk_last = gradpk_m0
        gradpkf_last = gradpkf_m0

        if m0 == 1
            Threads.atomic_add!(_dump_m0_one, 1)
            return no_collapse()
        end

        # Step inward from m0-1 to 1
        mupp_p = mupp
        for mp in m0-1:-1:1
            rupp_mp = rad[mp] + 2.0
            rlow_mp = max(rad[mp] - 2.0, 0.0)

            if rad[mupp] > rupp_mp
                mupp_new = mp
                for m1 in mp+1:mupp
                    rad[m1] < rupp_mp && (mupp_new = m1)
                end
                mupp = mupp_new
            end

            if mlow > 1 && rad[mlow-1] > rlow_mp
                mlownew = mlow
                for m1 in mlow-1:-1:1
                    rad[m1] > rlow_mp && (mlownew = m1)
                end
                mlow = mlownew
            end

            Ebar_mp, gradpk_mp, gradpkf_mp = kernel_strain(rad, Gshell, Gshellf, SRshell,
                                                   mp, mlow, mupp, akk_tab,
                                                   wRnor, aRnor, 1.0, 1.0)

            Lam_mp, _, iflag_mp = get_evals(Ebar_mp)
            Frho_mp = Lam_mp[1] + Lam_mp[2] + Lam_mp[3]
            e_v_mp = 0.0; p_v_mp = 0.0
            if Frho_mp > 0.0
                e_v_mp = 0.5 * (Lam_mp[3] - Lam_mp[1]) / Frho_mp
                p_v_mp = 0.5 * (Lam_mp[3] + Lam_mp[1] - 2.0 * Lam_mp[2]) / Frho_mp
            end

            Frhoh_mp = Frho_mp
            Frho_mp = Fbar[mp]

            zvir1p_mp = -1.0
            if iflag_mp == 0 && Frho_mp > 0.0
                poe_mp = e_v_mp < 1e-5 ? 0.0 : p_v_mp / e_v_mp
                zvir1p_mp = interpolate(ct, log10(Frho_mp), e_v_mp, poe_mp)
            end

            m0_next = mp + 1
            if zvir1p_mp >= ZZon
                m0 = m0_next
                zvir1p = zvir1p_mp
                Ebar_last = Ebar_mp
                gradpk_last = gradpk_mp
                gradpkf_last = gradpkf_mp
                # Fortran updates these BEFORE the zvir1p check (lines 773-774),
                # so they reflect the current mp, not the previous iteration
                Frhoh = Frhoh_mp
                Frho = Frho_mp
                e_v = e_v_mp; p_v = p_v_mp
                # Keep zvir1, Frhpk, Fnupk, etc. from previous iteration
                collapsed = true
                break
            end
            mupp_p = mupp
            zvir1 = zvir1p_mp
            Frhpk = Frhoh_mp
            Fnupk = Frho_mp
            Fevpk = e_v_mp
            Fpvpk = p_v_mp
            gradpk_prev = gradpk_mp
            gradpkf_prev = gradpkf_mp
            strain_mat = Ebar_mp
            Frhoh = Frhoh_mp
            Frho = Frho_mp
            e_v = e_v_mp; p_v = p_v_mp
            Ebar_last = Ebar_mp
            gradpk_last = gradpk_mp
            gradpkf_last = gradpkf_mp
        end
    end

    if !collapsed
        Threads.atomic_add!(_dump_no_collapse, 1)
        return no_collapse()
    end

    # ---- Phase 6: Interpolate to exact collapse radius ----
    dZvir = zvir1p - zvir1
    if zvir1 > 0.0 && dZvir != 0.0
        RTHL3 = rad[m0-1]^3 + (rad[m0]^3 - rad[m0-1]^3) * (zvir1p - ZZon) / dZvir
        RTHL = RTHL3^(1.0/3.0)
    else
        RTHL = rad[m0-1]
    end

    if RTHL <= 0.0
        Threads.atomic_add!(_dump_rthl_neg, 1)
        return no_collapse()
    end

    RTHL3 = RTHL^3
    RTHL5 = RTHL3 * RTHL * RTHL

    rad3p = rad[m0-1]^3
    rad3  = rad[m0]^3
    drad3 = rad3 - rad3p
    dFbar = Fbar[m0] - Fbar[m0-1]

    if zvir1 > 0.0 && drad3 != 0.0
        frac = (RTHL3 - rad3p) / drad3
        Fbarx    = Fbar[m0-1] + frac * dFbar
        Frhpk    = Frhoh + frac * (Frhpk - Frhoh)
        Fnupk    = Frho + frac * (Fnupk - Frho)
        Fevpk    = e_v + frac * (Fevpk - e_v)
        Fpvpk    = p_v + frac * (Fpvpk - p_v)
        strain_mat = Ebar_last + frac * (strain_mat - Ebar_last)
        gradpk_final = gradpk_last + frac * (gradpk_prev - gradpk_last)
        gradpkf_final = gradpkf_last + frac * (gradpkf_prev - gradpkf_last)
    else
        Fbarx = Fbar[m0-1]
        Fnupk = Frho
        Fevpk = e_v
        Fpvpk = p_v
        strain_mat = Ebar_last
        gradpk_final = gradpk_prev
        gradpkf_final = gradpkf_prev
    end

    strain_mat = normalize_strain(strain_mat, Fbarx)
    Lam_f, _, _ = get_evals(strain_mat)
    Frhoc = Lam_f[1] + Lam_f[2] + Lam_f[3]
    e_v_f = Frhoc > 0.0 ? 0.5 * (Lam_f[3] - Lam_f[1]) / Frhoc : 0.0
    p_v_f = Frhoc > 0.0 ? 0.5 * (Lam_f[3] + Lam_f[1] - 2.0 * Lam_f[2]) / Frhoc : 0.0

    strain_mat = Ebar_last * (Fbarx / Frhoh)

    # ---- Phase 7: Final properties ----
    # Average displacement over collapsed shells
    Sbar  = @MVector zeros(3)
    Sbar2 = @MVector zeros(3)
    nSbar = 0
    for m1 in 1:m0-1
        if nshell[m1] > 0
            for L in 1:3
                Sbar[L]  += Sshell[L][m1]
                Sbar2[L] += S2shell[L][m1]
            end
            nSbar += nshell[m1]
        end
    end
    if nSbar > 0
        Sbar  ./= nSbar
        Sbar2 ./= nSbar
    end

    # Energy factor Srb
    Srb = 0.0
    rad5p = 0.0
    if m0 > 2
        for mp in 2:m0-1
            rad5 = rad[mp]^5
            Srb += 0.5 * (Fbar[mp-1] + Fbar[mp]) * (rad5 - rad5p)
            rad5p = rad5
        end
    end
    if zvir1 > 0.0 && dZvir != 0.0
        Srb += 0.5 * (Fbar[m0-1] + Fbarx) * (RTHL5 - rad5p)
    end
    if Fbarx * RTHL5 > 0.0
        Srb /= (Fbarx * RTHL5)
    end

    # Mask particles within RTHL
    if pg.mask !== nothing
        for jp in 1:npart
            sqrt(Float64(shells[jp].r2)) > RTHL && continue
            iv = (ipk[1] + shells[jp].di, ipk[2] + shells[jp].dj, ipk[3] + shells[jp].dk)
            (iv[1] < nbuff + 1 || iv[1] > n1 - nbuff ||
             iv[2] < nbuff + 1 || iv[2] > n2 - nbuff ||
             iv[3] < nbuff + 1 || iv[3] > n3 - nbuff) && continue
            pg.mask[iv[1], iv[2], iv[3]] = 1
        end
    end

    # Compute d2F: average Laplacian over shells within RTHL (Fortran lines 909-927)
    d2F = 0.0
    if pg.lapd !== nothing
        navg = 0
        for jp in 2:npart
            sqrt(Float64(shells[jp].r2)) > RTHL && continue
            iv = (ipk[1] + shells[jp].di, ipk[2] + shells[jp].dj, ipk[3] + shells[jp].dk)
            (iv[1] < 1 || iv[1] > n1 || iv[2] < 1 || iv[2] > n2 ||
             iv[3] < 1 || iv[3] > n3) && break
            navg += 1
            d2F += Float64(pg.lapd[iv[1], iv[2], iv[3]])
        end
        if navg > 0
            d2F /= navg
        end
    end

    # ---- Formation redshift (zvir_half) ----
    # Collapse redshift at RTHL/2, computed by evaluating strain tensor at
    # 10 radii between RTHL/2 and RTHL and taking max(zvir).
    # Matches Fortran lines 920-1000 of peakvoidsubs.f90.
    zvir_half = -1.0
    if RTHL / hlatt >= 3.0
        for jj in 1:10
            rcur = RTHL * (((jj - 1) / 9.0 * 0.5 + 0.5)^(1.0/3.0))

            # Find mlow, mupp, m0 for rcur
            rupp_rc = rcur + 2.0 * hlatt
            mupp_rc = 1
            for m1 in 2:m-1
                if rad[m1] > rupp_rc
                    mupp_rc = m1
                    break
                elseif rad[m1] == 0.0
                    mupp_rc = m1 - 1
                    break
                end
                mupp_rc = m1
            end

            rlow_rc = max(rcur - 2.0 * hlatt, 0.0)
            mlow_rc = 1
            for m1 in 1:m-1
                if rad[m1] > rlow_rc
                    mlow_rc = m1
                    break
                end
            end

            m0_rc = 1
            for m1 in 1:m-1
                if rad[m1] > rcur
                    m0_rc = m1
                    break
                end
                m0_rc = m1
            end

            # Kernel-weighted strain at rcur
            # kernel_strain already normalizes by wRnor/rad[mp] internally
            Ebar_rc, _, _ = kernel_strain(rad, Gshell, Gshellf, SRshell,
                                          m0_rc, mlow_rc, mupp_rc, akk_tab,
                                          wRnor, aRnor, 1.0, 1.0)

            # Eigenvalue decomposition
            Lam_rc, _, iflag_rc = get_evals(Ebar_rc)
            Frho_rc = Lam_rc[1] + Lam_rc[2] + Lam_rc[3]
            e_v_rc = 0.0; p_v_rc = 0.0
            if Frho_rc > 0.0
                e_v_rc = 0.5 * (Lam_rc[3] - Lam_rc[1]) / Frho_rc
                p_v_rc = 0.5 * (Lam_rc[3] + Lam_rc[1] - 2.0 * Lam_rc[2]) / Frho_rc
            end

            Frhoc_rc = Fbar[m0_rc]
            poe_rc = e_v_rc < 1e-5 ? 0.0 : p_v_rc / e_v_rc
            zvir_rc = -1.0
            if iflag_rc == 0 && Frhoc_rc > 0.0
                zvir_rc = interpolate(ct, log10(Frhoc_rc), e_v_rc, poe_rc)
            end
            zvir_half = max(zvir_rc - 1.0, zvir_half)
        end
    end

    return PeakResult(RTHL, Srb, SVector{3}(Sbar), SVector{3}(Sbar2),
                      strain_mat, gradpk_final,
                      gradpkrf, Fbarx, e_v_f, p_v_f, zvir_half,
                      gradpkf_final, d2F)
end

end # module RadialShell
