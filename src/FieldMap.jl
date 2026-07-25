# FieldMap.jl — field-matter lightcone painting (included inside module MultiResolution).
#
# Websky-style field component (Stein+2020 §3.1.3): every Lagrangian lattice cell is a
# mass element a_latt³·ρ̄ that gets 2LPT-displaced to its lightcone position and painted
# onto a sky map with a redshift-dependent kernel W_F(z):
#
#     δF = (a_latt³/Ω_pix) · W_F(z_i),   z_i = z(|q_i − obs|)  (Lagrangian distance),
#     x_i = q_i + D(z_i)·ψ₁ + (3/7)·Ω_m(a)^(−1/143)·D(z_i)²·ψ₂
#
# The displacement convention (and the +3/7 sign carried by the pipeline's raw ψ₂) is
# copied from the halo record packing / Merger.finalize_eulerian — do NOT re-derive.
# Cells subtending more than one pixel are split into n³ sub-volumes (n ≤ subdiv_max),
# each painted independently with weight/n³ (suppresses near-observer shot noise).
#
# This is a SEPARATE pass from the halo pipeline: the Threefry RNG is deterministic, so
# the same seed regenerates bit-identical tile fields; only Phase-1 of the tile loop
# (residual → isolated FFT → ψ₁/ψ₂) is executed — no filters/peaks/shells/merge.
# See docs/field_lightcone_plan.md.

using ..Cosmology: chi_to_z

export run_multitile_fieldmap

"""
    ang2pix_ring(nside, x, y, z) -> pixel index (1-based, RING ordering)

HEALPix RING ang2pix for a direction vector (need not be normalized). Implemented from
the published HEALPix algorithm (Górski et al. 2005); pure scalar code so it runs both
on CPU and inside CUDA kernels (Phase-B GPU painting). Validated against Healpix.jl in
`test_fieldmap_smoke.jl`.
"""
@inline function ang2pix_ring(nside::Integer, x::Real, y::Real, z::Real)
    r = sqrt(x * x + y * y + z * z)
    zn = z / r
    za = abs(zn)
    tt = mod(atan(y, x), 2π) * (2 / π)          # in [0,4)
    if za <= 2 / 3
        temp1 = nside * (0.5 + tt)
        temp2 = nside * zn * 0.75
        jp = floor(Int, temp1 - temp2)          # ascending edge line
        jm = floor(Int, temp1 + temp2)          # descending edge line
        ir = nside + 1 + jp - jm                # ring counted from z=2/3
        kshift = 1 - (ir & 1)
        ip = mod(div(jp + jm - nside + kshift + 1, 2), 4nside)
        return 2nside * (nside - 1) + 4nside * (ir - 1) + ip + 1
    else
        tp = tt - floor(tt)
        tmp = nside * sqrt(3 * (1 - za))
        jp = floor(Int, tp * tmp)
        jm = floor(Int, (1 - tp) * tmp)
        ir = jp + jm + 1                        # ring counted from the nearest pole
        ip = mod(floor(Int, tt * ir), 4ir)
        return zn > 0 ? 2ir * (ir - 1) + ip + 1 : 12nside^2 - 2ir * (ir + 1) + ip + 1
    end
end

@inline function _rt_lerp(tab::Vector{Float64}, r::Float64, inv_dr::Float64, n::Int)
    x = r * inv_dr + 1.0
    i = unsafe_trunc(Int, x)
    i = ifelse(i < 1, 1, ifelse(i > n - 1, n - 1, i))
    t = x - i
    @inbounds return tab[i] * (1.0 - t) + tab[i+1] * t
end

# Paint the core cells of one tile into per-kernel maps. Function barrier: specializes
# on the vec2pix callable. ψ arrays are host Arrays (downloaded from device if needed).
function _paint_tile_field!(maps::Vector{Vector{Float64}},
                            p1x, p1y, p1z, p2x, p2y, p2z, pot,
                            nmesh::Int, nbuff::Int, alatt::Float64,
                            xbx::Float64, ybx::Float64, zbx::Float64,
                            obs::NTuple{3,Float64}, rmin::Float64, chi_max::Float64,
                            inv_dr::Float64, nrt::Int,
                            rt_D::Vector{Float64}, rt_c2::Vector{Float64},
                            rt_vf::Vector{Float64}, rt_w::Vector{Vector{Float64}},
                            vw::Vector{Bool}, pw::Vector{Bool},
                            excl::Union{Nothing,Array{Bool,3}},
                            theta_pix::Float64, subdiv_max::Int, vec2pix::F) where {F}
    cen = 0.5 * (nmesh + 1)
    has2 = p2x !== nothing
    nk = length(maps)
    need_v = any(vw)
    need_p = any(pw)
    @inbounds for k in (nbuff+1):(nmesh-nbuff)
        qz = zbx + alatt * (k - cen)
        for j in (nbuff+1):(nmesh-nbuff)
            qy = ybx + alatt * (j - cen)
            for i in (nbuff+1):(nmesh-nbuff)
                excl === nothing || !excl[i-nbuff, j-nbuff, k-nbuff] || continue
                qx = xbx + alatt * (i - cen)
                dqx = qx - obs[1]; dqy = qy - obs[2]; dqz = qz - obs[3]
                rq = sqrt(dqx*dqx + dqy*dqy + dqz*dqz)
                (rmin <= rq <= chi_max) || continue
                s1x = Float64(p1x[i,j,k]); s1y = Float64(p1y[i,j,k]); s1z = Float64(p1z[i,j,k])
                s2x = 0.0; s2y = 0.0; s2z = 0.0
                if has2
                    s2x = Float64(p2x[i,j,k]); s2y = Float64(p2y[i,j,k]); s2z = Float64(p2z[i,j,k])
                end
                pv = need_p ? Float64(pot[i,j,k]) : 0.0
                ns = min(subdiv_max, max(1, ceil(Int, (alatt / rq) / theta_pix)))
                if ns == 1
                    D  = _rt_lerp(rt_D,  rq, inv_dr, nrt)
                    c2 = _rt_lerp(rt_c2, rq, inv_dr, nrt)
                    ex = qx + D*s1x + c2*s2x - obs[1]
                    ey = qy + D*s1y + c2*s2y - obs[2]
                    ez = qz + D*s1z + c2*s2z - obs[3]
                    pix = vec2pix(ex, ey, ez)
                    pixq = need_p ? vec2pix(dqx, dqy, dqz) : pix
                    vr = 0.0
                    if need_v
                        vf = _rt_lerp(rt_vf, rq, inv_dr, nrt)
                        vr = vf * ((D*s1x + 2*c2*s2x)*dqx + (D*s1y + 2*c2*s2y)*dqy +
                                   (D*s1z + 2*c2*s2z)*dqz) / rq
                    end
                    for ik in 1:nk
                        w = _rt_lerp(rt_w[ik], rq, inv_dr, nrt)
                        vw[ik] && (w *= vr)
                        # pw kernels: LAGRANGIAN deposit (uniform lattice) — Eulerian
                        # painting imprints spurious delta*phi power on the smooth value
                        pw[ik] ? (maps[ik][pixq] += w * pv) : (maps[ik][pix] += w)
                    end
                else
                    wsub = 1.0 / (ns * ns * ns)
                    for c3 in 1:ns, c2i in 1:ns, c1 in 1:ns
                        qsx = qx + ((c1 - 0.5) / ns - 0.5) * alatt
                        qsy = qy + ((c2i - 0.5) / ns - 0.5) * alatt
                        qsz = qz + ((c3 - 0.5) / ns - 0.5) * alatt
                        dsx = qsx - obs[1]; dsy = qsy - obs[2]; dsz = qsz - obs[3]
                        rqs = sqrt(dsx*dsx + dsy*dsy + dsz*dsz)
                        (rmin <= rqs <= chi_max) || continue
                        D  = _rt_lerp(rt_D,  rqs, inv_dr, nrt)
                        c2 = _rt_lerp(rt_c2, rqs, inv_dr, nrt)
                        ex = qsx + D*s1x + c2*s2x - obs[1]
                        ey = qsy + D*s1y + c2*s2y - obs[2]
                        ez = qsz + D*s1z + c2*s2z - obs[3]
                        pix = vec2pix(ex, ey, ez)
                        pixq = need_p ? vec2pix(dsx, dsy, dsz) : pix
                        vr = 0.0
                        if need_v
                            vf = _rt_lerp(rt_vf, rqs, inv_dr, nrt)
                            vr = vf * ((D*s1x + 2*c2*s2x)*dsx + (D*s1y + 2*c2*s2y)*dsy +
                                       (D*s1z + 2*c2*s2z)*dsz) / rqs
                        end
                        for ik in 1:nk
                            w = wsub * _rt_lerp(rt_w[ik], rqs, inv_dr, nrt)
                            vw[ik] && (w *= vr)
                            pw[ik] ? (maps[ik][pixq] += w * pv) : (maps[ik][pix] += w)
                        end
                    end
                end
            end
        end
    end
    return nothing
end

# Rasterize the Lagrangian exclusion spheres of nearby halos into a core-sized Bool
# mask for one tile (true = cell inside a halo, skip). `bins` maps tile index →
# halo indices; ±1-neighbor search assumes max(R) < dcore (checked by the caller).
function _build_exclusion_mask(halos, bins::Dict{NTuple{3,Int},Vector{Int}},
                               it::Int, jt::Int, kt::Int, ntile::Int, dcore::Float64,
                               nmesh::Int, nbuff::Int, alatt::Float64,
                               xbx::Float64, ybx::Float64, zbx::Float64)
    ncore = nmesh - 2 * nbuff
    excl = zeros(Bool, ncore, ncore, ncore)
    cen = 0.5 * (nmesh + 1)
    half = dcore / 2
    @inbounds for dk in -1:1, dj in -1:1, di in -1:1
        bt = (it + di, jt + dj, kt + dk)
        haskey(bins, bt) || continue
        for n in bins[bt]
            hx = Float64(halos.x[n]); hy = Float64(halos.y[n]); hz = Float64(halos.z[n])
            R = Float64(halos.R[n])
            (abs(hx - xbx) <= half + R && abs(hy - ybx) <= half + R &&
             abs(hz - zbx) <= half + R) || continue
            R2 = R * R
            ilo = max(nbuff + 1, ceil(Int,  (hx - R - xbx) / alatt + cen))
            ihi = min(nmesh - nbuff, floor(Int, (hx + R - xbx) / alatt + cen))
            jlo = max(nbuff + 1, ceil(Int,  (hy - R - ybx) / alatt + cen))
            jhi = min(nmesh - nbuff, floor(Int, (hy + R - ybx) / alatt + cen))
            klo = max(nbuff + 1, ceil(Int,  (hz - R - zbx) / alatt + cen))
            khi = min(nmesh - nbuff, floor(Int, (hz + R - zbx) / alatt + cen))
            for k in klo:khi
                dz2 = (zbx + alatt * (k - cen) - hz)^2
                for j in jlo:jhi
                    dyz2 = dz2 + (ybx + alatt * (j - cen) - hy)^2
                    dyz2 > R2 && continue
                    for i in ilo:ihi
                        dx = xbx + alatt * (i - cen) - hx
                        if dx * dx + dyz2 <= R2
                            excl[i-nbuff, j-nbuff, k-nbuff] = true
                        end
                    end
                end
            end
        end
    end
    return excl
end

"""
    run_multitile_fieldmap(cfg; ntile, npix, vec2pix, kernels=[:kappa, :mass], ...)
        -> Dict{Symbol, Vector{Float64}}

Paint the LPT-displaced FIELD matter (every lattice cell) of a lightcone run onto sky
maps — the Websky "field component" (Stein+2020 §3.1.3, eq 3.11) generalized to a set of
kernels. Regenerates the tile fields deterministically from `seed` (same as
`run_multitile_split`) and executes only Phase-1 per tile (δ, ψ₁, ψ₂); no halo finding.

Required kwargs:
- `npix`: number of map pixels; `Ω_pix = 4π/npix` unless `omega_pix` is given.
- `vec2pix`: thread-safe callable `(x, y, z) -> pixel index in 1:npix` for an
  observer-centred direction vector (need not be normalized), e.g.
  `(x,y,z) -> Healpix.vec2pixRing(res, x, y, z)`.

Kernels (per-cell pixel contribution, lengths in Mpc/h):
- `:mass`  — ρ̄_m·a_latt³ [Msun/h]: total mass per pixel (map sums to the lightcone mass).
- `:kappa` — (3/2)·Ω_m·(H0/c)²·(1+z)·(1−χ/χ*)/χ · a_latt³/Ω_pix : Born CMB-lensing
  convergence of the full matter field (subtract the map mean for fluctuations).
  `chi_star` [Mpc/h] defaults to chi(z=1089) of the run cosmology (no radiation);
  pass Websky's 9656 (=14.2 Gpc × h) to match their hardwired value.
- `:tau`   — Thomson depth (Stein+2020 eq 3.22): σ_T·n_e,0·(1+z)²/χ² · a_latt³/Ω_pix
  with n_e,0 = f_e·ρ_b,0·x_e/m_p, f_e=0.9, He once-ionized above z=3, doubly below.
- `:ksz`   — −(v_r/c)·W_τ (eq 3.23): kSZ ΔT/T_CMB. v_r is the cell's LOS peculiar
  velocity v·q̂ with v = a·H·f·(D·ψ₁ + 2·D₂·ψ₂) (the `Merger.finalize_eulerian`
  convention, f₂≈2f) evaluated at the cell's Lagrangian distance.
- `:isw`   — linear ISW ΔT/T_CMB = 2∫dχ (aH/c)(f−1)·φ, φ = (3/2)Ωm(H0/c)²(D/a)·∇⁻²δ₀.
  The potential ∇⁻²δ₀ is generated per tile (coarse + fine, kernel −1/k²) and painted
  as the cell VALUE. Linear only (Websky isw.fits construction: halo catalogs unused);
  cells are painted at their (2LPT-displaced) positions — a second-order detail at ISW
  scales. NOTE: the potential is dominated by k ≲ 0.01 h/Mpc; use a LARGE coarse grid
  (coarse_factor ≥ 32 — same velocity-coherence argument as kSZ, but stronger).

Halo exclusion (Websky §3.1.3 "field component"): pass
`exclude_halos = (x=…, y=…, z=…, R=…)` (equal-length vectors; Lagrangian halo centers
in the same centered coordinates as the tiles, top-hat radii in Mpc/h, R < tile core
size). Cells whose Lagrangian centers fall inside any halo sphere are skipped.

Other kwargs mirror `run_multitile_split` (`ntile`, `seed`, `coarse_factor`,
`coarse_grid`, `use_gpu`, `devices`, `verbose`). `cpu_workers` spreads the tile loop
over that many CPU tasks when not using GPU devices (also lets the smoke test cover
the multi-worker dispatch without a GPU). `subdiv_max` caps the per-dimension
sub-cell splitting (Websky uses 5); `rmin` [Mpc/h] drops cells closer than this to the
observer (their splitting would be hopeless anyway; default 2 cells).

GPU painting (Phase B): `gpu_paint=true` with `nside` (and `npix == 12*nside^2`)
pixelizes on-device — own RING `ang2pix_ring` + Float64 atomic adds into a
device-resident map, no per-tile device→host ψ transfers, `vec2pix` not needed.
Requires `use_gpu=true`; maps match the CPU path up to float rounding
(cross-checked in `test_fieldmap_smoke.jl --gpu`).

Requires `ievol == 1` (lightcone). Returns one accumulated Float64 map per kernel.
"""
function run_multitile_fieldmap(cfg::PipelineConfig; ntile::Int, seed::Integer=42,
                                npix::Int, vec2pix=nothing, omega_pix::Float64=4π/npix,
                                kernels::Vector{Symbol}=[:kappa, :mass],
                                chi_star::Float64=0.0, subdiv_max::Int=3,
                                rmin::Float64=0.0,
                                exclude_halos::Union{Nothing,NamedTuple}=nothing,
                                coarse_factor::Int=0, coarse_grid::Int=0,
                                use_gpu::Bool=false,
                                devices::Union{Nothing,AbstractVector{Int}}=nothing,
                                gpu_paint::Bool=false, nside::Int=0,
                                cpu_workers::Int=1,
                                verbose::Bool=false)
    cfg.ievol == 1 || error("run_multitile_fieldmap requires ievol=1 (lightcone mode)")
    if use_gpu
        isdefined(_pp_parent(), :isolated_convolve_gpu) ||
            error("use_gpu=true requires CUDA.jl to be loaded (ext/CUDAExt.jl)")
    end
    if devices !== nothing && !isempty(devices)
        use_gpu || error("devices=$devices requires use_gpu=true")
    end
    if gpu_paint
        # Phase B: on-device RING pixelization + atomic accumulation (ext/CUDAExt.jl).
        use_gpu || error("gpu_paint=true requires use_gpu=true")
        nside > 0 || error("gpu_paint=true requires nside")
        npix == 12 * nside^2 || error("gpu_paint: npix must be 12*nside^2 (full-sky RING)")
    else
        vec2pix === nothing && error("vec2pix is required unless gpu_paint=true")
    end
    n_workers = if use_gpu && devices !== nothing && length(devices) >= 1
        length(devices)
    else
        gpu_paint && cpu_workers > 1 && error("cpu_workers > 1 requires gpu_paint=false")
        max(1, cpu_workers)
    end

    # ---- Geometry (identical to run_multitile_split) ----
    nmesh = cfg.n
    nbuff = cfg.nbuff
    nsub = nmesh - 2 * nbuff
    N = nsub * ntile + 2 * nbuff
    alatt = cfg.boxsize / nmesh
    boxsize_full = N * alatt
    dcore_box = nsub * alatt
    boxsize_local = nmesh * alatt

    # ---- Cosmology / lightcone ----
    Om_total = cfg.Omx + cfg.OmB
    cosmo = CosmologyParams(Om_total, cfg.OmB, cfg.Omvac, cfg.h, 0.965, 0.808)
    growth_tables = Dlinear_tables(cosmo)
    pk = load_pk(cfg.pkfile)
    ilpt = cfg.ilpt
    obs = (Float64(cfg.cenx), Float64(cfg.ceny), Float64(cfg.cenz))
    z_max = cfg.z_max
    chi2z = build_chi_to_z(cosmo; z_max=z_max + 1.0)
    chi_max = chi(z_max, cosmo)
    chistar = chi_star > 0 ? chi_star : chi(1089.0, cosmo)
    rho_m = 2.775e11 * cosmo.Om              # Msun/h per (Mpc/h)^3
    rmin_eff = rmin > 0 ? rmin : 2 * alatt
    theta_pix = sqrt(omega_pix)

    # Thomson-depth constant (Stein+2020 eq 3.22): σ_T·n_e,0 in (Mpc/h)⁻¹, comoving,
    # BEFORE the He-ionization factor x_e(z). n_e,0 = f_e·ρ_b,0·x_e/m_p with f_e=0.9;
    # ρ_crit,0/m_p = 11.2299·h² m⁻³, σ_T = 6.65246e-29 m², 1 Mpc = 3.0857e22 m.
    f_e = 0.9; Y_He = 0.245
    sigT_ne0 = 6.65246e-29 * 11.2299 * f_e * cfg.OmB * cosmo.h^2 * 3.0857e22 / cosmo.h
    c_kms = 299792.458

    # ---- Radial factor tables: everything per-cell depends only on r_Lagrangian ----
    nrt = 4096
    rt_dr = chi_max / (nrt - 1)
    inv_dr = 1.0 / rt_dr
    rt_D  = zeros(nrt); rt_c2 = zeros(nrt); rt_vf = zeros(nrt)
    rt_w  = [zeros(nrt) for _ in kernels]
    vw    = Bool[kern === :ksz for kern in kernels]   # velocity-weighted kernels (×v_r)
    pw    = Bool[kern === :isw for kern in kernels]   # potential-weighted kernels (×∇⁻²δ₀)
    needs_pot = any(pw)
    for irt in 2:nrt
        r = (irt - 1) * rt_dr
        z = chi_to_z(chi2z, r)
        a = 1.0 / (1.0 + z)
        D, f, _ = Dlinear_ab(a, growth_tables)   # 1st return = D (growth factor); 3rd is D/a
        # 2LPT coefficient: FIXED convention (MultiTile.jl post-95f04a1): -3/7, true D.
        # (MultiResolution's packing carried the pre-fix +3/7 and D/a — do not copy it.)
        Om_a = cosmo.Om * a^3 / (cosmo.Om * a^3 + cosmo.OL)
        rt_D[irt]  = D
        rt_c2[irt] = ilpt >= 2 ? (-3.0 / 7.0) * Om_a^(-1.0 / 143) * D^2 : 0.0
        # velocity factor km/s per Mpc/h of displacement (finalize_eulerian convention)
        rt_vf[irt] = a * 100.0 * sqrt(cosmo.Om * a^-3 + cosmo.OL) * f
        x_e = z < 3 ? (1 - Y_He / 2) : (1 - 3 * Y_He / 4)   # He doubly/once ionized
        w_tau = sigT_ne0 * x_e * (1 + z)^2 / r^2 * alatt^3 / omega_pix
        for (ik, kern) in enumerate(kernels)
            rt_w[ik][irt] = if kern === :mass
                rho_m * alatt^3
            elseif kern === :kappa
                1.5 * cosmo.Om * (1.0 / 2997.92458)^2 * (1 + z) * (1 - r / chistar) / r *
                    alatt^3 / omega_pix
            elseif kern === :tau
                w_tau
            elseif kern === :ksz
                -w_tau / c_kms          # painted weight × v_r [km/s] → −(v_r/c)·W_τ
            elseif kern === :isw
                # linear ISW ΔT/T = 2∫dχ (aH/c)(f−1)·φ with φ = (3/2)Ωm(H0/c)²(D/a)·pot0,
                # pot0 = ∇⁻²δ₀ painted as the cell VALUE (pw mask) — weight carries the rest
                3.0 * cosmo.Om * (1.0 / 2997.92458)^2 *
                    (a * 100.0 * sqrt(cosmo.Om * a^-3 + cosmo.OL) / c_kms) * (f - 1.0) *
                    (D / a) * alatt^3 / omega_pix / r^2
            else
                error("unknown field-map kernel: $kern (supported: :mass, :kappa, :tau, :ksz, :isw)")
            end
        end
    end

    # Radial table as a matrix for the GPU painter: columns D, coef2, vfac, then weights
    rtM = Matrix{Float64}(undef, 0, 0)
    if gpu_paint
        rtM = Matrix{Float64}(undef, nrt, 3 + length(kernels))
        rtM[:, 1] = rt_D; rtM[:, 2] = rt_c2; rtM[:, 3] = rt_vf
        for ik in eachindex(kernels)
            rtM[:, 3+ik] = rt_w[ik]
        end
    end
    vwmask = UInt32(0)
    for ik in eachindex(kernels)
        vw[ik] && (vwmask |= UInt32(1) << (ik - 1))
    end
    pwmask = UInt32(0)
    for ik in eachindex(kernels)
        pw[ik] && (pwmask |= UInt32(1) << (ik - 1))
    end

    # ---- Halo exclusion: bin halos by tile for fast per-tile mask rasterization ----
    excl_bins = nothing
    if exclude_halos !== nothing
        length(exclude_halos.x) == length(exclude_halos.y) == length(exclude_halos.z) ==
            length(exclude_halos.R) || error("exclude_halos: x/y/z/R length mismatch")
        isempty(exclude_halos.R) || maximum(exclude_halos.R) < dcore_box ||
            error("exclude_halos: max R must be smaller than the tile core size")
        x0 = -(ntile / 2) * dcore_box     # left edge of the tile grid, centered coords
        excl_bins = Dict{NTuple{3,Int},Vector{Int}}()
        for n in eachindex(exclude_halos.x)
            bt = (clamp(floor(Int, (exclude_halos.x[n] - x0) / dcore_box) + 1, 1, ntile),
                  clamp(floor(Int, (exclude_halos.y[n] - x0) / dcore_box) + 1, 1, ntile),
                  clamp(floor(Int, (exclude_halos.z[n] - x0) / dcore_box) + 1, 1, ntile))
            push!(get!(excl_bins, bt, Int[]), n)
        end
    end

    # ---- Coarse fields (Phase 1a; same as run_multitile_split, minus 2LPT/laplacian) ----
    if coarse_grid > 0
        M = coarse_grid
    elseif coarse_factor > 0
        M = ntile * coarse_factor
    else
        target_block = nmesh ÷ 3
        best_M = 0; best_dist = N
        for m in 2:N÷2
            N % m == 0 || continue
            dist = abs(N ÷ m - target_block)
            if dist < best_dist
                best_dist = dist; best_M = m
            end
        end
        M = best_M
    end
    @assert N % M == 0 "N=$N must be divisible by M=$M"

    coarse_noise = _downsample_noise(N, M, seed)
    coarse_k = rfft(coarse_noise)
    delta_coarse_k = copy(coarse_k)
    _periodic_convolve!(delta_coarse_k, pk, M, boxsize_full)
    delta_coarse = irfft(delta_coarse_k, M)
    psi_coarse = Vector{Array{Float32,3}}(undef, 3)
    for dim in 1:3
        psi_k = copy(coarse_k)
        _periodic_convolve!(psi_k, pk, M, boxsize_full; kernel_fn=_kernel_1lpt(dim))
        psi_coarse[dim] = irfft(psi_k, M)
    end
    pot_coarse = nothing
    if needs_pot
        pot_k = copy(coarse_k)
        _periodic_convolve!(pot_k, pk, M, boxsize_full; kernel_fn=_kernel_pot())
        pot_coarse = irfft(pot_k, M)
    end
    coarse_k = nothing
    verbose && @info "fieldmap Phase 1a: coarse fields done (N=$N, M=$M)"

    # ---- Tile list with horizon pruning (identical to run_multitile_split) ----
    tile_ids = NTuple{3,Int}[]
    for kt in 1:ntile, jt in 1:ntile, it in 1:ntile
        push!(tile_ids, (it, jt, kt))
    end
    half = dcore_box / 2.0
    filter!(tile_ids) do tid
        it, jt, kt = tid
        xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)
        dx = max(abs(xbx - obs[1]) - half, 0.0)
        dy = max(abs(ybx - obs[2]) - half, 0.0)
        dz = max(abs(zbx - obs[3]) - half, 0.0)
        sqrt(dx^2 + dy^2 + dz^2) <= chi_max
    end
    verbose && @info "fieldmap: $(length(tile_ids)) tiles inside z_max horizon"

    # ---- Per-tile work: Phase-1 fields, then paint core cells (host or device) ----
    # `acc` is a Vector{Vector{Float64}} of per-kernel host maps (CPU pixelization), or a
    # NamedTuple (maps_d, rt_d) of device arrays (gpu_paint).
    process_tile_field! = function (ti::Int, tid::NTuple{3,Int}, acc)
        it, jt, kt = tid
        residual = _generate_extended_residual(it, jt, kt, nsub, nmesh, N, seed,
                                                coarse_noise, M, 0)
        delta_tile = nothing
        psi_host = Vector{Array{Float32,3}}(undef, 3)
        psi2_host = nothing
        psi_dev = nothing
        psi2_dev = nothing
        pot_host = nothing
        pot_dev = nothing
        if use_gpu
            fn_multi = getglobal(_pp_parent(), :isolated_convolve_gpu_multi)
            klist = needs_pot ?
                [(0, 0, 0), (1, 1, 0), (1, 2, 0), (1, 3, 0), (5, 0, 0)] :
                [(0, 0, 0), (1, 1, 0), (1, 2, 0), (1, 3, 0)]
            outs = fn_multi(residual, pk, boxsize_local, nmesh;
                             kernels=klist, nshell=0, return_device=true)
            delta_self = outs[1]
            fn_interp = getglobal(_pp_parent(), :interpolate_to_tile_gpu)
            delta_tile = delta_self .+ fn_interp(delta_coarse, it, jt, kt, nsub, nmesh, N, M;
                                                  return_device=true)
            delta_self = nothing
            psi_dev = Vector{Any}(undef, 3)
            for dim in 1:3
                psi_long = fn_interp(psi_coarse[dim], it, jt, kt, nsub, nmesh, N, M;
                                      return_device=true)
                psi_dev[dim] = outs[dim+1] .+ psi_long
                psi_long = nothing
            end
            if needs_pot
                pot_long = fn_interp(pot_coarse, it, jt, kt, nsub, nmesh, N, M;
                                      return_device=true)
                pot_dev = outs[5] .+ pot_long
                pot_long = nothing
            end
            if ilpt >= 2
                fn2 = getglobal(_pp_parent(), :compute_2lpt_gpu)
                psi2_dev = fn2(delta_tile, nmesh, boxsize_local; return_device=true)
            end
            delta_tile = nothing
            if !gpu_paint
                for dim in 1:3
                    psi_host[dim] = Array(psi_dev[dim])
                    psi_dev[dim] = nothing
                end
                psi_dev = nothing
                if psi2_dev !== nothing
                    psi2_host = Vector{Array{Float32,3}}(undef, 3)
                    for dim in 1:3
                        psi2_host[dim] = Array(psi2_dev[dim])
                    end
                    psi2_dev = nothing
                end
                if pot_dev !== nothing
                    pot_host = Array(pot_dev)
                    pot_dev = nothing
                end
            end
        else
            delta_self = _isolated_convolve_dispatch(false, residual, pk,
                                                      boxsize_local, nmesh, 0, 0, 0)
            delta_tile = delta_self .+ _interpolate_to_tile(delta_coarse, it, jt, kt,
                                                             nsub, nmesh, N, M)
            delta_self = nothing
            for dim in 1:3
                psi_self = _isolated_convolve_dispatch(false, residual, pk,
                                                        boxsize_local, nmesh, 1, dim, 0)
                psi_host[dim] = psi_self .+ _interpolate_to_tile(psi_coarse[dim], it, jt, kt,
                                                                  nsub, nmesh, N, M)
            end
            if needs_pot
                pot_self = _isolated_convolve_dispatch(false, residual, pk,
                                                        boxsize_local, nmesh, 5, 0, 0)
                pot_host = pot_self .+ _interpolate_to_tile(pot_coarse, it, jt, kt,
                                                             nsub, nmesh, N, M)
            end
            if ilpt >= 2
                delta_tile_k = rfft(delta_tile)
                src2_local = delta_tile .^ 2 .* 0.5f0
                for d in 1:3
                    phi_k = copy(delta_tile_k)
                    _apply_kernel_inplace!(phi_k, nmesh, boxsize_local, _kernel_phi_ij(d, d))
                    src2_local .-= irfft(phi_k, nmesh) .^ 2 .* 0.5f0
                end
                for (di, dj) in ((1,2), (1,3), (2,3))
                    phi_k = copy(delta_tile_k)
                    _apply_kernel_inplace!(phi_k, nmesh, boxsize_local, _kernel_phi_ij(di, dj))
                    src2_local .-= irfft(phi_k, nmesh) .^ 2
                end
                src2_local_k = rfft(src2_local)
                psi2_host = Vector{Array{Float32,3}}(undef, 3)
                for dim in 1:3
                    psi2_k = copy(src2_local_k)
                    _apply_kernel_inplace!(psi2_k, nmesh, boxsize_local, _kernel_2lpt(dim))
                    psi2_host[dim] = irfft(psi2_k, nmesh)
                end
            end
            delta_tile = nothing
            GC.gc()
        end
        residual = nothing

        xbx, ybx, zbx = tile_center(it, jt, kt, ntile, dcore_box)
        excl = excl_bins === nothing ? nothing :
            _build_exclusion_mask(exclude_halos, excl_bins, it, jt, kt, ntile, dcore_box,
                                  nmesh, nbuff, alatt, xbx, ybx, zbx)
        if gpu_paint
            fnp = getglobal(_pp_parent(), :paint_tile_field_gpu!)
            has2 = psi2_dev !== nothing
            p2 = has2 ? psi2_dev : psi_dev
            fnp(acc.maps_d, psi_dev[1], psi_dev[2], psi_dev[3], p2[1], p2[2], p2[3],
                acc.rt_d, nmesh, nbuff, alatt, xbx, ybx, zbx, obs, rmin_eff, chi_max,
                inv_dr, nrt, theta_pix, subdiv_max, nside, length(kernels), has2,
                vwmask, excl, pot_dev, pwmask)
            psi_dev = nothing; psi2_dev = nothing; pot_dev = nothing
        else
            _paint_tile_field!(acc, psi_host[1], psi_host[2], psi_host[3],
                               psi2_host === nothing ? nothing : psi2_host[1],
                               psi2_host === nothing ? nothing : psi2_host[2],
                               psi2_host === nothing ? nothing : psi2_host[3],
                               pot_host,
                               nmesh, nbuff, alatt, xbx, ybx, zbx, obs, rmin_eff, chi_max,
                               inv_dr, nrt, rt_D, rt_c2, rt_vf, rt_w, vw, pw, excl,
                               theta_pix, subdiv_max, vec2pix)
        end
        verbose && @info "  fieldmap tile $ti/$(length(tile_ids)) ($it,$jt,$kt) painted"
        return nothing
    end

    _make_acc = function ()
        if gpu_paint
            alloc = getglobal(_pp_parent(), :fieldmap_gpu_alloc)
            maps_d, rt_d = alloc(npix, length(kernels), rtM)
            (maps_d=maps_d, rt_d=rt_d)
        else
            [zeros(Float64, npix) for _ in kernels]
        end
    end
    _collect_acc = function (acc)
        if gpu_paint
            coll = getglobal(_pp_parent(), :fieldmap_gpu_collect)
            Mh = coll(acc.maps_d)
            [Mh[:, ik] for ik in 1:length(kernels)]
        else
            acc
        end
    end

    # ---- Dispatch (mirrors run_multitile_split worker pattern) ----
    worker_maps = Vector{Any}(undef, n_workers)
    if n_workers > 1
        wid_of = Int[mod1(i, n_workers) for i in 1:length(tile_ids)]
        tasks = Task[]
        for wid in 1:n_workers
            my_indices = [i for i in 1:length(tile_ids) if wid_of[i] == wid]
            my_device = (use_gpu && devices !== nothing) ? devices[wid] : -1
            t = Threads.@spawn begin
                if my_device >= 0
                    set_dev = getglobal(_pp_parent(), :set_cuda_device!)
                    set_dev(my_device)
                end
                # `wacc` must not be assigned anywhere else in this function: a name
                # shared with the enclosing scope is captured by ALL worker closures
                # as one boxed variable, so every worker paints into the same map and
                # the reduction below self-adds it (2^(n_workers-1)× inflation).
                local wacc = _make_acc()
                for idx in my_indices
                    process_tile_field!(idx, tile_ids[idx], wacc)
                end
                worker_maps[wid] = _collect_acc(wacc)
            end
            push!(tasks, t)
        end
        foreach(wait, tasks)
    else
        my_acc = _make_acc()
        for (ti, tid) in enumerate(tile_ids)
            process_tile_field!(ti, tid, my_acc)
        end
        worker_maps[1] = _collect_acc(my_acc)
    end

    out = Dict{Symbol,Vector{Float64}}()
    for (ik, kern) in enumerate(kernels)
        total = copy(worker_maps[1][ik])
        for wid in 2:n_workers
            worker_maps[wid][ik] === worker_maps[1][ik] &&
                error("fieldmap: worker accumulators alias each other")
            total .+= worker_maps[wid][ik]
        end
        out[kern] = total
    end
    return out
end
