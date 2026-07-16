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
                            p1x, p1y, p1z, p2x, p2y, p2z,
                            nmesh::Int, nbuff::Int, alatt::Float64,
                            xbx::Float64, ybx::Float64, zbx::Float64,
                            obs::NTuple{3,Float64}, rmin::Float64, chi_max::Float64,
                            inv_dr::Float64, nrt::Int,
                            rt_D::Vector{Float64}, rt_c2::Vector{Float64},
                            rt_w::Vector{Vector{Float64}},
                            theta_pix::Float64, subdiv_max::Int, vec2pix::F) where {F}
    cen = 0.5 * (nmesh + 1)
    has2 = p2x !== nothing
    nk = length(maps)
    @inbounds for k in (nbuff+1):(nmesh-nbuff)
        qz = zbx + alatt * (k - cen)
        for j in (nbuff+1):(nmesh-nbuff)
            qy = ybx + alatt * (j - cen)
            for i in (nbuff+1):(nmesh-nbuff)
                qx = xbx + alatt * (i - cen)
                dqx = qx - obs[1]; dqy = qy - obs[2]; dqz = qz - obs[3]
                rq = sqrt(dqx*dqx + dqy*dqy + dqz*dqz)
                (rmin <= rq <= chi_max) || continue
                s1x = Float64(p1x[i,j,k]); s1y = Float64(p1y[i,j,k]); s1z = Float64(p1z[i,j,k])
                s2x = 0.0; s2y = 0.0; s2z = 0.0
                if has2
                    s2x = Float64(p2x[i,j,k]); s2y = Float64(p2y[i,j,k]); s2z = Float64(p2z[i,j,k])
                end
                ns = min(subdiv_max, max(1, ceil(Int, (alatt / rq) / theta_pix)))
                if ns == 1
                    D  = _rt_lerp(rt_D,  rq, inv_dr, nrt)
                    c2 = _rt_lerp(rt_c2, rq, inv_dr, nrt)
                    ex = qx + D*s1x + c2*s2x - obs[1]
                    ey = qy + D*s1y + c2*s2y - obs[2]
                    ez = qz + D*s1z + c2*s2z - obs[3]
                    pix = vec2pix(ex, ey, ez)
                    for ik in 1:nk
                        maps[ik][pix] += _rt_lerp(rt_w[ik], rq, inv_dr, nrt)
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
                        for ik in 1:nk
                            maps[ik][pix] += wsub * _rt_lerp(rt_w[ik], rqs, inv_dr, nrt)
                        end
                    end
                end
            end
        end
    end
    return nothing
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

Other kwargs mirror `run_multitile_split` (`ntile`, `seed`, `coarse_factor`,
`coarse_grid`, `use_gpu`, `devices`, `verbose`). `subdiv_max` caps the per-dimension
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
                                coarse_factor::Int=0, coarse_grid::Int=0,
                                use_gpu::Bool=false,
                                devices::Union{Nothing,AbstractVector{Int}}=nothing,
                                gpu_paint::Bool=false, nside::Int=0,
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
    n_workers = (use_gpu && devices !== nothing && length(devices) >= 1) ? length(devices) : 1

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

    # ---- Radial factor tables: everything per-cell depends only on r_Lagrangian ----
    nrt = 4096
    rt_dr = chi_max / (nrt - 1)
    inv_dr = 1.0 / rt_dr
    rt_D  = zeros(nrt); rt_c2 = zeros(nrt)
    rt_w  = [zeros(nrt) for _ in kernels]
    for irt in 2:nrt
        r = (irt - 1) * rt_dr
        z = chi_to_z(chi2z, r)
        a = 1.0 / (1.0 + z)
        D, _, _ = Dlinear_ab(a, growth_tables)   # 1st return = D (growth factor); 3rd is D/a
        # 2LPT coefficient: FIXED convention (MultiTile.jl post-95f04a1): -3/7, true D.
        # (MultiResolution's packing carried the pre-fix +3/7 and D/a — do not copy it.)
        Om_a = cosmo.Om * a^3 / (cosmo.Om * a^3 + cosmo.OL)
        rt_D[irt]  = D
        rt_c2[irt] = ilpt >= 2 ? (-3.0 / 7.0) * Om_a^(-1.0 / 143) * D^2 : 0.0
        for (ik, kern) in enumerate(kernels)
            rt_w[ik][irt] = if kern === :mass
                rho_m * alatt^3
            elseif kern === :kappa
                1.5 * cosmo.Om * (1.0 / 2997.92458)^2 * (1 + z) * (1 - r / chistar) / r *
                    alatt^3 / omega_pix
            else
                error("unknown field-map kernel: $kern (supported: :mass, :kappa)")
            end
        end
    end

    # Radial table as a matrix for the GPU painter: columns D, coef2, then weights
    rtM = Matrix{Float64}(undef, 0, 0)
    if gpu_paint
        rtM = Matrix{Float64}(undef, nrt, 2 + length(kernels))
        rtM[:, 1] = rt_D; rtM[:, 2] = rt_c2
        for ik in eachindex(kernels)
            rtM[:, 2+ik] = rt_w[ik]
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
        if use_gpu
            fn_multi = getglobal(_pp_parent(), :isolated_convolve_gpu_multi)
            outs = fn_multi(residual, pk, boxsize_local, nmesh;
                             kernels=[(0, 0, 0), (1, 1, 0), (1, 2, 0), (1, 3, 0)],
                             nshell=0, return_device=true)
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
        if gpu_paint
            fnp = getglobal(_pp_parent(), :paint_tile_field_gpu!)
            has2 = psi2_dev !== nothing
            p2 = has2 ? psi2_dev : psi_dev
            fnp(acc.maps_d, psi_dev[1], psi_dev[2], psi_dev[3], p2[1], p2[2], p2[3],
                acc.rt_d, nmesh, nbuff, alatt, xbx, ybx, zbx, obs, rmin_eff, chi_max,
                inv_dr, nrt, theta_pix, subdiv_max, nside, length(kernels), has2)
            psi_dev = nothing; psi2_dev = nothing
        else
            _paint_tile_field!(acc, psi_host[1], psi_host[2], psi_host[3],
                               psi2_host === nothing ? nothing : psi2_host[1],
                               psi2_host === nothing ? nothing : psi2_host[2],
                               psi2_host === nothing ? nothing : psi2_host[3],
                               nmesh, nbuff, alatt, xbx, ybx, zbx, obs, rmin_eff, chi_max,
                               inv_dr, nrt, rt_D, rt_c2, rt_w, theta_pix, subdiv_max, vec2pix)
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
            my_device = devices[wid]
            t = Threads.@spawn begin
                set_dev = getglobal(_pp_parent(), :set_cuda_device!)
                set_dev(my_device)
                my_acc = _make_acc()
                for idx in my_indices
                    process_tile_field!(idx, tile_ids[idx], my_acc)
                end
                worker_maps[wid] = _collect_acc(my_acc)
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
        total = worker_maps[1][ik]
        for wid in 2:n_workers
            total .+= worker_maps[wid][ik]
        end
        out[kern] = total
    end
    return out
end
