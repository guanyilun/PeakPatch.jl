module Merger

import ..Catalog: HaloRecord, ExtHaloRecord
import ..Exclusion: SpatialHash, build_hash, lagrangian_exclusion!, volume_reduction!
import ..Cosmology: CosmologyParams, Dlinear_tables, Dlinear_ab, build_chi_to_z, chi_to_z

export merge_catalog, finalize_eulerian

"""
    merge_catalog(halos; verbose=false) -> Vector{HaloRecord} or Vector{ExtHaloRecord}

Apply Lagrangian exclusion and volume reduction to a raw halo catalog.

Two-pass algorithm (matching Fortran `merge_pkvd`):
1. **Exclusion**: Sort halos by radius (largest first). Remove halos whose
   Lagrangian center falls inside a larger halo's sphere.
2. **Reduction**: For surviving halos with overlapping spheres, subtract the
   intersection volume and shrink radii accordingly.

Works with both `HaloRecord` (11-field) and `ExtHaloRecord` (33-field) catalogs.
"""
function merge_catalog(halos::Vector{HaloRecord}; verbose::Bool=false)
    isempty(halos) && return HaloRecord[]
    _merge_impl(halos, verbose)
end

function merge_catalog(halos::Vector{ExtHaloRecord}; verbose::Bool=false)
    isempty(halos) && return ExtHaloRecord[]
    _merge_impl(halos, verbose)
end

function _merge_impl(halos::Vector{T}, verbose::Bool) where T <: Union{HaloRecord, ExtHaloRecord}
    nhalo = length(halos)

    # Extract position and radius arrays
    x = Float64[h.x for h in halos]
    y = Float64[h.y for h in halos]
    z = Float64[h.z for h in halos]
    r = Float64[h.RTHL for h in halos]

    # Sort order: descending by radius (largest first)
    order = sortperm(r; rev=true)

    # Domain bounds (with padding of max radius)
    rmax = maximum(r)
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    zmin, zmax = extrema(z)
    pad = 2.0 * rmax
    domain_min = (xmin - pad, ymin - pad, zmin - pad)
    domain_max = (xmax + pad, ymax + pad, zmax + pad)

    # Build spatial hash
    sh = build_hash(x, y, z, nhalo; domain_min=domain_min, domain_max=domain_max)

    survived = fill(true, nhalo)

    # Pass 1: Lagrangian exclusion
    lagrangian_exclusion!(survived, x, y, z, r, order, sh)
    n_after_exc = count(survived)
    verbose && @info "Exclusion: $nhalo → $n_after_exc halos (removed $(nhalo - n_after_exc))"

    # Pass 2: Volume reduction (disabled — Fortran has this commented out)
    new_r = copy(r)
    n_after_red = n_after_exc

    # Build output catalog with updated radii
    out = T[]
    sizehint!(out, n_after_red)
    for i in 1:nhalo
        !survived[i] && continue
        h = halos[i]
        push!(out, _replace_rthl(h, Float32(new_r[i])))
    end

    verbose && @info "Merged catalog: $nhalo → $(length(out)) halos"
    return out
end

"""Replace RTHL in a HaloRecord, keeping all other fields."""
function _replace_rthl(h::HaloRecord, new_rthl::Float32)
    HaloRecord(h.x, h.y, h.z, h.vx, h.vy, h.vz, new_rthl,
               h.vx2, h.vy2, h.vz2, h.overdensity)
end

function _replace_rthl(h::ExtHaloRecord, new_rthl::Float32)
    ExtHaloRecord(h.x, h.y, h.z, h.vx, h.vy, h.vz, new_rthl,
                  h.vx2, h.vy2, h.vz2, h.overdensity,
                  h.e_v, h.p_v,
                  h.strain_11, h.strain_22, h.strain_33,
                  h.strain_23, h.strain_13, h.strain_12,
                  h.d2F, h.zform,
                  h.grad_x, h.grad_y, h.grad_z,
                  h.gradf_x, h.gradf_y, h.gradf_z,
                  h.Rf, h.FcollvRf, h.d2FRf,
                  h.gradrf_x, h.gradrf_y, h.gradrf_z)
end

"""
    finalize_eulerian(halos, cosmo, obs; ievol=1, z_out=0.0) -> same-type catalog

Convert a RAW catalog (Lagrangian positions in x/y/z, 1LPT displacement in vx/vy/vz, and
2LPT displacement in vx2/vy2/vz2 — all Mpc/h, as produced by the fixed pipeline) into the final
Websky-style catalog: **Eulerian positions + km/s peculiar velocities**. Port of Fortran
`merge_pkvd` (lines 145-172):

    x_Eulerian = q + disp1 + disp2
    v[km/s]    = a·H(a)·f(a)·(disp1 + 2·disp2) = a·100·E(a)·f(a)·(disp1 + 2·disp2)

(H = 100·h·E; the h cancels for displacements in Mpc/h.) `f` is the exact growth rate
(`Dlinear_ab[2]`). Run AFTER `merge_catalog` (which does the Lagrangian-space exclusion).
`obs` is the observer position (Mpc/h); for ievol==1 the redshift comes from the Eulerian
distance to `obs`, for ievol==0 from `z_out`. The 2LPT fields (vx2/vy2/vz2) are consumed → set 0.
"""
function finalize_eulerian(halos::Vector{T}, cosmo::CosmologyParams, obs::NTuple{3,<:Real};
                           ievol::Integer=1, z_out::Real=0.0) where {T<:Union{HaloRecord,ExtHaloRecord}}
    isempty(halos) && return halos
    gt = Dlinear_tables(cosmo)
    Om, OL = cosmo.Om, cosmo.OL
    ox, oy, oz = Float64(obs[1]), Float64(obs[2]), Float64(obs[3])
    chi2z = ievol == 1 ? build_chi_to_z(cosmo; z_max=10.0) : nothing
    a_snap = 1.0 / (1.0 + z_out)
    out = Vector{T}(undef, length(halos))
    Threads.@threads for i in eachindex(halos)
        h = halos[i]
        d1x=Float64(h.vx); d1y=Float64(h.vy); d1z=Float64(h.vz)       # 1LPT displacement (Mpc/h)
        d2x=Float64(h.vx2); d2y=Float64(h.vy2); d2z=Float64(h.vz2)    # 2LPT displacement (Mpc/h)
        ex=Float64(h.x)+d1x+d2x; ey=Float64(h.y)+d1y+d2y; ez=Float64(h.z)+d1z+d2z   # Eulerian
        a = if ievol == 1
            r = sqrt((ex-ox)^2+(ey-oy)^2+(ez-oz)^2)
            r > 0 ? 1.0/(1.0+chi_to_z(chi2z, r)) : 1.0
        else
            a_snap
        end
        E = sqrt(Om*a^-3 + OL); _, f, _ = Dlinear_ab(a, gt); vfac = a*100.0*E*f   # km/s per Mpc/h
        vx=vfac*(d1x+2*d2x); vy=vfac*(d1y+2*d2y); vz=vfac*(d1z+2*d2z)
        out[i] = _replace_pos_vel(h, ex, ey, ez, vx, vy, vz)
    end
    return out
end

"""Replace position (x,y,z) and velocity (vx,vy,vz); zero the consumed 2LPT fields; keep the rest."""
function _replace_pos_vel(h::HaloRecord, x, y, z, vx, vy, vz)
    HaloRecord(Float32(x), Float32(y), Float32(z), Float32(vx), Float32(vy), Float32(vz),
               h.RTHL, 0f0, 0f0, 0f0, h.overdensity)
end
function _replace_pos_vel(h::ExtHaloRecord, x, y, z, vx, vy, vz)
    ExtHaloRecord(Float32(x), Float32(y), Float32(z), Float32(vx), Float32(vy), Float32(vz),
                  h.RTHL, 0f0, 0f0, 0f0, h.overdensity,
                  h.e_v, h.p_v,
                  h.strain_11, h.strain_22, h.strain_33,
                  h.strain_23, h.strain_13, h.strain_12,
                  h.d2F, h.zform,
                  h.grad_x, h.grad_y, h.grad_z,
                  h.gradf_x, h.gradf_y, h.gradf_z,
                  h.Rf, h.FcollvRf, h.d2FRf,
                  h.gradrf_x, h.gradrf_y, h.gradrf_z)
end

end # module Merger
