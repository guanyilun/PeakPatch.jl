module Merger

import ..Catalog: HaloRecord, ExtHaloRecord
import ..Exclusion: SpatialHash, build_hash, lagrangian_exclusion!, volume_reduction!

export merge_catalog

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

    # Pass 2: Volume reduction
    new_r = volume_reduction!(survived, x, y, z, r, order, sh)
    n_after_red = count(survived)
    verbose && @info "Reduction: $n_after_exc → $n_after_red halos (removed $(n_after_exc - n_after_red) by volume)"

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

end # module Merger
