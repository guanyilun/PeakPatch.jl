module Exclusion

export SpatialHash, build_hash, auto_hash_nc, sphere_overlap,
       lagrangian_exclusion!, volume_reduction!

"""
    SpatialHash

256³ linked-list spatial hash for O(1) neighbor lookup.
Halos are assigned to cells based on their (x,y,z) positions.
`hoc[i,j,k]` points to the first halo in cell (i,j,k), and
`ll[m]` points to the next halo in the same cell (0 = end).
"""
struct SpatialHash
    hoc::Array{Int32,3}     # head-of-chain: nc×nc×nc
    ll::Vector{Int32}       # linked list: nhalo
    nc::Int                 # number of cells per dimension
    cell_size::Float64      # physical size of each cell
    origin::NTuple{3,Float64}  # minimum corner of domain
end

const NC = 256  # cells per dimension, matching Fortran

"""
    build_hash(x, y, z, nhalo; domain_min, domain_max, nc=NC) -> SpatialHash

Build a spatial hash for `nhalo` halos with positions `x[1:nhalo]`, etc.
Domain is divided into nc³ cells spanning [domain_min, domain_max] per axis.
The exclusion searches are exhaustive within their radius, so `nc` changes only
the cost, never the survivors (see `auto_hash_nc`).
"""
function build_hash(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                    z::AbstractVector{<:Real}, nhalo::Int;
                    domain_min::NTuple{3,Float64}, domain_max::NTuple{3,Float64},
                    nc::Int=NC)
    Lx = domain_max[1] - domain_min[1]
    Ly = domain_max[2] - domain_min[2]
    Lz = domain_max[3] - domain_min[3]
    L = max(Lx, Ly, Lz)
    cell_size = L / nc

    hoc = zeros(Int32, nc, nc, nc)
    ll = zeros(Int32, nhalo)

    for m in 1:nhalo
        ix = clamp(floor(Int, (Float64(x[m]) - domain_min[1]) / cell_size) + 1, 1, nc)
        iy = clamp(floor(Int, (Float64(y[m]) - domain_min[2]) / cell_size) + 1, 1, nc)
        iz = clamp(floor(Int, (Float64(z[m]) - domain_min[3]) / cell_size) + 1, 1, nc)
        ll[m] = hoc[ix, iy, iz]
        hoc[ix, iy, iz] = Int32(m)
    end

    return SpatialHash(hoc, ll, nc, cell_size, domain_min)
end

"""
    auto_hash_nc(nhalo) -> Int

Hash resolution of about one halo per cell, clamped to 16..1024 per axis (≤ 4.3 GB). The
fixed Fortran NC=256 gives ~20 halos per 20.9 Mpc/h cell at Websky scale, and the
linked-list walks dominate (validation/performance/PERFORMANCE_2026-09-26.md: 346 s →
13 s on 9.8M halos with 5 Mpc/h cells, identical survivors).
"""
auto_hash_nc(nhalo::Integer) = clamp(round(Int, cbrt(nhalo)), 16, 1024)

"""Cell indices for a position."""
@inline function _cell_idx(sh::SpatialHash, px::Real, py::Real, pz::Real)
    ix = clamp(floor(Int, (Float64(px) - sh.origin[1]) / sh.cell_size) + 1, 1, sh.nc)
    iy = clamp(floor(Int, (Float64(py) - sh.origin[2]) / sh.cell_size) + 1, 1, sh.nc)
    iz = clamp(floor(Int, (Float64(pz) - sh.origin[3]) / sh.cell_size) + 1, 1, sh.nc)
    return (ix, iy, iz)
end

"""
    sphere_overlap(d, r1, r2) -> (v1, v2)

Analytic volume of intersection caps when two spheres of radii r1, r2
are separated by distance d.  Returns (v1, v2) where v1 is the cap
volume subtracted from sphere 1, v2 from sphere 2.
Uses the Wolfram Sphere-Sphere Intersection formula.
Returns (0,0) if spheres don't overlap or one is inside the other.
"""
function sphere_overlap(d::Float64, r1::Float64, r2::Float64)
    if d >= r1 + r2
        return (0.0, 0.0)  # no overlap
    end
    if d <= abs(r1 - r2)
        # one sphere inside the other — handled by exclusion, not reduction
        return (0.0, 0.0)
    end

    # Cap heights (Wolfram formula)
    h1 = (r2 - r1 + d) * (r2 + r1 - d) / (2.0 * d)
    h2 = (r1 - r2 + d) * (r1 + r2 - d) / (2.0 * d)

    # Cap volumes: V = (π/3) h² (3r - h)
    v1 = (π / 3.0) * h1^2 * (3.0 * r1 - h1)
    v2 = (π / 3.0) * h2^2 * (3.0 * r2 - h2)

    return (v1, v2)
end

"""
    lagrangian_exclusion!(survived, x, y, z, r, order, sh)

Pass 1: Lagrangian exclusion.  Process halos in decreasing radius order.
For each surviving halo i, find neighbors with centers inside r_i and
mark them as dead (survived[j] = false).

`order` is a permutation such that `r[order[1]] >= r[order[2]] >= ...`.
"""
function lagrangian_exclusion!(survived::Vector{Bool},
                               x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                               z::AbstractVector{<:Real}, r::AbstractVector{<:Real},
                               order::Vector{Int}, sh::SpatialHash)
    nhalo = length(order)

    for rank in 1:nhalo
        i = order[rank]
        !survived[i] && continue

        ri = Float64(r[i])
        xi, yi, zi = Float64(x[i]), Float64(y[i]), Float64(z[i])

        # Search radius in cells
        dcell = ceil(Int, ri / sh.cell_size)
        ci, cj, ck = _cell_idx(sh, xi, yi, zi)

        for diz in -dcell:dcell, diy in -dcell:dcell, dix in -dcell:dcell
            jx = dix + ci
            jy = diy + cj
            jz = diz + ck
            (jx < 1 || jx > sh.nc || jy < 1 || jy > sh.nc || jz < 1 || jz > sh.nc) && continue

            j = sh.hoc[jx, jy, jz]
            while j > 0
                if j != i && survived[j]
                    rj = Float64(r[j])
                    if ri >= rj
                        dx = Float64(x[j]) - xi
                        dy = Float64(y[j]) - yi
                        dz = Float64(z[j]) - zi
                        dist = sqrt(dx^2 + dy^2 + dz^2)
                        if dist < ri
                            survived[j] = false
                        end
                    end
                end
                j = sh.ll[j]
            end
        end
    end
end

"""
    volume_reduction!(survived, x, y, z, r, order, sh) -> new_r

Pass 2: Volume reduction.  For each pair of surviving halos whose spheres
overlap (d < r_i + r_j), compute the intersection volume and accumulate
a volume deficit dV.  After the sweep, reduce each halo's radius:
  r_new = (r_old³ - 3 dV / 4π)^(1/3)
Halos with negative remaining volume are killed.

Returns a new radius vector (original `r` is not modified).
"""
function volume_reduction!(survived::Vector{Bool},
                           x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                           z::AbstractVector{<:Real}, r::AbstractVector{<:Real},
                           order::Vector{Int}, sh::SpatialHash; fortran_ties::Bool=false)
    # Fortran merge_pkvd 'shared' reduction (exclusion.f90:120-204): every surviving pair
    # with d < r_i + r_j adds to EACH halo its own cap beyond the mid-plane (sphere_overlap),
    # all overlaps are accumulated first, then r_new = (r³ − 3ΔV/4π)^(1/3).
    # Each pair is visited once, from its larger member (`order` rank breaks radius ties), so
    # neighbours lie within 2 r_i. Fortran's neighbour test `ri >= rj` (exclusion.f90:122) puts an
    # equal-radius pair in both lists, so it is reduced twice; `fortran_ties=true` reproduces that.
    rank = similar(order); rank[order] = eachindex(order)
    dV = zeros(Float64, length(r))
    for i in order
        survived[i] || continue
        ri = Float64(r[i]); xi, yi, zi = Float64(x[i]), Float64(y[i]), Float64(z[i])
        dcell = ceil(Int, 2ri / sh.cell_size)
        ci, cj, ck = _cell_idx(sh, xi, yi, zi)
        for diz in -dcell:dcell, diy in -dcell:dcell, dix in -dcell:dcell
            jx = dix + ci; jy = diy + cj; jz = diz + ck
            (jx < 1 || jx > sh.nc || jy < 1 || jy > sh.nc || jz < 1 || jz > sh.nc) && continue
            j = sh.hoc[jx, jy, jz]
            while j > 0
                if j != i && survived[j] && rank[j] > rank[i]        # j is the smaller of the pair
                    rj = Float64(r[j])
                    d = sqrt((Float64(x[j]) - xi)^2 + (Float64(y[j]) - yi)^2 + (Float64(z[j]) - zi)^2)
                    if d < ri + rj
                        v1, v2 = sphere_overlap(d, ri, rj)
                        w = (fortran_ties && rj == ri) ? 2.0 : 1.0
                        dV[i] += w * v1; dV[j] += w * v2
                    end
                end
                j = sh.ll[j]
            end
        end
    end
    new_r = copy(Vector{Float64}(r))
    for i in eachindex(r)
        survived[i] || continue
        vol_new = 4π / 3 * Float64(r[i])^3 - dV[i]
        if vol_new <= 0.0
            survived[i] = false; new_r[i] = 0.0          # Fortran: NaN radius (r³−ΔV<0), never flagged
        else
            new_r[i] = cbrt(vol_new * 3 / (4π))
        end
    end
    return new_r
end

end # module Exclusion
