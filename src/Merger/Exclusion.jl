module Exclusion

export SpatialHash, build_hash, sphere_overlap,
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
    build_hash(x, y, z, nhalo; domain_min, domain_max) -> SpatialHash

Build a spatial hash for `nhalo` halos with positions `x[1:nhalo]`, etc.
Domain is divided into NC³ cells spanning [domain_min, domain_max] per axis.
"""
function build_hash(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                    z::AbstractVector{<:Real}, nhalo::Int;
                    domain_min::NTuple{3,Float64}, domain_max::NTuple{3,Float64})
    Lx = domain_max[1] - domain_min[1]
    Ly = domain_max[2] - domain_min[2]
    Lz = domain_max[3] - domain_min[3]
    L = max(Lx, Ly, Lz)
    cell_size = L / NC

    hoc = zeros(Int32, NC, NC, NC)
    ll = zeros(Int32, nhalo)

    for m in 1:nhalo
        ix = clamp(floor(Int, (Float64(x[m]) - domain_min[1]) / cell_size) + 1, 1, NC)
        iy = clamp(floor(Int, (Float64(y[m]) - domain_min[2]) / cell_size) + 1, 1, NC)
        iz = clamp(floor(Int, (Float64(z[m]) - domain_min[3]) / cell_size) + 1, 1, NC)
        ll[m] = hoc[ix, iy, iz]
        hoc[ix, iy, iz] = Int32(m)
    end

    return SpatialHash(hoc, ll, NC, cell_size, domain_min)
end

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
                           order::Vector{Int}, sh::SpatialHash)
    nhalo = length(order)
    dV = zeros(Float64, length(r))

    for rank in 1:nhalo
        i = order[rank]
        !survived[i] && continue

        ri = Float64(r[i])
        xi, yi, zi = Float64(x[i]), Float64(y[i]), Float64(z[i])

        # Search radius: r_i + r_max_possible, but conservatively use 2*r_i
        # (matches Fortran which uses 2*int(ri) cells)
        search_r = ri + Float64(r[order[1]])  # ri + largest surviving radius
        dcell = ceil(Int, search_r / sh.cell_size)
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
                    dx = Float64(x[j]) - xi
                    dy = Float64(y[j]) - yi
                    dz = Float64(z[j]) - zi
                    dist = sqrt(dx^2 + dy^2 + dz^2)

                    if dist < ri + rj && dist > abs(ri - rj)
                        v1, v2 = sphere_overlap(dist, ri, rj)
                        # Only accumulate for the smaller halo to avoid
                        # double-counting: each pair adds to the smaller one.
                        # Following Fortran: process in order, accumulate both.
                        dV[i] += v1
                        dV[j] += v2
                    end
                end
                j = sh.ll[j]
            end
        end
    end

    # Correct for double-counting: each pair (i,j) was visited twice
    # (once when processing i, once when processing j)
    dV .*= 0.5

    # Reduce radii
    new_r = copy(Vector{Float64}(r))
    four_pi_thirds = 4.0 * π / 3.0
    for i in 1:length(r)
        !survived[i] && continue
        vol_old = four_pi_thirds * Float64(r[i])^3
        vol_new = vol_old - dV[i]
        if vol_new <= 0.0
            survived[i] = false
            new_r[i] = 0.0
        else
            new_r[i] = (vol_new / four_pi_thirds)^(1.0/3.0)
        end
    end

    return new_r
end

end # module Exclusion
