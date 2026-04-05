module HDF5Ext

using PeakPatch
using HDF5

import PeakPatch: write_catalog_hdf5

import PeakPatch.Cosmology: CosmologyParams, chi, E2,
    ChiToZTable, build_chi_to_z, chi_to_z
import PeakPatch.Catalog: HaloRecord, ExtHaloRecord

"""
    _rect_to_radec(x, y, z)

Convert Cartesian (x, y, z) in Mpc/h to (ra, dec) in radians.
Follows the same convention as pixell's `rect2ang`.
"""
function _rect_to_radec(x::Float64, y::Float64, z::Float64)
    r = sqrt(x^2 + y^2 + z^2)
    if r == 0.0
        return (0.0, 0.0)
    end
    dec = asin(clamp(z / r, -1.0, 1.0))
    ra  = atan(y, x)
    return (ra, dec)
end

"""
    write_catalog_hdf5(path, halos, cosmo; z_max=10.0)

Write a halo catalog to HDF5 in the format expected by XGPaint.jl and other
downstream tools.

# Output datasets
- `x`, `y`, `z`: comoving position [Mpc/h] (Float32)
- `ra`, `dec`: sky coordinates [radians] (Float32)
- `redshift`: cosmological redshift (Float32)
- `chi`: comoving distance [Mpc/h] (Float32)
- `M200m`: halo mass M₂₀₀ₘ [Msun/h] (Float32)
- `RTHL`: Lagrangian radius [Mpc/h] (Float32)
- `vx`, `vy`, `vz`: 1LPT velocities (Float32)

If `halos` are `ExtHaloRecord`, additional datasets are written:
- `vx2`, `vy2`, `vz2`: 2LPT velocities
- `overdensity`, `e_v`, `p_v`, `Rf`, `zform`
"""
function PeakPatch.write_catalog_hdf5(path::String,
                                       halos::Vector{<:Union{HaloRecord, ExtHaloRecord}},
                                       cosmo::CosmologyParams;
                                       z_max::Float64=10.0)
    isempty(halos) && (h5open(path, "w") do f; end; return)

    chi2z_table = build_chi_to_z(cosmo; z_max=z_max)

    nhalo = length(halos)

    # Mean matter density rho_m = 2.775e11 * Om * h² [Msun/h / (Mpc/h)³]
    rho_m = 2.775e11 * cosmo.Om * cosmo.h^2

    # Extract arrays
    pos_x = Float32[h.x for h in halos]
    pos_y = Float32[h.y for h in halos]
    pos_z = Float32[h.z for h in halos]
    rthl  = Float32[h.RTHL for h in halos]

    chi_arr = Float32[sqrt(Float64(h.x)^2 + Float64(h.y)^2 + Float64(h.z)^2) for h in halos]
    z_arr   = Float32[chi_to_z(chi2z_table, Float64(c)) for c in chi_arr]
    mass    = Float32[Float32(4.0/3.0 * π * rho_m * Float64(h.RTHL)^3) for h in halos]

    ra_arr  = Vector{Float32}(undef, nhalo)
    dec_arr = Vector{Float32}(undef, nhalo)
    for i in 1:nhalo
        ra, dec = _rect_to_radec(Float64(pos_x[i]), Float64(pos_y[i]), Float64(pos_z[i]))
        ra_arr[i]  = Float32(ra)
        dec_arr[i] = Float32(dec)
    end

    h5open(path, "w") do f
        f["x"]    = pos_x
        f["y"]    = pos_y
        f["z"]    = pos_z
        f["ra"]   = ra_arr
        f["dec"]  = dec_arr
        f["redshift"] = z_arr
        f["chi"]  = chi_arr
        f["M200m"] = mass
        f["RTHL"] = rthl
        f["vx"]   = Float32[h.vx for h in halos]
        f["vy"]   = Float32[h.vy for h in halos]
        f["vz"]   = Float32[h.vz for h in halos]

        # Extended fields
        if eltype(halos) <: ExtHaloRecord
            f["vx2"] = Float32[h.vx2 for h in halos]
            f["vy2"] = Float32[h.vy2 for h in halos]
            f["vz2"] = Float32[h.vz2 for h in halos]
            f["overdensity"] = Float32[h.overdensity for h in halos]
            f["e_v"]  = Float32[h.e_v for h in halos]
            f["p_v"]  = Float32[h.p_v for h in halos]
            f["Rf"]   = Float32[h.Rf for h in halos]
            f["zform"] = Float32[h.zform for h in halos]
        end
    end
end

end # module HDF5Ext
