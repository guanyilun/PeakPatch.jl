module HDF5Ext

using PeakPatch
using HDF5

import PeakPatch: write_catalog_hdf5, write_compact_catalog, read_compact_catalog
import TOML

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

# ---------------------------------------------------------------------------------------------------------
# Compact, self-contained catalog: quantised columns + full provenance (docs in PeakPatch.write_compact_catalog)
# ---------------------------------------------------------------------------------------------------------
const COMPACT_FORMAT = "peakpatch-compact-catalog/1"

function _write_q(f, name, vals::AbstractVector{<:Real}, ::Type{T}, step, offset; units="", desc="",
                  compress=4) where {T<:Integer}
    lo, hi = typemin(T), typemax(T)
    q = Vector{T}(undef, length(vals)); nclip = 0
    @inbounds for i in eachindex(vals)
        v = round(Int64, (Float64(vals[i]) - offset) / step)
        (v < lo || v > hi) && (nclip += 1; v = clamp(v, lo, hi))
        q[i] = T(v)
    end
    n = length(q)
    d = create_dataset(f, name, datatype(T), dataspace((n,)); chunk=(max(1, min(n, 1 << 20)),), shuffle=true, deflate=compress)
    write(d, q)
    attrs(d)["step"] = Float64(step); attrs(d)["offset"] = Float64(offset)
    attrs(d)["units"] = units; attrs(d)["description"] = desc; attrs(d)["n_clipped"] = nclip
    nclip > 0 && @warn "write_compact_catalog: $nclip values of $name clipped to the $(T) range"
    return nothing
end

function _read_q(f, name)
    d = f[name]; q = read(d)
    st = Float64(attrs(d)["step"]); of = Float64(attrs(d)["offset"])
    return of .+ st .* Float64.(q)
end

function PeakPatch.write_compact_catalog(path::AbstractString, halos::Vector{<:Union{HaloRecord, ExtHaloRecord}},
                                          prov::AbstractDict; extra_fields::Bool=false, compress::Integer=4)
    n = length(halos)
    cfg = TOML.parse(get(prov, "recipe.config_toml", ""))
    Om = Float64(get(get(cfg, "cosmology", Dict()), "Om", NaN))
    isfinite(Om) || error("write_compact_catalog: provenance has no [cosmology] Om (needed for the mass definition)")
    rho_m = 2.775e11 * Om                                  # Msun/h per (Mpc/h)^3
    h5open(path, "w") do f
        attrs(f)["format"] = COMPACT_FORMAT
        attrs(f)["n_halos"] = n
        attrs(f)["mass_definition"] = "M = 4π/3 ρ_m R_TH³, ρ_m = 2.775e11 Ω_m Msun/h (Mpc/h)^-3, Ω_m = $(Om); R_TH as in the source catalog (abundance-matched if AM was applied)"
        attrs(f)["coordinates"] = "Eulerian comoving positions [Mpc/h] in the PeakPatch centred frame; observer in provenance recipe.config_toml [grid] cen*"
        for (c, nm) in ((h -> h.x, "x"), (h -> h.y, "y"), (h -> h.z, "z"))
            v = Float64[c(h) for h in halos]
            off = n > 0 ? floor(minimum(v)) - 1.0 : 0.0
            _write_q(f, nm, v, Int32, 0.01, off; units="Mpc/h", desc="comoving position", compress=compress)
        end
        for (c, nm) in ((h -> h.vx, "vx"), (h -> h.vy, "vy"), (h -> h.vz, "vz"))
            _write_q(f, nm, Float64[c(h) for h in halos], Int16, 1.0, 0.0; units="km/s", desc="peculiar velocity", compress=compress)
        end
        lm = Float64[log10(4π / 3 * rho_m * Float64(h.RTHL)^3) for h in halos]
        _write_q(f, "log10M", lm, UInt16, 1e-4, 10.0; units="log10(Msun/h)", desc="halo mass (see mass_definition)", compress=compress)
        if extra_fields
            eltype(halos) <: ExtHaloRecord || error("extra_fields=true needs an extended (ioutshear ≥ 1) catalog")
            _write_q(f, "zform", Float64[h.zform for h in halos], UInt16, 2e-4, -1.0; desc="formation redshift (collapse of R_TH/2)", compress=compress)
            _write_q(f, "e_v", Float64[h.e_v for h in halos], UInt16, 1e-5, 0.0; desc="ellipticity of the strain at R_TH", compress=compress)
            _write_q(f, "p_v", Float64[h.p_v for h in halos], UInt16, 2e-5, -0.6; desc="prolateness of the strain at R_TH", compress=compress)
        end
        g = create_group(f, "provenance")
        for (k, v) in prov
            attrs(g)[string(k)] = string(v)
        end
    end
    return path
end

function PeakPatch.read_compact_catalog(path::AbstractString)
    h5open(path, "r") do f
        attrs(f)["format"] == COMPACT_FORMAT || error("not a $(COMPACT_FORMAT) file: $path")
        g = f["provenance"]
        prov = Dict{String,String}(k => string(attrs(g)[k]) for k in keys(attrs(g)))
        cols = Dict{Symbol,Any}()
        for nm in ("x", "y", "z", "vx", "vy", "vz")
            cols[Symbol(nm)] = Float32.(_read_q(f, nm))
        end
        cols[:M] = 10 .^ _read_q(f, "log10M")
        for nm in ("zform", "e_v", "p_v")
            haskey(f, nm) && (cols[Symbol(nm)] = Float32.(_read_q(f, nm)))
        end
        cols[:provenance] = prov
        return (; cols...)
    end
end

end # module HDF5Ext
