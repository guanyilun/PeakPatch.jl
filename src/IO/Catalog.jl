module Catalog

"""Basic 11-field halo record (ioutshear==0)."""
struct HaloRecord
    x::Float32;  y::Float32;  z::Float32
    vx::Float32; vy::Float32; vz::Float32
    RTHL::Float32
    vx2::Float32; vy2::Float32; vz2::Float32
    overdensity::Float32
end

"""Extended 33-field halo record (ioutshear≥1).

Additional fields beyond the basic 11:
  e_v, p_v           — ellipticity and prolateness at virial radius
  strain(6)          — strain tensor (1,1),(2,2),(3,3),(2,3),(1,3),(1,2)
  d2F                — Laplacian of overdensity at RTHL
  zform              — formation redshift (collapse at RTHL/2)
  grad(3)            — gradient of overdensity at RTHL
  gradf(3)           — gradient of filtered overdensity at RTHL
  Rf                 — filter radius (Mpc/h)
  FcollvRf           — collapse factor at filter scale
  d2FRf              — Laplacian at filter scale
  gradrf(3)          — gradient at filter scale
"""
struct ExtHaloRecord
    # Basic 11 fields (same as HaloRecord)
    x::Float32;  y::Float32;  z::Float32
    vx::Float32; vy::Float32; vz::Float32
    RTHL::Float32
    vx2::Float32; vy2::Float32; vz2::Float32
    overdensity::Float32
    # Extended 22 fields
    e_v::Float32
    p_v::Float32
    strain_11::Float32; strain_22::Float32; strain_33::Float32
    strain_23::Float32; strain_13::Float32; strain_12::Float32
    d2F::Float32
    zform::Float32
    grad_x::Float32; grad_y::Float32; grad_z::Float32
    gradf_x::Float32; gradf_y::Float32; gradf_z::Float32
    Rf::Float32
    FcollvRf::Float32
    d2FRf::Float32
    gradrf_x::Float32; gradrf_y::Float32; gradrf_z::Float32
end

"""Write a halo catalog in .pksc binary format (Fortran STREAM compatible).

Header: Int32 Nhalo, Float32 RTHLmax, Float32 z_out
Data: Nhalo × 11 Float32 records (basic) or Nhalo × 33 Float32 records (extended)
"""
function write_pksc(path::String, halos::Vector{HaloRecord}, RTHLmax::Float32, z_out::Float32)
    open(path, "w") do io
        write(io, Int32(length(halos)))
        write(io, RTHLmax)
        write(io, z_out)
        for h in halos
            write(io, h.x,  h.y,  h.z)
            write(io, h.vx, h.vy, h.vz)
            write(io, h.RTHL)
            write(io, h.vx2, h.vy2, h.vz2)
            write(io, h.overdensity)
        end
    end
end

function write_pksc(path::String, halos::Vector{ExtHaloRecord}, RTHLmax::Float32, z_out::Float32)
    open(path, "w") do io
        write(io, Int32(length(halos)))
        write(io, RTHLmax)
        write(io, z_out)
        for h in halos
            # Basic 11
            write(io, h.x,  h.y,  h.z)
            write(io, h.vx, h.vy, h.vz)
            write(io, h.RTHL)
            write(io, h.vx2, h.vy2, h.vz2)
            write(io, h.overdensity)
            # Extended 22
            write(io, h.e_v, h.p_v)
            write(io, h.strain_11, h.strain_22, h.strain_33)
            write(io, h.strain_23, h.strain_13, h.strain_12)
            write(io, h.d2F)
            write(io, h.zform)
            write(io, h.grad_x, h.grad_y, h.grad_z)
            write(io, h.gradf_x, h.gradf_y, h.gradf_z)
            write(io, h.Rf)
            write(io, h.FcollvRf)
            write(io, h.d2FRf)
            write(io, h.gradrf_x, h.gradrf_y, h.gradrf_z)
        end
    end
end

"""Read a halo catalog from .pksc binary format.

Returns (halos::Vector{HaloRecord}, RTHLmax::Float32, z_out::Float32).
"""
function read_pksc(path::String)
    open(path, "r") do io
        nhalo   = read(io, Int32)
        RTHLmax = read(io, Float32)
        z_out   = read(io, Float32)

        # Determine record size from file: (filesize - 12 header bytes) / nhalo / 4
        data_bytes = filesize(path) - 12
        if nhalo > 0
            nfields = data_bytes ÷ (nhalo * 4)
        else
            nfields = 11
        end

        if nfields == 33
            return _read_ext_records(io, nhalo), RTHLmax, z_out
        else
            return _read_basic_records(io, nhalo), RTHLmax, z_out
        end
    end
end

function _read_basic_records(io::IO, nhalo::Int32)
    halos = Vector{HaloRecord}(undef, nhalo)
    for i in 1:nhalo
        x  = read(io, Float32); y  = read(io, Float32); z  = read(io, Float32)
        vx = read(io, Float32); vy = read(io, Float32); vz = read(io, Float32)
        RTHL = read(io, Float32)
        vx2 = read(io, Float32); vy2 = read(io, Float32); vz2 = read(io, Float32)
        overdensity = read(io, Float32)
        halos[i] = HaloRecord(x, y, z, vx, vy, vz, RTHL, vx2, vy2, vz2, overdensity)
    end
    return halos
end

function _read_ext_records(io::IO, nhalo::Int32)
    halos = Vector{ExtHaloRecord}(undef, nhalo)
    for i in 1:nhalo
        x  = read(io, Float32); y  = read(io, Float32); z  = read(io, Float32)
        vx = read(io, Float32); vy = read(io, Float32); vz = read(io, Float32)
        RTHL = read(io, Float32)
        vx2 = read(io, Float32); vy2 = read(io, Float32); vz2 = read(io, Float32)
        overdensity = read(io, Float32)
        e_v = read(io, Float32); p_v = read(io, Float32)
        s11 = read(io, Float32); s22 = read(io, Float32); s33 = read(io, Float32)
        s23 = read(io, Float32); s13 = read(io, Float32); s12 = read(io, Float32)
        d2F = read(io, Float32)
        zform = read(io, Float32)
        gx = read(io, Float32); gy = read(io, Float32); gz = read(io, Float32)
        gfx = read(io, Float32); gfy = read(io, Float32); gfz = read(io, Float32)
        Rf = read(io, Float32)
        FcollvRf = read(io, Float32)
        d2FRf = read(io, Float32)
        grx = read(io, Float32); gry = read(io, Float32); grz = read(io, Float32)
        halos[i] = ExtHaloRecord(x, y, z, vx, vy, vz, RTHL, vx2, vy2, vz2, overdensity,
                                  e_v, p_v, s11, s22, s33, s23, s13, s12,
                                  d2F, zform, gx, gy, gz, gfx, gfy, gfz,
                                  Rf, FcollvRf, d2FRf, grx, gry, grz)
    end
    return halos
end

end # module Catalog
