module CollapseTable

import ..Cosmology: CosmologyParams
import ..EllipsoidalCollapse: EllipsoidParams, evolve_ellipse_full

# Table grid parameters
struct CollapseTableParams
    Nx::Int; Ny::Int; Nz::Int
    X1::Float64; X2::Float64   # log10(Frho) range [log10(1.5), log10(8.0)] after conversion
    Y1::Float64; Y2::Float64   # e_v range [0.0, 0.5]
    Z1::Float64; Z2::Float64   # p_v/e_v range [-1+1e-4, 1-1e-4]
end

function CollapseTableParams(;
    Nx::Int = 50, Ny::Int = 20, Nz::Int = 20,
    X1::Float64 = log10(1.5), X2::Float64 = log10(8.0),
    Y1::Float64 = 0.0, Y2::Float64 = 0.5,
    Z1::Float64 = -1.0 + 1e-4, Z2::Float64 = 1.0 - 1e-4)
    CollapseTableParams(Nx, Ny, Nz, X1, X2, Y1, Y2, Z1, Z2)
end

# Wrapper that matches Fortran's evolve_ellipse_function
function evolve_ellipse_function(logF, e, poe, ep::EllipsoidParams)
    Frho = 10.0^logF
    p_v = poe * e
    zvir, _, _ = evolve_ellipse_full(Frho, e, p_v, ep; idynax = 1)
    return Float32(zvir)
end

# Generate the collapse table
function make_table(ep::EllipsoidParams, tp::CollapseTableParams; verbose::Bool = true)
    table = Array{Float32}(undef, tp.Nx, tp.Ny, tp.Nz)
    dX = (tp.X2 - tp.X1) / (tp.Nx - 1)
    dY = (tp.Y2 - tp.Y1) / (tp.Ny - 1)
    dZ = (tp.Z2 - tp.Z1) / (tp.Nz - 1)

    for k in 1:tp.Nz
        if verbose
            println("  k=$k/$(tp.Nz)")
        end
        z = tp.Z1 + (k - 1) * dZ
        for j in 1:tp.Ny
            y = tp.Y1 + (j - 1) * dY
            for i in 1:tp.Nx
                x = tp.X1 + (i - 1) * dX
                table[i, j, k] = evolve_ellipse_function(x, y, z, ep)
            end
        end
    end

    return table
end

# Threaded table generation
function make_table_threaded(ep::EllipsoidParams, tp::CollapseTableParams; verbose::Bool = true)
    table = Array{Float32}(undef, tp.Nx, tp.Ny, tp.Nz)
    dX = (tp.X2 - tp.X1) / (tp.Nx - 1)
    dY = (tp.Y2 - tp.Y1) / (tp.Ny - 1)
    dZ = (tp.Z2 - tp.Z1) / (tp.Nz - 1)

    for k in 1:tp.Nz
        if verbose
            println("  k=$k/$(tp.Nz)")
        end
        z = tp.Z1 + (k - 1) * dZ
        Threads.@threads for idx in 1:(tp.Nx * tp.Ny)
            j = ((idx - 1) ÷ tp.Nx) + 1
            i = ((idx - 1) % tp.Nx) + 1
            x = tp.X1 + (i - 1) * dX
            y = tp.Y1 + (j - 1) * dY
            table[i, j, k] = evolve_ellipse_function(x, y, z, ep)
        end
    end

    return table
end

# -------------------------------------------------------------------------
# Binary I/O — Fortran stream-compatible format
# Header: 3×Int32 + 6×Float32 = 36 bytes
# Data: Array{Float32}(Nx, Ny, Nz) in column-major order (matches Fortran)
# -------------------------------------------------------------------------
function write_homeltab(path::String, table::Array{Float32,3}, tp::CollapseTableParams)
    open(path, "w") do f
        # Header
        write(f, Int32(tp.Nx))
        write(f, Int32(tp.Ny))
        write(f, Int32(tp.Nz))
        write(f, Float32(tp.X1))
        write(f, Float32(tp.X2))
        write(f, Float32(tp.Y1))
        write(f, Float32(tp.Y2))
        write(f, Float32(tp.Z1))
        write(f, Float32(tp.Z2))
        # Data
        write(f, table)
    end
end

function read_homeltab(path::String)
    open(path, "r") do f
        # Header
        Nx = read(f, Int32)
        Ny = read(f, Int32)
        Nz = read(f, Int32)
        X1 = Float64(read(f, Float32))
        X2 = Float64(read(f, Float32))
        Y1 = Float64(read(f, Float32))
        Y2 = Float64(read(f, Float32))
        Z1 = Float64(read(f, Float32))
        Z2 = Float64(read(f, Float32))

        tp = CollapseTableParams(Nx=Int(Nx), Ny=Int(Ny), Nz=Int(Nz),
            X1=X1, X2=X2, Y1=Y1, Y2=Y2, Z1=Z1, Z2=Z2)

        table = Array{Float32,3}(undef, Int(Nx), Int(Ny), Int(Nz))
        read!(f, table)
        return (table, tp)
    end
end

# -------------------------------------------------------------------------
# Trilinear interpolation (matches TabInterpInterpolate from TabInterp.f90)
# -------------------------------------------------------------------------
struct CollapseTableInterp
    array::Array{Float32,3}
    X1::Float64; X2::Float64
    Y1::Float64; Y2::Float64
    Z1::Float64; Z2::Float64
    dX::Float64; dY::Float64; dZ::Float64
    Nx::Int; Ny::Int; Nz::Int
    out_val::Float32  # value to return for out-of-bounds queries
end

function CollapseTableInterp(table::Array{Float32,3}, tp::CollapseTableParams)
    dX = (tp.X2 - tp.X1) / (tp.Nx - 1)
    dY = (tp.Y2 - tp.Y1) / (tp.Ny - 1)
    dZ = (tp.Z2 - tp.Z1) / (tp.Nz - 1)
    CollapseTableInterp(table, tp.X1, tp.X2, tp.Y1, tp.Y2, tp.Z1, tp.Z2,
        dX, dY, dZ, tp.Nx, tp.Ny, tp.Nz, Float32(-1.0))
end

function interpolate(ct::CollapseTableInterp, x, y, z)
    if x > ct.X2 || x < ct.X1 || y > ct.Y2 || y < ct.Y1 || z > ct.Z2 || z < ct.Z1
        return ct.out_val
    end

    i1 = floor(Int, (x - ct.X1) / ct.dX) + 1
    j1 = floor(Int, (y - ct.Y1) / ct.dY) + 1
    k1 = floor(Int, (z - ct.Z1) / ct.dZ) + 1

    # Clamp to valid range
    i1 = min(i1, ct.Nx - 1)
    j1 = min(j1, ct.Ny - 1)
    k1 = min(k1, ct.Nz - 1)
    i1 = max(i1, 1)
    j1 = max(j1, 1)
    k1 = max(k1, 1)

    i2 = i1 + 1
    j2 = j1 + 1
    k2 = k1 + 1

    fx = (x - (i1 - 1) * ct.dX - ct.X1) / ct.dX
    fy = (y - (j1 - 1) * ct.dY - ct.Y1) / ct.dY
    fz = (z - (k1 - 1) * ct.dZ - ct.Z1) / ct.dZ

    return Float64(
        ct.array[i1, j1, k1] * (1 - fx) * (1 - fy) * (1 - fz) +
        ct.array[i2, j1, k1] * fx * (1 - fy) * (1 - fz) +
        ct.array[i1, j2, k1] * (1 - fx) * fy * (1 - fz) +
        ct.array[i2, j2, k1] * fx * fy * (1 - fz) +
        ct.array[i1, j1, k2] * (1 - fx) * (1 - fy) * fz +
        ct.array[i2, j1, k2] * fx * (1 - fy) * fz +
        ct.array[i1, j2, k2] * (1 - fx) * fy * fz +
        ct.array[i2, j2, k2] * fx * fy * fz
    )
end

end # module CollapseTable
