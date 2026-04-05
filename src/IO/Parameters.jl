module Parameters

using TOML

"""
    PipelineConfig

Clean configuration struct for the Julia pipeline. Uses keyword arguments with
sensible Planck 2018 defaults. This is the primary config type for `run_tile`,
`run_multitile`, and `run_multitile_mpi`.

Construct from:
- Keyword arguments: `PipelineConfig(n=128, boxsize=200.0, ...)`
- TOML config: `PipelineConfig(TOML.parsefile("config.toml"))`
- Legacy binary: `PipelineConfig(read_params_bin("hpkvd_params.bin"))`
"""
Base.@kwdef struct PipelineConfig
    # Cosmology (Omx is CDM-only; total Om = Omx + OmB)
    Omx::Float64 = 0.266
    OmB::Float64 = 0.049
    Omvac::Float64 = 0.685
    h::Float64 = 0.674
    # Grid
    n::Int = 142
    boxsize::Float64 = 200.0
    nbuff::Int = 4
    # Redshift & lightcone
    z_out::Float64 = 0.0
    z_max::Float64 = 0.0
    ievol::Int = 0
    cenx::Float64 = 0.0
    ceny::Float64 = 0.0
    cenz::Float64 = 0.0
    # Physics
    ilpt::Int = 2
    ioutshear::Int = 0
    wsmooth::Int = 0
    rmax2rs::Float64 = 1.0
    NonGauss::Int = 0
    fNL::Float64 = 0.0
    # Files
    pkfile::String = "pk.dat"
    filterfile::String = "filters.dat"
    tabfile::String = "HomelTab.dat"
    fileout::String = "catalog.pksc"
end

"""
    PipelineConfig(config::Dict{String,Any})

Construct from a parsed TOML configuration dictionary.

# TOML sections and keys

```toml
[cosmology]
Om   = 0.315    # total Omega_matter (CDM + baryon)
OB   = 0.049    # Omega_baryon
OL   = 0.685    # Omega_Lambda
h    = 0.674

[grid]
n      = 142    # grid cells per dimension
boxsize = 200.0 # box size [Mpc/h]
nbuff  = 4      # buffer cells
cenx   = 0.0    # observer x position [Mpc/h]
ceny   = 0.0    # observer y position [Mpc/h]
cenz   = 0.0    # observer z position [Mpc/h]

[run]
z_out     = 0.0   # output redshift
z_max     = 0.0   # maximum redshift (lightcone horizon)
ievol     = 0     # 0: global z, 1: lightcone (per-peak z)
ilpt      = 2     # LPT order (1 or 2)
ioutshear = 0     # 0: basic catalog, ≥1: extended catalog
wsmooth   = 0     # smoothing window (0=Gaussian, 1=tophat)
rmax2rs   = 1.0   # max shell radius / filter radius
NonGauss  = 0     # non-Gaussian mode (0=none, 1=correlated fNL, 2=uncorrelated)
fNL       = 0.0

[files]
pk         = "pk.dat"
filterbank = "filters.dat"
homeltab   = "HomelTab.dat"
output     = "catalog.pksc"
```
"""
function PipelineConfig(config::Dict{String,Any})
    cosmo = get(config, "cosmology", Dict{String,Any}())
    grid  = get(config, "grid",      Dict{String,Any}())
    run   = get(config, "run",       Dict{String,Any}())
    files = get(config, "files",     Dict{String,Any}())

    Om_total = Float64(get(cosmo, "Om", 0.315))
    OmB      = Float64(get(cosmo, "OB", 0.049))
    Omx      = Om_total - OmB

    z_out = Float64(get(run, "z_out", 0.0))

    PipelineConfig(
        Omx      = Omx,
        OmB      = OmB,
        Omvac    = Float64(get(cosmo, "OL", 0.685)),
        h        = Float64(get(cosmo, "h",  0.674)),
        n        = Int(get(grid, "n", 142)),
        boxsize  = Float64(get(grid, "boxsize", get(grid, "dL_box", 200.0))),
        nbuff    = Int(get(grid, "nbuff", 4)),
        z_out    = z_out,
        z_max    = Float64(get(run, "z_max", z_out)),
        ievol    = Int(get(run, "ievol", 0)),
        cenx     = Float64(get(grid, "cenx", 0.0)),
        ceny     = Float64(get(grid, "ceny", 0.0)),
        cenz     = Float64(get(grid, "cenz", 0.0)),
        ilpt      = Int(get(run, "ilpt", 2)),
        ioutshear = Int(get(run, "ioutshear", 0)),
        wsmooth   = Int(get(run, "wsmooth", 0)),
        rmax2rs   = Float64(get(run, "rmax2rs", 1.0)),
        NonGauss  = Int(get(run, "NonGauss", 0)),
        fNL       = Float64(get(run, "fNL", 0.0)),
        pkfile        = get(files, "pk", "pk.dat"),
        filterfile    = get(files, "filterbank", "filters.dat"),
        tabfile       = get(files, "homeltab", "HomelTab.dat"),
        fileout       = get(files, "output", "catalog.pksc"),
    )
end

# ============================================================
# Legacy Fortran binary format (FortranParams)
# ============================================================

# All parameters from hpkvd_params.bin (Fortran ACCESS='STREAM')
struct FortranParams
    # Integer-ish flags (stored as Float32 in binary for uniformity)
    ireadfield::Int32
    ioutshear::Int32
    global_redshift::Float32
    maximum_redshift::Float32
    num_redshifts::Int32
    # Cosmology
    Omx::Float32
    OmB::Float32
    Omvac::Float32
    h::Float32
    # Grid layout
    nlx::Int32
    nly::Int32
    nlz::Int32
    dcore_box::Float32
    dL_box::Float32
    cenx::Float32
    ceny::Float32
    cenz::Float32
    nbuff::Int32
    next::Int32
    ievol::Int32
    # Collapse strategy
    ivir_strat::Int32
    fcoll_3::Float32
    fcoll_2::Float32
    fcoll_1::Float32
    dcrit::Float32
    iforce_strat::Int32
    # Interpolation table
    TabInterpNx::Int32
    TabInterpNy::Int32
    TabInterpNz::Int32
    TabInterpX1::Float32
    TabInterpX2::Float32
    TabInterpY1::Float32
    TabInterpY2::Float32
    TabInterpZ1::Float32
    TabInterpZ2::Float32
    # Smoothing
    wsmooth::Int32
    rmax2rs::Float32
    ioutfield::Int32
    # Non-Gaussian
    NonGauss::Int32
    fNL::Float32
    A_nG::Float32
    B_nG::Float32
    R_nG::Float32
    # LPT and field
    ilpt::Int32
    iwant_field_part::Int32
    largerun::Int32
    # Strings
    fielddir::String
    densfilein::String
    filein::String
    pkfile::String
    filterfile::String
    fileout::String
    TabInterpFile::String
end

"""Read 46 mixed Int32/Float32 values + 7 length-prefixed strings from Fortran binary.

The binary is written by peak-patch.py with np.int32 for integer fields and
np.float32 for real fields.  Each field occupies exactly 4 bytes.  The field
types must match the Python writer exactly — reading Int32 as Float32 (or
vice versa) reinterprets the bit pattern and gives wrong values.
"""
function read_params_bin(path::String)
    open(path, "r") do io
        # Part 1: 46 fields, each 4 bytes, in the exact order and type
        # written by peak-patch.py (lines 299-321)
        ireadfield       = read(io, Int32)
        ioutshear        = read(io, Int32)
        global_redshift  = read(io, Float32)
        maximum_redshift = read(io, Float32)
        num_redshifts    = read(io, Int32)
        Omx              = read(io, Float32)
        OmB              = read(io, Float32)
        Omvac            = read(io, Float32)
        h                = read(io, Float32)
        nlx              = read(io, Int32)
        nly              = read(io, Int32)
        nlz              = read(io, Int32)
        dcore_box        = read(io, Float32)
        dL_box           = read(io, Float32)
        cenx             = read(io, Float32)
        ceny             = read(io, Float32)
        cenz             = read(io, Float32)
        nbuff            = read(io, Int32)
        next             = read(io, Int32)
        ievol            = read(io, Int32)
        ivir_strat       = read(io, Int32)
        fcoll_3          = read(io, Float32)
        fcoll_2          = read(io, Float32)
        fcoll_1          = read(io, Float32)
        dcrit            = read(io, Float32)
        iforce_strat     = read(io, Int32)
        TabInterpNx      = read(io, Int32)
        TabInterpNy      = read(io, Int32)
        TabInterpNz      = read(io, Int32)
        TabInterpX1      = read(io, Float32)
        TabInterpX2      = read(io, Float32)
        TabInterpY1      = read(io, Float32)
        TabInterpY2      = read(io, Float32)
        TabInterpZ1      = read(io, Float32)
        TabInterpZ2      = read(io, Float32)
        wsmooth          = read(io, Int32)
        rmax2rs          = read(io, Float32)
        ioutfield        = read(io, Int32)
        NonGauss         = read(io, Int32)
        fNL              = read(io, Float32)
        A_nG             = read(io, Float32)
        B_nG             = read(io, Float32)
        R_nG             = read(io, Float32)
        ilpt             = read(io, Int32)
        iwant_field_part = read(io, Int32)
        largerun         = read(io, Int32)

        # Part 2: 7 length-prefixed strings (Int32 length + chars)
        fielddir      = _read_string(io)
        densfilein    = _read_string(io)
        filein        = _read_string(io)
        pkfile        = _read_string(io)
        filterfile    = _read_string(io)
        fileout       = _read_string(io)
        TabInterpFile = _read_string(io)

        FortranParams(
            ireadfield, ioutshear, global_redshift, maximum_redshift, num_redshifts,
            Omx, OmB, Omvac, h, nlx, nly, nlz,
            dcore_box, dL_box, cenx, ceny, cenz, nbuff, next, ievol,
            ivir_strat, fcoll_3, fcoll_2, fcoll_1, dcrit, iforce_strat,
            TabInterpNx, TabInterpNy, TabInterpNz,
            TabInterpX1, TabInterpX2, TabInterpY1, TabInterpY2, TabInterpZ1, TabInterpZ2,
            wsmooth, rmax2rs, ioutfield,
            NonGauss, fNL, A_nG, B_nG, R_nG,
            ilpt, iwant_field_part, largerun,
            fielddir, densfilein, filein, pkfile, filterfile, fileout, TabInterpFile
        )
    end
end

function _read_string(io::IO)
    slen = read(io, Int32)
    if slen == 0
        return ""
    end
    data = read(io, slen)
    return String(data)
end

function _write_string(io::IO, s::String)
    write(io, Int32(length(s)))
    write(io, s)
end

"""Write parameters in Fortran-compatible binary format (matching peak-patch.py layout)."""
function write_params_bin(path::String, p::FortranParams)
    open(path, "w") do io
        # 46 fields with correct Int32/Float32 types matching the Python writer
        write(io, Int32(p.ireadfield))
        write(io, Int32(p.ioutshear))
        write(io, Float32(p.global_redshift))
        write(io, Float32(p.maximum_redshift))
        write(io, Int32(p.num_redshifts))
        write(io, Float32(p.Omx))
        write(io, Float32(p.OmB))
        write(io, Float32(p.Omvac))
        write(io, Float32(p.h))
        write(io, Int32(p.nlx))
        write(io, Int32(p.nly))
        write(io, Int32(p.nlz))
        write(io, Float32(p.dcore_box))
        write(io, Float32(p.dL_box))
        write(io, Float32(p.cenx))
        write(io, Float32(p.ceny))
        write(io, Float32(p.cenz))
        write(io, Int32(p.nbuff))
        write(io, Int32(p.next))
        write(io, Int32(p.ievol))
        write(io, Int32(p.ivir_strat))
        write(io, Float32(p.fcoll_3))
        write(io, Float32(p.fcoll_2))
        write(io, Float32(p.fcoll_1))
        write(io, Float32(p.dcrit))
        write(io, Int32(p.iforce_strat))
        write(io, Int32(p.TabInterpNx))
        write(io, Int32(p.TabInterpNy))
        write(io, Int32(p.TabInterpNz))
        write(io, Float32(p.TabInterpX1))
        write(io, Float32(p.TabInterpX2))
        write(io, Float32(p.TabInterpY1))
        write(io, Float32(p.TabInterpY2))
        write(io, Float32(p.TabInterpZ1))
        write(io, Float32(p.TabInterpZ2))
        write(io, Int32(p.wsmooth))
        write(io, Float32(p.rmax2rs))
        write(io, Int32(p.ioutfield))
        write(io, Int32(p.NonGauss))
        write(io, Float32(p.fNL))
        write(io, Float32(p.A_nG))
        write(io, Float32(p.B_nG))
        write(io, Float32(p.R_nG))
        write(io, Int32(p.ilpt))
        write(io, Int32(p.iwant_field_part))
        write(io, Int32(p.largerun))
        # 7 length-prefixed strings
        _write_string(io, p.fielddir)
        _write_string(io, p.densfilein)
        _write_string(io, p.filein)
        _write_string(io, p.pkfile)
        _write_string(io, p.filterfile)
        _write_string(io, p.fileout)
        _write_string(io, p.TabInterpFile)
    end
end

"""
    PipelineConfig(sp::FortranParams)

Convert a legacy `FortranParams` (from Fortran binary) to `PipelineConfig`.
"""
function PipelineConfig(sp::FortranParams)
    PipelineConfig(
        Omx       = Float64(sp.Omx),
        OmB       = Float64(sp.OmB),
        Omvac     = Float64(sp.Omvac),
        h         = Float64(sp.h),
        n         = Int(sp.nlx),
        boxsize   = Float64(sp.dL_box),
        nbuff     = Int(sp.nbuff),
        z_out     = Float64(sp.global_redshift),
        z_max     = Float64(sp.maximum_redshift),
        ievol     = Int(sp.ievol),
        cenx      = Float64(sp.cenx),
        ceny      = Float64(sp.ceny),
        cenz      = Float64(sp.cenz),
        ilpt      = Int(sp.ilpt),
        ioutshear = Int(sp.ioutshear),
        wsmooth   = Int(sp.wsmooth),
        rmax2rs   = Float64(sp.rmax2rs),
        NonGauss  = Int(sp.NonGauss),
        fNL       = Float64(sp.fNL),
        pkfile    = sp.pkfile,
        filterfile = sp.filterfile,
        tabfile   = sp.TabInterpFile,
        fileout   = sp.fileout,
    )
end

end # module Parameters
