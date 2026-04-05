module Parameters

using TOML

# All parameters from hpkvd_params.bin (Fortran ACCESS='STREAM')
struct SimParams
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

        SimParams(
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
function write_params_bin(path::String, p::SimParams)
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
    SimParams(config::Dict{String,Any})

Construct `SimParams` from a parsed TOML configuration.

# TOML sections and keys

```toml
[cosmology]
Om   = 0.315    # total Omega_matter (CDM + baryon)
OB   = 0.049    # Omega_baryon
OL   = 0.685    # Omega_Lambda
h    = 0.674

[grid]
n      = 142    # grid cells per dimension (nlx = nly = nlz)
dL_box = 200.0  # box size [Mpc/h]
nbuff  = 4      # buffer cells

[run]
z_out     = 0.0   # output redshift
ilpt      = 2     # LPT order (1 or 2)
seed      = 42    # RNG seed (not stored in SimParams, used by driver)
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

[output]
format = "pksc"   # "pksc", "hdf5", or "both" (read by driver, not SimParams)
```
"""
function SimParams(config::Dict{String,Any})
    cosmo = get(config, "cosmology", Dict{String,Any}())
    grid  = get(config, "grid",      Dict{String,Any}())
    run   = get(config, "run",       Dict{String,Any}())
    files = get(config, "files",     Dict{String,Any}())

    # Cosmology: Om is total matter, Omx = Om - OB (CDM-only, as Fortran expects)
    Om_total = Float32(get(cosmo, "Om", 0.315))
    OmB      = Float32(get(cosmo, "OB", 0.049))
    Omx      = Om_total - OmB
    Omvac    = Float32(get(cosmo, "OL", 0.685))
    h        = Float32(get(cosmo, "h",  0.674))

    # Grid
    n     = Int32(get(grid, "n",      142))
    dL    = Float32(get(grid, "dL_box", 200.0))
    nbuff = Int32(get(grid, "nbuff",  4))

    # Run parameters
    z_out     = Float32(get(run, "z_out",     0.0))
    z_max     = Float32(get(run, "z_max",     z_out))
    ievol     = Int32(get(run, "ievol",       0))
    ilpt      = Int32(get(run, "ilpt",        2))
    ioutshear = Int32(get(run, "ioutshear",   0))
    wsmooth   = Int32(get(run, "wsmooth",     0))
    rmax2rs   = Float32(get(run, "rmax2rs",   1.0))
    NonGauss  = Int32(get(run, "NonGauss",    0))
    fNL       = Float32(get(run, "fNL",       0.0))

    # Files
    pkfile        = get(files, "pk",         "pk.dat")
    filterfile    = get(files, "filterbank", "filters.dat")
    TabInterpFile = get(files, "homeltab",   "HomelTab.dat")
    fileout       = get(files, "output",     "catalog.pksc")

    # Read collapse table to get its grid dimensions
    ct_nx = Int32(0); ct_ny = Int32(0); ct_nz = Int32(0)
    ct_x1 = Float32(0); ct_x2 = Float32(0)
    ct_y1 = Float32(0); ct_y2 = Float32(0)
    ct_z1 = Float32(0); ct_z2 = Float32(0)
    if isfile(TabInterpFile)
        open(TabInterpFile, "r") do io
            ct_nx = read(io, Int32)
            ct_ny = read(io, Int32)
            ct_nz = read(io, Int32)
            ct_x1 = read(io, Float32); ct_x2 = read(io, Float32)
            ct_y1 = read(io, Float32); ct_y2 = read(io, Float32)
            ct_z1 = read(io, Float32); ct_z2 = read(io, Float32)
        end
    end

    # Observer position (defaults to grid center)
    cenx = Float32(get(grid, "cenx", 0.0))
    ceny = Float32(get(grid, "ceny", 0.0))
    cenz = Float32(get(grid, "cenz", 0.0))

    SimParams(
        Int32(0),           # ireadfield
        ioutshear,
        z_out,
        z_max,              # maximum_redshift
        Int32(1),           # num_redshifts
        Omx, OmB, Omvac, h,
        n, n, n,            # nlx, nly, nlz
        Float32((n - 2*nbuff) * Float64(dL) / n),  # dcore_box
        dL,
        cenx, ceny, cenz,
        nbuff,
        Int32(0),           # next
        ievol,
        Int32(1),           # ivir_strat
        Float32(1.686), Float32(1.686), Float32(1.686),  # fcoll_3, _2, _1
        Float32(1.686),     # dcrit
        Int32(0),           # iforce_strat
        ct_nx, ct_ny, ct_nz,
        ct_x1, ct_x2, ct_y1, ct_y2, ct_z1, ct_z2,
        wsmooth,
        rmax2rs,
        Int32(0),           # ioutfield
        NonGauss, fNL,
        Float32(0), Float32(0), Float32(0),  # A_nG, B_nG, R_nG
        ilpt,
        Int32(0),           # iwant_field_part
        Int32(0),           # largerun
        "",                 # fielddir
        "",                 # densfilein
        "",                 # filein
        pkfile,
        filterfile,
        fileout,
        TabInterpFile
    )
end

end # module Parameters
