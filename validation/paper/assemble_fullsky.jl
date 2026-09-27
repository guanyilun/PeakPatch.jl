#!/usr/bin/env julia
# Assemble full-sky frozen-campaign maps from the 8 octants (paper_comparison_plan A1).
# Each octant map is a full-sky HEALPix array that is (almost) zero outside its octant;
# halo/field terms that spill across an octant boundary belong to that sky, so the full
# sky is the plain sum over octants.
#
#   CAMPAIGN=prod|v2 julia --project=validation -t 8 assemble_fullsky.jl [product ...]
# (prod = frozen 2026-07 campaign with AMv2 halos; v2 = periodic-core rerun 2026-09)
#
# products (default all): kappa ksz tsz isw cib
#   kappa_lt4.5  = Σ (field κ + κ_halo_comp)         [construction B, z<4.5]
#   kappa_field, kappa_halo_comp, kappa_halo_plain    [components]
#   ksz_total_uK = Σ (T_CMB·field kSZ + ksz_halo_Wc)  [μK, Wc construction]
#   ksz_field_uK, ksz_halo_Wc_uK, ksz_halo_W_uK       [components, μK]
#   tsz_y        = Σ y
#   isw_uK       = Σ T_CMB·ISW                         [μK]
#   cib_nuXXXX_wcut (MJy/sr) for each frequency present
using Healpix, Printf

const D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"
const CAMP = get(ENV, "CAMPAIGN", "prod")
CAMP in ("prod", "v2") || error("CAMPAIGN must be prod or v2")
const FM = joinpath(D, "fieldmaps_$(CAMP)"); const HM = joinpath(D, "halomaps_$(CAMP)")
const CM = joinpath(D, "cibmaps_$(CAMP)"); const OUT = joinpath(D, "fullsky_$(CAMP)")
const HTAG = CAMP == "prod" ? (oct -> "prod_oct$(oct)_AMv2") : (oct -> "v2_oct$(oct)_AM")
const OCTS = ["000", "001", "010", "011", "100", "101", "110", "111"]
const TCMB_UK = 2.7255e6
const NSIDE = 4096

readmap(p) = (isfile(p) || error("missing $p"); Healpix.readMapFromFITS(p, 1, Float32))

# sum Σ_oct Σ_(pattern,scale) scale·map(pattern(oct)) → write
function assemble(name, terms)
    acc = zeros(Float64, 12NSIDE^2)
    for oct in OCTS, (pat, scale) in terms
        m = readmap(pat(oct))
        m.resolution.nside == NSIDE || error("nside mismatch in $(pat(oct))")
        acc .+= scale .* m.pixels
    end
    hm = HealpixMap{Float32,RingOrder}(NSIDE); hm.pixels .= Float32.(acc)
    out = joinpath(OUT, "$(name)_$(CAMP)_fullsky_nside$(NSIDE).fits")
    Healpix.saveToFITS(hm, "!" * out, typechar="E")
    @printf("%-22s mean %+.4e  rms %.4e  zero-pix %d -> %s\n", name, sum(acc) / length(acc),
            sqrt(sum(abs2, acc) / length(acc)), count(iszero, acc), out)
end

fld(k) = oct -> joinpath(FM, "$(k)_$(CAMP)_oct$(oct)_nside$(NSIDE).fits")
hal(k) = oct -> joinpath(HM, "$(k)_$(HTAG(oct))_nside$(NSIDE).fits")

prods = isempty(ARGS) ? ["kappa", "ksz", "tsz", "isw", "cib"] : ARGS
mkpath(OUT)
if "kappa" in prods
    assemble("kappa_lt4.5", [(fld("kappa"), 1.0), (hal("kappa_halo_comp"), 1.0)])
    assemble("kappa_field", [(fld("kappa"), 1.0)])
    assemble("kappa_halo_comp", [(hal("kappa_halo_comp"), 1.0)])
    assemble("kappa_halo_plain", [(hal("kappa_halo_plain"), 1.0)])
end
if "ksz" in prods
    assemble("ksz_total_uK", [(fld("ksz"), TCMB_UK), (hal("ksz_halo_Wc"), 1.0)])
    assemble("ksz_field_uK", [(fld("ksz"), TCMB_UK)])
    assemble("ksz_halo_Wc_uK", [(hal("ksz_halo_Wc"), 1.0)])
    assemble("ksz_halo_W_uK", [(hal("ksz_halo_W"), 1.0)])
end
"tsz" in prods && assemble("tsz_y", [(hal("tsz_y"), 1.0)])
"isw" in prods && assemble("isw_uK", [(fld("isw"), TCMB_UK)])
if "cib" in prods
    for ν in (100, 143, 217, 353, 545, 857)
        tag = "cib_nu$(lpad(ν, 4, '0'))_wcut"
        assemble(tag, [(oct -> joinpath(CM, "$(tag)_$(HTAG(oct))_nside$(NSIDE).fits"), 1.0)])
    end
end
