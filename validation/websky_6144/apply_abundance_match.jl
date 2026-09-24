#!/usr/bin/env julia
# Apply abundance matching to the corrected production octant catalog and measure
# how much it closes the gap to Websky (N(>1.69e12): Websky ≈5.09e7/octant; the old "110M" conflated total vs >1.69e12).
# AM remaps raw peak-patch top-hat masses to Tinker M200 in z-bins (Websky's step).
#
# NOTE: build_abundance_table needs the PHYSICAL P(k) (σ8=0.81) for σ(M)/Tinker —
# use the RAW (un-normalized) pk file, NOT the ÷(2π)³ field-gen one.
#
# Usage: julia --project=. -t 8 validation/websky_6144/apply_abundance_match.jl

using PeakPatch, Printf
const D = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"

# ARGS: [1]=input catalog, [2]=output catalog, [3]=observer (Mpc/h): either a single
#        value applied to all axes ("-2618") or comma-separated per-axis
#        ("2618,-2618,2618" — REQUIRED for mixed-sign octants: the observer sets each
#        halo's chi→z for the AM z-binning), [4]=z_max (default 4.6).
# Defaults = old coarse pkfix octant (obs -3850). For finecell octant pass obs -2618.
cat     = length(ARGS) >= 1 ? ARGS[1] : joinpath(D, "catalog_websky_6144_oct000_pkfix.pksc")
out     = length(ARGS) >= 2 ? ARGS[2] : joinpath(D, "catalog_websky_6144_oct000_pkfix_AM.pksc")
obs_v   = length(ARGS) >= 3 ? parse.(Float64, split(ARGS[3], ",")) : [-3850.0]
length(obs_v) in (1, 3) || error("observer must be 1 or 3 comma-separated values")
z_max_am = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 4.6
@info "reading catalog..." cat
halos, RTHLmax, z_out = read_pksc(cat)
@info "read" n=length(halos)

# CosmologyParams first arg is Om_TOTAL (=Omx+OmB). Websky Om=0.31 total, OmB=0.049.
# (Was previously 0.31+0.049=0.359 — wrong; double-counted baryons, corrupting rho_mean/D(z)/volumes.)
cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.81)
pk = PeakPatch.PowerSpectrum.load_pk(joinpath(@__DIR__, "data", "pk_websky_RAW_unnormalized.dat"))
obs = length(obs_v) == 3 ? (obs_v[1], obs_v[2], obs_v[3]) : (obs_v[1], obs_v[1], obs_v[1])
rho_m = 2.775e11 * 0.31

NgtM(hs, M0) = count(h -> (4π/3)*rho_m*Float64(h.RTHL)^3 > M0, hs)
function report(label, hs)
    @printf("%-18s  N(>1.7e12)=%.3e  N(>1e13)=%.3e  N(>1e14)=%.3e\n",
            label, NgtM(hs,1.7e12), NgtM(hs,1e13), NgtM(hs,1e14))
end

report("RAW (pre-AM)", halos)

@info "building abundance table (Tinker)..." z_max_am obs
table = build_abundance_table(halos, cosmo, pk; hmf=:tinker, z_max=z_max_am, obs=obs,
                              fsky=1/8, verbose=true)   # single octant = 1/8 sky
@info "applying abundance match..."
halos_am = abundance_match(halos, table, cosmo; obs=obs)

report("AFTER AM (Tinker)", halos_am)
@printf("\nReference: Websky measured N(>1.69e12) ≈ 5.09e7/octant (finecell+AM validation job 4033130: 5.08e7; see LOW_MASS_COMPLETENESS_2026-06-15.md)\n")

# write the AM'd catalog
Rmax = isempty(halos_am) ? 0f0 : maximum(h.RTHL for h in halos_am)
write_pksc(out, halos_am, Float32(Rmax), Float32(z_out))
@info "wrote AM catalog" out
