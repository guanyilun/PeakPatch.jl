#!/usr/bin/env julia
# Smoking-gun: does the collapse-table threshold fsc_of_z(z) track the analytic
# δ_c/D(z)?  If fsc_of_z rises FASTER than δ_c/D(z) at high z, too few peaks
# collapse → the dN/dz high-z deficit. Also exposes the Frho=8 table cap.
using PeakPatch, Printf
import PeakPatch.Cosmology: growth_factor, delta_c, CosmologyParams
import PeakPatch: fsc_of_z, read_homeltab, CollapseTableInterp

cosmo = CosmologyParams(0.31, 0.049, 0.69, 0.68, 0.965, 0.808)   # Om=0.31 TOTAL (correct)
ct_array, ct_params = read_homeltab(joinpath(@__DIR__, "data", "HomelTab_websky.dat"))
ct = CollapseTableInterp(ct_array, ct_params)

dc0 = 1.686
println("ct x-range (log10 Frho): x1=$(ct_params.X1)  x2=$(ct_params.X2)  -> Frho in [$(10^ct_params.X1), $(10^ct_params.X2)]")
@printf "\n%-6s %-10s %-10s %-12s %-12s %-8s\n" "z" "D(z)" "dc/D(z)" "fsc_of_z" "fsc/(dc/D)" "capped?"
for z in [0.0,0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.4,4.6,5.0]
    D = growth_factor(z, cosmo)
    analytic = dc0 / D
    fsc = fsc_of_z(Float64(z), ct)
    capped = fsc >= 10^ct_params.X2 * 0.999 ? "YES" : ""
    @printf "%-6.2f %-10.4f %-10.4f %-12.4f %-12.4f %-8s\n" z D analytic fsc fsc/analytic capped
end
