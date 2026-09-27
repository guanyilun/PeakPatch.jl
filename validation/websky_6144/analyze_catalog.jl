using PeakPatch

halos, RTHLmax, z_out = read_pksc("/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000.pksc")
println("Total halos: ", length(halos))
println("RTHLmax: ", RTHLmax)
println("z_out: ", z_out)

rthls = [h.RTHL for h in halos]
println("\nRTHL stats (Mpc/h):")
println("  min: ", minimum(rthls))
println("  max: ", maximum(rthls))
sorted_rthls = sort(rthls)
println("  median: ", sorted_rthls[length(rthls)÷2])
println("  mean: ", sum(Float64.(rthls))/length(rthls))

Om = 0.31; hh = 0.68
rho_mean = 2.775e11 * Om * hh^2
masses = [(4pi/3) * rho_mean * Float64(r)^3 for r in rthls]
println("\nMass stats (M_sun/h):")
println("  min: ", minimum(masses), " = ", minimum(masses)/1e12, " x 10^12 M_sun/h")
println("  max: ", maximum(masses), " = ", maximum(masses)/1e12, " x 10^12 M_sun/h")
sorted_masses = sort(masses)
println("  median: ", sorted_masses[length(masses)÷2]/1e12, " x 10^12 M_sun/h")
println("  mean: ", (sum(masses)/length(masses))/1e12, " x 10^12 M_sun/h")

log_masses = log10.(masses)
edges = 12.0:0.5:16.5
println("\nMass function (log10 M [M_sun/h]):")
for i in 1:length(edges)-1
    n = count(m -> edges[i] <= m < edges[i+1], log_masses)
    n > 0 && println("  $(edges[i])-$(edges[i+1]): $n halos")
end

obs = (-3850.0, -3850.0, -3850.0)
dists = [sqrt((Float64(ha.x) - obs[1])^2 + (Float64(ha.y) - obs[2])^2 + (Float64(ha.z) - obs[3])^2) for ha in halos]
println("\nDistance from observer (Mpc/h):")
println("  min: ", minimum(dists))
println("  max: ", maximum(dists))
sorted_dists = sort(dists)
println("  median: ", sorted_dists[length(dists)÷2])

edges_d = 0:500:8000
println("\nHalos by distance shell:")
for i in 1:length(edges_d)-1
    n = count(d -> edges_d[i] <= d < edges_d[i+1], dists)
    n > 0 && println("  $(edges_d[i])-$(edges_d[i+1]) Mpc/h: $n halos")
end

println("\nPosition ranges (Mpc/h):")
println("  x: ", minimum(ha.x for ha in halos), " to ", maximum(ha.x for ha in halos))
println("  y: ", minimum(ha.y for ha in halos), " to ", maximum(ha.y for ha in halos))
println("  z: ", minimum(ha.z for ha in halos), " to ", maximum(ha.z for ha in halos))

ovds = [ha.overdensity for ha in halos]
println("\nOverdensity stats:")
println("  min: ", minimum(ovds))
println("  max: ", maximum(ovds))
println("  mean: ", sum(Float64.(ovds))/length(ovds))

if halos[1] isa PeakPatch.Catalog.ExtHaloRecord
    println("\nExtended catalog (33 fields)")
    rfs = [ha.Rf for ha in halos]
    println("  Rf range: ", minimum(rfs), " to ", maximum(rfs), " Mpc/h")
    zforms = [ha.zform for ha in halos]
    valid_zf = filter(z -> z >= 0, zforms)
    if !isempty(valid_zf)
        println("  zform range (valid): ", minimum(valid_zf), " to ", maximum(valid_zf))
    end
    println("  zform invalid count: ", count(z -> z < 0, zforms))
else
    println("\nBasic catalog (11 fields)")
end
