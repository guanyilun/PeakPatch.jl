module AbundanceMatch

using Interpolations
import ..Cosmology: CosmologyParams, growth_factor, build_chi_to_z, chi_to_z
import ..MassFunction: rho_mean, R_of_M, M_of_R, sigma_R, sigma_M,
    dlnsigma_dlnM, tinker_dndlnM, sheth_tormen_dndlnM,
    cumulative_ngtm, precompute_sigma
import ..Catalog: HaloRecord, ExtHaloRecord

"""
    AbundanceTable

Precomputed 2D lookup table mapping (M_TH, z) → M_target.

Fields:
- `Medge`: mass bin edges [M_sun/h] (nMbins+1)
- `zcent`: redshift bin centers (nzbins)
- `M_target`: 2D array M_target[i,j] for mass bin i, redshift bin j
- `interp`: bilinear interpolation function (log10M, z) → log10(M_target)
"""
struct AbundanceTable
    Medge::Vector{Float64}
    zcent::Vector{Float64}
    M_target::Matrix{Float64}
    interp  # interpolation object
end

"""
    build_abundance_table(halos, cosmo, pk;
        nMbins=10000, z_min=0.0, z_max=4.6, nzbins=46,
        Mmin=5e11, Mmax=1e16, hmf=:tinker,
        obs=(0.0, 0.0, 0.0), verbose=false)

Build the abundance matching lookup table M_target(M_TH, z).

For each redshift bin:
1. Count PeakPatch halos → cumulative N_PP(>M)
2. Compute target HMF (Tinker or Sheth-Tormen) → N_target(>M)
3. Invert: for each M_TH with N_PP(>M_TH) > 0, find M_target such that
   N_target(>M_target) = N_PP(>M_TH)

Arguments:
- `halos`: vector of HaloRecord or ExtHaloRecord
- `cosmo`: CosmologyParams
- `pk`: P(k) interpolator from load_pk
- `hmf`: mass function to match (:tinker or :sheth_tormen)
- `obs`: observer position (x,y,z) in Mpc/h for redshift from distance
"""
function build_abundance_table(halos::AbstractVector, cosmo::CosmologyParams, pk;
        nMbins::Int=10000, z_min::Real=0.0, z_max::Real=4.6, nzbins::Int=46,
        Mmin::Real=5e11, Mmax::Real=1e16, hmf::Symbol=:tinker,
        obs::Tuple=(0.0, 0.0, 0.0),
        nsub_integral::Int=10,
        verbose::Bool=false)

    Om = cosmo.Om

    # ---- Mass bins (log-spaced) ----
    Medge = 10.0 .^ range(log10(Mmin), log10(Mmax); length=nMbins+1)
    Mcent = @. 10.0^((log10(Medge[1:end-1]) + log10(Medge[2:end])) / 2.0)
    dlnM = diff(log.(Medge))  # bin widths in ln(M)

    # ---- Redshift bins ----
    zedge = range(z_min, z_max; length=nzbins+1)
    zcent = [(zedge[i] + zedge[i+1]) / 2.0 for i in 1:nzbins]

    # ---- Chi-to-z table for halo redshifts ----
    chi2z = build_chi_to_z(cosmo; z_max=z_max + 1.0)

    # ---- Compute per-halo masses and redshifts ----
    ρ_m = rho_mean(Om)
    halo_M = Vector{Float64}(undef, length(halos))
    halo_z = Vector{Float64}(undef, length(halos))

    for (i, h) in enumerate(halos)
        halo_M[i] = (4π/3.0) * ρ_m * Float64(h.RTHL)^3
        r = sqrt(Float64(h.x - obs[1])^2 + Float64(h.y - obs[2])^2 + Float64(h.z - obs[3])^2)
        halo_z[i] = r > 0 ? chi_to_z(chi2z, r) : 0.0
    end

    # ---- Comoving distance edges for redshift bins ----
    redge = Vector{Float64}(undef, length(zedge))
    for i in eachindex(zedge)
        # Simple integration for chi(z)
        z = zedge[i]
        if z == 0.0
            redge[i] = 0.0
        else
            npts = 1000
            zz = range(0.0, z; length=npts)
            dz = zz[2] - zz[1]
            chi_val = 0.0
            for zp in zz
                E = sqrt(cosmo.Om * (1 + zp)^3 + cosmo.OL)
                chi_val += (2.998e5 / 100.0) * cosmo.h / E * dz
            end
            redge[i] = chi_val
        end
    end

    # ---- Bin PeakPatch halos by mass and redshift → N(>M|z) ----
    N_pp = zeros(nMbins, nzbins)
    for i in eachindex(halo_M)
        M = halo_M[i]
        z = halo_z[i]
        (M < Mmin || M >= Mmax) && continue

        # Find mass bin
        iM = searchsortedfirst(Medge, M) - 1
        (iM < 1 || iM > nMbins) && continue

        # Find redshift bin
        iz = searchsortedfirst(collect(zedge), z) - 1
        (iz < 1 || iz > nzbins) && continue

        N_pp[iM, iz] += 1.0
    end

    # Cumulative N(>M) per z-bin
    NgtM_pp = similar(N_pp)
    for iz in 1:nzbins
        for iM in nMbins:-1:1
            NgtM_pp[iM, iz] = (iM == nMbins) ? N_pp[iM, iz] : NgtM_pp[iM+1, iz] + N_pp[iM, iz]
        end
    end

    # ---- Precompute σ(M) at z=0 ----
    verbose && @info "Precomputing σ(M) at $(nMbins+1) mass bin edges..."
    sigma_edge = precompute_sigma(Medge, pk, Om)
    sigma_cent = [(sigma_edge[i] + sigma_edge[i+1]) / 2.0 for i in 1:nMbins]

    # ---- Compute target N(>M|z) per redshift bin ----
    hmf_func = hmf == :sheth_tormen ? sheth_tormen_dndlnM : tinker_dndlnM

    NgtM_target = zeros(nMbins, nzbins)

    verbose && @info "Computing target HMF ($hmf) at $nzbins redshift bins..."

    for iz in 1:nzbins
        # Integrate N(>M|z) over the redshift shell volume
        dr = (redge[iz+1] - redge[iz]) / nsub_integral

        for jsub in 1:nsub_integral
            r_lo = redge[iz] + (jsub - 1) * dr
            r_hi = redge[iz] + jsub * dr
            r_mid = (r_lo + r_hi) / 2.0
            z_mid = chi_to_z(chi2z, r_mid)
            D_z = growth_factor(z_mid, cosmo)

            # Shell volume (full sky)
            dV = (4π/3.0) * (r_hi^3 - r_lo^3)

            # σ(M, z) = D(z) × σ(M, z=0)
            sigma_cent_z = sigma_cent .* D_z
            sigma_edge_z = sigma_edge .* D_z

            # dn/dlnM at each mass bin center
            dndlnM = Vector{Float64}(undef, nMbins)
            for iM in 1:nMbins
                dlnsig = (log(sigma_edge_z[iM+1]) - log(sigma_edge_z[iM])) /
                         (log(Medge[iM+1]) - log(Medge[iM]))
                dndlnM[iM] = hmf_func(Mcent[iM], sigma_cent_z[iM], dlnsig, z_mid, Om)
            end

            # N(>M) = ∫_M^∞ (dn/dlnM) dlnM, then multiply by shell volume
            ngtm = cumulative_ngtm(dndlnM, dlnM[1])  # approximate with uniform dlnM
            NgtM_target[:, iz] .+= ngtm .* dV
        end
    end

    # ---- Build abundance match table: M_target(M_TH, z) ----
    M_target = zeros(nMbins, nzbins)

    for iz in 1:nzbins
        ngtm_pp = NgtM_pp[:, iz]
        ngtm_tgt = NgtM_target[:, iz]

        # Default: M_target = M_TH (no correction if outside interpolation range)
        M_target[:, iz] .= Medge[1:end-1]

        # Build inverse: M_target(N) from the target HMF
        # N_target(>M) is monotonically decreasing in M, so invert
        # Find indices where N_target > 0
        valid_tgt = findall(ngtm_tgt .> 0)
        length(valid_tgt) < 2 && continue

        # Interpolator: given N_PP(>M_TH), find M_target with same N_target(>M_target)
        # N_target is decreasing → reverse for interpolation
        M_of_N = let M_rev = Medge[valid_tgt[end]:-1:valid_tgt[1]],
                     N_rev = ngtm_tgt[valid_tgt[end]:-1:valid_tgt[1]]
            # N_rev is now increasing (reversed from decreasing)
            # Remove duplicates
            unique_idx = [1]
            for i in 2:length(N_rev)
                N_rev[i] > N_rev[unique_idx[end]] && push!(unique_idx, i)
            end
            length(unique_idx) < 2 && continue
            itp = interpolate((N_rev[unique_idx],), M_rev[unique_idx], Gridded(Linear()))
            extrapolate(itp, Flat())
        end

        # For each mass bin where PP has halos, find the matched mass
        for iM in 1:nMbins
            npp = ngtm_pp[iM]
            npp > 0 || continue
            M_target[iM, iz] = M_of_N(npp)
        end
    end

    # ---- Build 2D interpolator ----
    log10M = log10.(Medge[1:end-1])
    log10M_target = log10.(clamp.(M_target, Mmin, Mmax))

    itp = interpolate((log10M, zcent), log10M_target, Gridded(Linear()))
    interp = extrapolate(itp, Flat())

    verbose && @info "Abundance table built: $(nMbins) mass bins × $(nzbins) z bins"

    return AbundanceTable(Medge, zcent, M_target, interp)
end

"""
    abundance_match(halos, table, Om; obs=(0.0, 0.0, 0.0))

Apply abundance matching to a halo catalog. Returns a new vector of halos
with remapped RTHL (and hence mass).
"""
function abundance_match(halos::Vector{HaloRecord}, table::AbundanceTable,
                         cosmo::CosmologyParams; obs::Tuple=(0.0, 0.0, 0.0))
    Om = cosmo.Om
    ρ_m = rho_mean(Om)
    chi2z = build_chi_to_z(cosmo; z_max=maximum(table.zcent) + 1.0)

    new_halos = Vector{HaloRecord}(undef, length(halos))
    for i in eachindex(halos)
        h = halos[i]

        # Get halo redshift from position
        r = sqrt(Float64(h.x - obs[1])^2 + Float64(h.y - obs[2])^2 + Float64(h.z - obs[3])^2)
        z = r > 0 ? chi_to_z(chi2z, r) : 0.0

        # Current mass
        M_TH = (4π/3.0) * ρ_m * Float64(h.RTHL)^3

        # Look up matched mass
        log10M_new = table.interp(log10(M_TH), z)
        M_new = 10.0^log10M_new

        # Convert back to RTHL
        R_new = R_of_M(M_new, Om)

        new_halos[i] = HaloRecord(
            h.x, h.y, h.z,
            h.vx, h.vy, h.vz,
            Float32(R_new),
            h.vx2, h.vy2, h.vz2,
            h.overdensity
        )
    end
    return new_halos
end

function abundance_match(halos::Vector{ExtHaloRecord}, table::AbundanceTable,
                         cosmo::CosmologyParams; obs::Tuple=(0.0, 0.0, 0.0))
    Om = cosmo.Om
    ρ_m = rho_mean(Om)
    chi2z = build_chi_to_z(cosmo; z_max=maximum(table.zcent) + 1.0)

    new_halos = Vector{ExtHaloRecord}(undef, length(halos))
    for i in eachindex(halos)
        h = halos[i]
        r = sqrt(Float64(h.x - obs[1])^2 + Float64(h.y - obs[2])^2 + Float64(h.z - obs[3])^2)
        z = r > 0 ? chi_to_z(chi2z, r) : 0.0

        M_TH = (4π/3.0) * ρ_m * Float64(h.RTHL)^3
        log10M_new = table.interp(log10(M_TH), z)
        M_new = 10.0^log10M_new
        R_new = R_of_M(M_new, Om)

        new_halos[i] = ExtHaloRecord(
            h.x, h.y, h.z,
            h.vx, h.vy, h.vz,
            Float32(R_new),
            h.vx2, h.vy2, h.vz2,
            h.overdensity,
            h.e_v, h.p_v,
            h.strain_11, h.strain_22, h.strain_33,
            h.strain_23, h.strain_13, h.strain_12,
            h.d2F, h.zform,
            h.grad_x, h.grad_y, h.grad_z,
            h.gradf_x, h.gradf_y, h.gradf_z,
            h.Rf, h.FcollvRf, h.d2FRf,
            h.gradrf_x, h.gradrf_y, h.gradrf_z
        )
    end
    return new_halos
end

"""
    save_abundance_table(path, table)

Save abundance matching table to a text file.
"""
function save_abundance_table(path::String, table::AbundanceTable)
    open(path, "w") do f
        nM = length(table.Medge) - 1
        nz = length(table.zcent)
        println(f, "# Abundance matching table")
        println(f, "# nMbins=$nM nzbins=$nz")
        println(f, "# Medge (nMbins+1 values):")
        println(f, join(table.Medge, " "))
        println(f, "# zcent (nzbins values):")
        println(f, join(table.zcent, " "))
        println(f, "# M_target[iM, iz] (nMbins × nzbins, row-major):")
        for iM in 1:nM
            println(f, join(table.M_target[iM, :], " "))
        end
    end
end

"""
    load_abundance_table(path)

Load a previously saved abundance matching table.
"""
function load_abundance_table(path::String)
    lines = readlines(path)
    idx = 1
    while startswith(lines[idx], "#"); idx += 1; end

    Medge = parse.(Float64, split(lines[idx]))
    idx += 1
    while startswith(lines[idx], "#"); idx += 1; end

    zcent = parse.(Float64, split(lines[idx]))
    idx += 1
    while startswith(lines[idx], "#"); idx += 1; end

    nM = length(Medge) - 1
    nz = length(zcent)
    M_target = zeros(nM, nz)
    for iM in 1:nM
        M_target[iM, :] .= parse.(Float64, split(lines[idx]))
        idx += 1
    end

    log10M = log10.(Medge[1:end-1])
    Mmin = Medge[1]; Mmax = Medge[end]
    log10M_target = log10.(clamp.(M_target, Mmin, Mmax))
    itp = interpolate((log10M, zcent), log10M_target, Gridded(Linear()))
    interp = extrapolate(itp, Flat())

    return AbundanceTable(Medge, zcent, M_target, interp)
end

end # module AbundanceMatch
