module Pipeline

import ..Cosmology: CosmologyParams, Dlinear_tables, Dlinear_ab
import ..PowerSpectrum: load_pk, load_pk_nongaussian
import ..NonGaussian: apply_fnl_correlated!, apply_fnl_uncorrelated!
import ..RandomField: generate_grf, generate_grf_lcg
import ..LPT: displacements_1lpt, displacements_2lpt
import ..Filters: read_filterbank, smooth_field
import ..PeakFind: PeakCandidate, find_peaks
import ..RadialShell: ShellCell, PeakGrid, PeakResult, no_collapse,
    precompute_shells, analyse_peak, fsc_of_z
import ..Parameters: SimParams
import ..Catalog: HaloRecord, ExtHaloRecord, write_pksc, read_pksc
import ..CollapseTable: CollapseTableInterp, read_homeltab

using FFTW
using StaticArrays

"""Convert 1-based flat index to (i,j,k) tuple for column-major n×n×n array."""
function _ipp_to_ijk(ipp::Int, n::Int)
    k = (ipp - 1) ÷ (n * n) + 1
    rem = (ipp - 1) % (n * n)
    j = rem ÷ n + 1
    i = rem % n + 1
    return (i, j, k)
end

"""Compute k²-weighted field (Fortran 'lapld'): ifft(k² × δ_k) / n³.

This matches the Fortran convention where the 'Laplacian' field is
`ctemp = ctemp * ak**2` (positive k², no sign flip), used for the
`d2F` diagnostic in the extended catalog.
"""
function _compute_laplacian(delta_k::Array{Complex{Float32},3}, n::Int, boxsize::Float64)
    dk = Float32(2π / boxsize)
    nk = n ÷ 2 + 1
    nyq = n ÷ 2

    lapd_k = similar(delta_k)
    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = Float32(ix - 1) * dk
        ky = Float32(iy <= nyq + 1 ? iy - 1 : iy - 1 - n) * dk
        kz = Float32(iz <= nyq + 1 ? iz - 1 : iz - 1 - n) * dk
        k2 = kx^2 + ky^2 + kz^2
        lapd_k[ix, iy, iz] = k2 * delta_k[ix, iy, iz]
    end
    lapd = irfft(lapd_k, n) ./ Float32(n)^3
    return Float32.(lapd)
end

"""
    run_tile(sp::SimParams; seed::Integer=42, verbose::Bool=false) -> Vector{HaloRecord}

Single-tile, single-redshift driver that runs the full PeakPatch pipeline:
  Phase 0: Initialize cosmology, growth tables, collapse table, filter bank
  Phase 1: Generate density field and LPT displacements
  Phase 2: Multi-scale peak finding (largest filter first)
  Phase 3: Reset mask for shell analysis
  Phase 4: Radial shell analysis → halo catalog
  Phase 5: Write catalog to disk

Returns the vector of `HaloRecord`.
"""
function run_tile(sp::SimParams; seed::Integer=42, verbose::Bool=false,
                     use_lcg::Bool=false, fortran_compat::Bool=false)
    # ---- Phase 0: Initialization ----
    Om_total = Float64(sp.Omx) + Float64(sp.OmB)   # Omx is CDM-only; total = CDM + baryon
    cosmo = CosmologyParams(Om_total, Float64(sp.OmB), Float64(sp.Omvac),
                            Float64(sp.h), 0.965, 0.808)
    growth_tables = Dlinear_tables(cosmo)

    ct_array, ct_params = read_homeltab(sp.TabInterpFile)
    ct = CollapseTableInterp(ct_array, ct_params)

    filters = read_filterbank(sp.filterfile)
    sort!(filters; by=f -> -f[3])  # sort descending by Rf (largest first)

    n = Int(sp.nlx)
    boxsize = Float64(sp.dL_box)
    alatt = boxsize / n

    z_out = Float64(sp.global_redshift)
    a_out = 1.0 / (1.0 + z_out)
    ZZon = 1.0 + z_out

    fcrit = Float32(fsc_of_z(z_out, growth_tables))

    _, _, D_out = Dlinear_ab(a_out, growth_tables)

    Rfclmax = filters[1][3]
    nhunt = min(Int(sp.nbuff) - 1, floor(Int, Rfclmax * 1.75 / alatt))
    shells = precompute_shells(nhunt)

    Omnr = cosmo.Om
    vTHvir0 = 100.0 * cosmo.h * sqrt(Omnr)

    verbose && @info "Phase 0 done: n=$n, box=$boxsize, z=$z_out, fcrit=$fcrit, $(length(filters)) filters"

    # ---- Phase 1: Field generation ----
    NonGauss = Int(sp.NonGauss)
    pk = load_pk(sp.pkfile)
    delta = if use_lcg
        generate_grf_lcg(n, pk, boxsize, seed; fortran_compat=fortran_compat)
    else
        generate_grf(n, pk, boxsize, seed; fortran_compat=fortran_compat)
    end
    delta_k = rfft(delta)

    # Apply non-Gaussian corrections (modes 1-2)
    if NonGauss in (1, 2)
        _, tf_interp = load_pk_nongaussian(sp.pkfile)
        if NonGauss == 1
            apply_fnl_correlated!(delta, delta_k, pk, tf_interp,
                                  n, boxsize, Float64(sp.fNL))
            delta_k = rfft(delta)
        else  # NonGauss == 2
            apply_fnl_uncorrelated!(delta, delta_k, pk, tf_interp,
                                    n, boxsize, Float64(sp.fNL), seed;
                                    use_lcg=use_lcg)
        end
        verbose && @info "  Applied fNL=$(sp.fNL) (mode $NonGauss)"
    end

    psi_x, psi_y, psi_z = displacements_1lpt(delta_k, n, boxsize)

    if Int(sp.ilpt) >= 2
        psi2_x, psi2_y, psi2_z = displacements_2lpt(delta_k, n, boxsize)
    else
        psi2_x = psi2_y = psi2_z = nothing
    end

    verbose && @info "Phase 1 done: field generated, $(sum(delta .> fcrit)) cells above fcrit"

    # ---- Phase 2: Multi-scale peak finding ----
    mask = zeros(Int8, n, n, n)
    all_peaks = PeakCandidate[]
    peak_Rf = Float64[]
    peak_FcollvRf = Float32[]   # smoothed delta at peak location (per filter)
    peak_d2FRf = Float32[]      # Laplacian of smoothed delta at peak location

    ioutshear = Int(sp.ioutshear)

    for ic in 1:length(filters)
        Rf = filters[ic][3]
        delta_s = smooth_field(delta_k, n, boxsize, Rf, Int(sp.wsmooth);
                               fortran_compat=fortran_compat)

        # Compute smoothed Laplacian at this filter scale if needed
        lapd_s = nothing
        if ioutshear >= 1
            lapd_s = smooth_field(delta_k, n, boxsize, Rf, 3;
                                 fortran_compat=fortran_compat)  # wsmooth=3 → SIGMA_2 (k²-weighted)
        end

        new_peaks = find_peaks(delta_s, mask, 0.0, 0.0, 0.0, alatt, Int(sp.nbuff),
                               fcrit, Rf)
        append!(all_peaks, new_peaks)
        append!(peak_Rf, fill(Rf, length(new_peaks)))

        # Store per-peak filter-scale values for extended catalog
        for pk in new_peaks
            i, j, k = _ipp_to_ijk(pk.ipp, n)
            push!(peak_FcollvRf, delta_s[i, j, k])
            if lapd_s !== nothing
                push!(peak_d2FRf, lapd_s[i, j, k])
            else
                push!(peak_d2FRf, 0.0f0)
            end
        end

        verbose && @info "  Filter $ic: Rf=$Rf, found $(length(new_peaks)) peaks"
    end

    verbose && @info "Phase 2 done: $(length(all_peaks)) total peaks"

    # ---- Phase 3: Reset mask ----
    fill!(mask, 0)

    # ---- Phase 4: Radial shell analysis → catalog ----
    # Compute Laplacian field if extended output is needed
    lapd = nothing
    if Int(sp.ioutshear) >= 1
        lapd = _compute_laplacian(delta_k, n, boxsize)
    end

    pg = PeakGrid(delta, psi_x, psi_y, psi_z,
                  psi2_x, psi2_y, psi2_z,
                  mask, (n, n, n), lapd)

    # Analyse peaks in parallel — results stored per-peak, then filtered
    npeaks = length(all_peaks)
    results = Vector{PeakResult}(undef, npeaks)

    # Pre-compute per-peak ir2min
    ir2min_vec = Vector{Int}(undef, npeaks)
    nbuff_int = Int(sp.nbuff)
    rmax2rs_f = Float64(sp.rmax2rs)
    for idx in 1:npeaks
        Rf = peak_Rf[idx]
        ir2min_vec[idx] = min(floor(Int, (1.75 * Rf / alatt)^2),
                              floor(Int, (40.0 / alatt - 1)^2))
    end

    Threads.@threads for idx in 1:npeaks
        results[idx] = analyse_peak(pg, all_peaks[idx].ipp, alatt,
                                    ir2min_vec[idx], ZZon, peak_Rf[idx],
                                    ct, shells;
                                    nbuff=nbuff_int,
                                    growth_tables=growth_tables,
                                    rmax2rs=rmax2rs_f,
                                    fortran_compat=fortran_compat)
    end

    # Collect halos sequentially (deterministic ordering)
    halos_basic = HaloRecord[]
    halos_ext = ExtHaloRecord[]

    # Pre-compute 2LPT velocity factor
    Om_a_factor = if psi2_x !== nothing
        Om_a = Omnr * a_out^3 / (Omnr * a_out^3 + cosmo.OL)
        -(-3.0/7.0 * Om_a^(-1.0/143) * D_out^2)
    else
        0.0
    end

    for idx in 1:npeaks
        result = results[idx]
        result.RTHL <= 0 && continue

        peak = all_peaks[idx]
        Rf = peak_Rf[idx]
        RTHL_phys = Float32(result.RTHL * alatt)

        # 1LPT velocity: Sbar × D(z_out)
        Sbar_vel = result.Sbar .* D_out

        # 2LPT velocity
        Sbar2_vel = psi2_x !== nothing ? result.Sbar2 .* Om_a_factor : @SVector zeros(3)

        # Virial velocity: v² = vTHvir0² × Fbarx × Srb × RTHL²
        vE2 = vTHvir0^2 * result.Fbarx * result.Srb * Float64(RTHL_phys)^2

        if ioutshear >= 1
            sm = result.strain_mat
            push!(halos_ext, ExtHaloRecord(
                Float32(peak.x), Float32(peak.y), Float32(peak.z),
                Float32(Sbar_vel[1]), Float32(Sbar_vel[2]), Float32(Sbar_vel[3]),
                RTHL_phys,
                Float32(Sbar2_vel[1]), Float32(Sbar2_vel[2]), Float32(Sbar2_vel[3]),
                Float32(result.Fbarx),
                Float32(result.e_v), Float32(result.p_v),
                Float32(sm[1,1]), Float32(sm[2,2]), Float32(sm[3,3]),
                Float32(sm[2,3]), Float32(sm[1,3]), Float32(sm[1,2]),
                Float32(result.d2F),
                Float32(result.zvir_half),
                Float32(result.gradpk[1]), Float32(result.gradpk[2]), Float32(result.gradpk[3]),
                Float32(result.gradpkf[1]), Float32(result.gradpkf[2]), Float32(result.gradpkf[3]),
                Float32(Rf),
                peak_FcollvRf[idx],
                peak_d2FRf[idx],
                Float32(result.gradpkrf[1]), Float32(result.gradpkrf[2]), Float32(result.gradpkrf[3])
            ))
        else
            push!(halos_basic, HaloRecord(
                Float32(peak.x), Float32(peak.y), Float32(peak.z),
                Float32(Sbar_vel[1]), Float32(Sbar_vel[2]), Float32(Sbar_vel[3]),
                RTHL_phys,
                Float32(Sbar2_vel[1]), Float32(Sbar2_vel[2]), Float32(Sbar2_vel[3]),
                Float32(result.Fbarx)
            ))
        end
    end

    verbose && @info "Phase 4 done: $(ioutshear >= 1 ? length(halos_ext) : length(halos_basic)) halos"

    # ---- Phase 5: Write ----
    if ioutshear >= 1
        RTHLmax = isempty(halos_ext) ? Float32(0) : maximum(h.RTHL for h in halos_ext)
        write_pksc(sp.fileout, halos_ext, RTHLmax, Float32(z_out))
    else
        RTHLmax = isempty(halos_basic) ? Float32(0) : maximum(h.RTHL for h in halos_basic)
        write_pksc(sp.fileout, halos_basic, RTHLmax, Float32(z_out))
    end

    halos = ioutshear >= 1 ? halos_ext : halos_basic
    verbose && @info "Phase 5 done: wrote $(length(halos)) halos to $(sp.fileout)"

    return halos
end

end # module Pipeline
