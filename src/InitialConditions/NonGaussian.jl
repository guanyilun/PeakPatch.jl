module NonGaussian

using FFTW

import ..PowerSpectrum: load_pk_nongaussian
import ..RandomField: generate_grf, generate_grf_lcg

"""
    apply_fnl_correlated!(delta::Array{Float32,3}, delta_k, pk_interp, tf_interp,
                          n::Int, boxsize::Float64, fNL::Float64)

Local-type non-Gaussianity (mode 1, correlated):
    ζ_NL(x) = ζ_G(x) + fNL × [ζ_G(x)² - ⟨ζ_G²⟩]
    δ_NL(k) = T(k)² × ζ_NL(k)

Modifies `delta` in place. `delta_k` is the rfft of the Gaussian delta field.
`pk_interp` returns sqrt(P(k)) and `tf_interp` returns sqrt(k² T(k)).
"""
function apply_fnl_correlated!(delta::Array{Float32,3}, delta_k,
                                pk_interp, tf_interp,
                                n::Int, boxsize::Float64, fNL::Float64)
    dk = 2π / boxsize
    nk = n ÷ 2 + 1
    nyq = n ÷ 2
    n3 = Float64(n)^3

    # Step 1: Convert delta_k → zeta_k by dividing by T(k)²
    # delta_k = amplitude × pk(k), so zeta_k = delta_k × pk/tf² / pk = delta_k / tf²
    # But actually: delta = noise * pk, and zeta = noise * pk/tf²
    # So zeta_k = delta_k / tf² ... but tf = sqrt(k²T(k)), so tf² = k²T(k)
    # This means we need to divide each k-mode by tf(k)²
    zeta_k = similar(delta_k)
    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = Float64(ix - 1) * dk
        ky = Float64(iy <= nyq + 1 ? iy - 1 : iy - 1 - n) * dk
        kz = Float64(iz <= nyq + 1 ? iz - 1 : iz - 1 - n) * dk
        k = sqrt(kx^2 + ky^2 + kz^2)
        if k == 0.0
            zeta_k[ix, iy, iz] = 0.0f0
        else
            tf2 = tf_interp(k)^2
            if tf2 > 0.0
                zeta_k[ix, iy, iz] = delta_k[ix, iy, iz] / Float32(tf2)
            else
                zeta_k[ix, iy, iz] = 0.0f0
            end
        end
    end

    # Step 2: IFFT to get zeta_G(x)
    zeta_G = irfft(zeta_k, n) ./ Float32(n3)

    # Step 3: Apply fNL: zeta_NL = zeta_G + fNL * (zeta_G² - ⟨zeta_G²⟩)
    avg_zeta2 = sum(zeta_G .^ 2) / n3
    zeta_NL = zeta_G .+ Float32(fNL) .* (zeta_G .^ 2 .- Float32(avg_zeta2))

    # Step 4: FFT back to k-space and multiply by T(k)² to get delta_NL
    zeta_NL_k = rfft(Float32.(zeta_NL))
    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = Float64(ix - 1) * dk
        ky = Float64(iy <= nyq + 1 ? iy - 1 : iy - 1 - n) * dk
        kz = Float64(iz <= nyq + 1 ? iz - 1 : iz - 1 - n) * dk
        k = sqrt(kx^2 + ky^2 + kz^2)
        if k > 0.0
            tf2 = tf_interp(k)^2
            delta_k[ix, iy, iz] = zeta_NL_k[ix, iy, iz] * Float32(tf2)
        else
            delta_k[ix, iy, iz] = 0.0f0
        end
    end

    # Step 5: IFFT to get delta_NL(x)
    delta .= Float32.(irfft(delta_k, n) ./ Float32(n3))
    return delta
end

"""
    apply_fnl_uncorrelated!(delta::Array{Float32,3}, delta_k,
                            pk_interp, tf_interp,
                            n::Int, boxsize::Float64, fNL::Float64, seed::Integer;
                            use_lcg::Bool=false)

Local-type non-Gaussianity (mode 2, uncorrelated):
    δ_NL = δ_G + T(k)² × fNL × [χ² - ⟨χ²⟩]

where χ is an independent Gaussian random field with the same power spectrum,
generated with seed×5.  Modifies `delta` in place.
"""
function apply_fnl_uncorrelated!(delta::Array{Float32,3}, delta_k,
                                  pk_interp, tf_interp,
                                  n::Int, boxsize::Float64, fNL::Float64,
                                  seed::Integer; use_lcg::Bool=false)
    dk = 2π / boxsize
    nk = n ÷ 2 + 1
    nyq = n ÷ 2
    n3 = Float64(n)^3

    # Step 1: delta_G is already in delta (the Gaussian overdensity field)
    # Apply Fortran's 0.9925 correction factor
    delta .*= 0.9925f0

    # Step 2: Generate independent chi field (zeta_G) with seed*5
    chi_seed = seed * 5
    chi_pk = pk_interp  # same P(k)
    chi = if use_lcg
        generate_grf_lcg(n, chi_pk, boxsize, chi_seed)
    else
        generate_grf(n, chi_pk, boxsize, chi_seed)
    end

    # Step 3: Convert chi to zeta space: multiply by pk/tf² in k-space
    chi_k = rfft(chi)
    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = Float64(ix - 1) * dk
        ky = Float64(iy <= nyq + 1 ? iy - 1 : iy - 1 - n) * dk
        kz = Float64(iz <= nyq + 1 ? iz - 1 : iz - 1 - n) * dk
        k = sqrt(kx^2 + ky^2 + kz^2)
        if k == 0.0
            chi_k[ix, iy, iz] = 0.0f0
        else
            tf2 = tf_interp(k)^2
            if tf2 > 0.0
                chi_k[ix, iy, iz] /= Float32(tf2)
            else
                chi_k[ix, iy, iz] = 0.0f0
            end
        end
    end
    zeta_chi = irfft(chi_k, n) ./ Float32(n3)

    # Step 4: fNL contribution: fNL * (zeta_chi² - ⟨zeta_chi²⟩)
    avg_chi2 = sum(zeta_chi .^ 2) / n3
    delta_ng = Float32(fNL) .* (zeta_chi .^ 2 .- Float32(avg_chi2))

    # Step 5: Convert delta_ng from zeta space to delta space via T(k)²
    delta_ng_k = rfft(Float32.(delta_ng))
    for iz in 1:n, iy in 1:n, ix in 1:nk
        kx = Float64(ix - 1) * dk
        ky = Float64(iy <= nyq + 1 ? iy - 1 : iy - 1 - n) * dk
        kz = Float64(iz <= nyq + 1 ? iz - 1 : iz - 1 - n) * dk
        k = sqrt(kx^2 + ky^2 + kz^2)
        if k > 0.0
            tf2 = tf_interp(k)^2
            delta_ng_k[ix, iy, iz] *= Float32(tf2)
        else
            delta_ng_k[ix, iy, iz] = 0.0f0
        end
    end
    delta_ng_real = irfft(delta_ng_k, n) ./ Float32(n3)

    # Step 6: Add NG contribution to Gaussian delta
    delta .+= Float32.(delta_ng_real)

    # Update delta_k
    delta_k .= rfft(delta)
    return delta
end

end # module NonGaussian
