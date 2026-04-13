module PeakFind

struct PeakCandidate{T<:AbstractFloat}
    i::Int           # 1-based grid x-index
    j::Int           # 1-based grid y-index
    k::Int           # 1-based grid z-index
    ipp::Int         # column-major flat index: i + (j-1)*n1 + (k-1)*n1*n2
    x::Float64       # physical x coordinate (Mpc/h)
    y::Float64       # physical y coordinate (Mpc/h)
    z::Float64       # physical z coordinate (Mpc/h)
    delta::T         # overdensity at peak
    Rsmooth::Float64 # smoothing scale that found this peak
end

function find_peaks(delta::Array{T,3}, mask::Array{Int8,3},
                    xbx::Float64, ybx::Float64, zbx::Float64,
                    alatt::Float64, nbuff::Int, fcrit,
                    Rsmooth::Float64;
                    max_peaks::Int=2000000) where {T<:AbstractFloat}
    n1, n2, n3 = size(delta)
    cen1 = 0.5 * (n1 + 1)
    cen2 = 0.5 * (n2 + 1)
    cen3 = 0.5 * (n3 + 1)

    peaks = PeakCandidate[]

    for k in (nbuff+1):(n3-nbuff)
        zc = zbx + alatt * (k - cen3)
        for j in (nbuff+1):(n2-nbuff)
            yc = ybx + alatt * (j - cen2)
            for i in (nbuff+1):(n1-nbuff)
                xc = xbx + alatt * (i - cen1)

                if mask[i, j, k] == 1
                    continue
                end

                ff = delta[i, j, k]

                if ff < fcrit
                    continue
                end

                # Check 3x3x3 neighbourhood for strict local maximum
                is_max = true
                for kk in -1:1
                    kkk = k + kk
                    for jj in -1:1
                        jjj = j + jj
                        for ii in -1:1
                            iii = i + ii
                            if ff < delta[iii, jjj, kkk]
                                is_max = false
                                break
                            end
                        end
                        if !is_max break end
                    end
                    if !is_max break end
                end

                if !is_max
                    continue
                end

                ipp = i + (j - 1) * n1 + (k - 1) * n1 * n2
                push!(peaks, PeakCandidate(i, j, k, ipp, xc, yc, zc, ff, Rsmooth))
                mask[i, j, k] = 1

                if length(peaks) >= max_peaks
                    return peaks
                end
            end
        end
    end

    return peaks
end

end # module PeakFind
