module PowerSpectrum

using DelimitedFiles: readdlm
using Interpolations

"""Load a P(k) table from ASCII file and return a log-log interpolator.

The returned function `pk(k)` gives P(k) for any k, with log-linear
extrapolation beyond the table range.
"""
function load_pk(path::String)
    data = readdlm(path, comments=true)
    ks  = Float64.(data[:, 1])
    pks = Float64.(data[:, 2])

    logk  = log.(ks)
    logpk = log.(pks)

    itp = interpolate((logk,), logpk, Gridded(Linear()))
    extrap = extrapolate(itp, Line())

    return k -> exp(extrap(log(k)))
end

"""Load a P(k) table with transfer function for non-Gaussian runs.

Reads 4-column format: k, sqrt(P(k)), sqrt(k²T(k)), sqrt(P_chi(k)).
Returns (pk_interp, tf_interp) where:
- pk_interp(k) gives sqrt(P(k))  (the amplitude, not P(k) itself)
- tf_interp(k) gives sqrt(k²T(k)) (the transfer function factor)
"""
function load_pk_nongaussian(path::String)
    data = readdlm(path, comments=true)
    ks  = Float64.(data[:, 1])
    pks = Float64.(data[:, 2])  # sqrt(P(k))
    tfs = Float64.(data[:, 3])  # sqrt(k² T(k))

    # For pk: log-log interpolation (same as load_pk but on sqrt(P))
    logk  = log.(ks)
    logpk = log.(pks)
    itp_pk = interpolate((logk,), logpk, Gridded(Linear()))
    ext_pk = extrapolate(itp_pk, Line())
    pk_interp = k -> exp(ext_pk(log(k)))

    # For tf: log-log interpolation
    # Filter out zero values
    mask = tfs .> 0
    if any(mask)
        logtf = log.(tfs[mask])
        logk_tf = logk[mask]
        itp_tf = interpolate((logk_tf,), logtf, Gridded(Linear()))
        ext_tf = extrapolate(itp_tf, Line())
        tf_interp = k -> exp(ext_tf(log(k)))
    else
        tf_interp = k -> 0.0
    end

    return pk_interp, tf_interp
end

end # module PowerSpectrum
