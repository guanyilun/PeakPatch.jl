module LCG

# 48-bit multiplicative LCG matching the Fortran random.f90 exactly.
# The 48-bit state is stored as 4 base-4096 digits (12 bits each).
# Seed(1) is the most-significant digit, Seed(4) the least-significant.

const NDIGIT = 4
const MULTIPLIER = (373, 3707, 1442, 647)
const DEFAULT_SEED = (3281, 4041, 595, 2376)
const DIVISOR = (2.81474976710656e14, 6.8719476736e10, 1.6777216e7, 4096.0)

# ---------- core modular multiplication mod 2^48 ----------
function modmult(A, B)
    j1 = A[1] * B[1]
    j2 = A[1] * B[2] + A[2] * B[1]
    j3 = A[1] * B[3] + A[2] * B[2] + A[3] * B[1]
    j4 = A[1] * B[4] + A[2] * B[3] + A[3] * B[2] + A[4] * B[1]

    k1 = j1
    k2 = j2 + k1 ÷ 4096
    k3 = j3 + k2 ÷ 4096
    k4 = j4 + k3 ÷ 4096

    return (k1 % 4096, k2 % 4096, k3 % 4096, k4 % 4096)
end

# ---------- convert seed to Float64 in [0, 1) and advance ----------
function ranf(seed)
    r = seed[4] / DIVISOR[4] +
        seed[3] / DIVISOR[3] +
        seed[2] / DIVISOR[2] +
        seed[1] / DIVISOR[1]
    newseed = modmult(MULTIPLIER, seed)
    return r, newseed
end

# ---------- seed initialisation for parallel streams ----------
# Port of rans(): returns N seeds that partition the 2^46 cycle.
function rans(N::Int, startval::Int = 0)
    nn = isodd(N) ? N : N + 1

    if startval == 0
        seed1 = DEFAULT_SEED
    else
        seed1 = (abs(startval), 0, 0, 0)
    end

    seeds = Vector{NTuple{4,Int}}(undef, nn)

    if nn == 1
        seeds[1] = seed1
        return seeds
    end

    K = _ranfk(nn)
    Kbin = _ranfkbinary(K)
    atothek = _ranfatok(MULTIPLIER, Kbin)

    seeds[1] = seed1
    for i in 2:nn
        seeds[i] = modmult(seeds[i-1], atothek)
    end
    return seeds
end

# ---------- helper: compute K = floor(2^46 / nn) ----------
function _ranfk(N::Int)
    nn = isodd(N) ? N : N + 1
    q4 = 1024 ÷ nn
    r4 = 1024 - nn * q4
    q3 = (r4 * 4096) ÷ nn
    r3 = (r4 * 4096) - nn * q3
    q2 = (r3 * 4096) ÷ nn
    r2 = (r3 * 4096) - nn * q2
    q1 = (r2 * 4096) ÷ nn
    return (q1, q2, q3, q4)
end

# ---------- helper: expand K to 48-bit binary representation ----------
function _ranfkbinary(K)
    Kbin = Vector{Int}(undef, 48)
    for i in 1:4
        X = K[i] ÷ 2
        Kbin[(i-1)*12 + 1] = K[i] % 2
        for j in 2:12
            Kbin[(i-1)*12 + j] = X % 2
            X = X ÷ 2
        end
    end
    return Kbin
end

# ---------- helper: exponentiation by squaring: a^K mod 2^48 ----------
function _ranfatok(a, Kbin)
    asubi = a
    atothek = (1, 0, 0, 0)  # identity

    for i in 1:45
        if Kbin[i] != 0
            atothek = modmult(atothek, asubi)
        end
        asubi = modmult(asubi, asubi)
    end
    return atothek
end

# ---------- Box-Muller Gaussian deviates ----------
# Stateful: each pair of uniforms produces two Gaussians.
# We return one at a time, caching the second.
mutable struct GaussState
    have_spare::Bool
    spare::Float64
    GaussState() = new(false, 0.0)
end

function gaussdev(seed, state::GaussState = GaussState())
    if state.have_spare
        state.have_spare = false
        return state.spare, seed
    end

    local v1::Float64, v2::Float64, rsq::Float64
    while true
        v1, seed = ranf(seed)
        v2, seed = ranf(seed)
        v1 = 2.0 * v1 - 1.0
        v2 = 2.0 * v2 - 1.0
        rsq = v1^2 + v2^2
        (rsq >= 1.0 || rsq <= 0.0) || break
    end
    fac = sqrt(-2.0 * log(rsq) / rsq)
    state.spare = v1 * fac
    state.have_spare = true
    return v2 * fac, seed
end

end # module LCG
