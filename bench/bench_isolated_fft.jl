# Benchmark CPU FFTW vs GPU cuFFT for the isolated convolution that
# powers run_multitile_split.
#
# Run: CUDA_VISIBLE_DEVICES=1 julia --project=/tmp/pp_gpu_sandbox \
#      --threads=64 bench/bench_isolated_fft.jl
#
# Override grid via env: N_GRID (default 384), NSHELL (default 0).

using PeakPatch
using CUDA
using Random
using Printf
import PeakPatch.MultiResolution: _isolated_convolve, _kernel_1lpt, _kernel_2lpt

const N      = parse(Int, get(ENV, "N_GRID", "384"))
const NSHELL = parse(Int, get(ENV, "NSHELL", "0"))

println("\n" * "="^60)
println("Isolated FFT convolution: CPU vs GPU")
println("="^60)
n_input = N + 2 * NSHELL
n2 = 2 * n_input
@printf "Tile output size  : %d³\n" N
@printf "Noise input size  : %d³ (NSHELL=%d)\n" n_input NSHELL
@printf "Padded FFT size   : %d³ (%.2f GB Float32)\n" n2 (Float64(n2)^3 * 4 / 1e9)
@printf "Threads (CPU FFTW): %d\n" Threads.nthreads()
@printf "GPU device        : %s\n" CUDA.functional() ? CUDA.name(CUDA.device()) : "n/a"
println()

Random.seed!(42)
noise = randn(Float32, n_input, n_input, n_input)
pk_test = k -> 100.0 / (1.0 + (k * 50.0)^4)
boxsize = 100.0

# ============================================================
# CPU baseline (uses FFTW with whatever Threads.nthreads() is set to)
# ============================================================
# FFTW uses its own thread count; we set it to match Julia threads.
import FFTW
FFTW.set_num_threads(Threads.nthreads())

# Warm up CPU
_ = _isolated_convolve(noise, pk_test, boxsize, N; nshell=NSHELL)

println("CPU _isolated_convolve (FFTW):")
for (label, kernel_fn) in [
    ("delta",  nothing),
    ("psi1_x", _kernel_1lpt(1)),
    ("psi2_x", _kernel_2lpt(1)),
]
    t = @elapsed begin
        for _ in 1:3
            _ = _isolated_convolve(noise, pk_test, boxsize, N;
                                    kernel_fn=kernel_fn, nshell=NSHELL)
        end
    end
    @printf "  %-7s : %6.3f s/call (avg of 3)\n" label (t / 3)
end

# ============================================================
# GPU
# ============================================================
println("\nGPU isolated_convolve_gpu (cuFFT):")
# Warm up GPU
CUDA.@sync _ = isolated_convolve_gpu(noise, pk_test, boxsize, N; nshell=NSHELL)

for (label, id, dim1) in [
    ("delta",  0, 0),
    ("psi1_x", 1, 1),
    ("psi2_x", 2, 1),
]
    t = CUDA.@elapsed begin
        for _ in 1:3
            _ = isolated_convolve_gpu(noise, pk_test, boxsize, N;
                                       kernel_fn_id=id, dim1=dim1, nshell=NSHELL)
            CUDA.synchronize()
        end
    end
    @printf "  %-7s : %6.3f s/call (avg of 3)\n" label (t / 3)
end

# ============================================================
# Eight-call sequence (typical run_multitile_split per-tile workload):
# δ, ψ1×3, ψ2×3, lap = 8 isolated FFTs per tile
# ============================================================
println("\nFull tile-load (8 FFTs: δ + 3×ψ1 + 3×ψ2 + lap):")

t_cpu = @elapsed begin
    _ = _isolated_convolve(noise, pk_test, boxsize, N; nshell=NSHELL)
    for dim in 1:3
        _ = _isolated_convolve(noise, pk_test, boxsize, N; kernel_fn=_kernel_1lpt(dim), nshell=NSHELL)
    end
    for dim in 1:3
        _ = _isolated_convolve(noise, pk_test, boxsize, N; kernel_fn=_kernel_2lpt(dim), nshell=NSHELL)
    end
    _ = _isolated_convolve(noise, pk_test, boxsize, N; kernel_fn=PeakPatch.MultiResolution._kernel_laplacian(), nshell=NSHELL)
end
@printf "  CPU FFTW (%2d threads): %6.3f s\n" Threads.nthreads() t_cpu

t_gpu = CUDA.@elapsed begin
    _ = isolated_convolve_gpu(noise, pk_test, boxsize, N; kernel_fn_id=0, nshell=NSHELL)
    for dim in 1:3
        _ = isolated_convolve_gpu(noise, pk_test, boxsize, N; kernel_fn_id=1, dim1=dim, nshell=NSHELL)
    end
    for dim in 1:3
        _ = isolated_convolve_gpu(noise, pk_test, boxsize, N; kernel_fn_id=2, dim1=dim, nshell=NSHELL)
    end
    _ = isolated_convolve_gpu(noise, pk_test, boxsize, N; kernel_fn_id=4, nshell=NSHELL)
    CUDA.synchronize()
end
@printf "  GPU cuFFT             : %6.3f s\n" t_gpu
@printf "  GPU speedup           : %.1fx\n" t_cpu / t_gpu
println("="^60)
