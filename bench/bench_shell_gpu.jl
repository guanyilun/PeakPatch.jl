# Large-scale CPU vs GPU comparison for PeakPatch shell analysis.
#
# Run: CUDA_VISIBLE_DEVICES=1 julia --project=/tmp/pp_gpu_sandbox bench/bench_shell_gpu.jl
#
# Generates a realistic 3D field on an N³ grid, runs `find_peaks` to locate
# peak candidates, then times CPU `analyse_peak_gpu` against GPU
# `analyse_peak_gpu_cuda`. Validates a random sample of peaks field-by-field
# and compares the full mask grid bit-by-bit.
#
# Defaults: N=512, rmax=12, ~200 Gaussian bumps. Override via CLI args
# or env vars: N_GRID, RMAX, NBUMPS, NSAMPLE.

using PeakPatch
using CUDA
using Random
using Printf

const N       = parse(Int, get(ENV, "N_GRID",   "512"))
const RMAX    = parse(Int, get(ENV, "RMAX",     "12"))
const NBUMPS  = parse(Int, get(ENV, "NBUMPS",   "200"))
const NSAMPLE = parse(Int, get(ENV, "NSAMPLE",  "100"))

println("\n" * "="^60)
println("PeakPatch CPU vs GPU — Shell analysis benchmark")
println("="^60)
@printf "Grid           : %d³ (%.2f GB per Float32 field)\n"  N (N^3 * 4 / 1e9)
@printf "Shell rmax     : %d\n" RMAX
@printf "Bump count     : %d\n" NBUMPS
@printf "Sample to val  : %d peaks (random subset)\n" NSAMPLE
@printf "CUDA functional: %s\n" CUDA.functional()
if CUDA.functional()
    @printf "CUDA device    : %s\n" CUDA.name(CUDA.device())
end
println()

# ============================================================
# Build shell tables (shared by CPU and GPU)
# ============================================================
println("[1/7] Building shell tables (rmax=$RMAX)…")
t_stab = @elapsed stab = build_shell_tables(RMAX)
@printf "       nshells=%d, time=%.3fs\n" stab.nshells t_stab

# ============================================================
# Generate field on host (this is the bulk of wall-clock setup)
# ============================================================
println("[2/7] Generating $(N)³ field (delta + 3η + 3η² + lapd)…")
Random.seed!(31415)
t_gen = @elapsed begin
    delta = 0.05f0 .* randn(Float32, N, N, N)
    for _ in 1:NBUMPS
        ci = rand(RMAX+2:N-RMAX-1); cj = rand(RMAX+2:N-RMAX-1); ck = rand(RMAX+2:N-RMAX-1)
        amp = 1.0f0 + 2.0f0 * rand(Float32)
        sigma2 = 15.0f0 + 10.0f0 * rand(Float32)
        # Only fill within a bounding box around the bump (8σ ~ 14 cells)
        half = 14
        i0 = max(1, ci - half); i1 = min(N, ci + half)
        j0 = max(1, cj - half); j1 = min(N, cj + half)
        k0 = max(1, ck - half); k1 = min(N, ck + half)
        for k in k0:k1, j in j0:j1, i in i0:i1
            d2 = (i-ci)^2 + (j-cj)^2 + (k-ck)^2
            @inbounds delta[i, j, k] += amp * exp(-d2 / sigma2)
        end
    end
    etax  = 0.3f0  .* randn(Float32, N, N, N)
    etay  = 0.3f0  .* randn(Float32, N, N, N)
    etaz  = 0.3f0  .* randn(Float32, N, N, N)
    eta2x = 0.15f0 .* randn(Float32, N, N, N)
    eta2y = 0.15f0 .* randn(Float32, N, N, N)
    eta2z = 0.15f0 .* randn(Float32, N, N, N)
    lapd  = 0.05f0 .* randn(Float32, N, N, N)
end
@printf "       time=%.2fs\n" t_gen

# ============================================================
# Find peaks via production path
# ============================================================
println("[3/7] Finding peaks (threshold δ≥0.8, Rsmooth=3)…")
t_find = @elapsed begin
    mask_find = zeros(Int8, N, N, N)
    peak_candidates = PeakPatch.find_peaks(delta, mask_find,
                                            0.0, 0.0, 0.0, 1.0, 1,
                                            0.8f0, 3.0)
end
npeaks = length(peak_candidates)
@printf "       npeaks=%d, time=%.2fs\n" npeaks t_find
if npeaks == 0
    error("No peaks found — adjust threshold or bump parameters")
end

peaks_i = Int32[p.i for p in peak_candidates]
peaks_j = Int32[p.j for p in peak_candidates]
peaks_k = Int32[p.k for p in peak_candidates]

# Collapse table with structure (not trivial) so interpolation is exercised
ct_Nx, ct_Ny, ct_Nz = 11, 11, 11
ct_X1, ct_X2 = -2.0, 2.0
ct_Y1, ct_Y2 =  0.0, 1.0
ct_Z1, ct_Z2 =  0.0, 2.0
ct_table = Float32[
    0.5f0 + 0.3f0*sin(Float32(i-1)*0.3f0) + 0.2f0*cos(Float32(j-1)*0.4f0) + 0.1f0*Float32(k-1)*0.1f0
    for i in 1:ct_Nx, j in 1:ct_Ny, k in 1:ct_Nz]
tp = PeakPatch.CollapseTableParams(Nx=ct_Nx, Ny=ct_Ny, Nz=ct_Nz,
                                    X1=ct_X1, X2=ct_X2, Y1=ct_Y1, Y2=ct_Y2,
                                    Z1=ct_Z1, Z2=ct_Z2)
ct = PeakPatch.CollapseTableInterp(ct_table, tp)

alatt  = 1.0
ir2min = 4
ZZon   = 0.0
Rfclvi = 4.0
fcrit  = 1.0

# ============================================================
# CPU baseline — single-threaded (authoritative; writes mask_cpu)
# ============================================================
println("[4a/7] CPU analyse_peak_gpu baseline (single-threaded, writes mask)…")
mask_cpu = zeros(Int8, N, N, N)
pg = PeakGrid(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
               mask_cpu, (N, N, N), lapd)
cpu_results = Vector{Any}(undef, npeaks)
t_cpu_st = @elapsed begin
    for idx in 1:npeaks
        ipp = Int(peaks_i[idx]) + (Int(peaks_j[idx]) - 1) * N + (Int(peaks_k[idx]) - 1) * N * N
        cpu_results[idx] = analyse_peak_gpu(pg, ipp, alatt, ir2min, ZZon, Rfclvi, ct, stab;
                                             fcrit_override=fcrit)
    end
end
ncollapsed_cpu = count(r -> r.RTHL > 0.0, cpu_results)
@printf "       total=%.2fs (%.2fms/peak), collapsed=%d/%d\n" t_cpu_st (1000*t_cpu_st/npeaks) ncollapsed_cpu npeaks

# ============================================================
# CPU baseline — multi-threaded (no mask side-effect for timing purity)
# ============================================================
nthreads = Threads.nthreads()
println("[4b/7] CPU analyse_peak_gpu baseline (threaded, $nthreads threads, no mask)…")
# Shared PeakGrid WITHOUT mask (to avoid racy writes). Physics-only timing.
pg_nomask = PeakGrid(delta, etax, etay, etaz, eta2x, eta2y, eta2z,
                      nothing, (N, N, N), lapd)
t_cpu_mt = @elapsed begin
    Threads.@threads for idx in 1:npeaks
        ipp = Int(peaks_i[idx]) + (Int(peaks_j[idx]) - 1) * N + (Int(peaks_k[idx]) - 1) * N * N
        _ = analyse_peak_gpu(pg_nomask, ipp, alatt, ir2min, ZZon, Rfclvi, ct, stab;
                              fcrit_override=fcrit)
    end
end
@printf "       total=%.2fs (%.2fms/peak)\n" t_cpu_mt (1000*t_cpu_mt/npeaks)

# ============================================================
# GPU
# ============================================================
println("[5/7] GPU analyse_peak_gpu_cuda…")
mask_gpu = zeros(Int8, N, N, N)
# Warm up — first launch pays JIT cost
CUDA.@sync begin
    _warmup = analyse_peak_gpu_cuda(
        delta, etax, etay, etaz, eta2x, eta2y, eta2z,
        peaks_i[1:1], peaks_j[1:1], peaks_k[1:1], stab,
        ct_table, ct_X1, ct_X2, ct_Y1, ct_Y2, ct_Z1, ct_Z2,
        alatt, ir2min, ZZon, Rfclvi;
        fcrit_override=fcrit, lapd=lapd, mask=zeros(Int8, N, N, N), nbuff=0)
end
t_gpu = CUDA.@elapsed begin
    global gpu_res = analyse_peak_gpu_cuda(
        delta, etax, etay, etaz, eta2x, eta2y, eta2z,
        peaks_i, peaks_j, peaks_k, stab,
        ct_table, ct_X1, ct_X2, ct_Y1, ct_Y2, ct_Z1, ct_Z2,
        alatt, ir2min, ZZon, Rfclvi;
        fcrit_override=fcrit, lapd=lapd, mask=mask_gpu, nbuff=0)
    CUDA.synchronize()
end
ncollapsed_gpu = count(r -> r > 0.0f0, gpu_res.RTHL)
@printf "       total=%.3fs (%.2fms/peak), collapsed=%d/%d\n" t_gpu (1000*t_gpu/npeaks) ncollapsed_gpu npeaks

# Measure upload cost alone (for amortization analysis)
println("[5b/7] Timing breakdown: upload vs kernel work…")
t_upload = CUDA.@elapsed begin
    _d_d = CuArray(delta); _ex_d = CuArray(etax); _ey_d = CuArray(etay); _ez_d = CuArray(etaz)
    _e2x_d = CuArray(eta2x); _e2y_d = CuArray(eta2y); _e2z_d = CuArray(eta2z); _l_d = CuArray(lapd)
    CUDA.synchronize()
end
# Free the probe uploads so they don't linger
CUDA.unsafe_free!(_d_d); CUDA.unsafe_free!(_ex_d); CUDA.unsafe_free!(_ey_d); CUDA.unsafe_free!(_ez_d)
CUDA.unsafe_free!(_e2x_d); CUDA.unsafe_free!(_e2y_d); CUDA.unsafe_free!(_e2z_d); CUDA.unsafe_free!(_l_d)
t_kernel = max(t_gpu - t_upload, 0.0)
@printf "       upload time (8 fields, %.1f GB): %.3f s  (%.1f GB/s effective)\n" (8*N^3*4/1e9) t_upload (8*N^3*4/t_upload/1e9)
@printf "       kernel work (gather + post-proc): %.3f s  (%.2fms/peak)\n" t_kernel (1000*t_kernel/npeaks)

# ============================================================
# Validate sample of peaks
# ============================================================
println("[6/7] Validating random sample of $NSAMPLE peaks…")
Random.seed!(42)
sample_idxs = randperm(npeaks)[1:min(NSAMPLE, npeaks)]

global n_collapse_agree = 0
global n_collapse_both  = 0
scalar_errs = Dict(:RTHL => Float64[], :Fbarx => Float64[], :e_v => Float64[],
                   :p_v => Float64[], :Srb => Float64[], :d2F => Float64[])
vec_errs    = Dict(:Sbar => Float64[], :Sbar2 => Float64[],
                   :gradpk => Float64[], :gradpkf => Float64[], :gradpkrf => Float64[])

for idx in sample_idxs
    cpu_r = cpu_results[idx]
    cpu_c = cpu_r.RTHL > 0.0
    gpu_c = gpu_res.RTHL[idx] > 0.0f0
    cpu_c == gpu_c && (global n_collapse_agree += 1)
    if cpu_c && gpu_c
        global n_collapse_both += 1
        # Relative error for non-zero, absolute for near-zero
        relerr(a, b) = abs(a) > 1f-4 ? abs(a - b) / abs(a) : abs(a - b)
        push!(scalar_errs[:RTHL],  relerr(Float32(cpu_r.RTHL),  gpu_res.RTHL[idx]))
        push!(scalar_errs[:Fbarx], relerr(Float32(cpu_r.Fbarx), gpu_res.Fbarx[idx]))
        push!(scalar_errs[:e_v],   relerr(Float32(cpu_r.e_v),   gpu_res.e_v[idx]))
        push!(scalar_errs[:p_v],   relerr(Float32(cpu_r.p_v),   gpu_res.p_v[idx]))
        push!(scalar_errs[:Srb],   relerr(Float32(cpu_r.Srb),   gpu_res.Srb[idx]))
        push!(scalar_errs[:d2F],   relerr(Float32(cpu_r.d2F),   gpu_res.d2F[idx]))
        for L in 1:3
            push!(vec_errs[:Sbar],    relerr(Float32(cpu_r.Sbar[L]),     gpu_res.Sbar[L, idx]))
            push!(vec_errs[:Sbar2],   relerr(Float32(cpu_r.Sbar2[L]),    gpu_res.Sbar2[L, idx]))
            push!(vec_errs[:gradpk],  relerr(Float32(cpu_r.gradpk[L]),   gpu_res.gradpk[L, idx]))
            push!(vec_errs[:gradpkf], relerr(Float32(cpu_r.gradpkf[L]),  gpu_res.gradpkf[L, idx]))
            push!(vec_errs[:gradpkrf],relerr(Float32(cpu_r.gradpkrf[L]), gpu_res.gradpkrf[L, idx]))
        end
    end
end

@printf "       collapse match : %d/%d\n" n_collapse_agree length(sample_idxs)
@printf "       both collapsed : %d\n" n_collapse_both
println("       field | max relerr  | mean relerr")
println("       ------|-------------|------------")
for (k, v) in scalar_errs
    isempty(v) && continue
    @printf "       %-7s %.2e    %.2e\n" string(k) maximum(v) (sum(v)/length(v))
end
for (k, v) in vec_errs
    isempty(v) && continue
    @printf "       %-7s %.2e    %.2e\n" string(k) maximum(v) (sum(v)/length(v))
end

# ============================================================
# Validate mask side-effect
# ============================================================
println("[7/7] Validating mask side-effect…")
mask_match = gpu_res.mask == mask_cpu
cpu_mask_sum = sum(mask_cpu)
gpu_mask_sum = sum(gpu_res.mask)
@printf "       CPU mask cells set : %d\n" cpu_mask_sum
@printf "       GPU mask cells set : %d\n" gpu_mask_sum
@printf "       masks bit-identical: %s\n" mask_match
if !mask_match
    diffs = findall(mask_cpu .!= gpu_res.mask)
    @printf "       differing cells    : %d\n" length(diffs)
end

# ============================================================
# Summary
# ============================================================
println("\n" * "="^60)
println("SUMMARY")
println("="^60)
@printf "CPU single-threaded    : %.3f s (%.2fms/peak)\n" t_cpu_st (1000*t_cpu_st/npeaks)
@printf "CPU multi-threaded (%2d): %.3f s (%.2fms/peak)\n" nthreads t_cpu_mt (1000*t_cpu_mt/npeaks)
@printf "CPU threading speedup  : %.1fx\n" t_cpu_st / t_cpu_mt
@printf "GPU total              : %.3f s (%.2fms/peak)\n" t_gpu (1000*t_gpu/npeaks)
@printf "  upload (amortisable) : %.3f s\n" t_upload
@printf "  kernel work only     : %.3f s (%.2fms/peak)\n" t_kernel (1000*t_kernel/npeaks)
@printf "GPU vs CPU-1T speedup  : %.1fx (total) / %.1fx (kernel-only)\n" t_cpu_st / t_gpu t_cpu_st / max(t_kernel, 1e-9)
@printf "GPU vs CPU-%2dT speedup : %.1fx (total) / %.1fx (kernel-only)\n" nthreads t_cpu_mt / t_gpu t_cpu_mt / max(t_kernel, 1e-9)
@printf "Peaks processed        : %d (collapsed: CPU=%d, GPU=%d)\n" npeaks ncollapsed_cpu ncollapsed_gpu
@printf "Sample field agreement : %d/%d both-collapse, max scalar relerr=%.2e\n" n_collapse_both length(sample_idxs) maximum(maximum.(values(scalar_errs); init=0.0))
@printf "Mask bit-identical     : %s\n" mask_match
println("="^60)
