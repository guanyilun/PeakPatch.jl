# Shell-analysis optimisation round — notes

**Branch:** `gpu-shell-opt-explore` (commit `aab1177`)
**Baseline (main):** `cda3b00` — measured 9.03 s GPU wall at nmesh=192,
ntile=2, ilpt=2, ioutshear=1 on A5000; `09_shell_analysis` 6.73 s
(76.8 % of GPU wall).

After the four main-branch optimisations landed (plan cache, fused ψ₁
FFTs, GPU-resident tiles, vectorised record packing), shell analysis
is the clear next target. This doc captures the first exploration
round, what was tried, what worked, and what's deferred.

## 1. What was changed

Two changes, layered on top of `main`:

### 1.1 Persistent per-call scratch pool *(active on branch)*

`analyse_peaks_gpu_cuda_multirf` previously called `CUDA.zeros(...)`
for ~20 peak-indexed buffers (Fshell, nshell, Sshell, S2shell, Gshell,
Gfshell, SRshell, lapdshell, RTHL, Fbarx, e_v, p_v, strain_final,
eigs, Srb, Sbar, Sbar2, gradpk, gradpkf, gradpkrf, d2F, zvir_half —
plus per-peak ZZon_pp / fcrit_pp and peak coordinate arrays) *per
batch*, then `CUDA.unsafe_free!`'d them all at the end. With many
batches per tile and many tiles per run, that's many thousands of
device allocations.

The branch refactors this to:

- Scan all batches, compute `max_np = max(npeaks_per_batch)`.
- Allocate each scratch buffer **once** at size `max_np` (or
  `(3, max_np)` / `(3, 3, max_np)` etc. as appropriate).
- Per-batch work takes `view(scratch, 1:npeaks)` slices — zero new
  allocations per batch, no per-batch `unsafe_free!`s.
- Free the pool once at call-end.

**Parity:** exact on halo count and every per-peak quantity
(verified via `test/test_multiresolution_gpu.jl`,
`test/test_multiresolution_gpu_ievol1.jl`,
`test/test_multigpu_dispatch.jl`).

**Measured impact:** below the noise floor of the shared benchmarking
machine at the time the branch was developed (load avg 30+, shared GPU
pinned at 100 % by another user). No reliable A/B; need a quiet box
to quantify.

### 1.2 Fused gather + post-process kernel *(dormant on branch)*

`_fused_shell_analysis_kernel!` (see `ext/CUDAExt.jl`) combines
`_shell_gather_full_kernel!` + `_post_process_kernel!` into one
launch. The per-shell sums stream directly through shared memory
instead of global:

```
Before (two kernels):
  gather   : per-peak block reduces 24 accumulators per shell,
             writes to global buffers sized (nshells, npeaks, ...)
  post-proc: reads those buffers back into shmem, runs phases 3–8
             on thread 0 of each block

After (one kernel):
  gather phase   : same block reductions, writes to shared memory
                   of the same block (no global traffic)
  post-proc phase: same phases 3–8 on thread 0, reads shmem
```

**Savings (theoretical):** ~100 MB of per-batch global traffic at
nshells=200, npeaks=5k (eliminated for both write-by-gather and
read-by-post-process).

**Parity:** exact on halo count at nmesh=68 (2127/2127 1LPT,
2127/2127 2LPT) and nmesh=192 (20467/20467). Top-10 per-peak RTHL,
pos, strain match bit-for-bit.

**Why it's dormant:** on SM_86 (A5000) it ran **~30 % slower** on the
shell-analysis stage. Root cause is the shmem budget: the post-process
needs ~13.6 KB of per-peak shell profile staging (rad, Fbar, Gn, Gfn,
SRn), and fusing with the gather means that full 13.6 KB is live
during the memory-bound gather phase. Shmem per block:

- Two-kernel: gather 256 B, post-process 13.6 KB (but only 1 thread
  per block does work).
- Fused: 20.3 KB live throughout (including Sshell/S2shell/lapd/n
  which the post-process previously read from global in Phase 7).

On SM_86 the shmem/occupancy interaction costs us:

| Variant | Shmem/block | Blocks/SM | Warps/SM | Warp occupancy |
|---|---|---|---|---|
| Two-kernel (gather) | 256 B | 12 | 48 | 48 / 48 = 100 % |
| Two-kernel (post-proc) | 13.6 KB | 3 | 12 | 12 / 48 = 25 % |
| Fused | 20.3 KB | 2 | 8 | 8 / 48 = 17 % |

The gather phase is memory-bandwidth bound; the ~6× drop in active
warps kills latency hiding for DRAM reads. The ~100 MB saved in
global traffic is simply less than the bandwidth we lose by
stalling on uncoalesced reads.

The post-process is *compute-serial* (1 thread per block runs the
sequential Phase 3–8 walk), so occupancy doesn't matter there —
which is why the two-kernel layout is actually the right fit for
SM_86. It trades gather-phase high-occupancy for post-process
sequentialism, and pays only one round-trip through global in
between.

**Where the fused kernel could still win:**

1. **H100 (SM_90)**: shmem/SM is 228 KB (vs 100 KB on SM_86), so
   20.3 KB/block allows 4+ blocks/SM — the occupancy hit would be
   roughly halved. Re-measure there.
2. **Large `npeaks` per batch**: the global-traffic savings scale
   linearly with `nshells × npeaks` while the occupancy hit is
   ~constant. At production scale (nmesh=399, ntile=16, thousands
   of peaks/batch) the break-even may flip.
3. **Reduced shmem footprint**: drop the Sshell/S2shell/n/lapd
   staging (go back to global for those) — brings shmem down to
   ~14 KB/block, same as post-process alone, so no occupancy
   penalty vs the post-process phase. But then most of the
   global-traffic savings also disappear.

The source is kept in `ext/CUDAExt.jl` with a block comment
describing the shmem layout (`_FSA_OFF_*` constants) and how to
re-enable it. To wire it back in, swap the `@cuda ... _post_process_kernel!`
call in `analyse_peaks_gpu_cuda_multirf` for a single
`@cuda ... _fused_shell_analysis_kernel!` call; remove the
`_launch_shell_gather_full!` invocation; drop the per-shell
intermediate views (Fshell_d … lapdshell_d). ~30 LOC change.

## 2. Decision

`main` (`cda3b00`) stays the known-good measured state. The branch
retains both changes so it's one `git checkout` away to:

- Re-measure the persistent-scratch delta on a quieter machine.
- Re-enable the fused kernel on H100 / L40S to see if the occupancy
  tradeoff flips.

Promote to `main` only when a clean A/B shows a sustained positive
delta across nmesh ∈ {128, 192, 256} on A5000 or L40S.

## 3. What else was considered but not attempted

From the earlier exploration planning:

- **Warp-level parallelisation of post-process** — phases with real
  data dependencies (inward walk `zvir1p` chain, fcrit crossing) are
  inherently serial. The parallelisable phases (Sbar/S2bar sum, Srb
  energy integral, zvir_half's 10-iteration outer loop) together are
  <20 % of shell_analysis time. Net expected speedup 5–10 %, behind a
  substantial rewrite of the control flow. Not attempted.

- **Better block geometry** (cooperative-groups / persistent kernels
  / flatter cells × peaks decomposition) — the current block-per-peak
  uses all 128 threads cooperatively during gather; the inner-shell
  under-utilisation (nc < bdim) is real but small. A flatter layout
  would require a multi-bin reduction across both (cell, peak) axes;
  non-trivial. Not attempted.

- **Float16 storage for tile fields** — ruled out in pre-implementation
  analysis: 5e-4 relative quantisation of δ ≈ 1.7 near the `fcrit`
  threshold can flip collapse decisions, producing a 0.1–0.5 %
  halo-count drift that violates the parity-exact invariant we've
  relied on for every other GPU change. Pre-decision conversation
  captured in chat; not pursued.

## 4. Reproducibility

To re-benchmark this branch on a clean box:

```bash
git checkout main
# ... baseline numbers ...
CUDA_VISIBLE_DEVICES=0 julia --project=/tmp/pp_gpu_sandbox --threads=32 \
    bench/profile_multitile.jl > /tmp/pp_main.log 2>&1

git checkout gpu-shell-opt-explore
CUDA_VISIBLE_DEVICES=0 julia --project=/tmp/pp_gpu_sandbox --threads=32 \
    bench/profile_multitile.jl > /tmp/pp_branch.log 2>&1

diff /tmp/pp_main.log /tmp/pp_branch.log
```

Compare the `09_shell_analysis` row and the total GPU wall between the
two logs. If branch is ≥5 % faster on total wall, merge. If neutral or
worse, leave on branch as documented exploration.

Check box hygiene first (bench is useless on contaminated hardware):

```bash
uptime   # load avg should be < 10
nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv
# target GPU should show ~0 MiB used, 0 % util
```
