# Plan: GPU shell analysis for lightcone mode (`ievol=1`)

## 1. Problem statement

The current GPU shell-analysis batch kernel
(`_post_process_kernel!` in `ext/CUDAExt.jl`) assumes a single scalar
collapse redshift `ZZon` for all peaks in a tile. That assumption is
only valid for **fixed-redshift snapshots** (`cfg.ievol == 0`).

In **lightcone mode** (`cfg.ievol == 1`), each halo collapses at the
redshift corresponding to its own comoving distance from the observer:

```julia
# src/MultiResolution.jl:881 (CPU per-peak path)
z_pk = peak_redshift(obs[1], obs[2], obs[3], peak.x, peak.y, peak.z, chi2z)
ZZon_pk = 1.0 + z_pk
z_pk > z_max && continue            # skip peaks beyond max redshift
```

Because the GPU batch kernel currently accepts only a scalar ZZon, the
dispatcher in `src/MultiResolution.jl:553` gates it off:

```julia
use_gpu_shell = use_gpu && ievol == 0
```

With `ievol == 1` the pipeline falls back to the CPU `analyse_peak`
per-peak loop, which is 55% of wall time at production scale. The
Killarney projection shows lightcone octants taking ~9 h instead of
~4.5 h on 4× L40S. Enabling GPU shell analysis for `ievol=1` is the
single biggest wall-time lever before production lightcone runs.

## 2. What changes per peak in `ievol=1`

In the GPU kernel, the scalar `ZZon` and `fcrit` (derived from ZZon via
`fsc_of_z(ZZon - 1.0, growth_tables)`) are used at three call sites
inside `_post_process_kernel!`:

| Line (current) | Use |
|---|---|
| `ext/CUDAExt.jl:1556` | `if zvir1p_m0 >= ZZon` — outer collapse test |
| `ext/CUDAExt.jl:1620` | `if zvir1p_mp >= ZZon` — inner shell-walk collapse test |
| `ext/CUDAExt.jl:1674` | `(zvir1p - ZZon) / dZvir` — RTHL interpolation factor |
| Phase 4 body (`~1448`) | `fcrit` threshold in the Fbar shell search |

The CPU reference (`src/HaloFinder/RadialShell.jl:283`) shows the same
structure: `ZZon` and `fcrit` are scalar within `analyse_peak`, but
`analyse_peak` is called per-peak in the ievol=1 CPU path, so they
vary across peaks of the same tile.

**Everything else in the kernel is already per-peak** (peak indices,
smoothing indices, strain tensor, etc.). This is a localised change.

## 3. Peak filtering (`z_pk > z_max`)

On top of the per-peak scalars, `ievol=1` also filters peaks with
`z_pk > z_max` — they are never analysed. The GPU batch should honor
this by:

- Computing `z_pk` for every candidate peak on the CPU (cheap — it's
  distance + `chi_to_z` interpolation).
- Dropping peaks with `z_pk > z_max` from the batch entirely, so the
  kernel runs over a pre-filtered peak list.
- Re-aligning the per-peak output slots so downstream record-packing
  sees the same peak count as if CPU had been used.

This avoids launching useless kernel threads and preserves catalog
identity with the CPU path.

## 4. Implementation phases

### Phase 1 — Backward-compatible scalar→array plumbing

No behaviour change. Just make the kernel and batch entry point
accept per-peak arrays internally, while callers can still pass a
scalar that gets broadcast.

Touchpoints:
- `ext/CUDAExt.jl:1285` — `_post_process_kernel!` signature: replace
  `ZZon::Float32, fcrit::Float32` with
  `ZZon_pp::CuDeviceVector{Float32}, fcrit_pp::CuDeviceVector{Float32}`.
- Inside the kernel, at the top of the `peak_id` body, read
  `ZZon = ZZon_pp[peak_id]; fcrit = fcrit_pp[peak_id]` as locals. The
  rest of the kernel body is unchanged — it already uses `ZZon` and
  `fcrit` as locals.
- `ext/CUDAExt.jl:~2154, ~2324` — at the `@cuda` launches, swap the
  scalar args for per-peak `CuArray`s. Build them on the host as
  `fill(Float32(ZZon), npeaks)` / `fill(Float32(fcrit), npeaks)` when
  the caller passed scalars.
- `PeakPatch.analyse_peaks_gpu_cuda`, `analyse_peaks_gpu_cuda_multirf`:
  add kwargs `ZZon_pp::Union{Nothing,AbstractVector}=nothing`,
  `fcrit_pp::Union{Nothing,AbstractVector}=nothing`. If `nothing`,
  broadcast the scalar. If given, use directly.

**Validation**: all current `ievol=0` tests (`test_shell_gpu.jl`,
`test_multiresolution_gpu.jl`, `test_multigpu_dispatch.jl`) must
continue to pass **bit-identically**, because the kernel with a
broadcast scalar is numerically the same as the scalar kernel.

Estimated effort: **~2 h** (mostly mechanical, plus one compile-time
check that the new `CuDeviceVector` reads don't bloat register
pressure — the kernel is already shmem-heavy).

### Phase 2 — Host-side per-peak ZZon / fcrit assembly

In `_gpu_analyse_peaks_batch`
(`src/MultiResolution.jl:374`) and its one caller
(`src/MultiResolution.jl:855`), build the per-peak arrays when
`ievol == 1`:

```julia
if ievol == 1
    ZZon_pp  = Vector{Float32}(undef, length(tile_peaks))
    fcrit_pp = Vector{Float32}(undef, length(tile_peaks))
    @inbounds for idx in eachindex(tile_peaks)
        p = tile_peaks[idx]
        z_pk = peak_redshift(obs[1], obs[2], obs[3], p.x, p.y, p.z, chi2z)
        ZZon_pp[idx]  = Float32(1.0 + z_pk)
        fcrit_pp[idx] = Float32(fsc_of_z(z_pk, growth_tables))
    end
else
    ZZon_pp  = fill(Float32(ZZon),  length(tile_peaks))
    fcrit_pp = fill(Float32(fcrit), length(tile_peaks))
end
```

Then pass `ZZon_pp, fcrit_pp` through to the batch helper, which
forwards them to `analyse_peaks_gpu_cuda_multirf`.

Estimated effort: **~1 h**.

### Phase 3 — Peak filtering by `z_max`

Before constructing the batch, drop peaks with `z_pk > z_max`. Do this
with an index remap so the GPU output slots align with the surviving
peak list:

```julia
if ievol == 1
    # compute z_pk per peak (as above)
    keep = findall(z_pk_arr .<= z_max)
    tile_peaks = tile_peaks[keep]
    tile_Rf    = tile_Rf[keep]
    tile_FcRf  = tile_FcRf[keep]
    tile_d2Rf  = tile_d2Rf[keep]
    ZZon_pp    = ZZon_pp[keep]
    fcrit_pp   = fcrit_pp[keep]
end
```

The post-shell record-packing loop (`src/MultiResolution.jl:893`)
then iterates over the filtered arrays, identical to CPU behaviour.

Remove the `z_pk > z_max && continue` branch from the record-packing
loop when `use_gpu_shell` path is active (already filtered upstream).

Estimated effort: **~1 h**.

### Phase 4 — Remove the `ievol == 0` gate

Change `src/MultiResolution.jl:553` from:

```julia
use_gpu_shell = use_gpu && ievol == 0
```

to simply:

```julia
use_gpu_shell = use_gpu
```

The per-peak `ZZon_pp`/`fcrit_pp` arrays now make the batch valid
for both cases.

Keep a feature flag `cfg_use_gpu_shell=true` kwarg on
`run_multitile_split` in case we ever need to disable GPU shell for
debugging; default is enabled whenever `use_gpu=true`.

Estimated effort: **~30 min** (small code change, plus adjust
fall-through branches in `MultiResolution.jl:881-884` where CPU path
handles ievol=1 per-peak).

### Phase 5 — Parity test and benchmark

Add a new test analogous to `test_multiresolution_gpu.jl` but with
`ievol=1`:

- Observer at `(0,0,0)`, box sized so the far corner has a
  meaningful `z > 0` (e.g. N=68 at the existing test config already
  yields box ≈ 227 Mpc/h, z at corner ~0.08).
- Run CPU and GPU paths; compare halo counts, top-10 positions,
  RTHL, and the zvir field (which depends on per-peak ZZon).
- Additionally test `z_max` cutoff: set `z_max` small enough that
  some peaks are filtered; verify both paths produce the same
  post-filter catalog.

For the bench, extend `bench/profile_multitile.jl` with an
`IEVOL=1` knob and compare stage breakdowns. Expected:
- Before this change: `09_shell_analysis` (CPU fallback) dominates
- After: comparable to ievol=0 GPU timings

Estimated effort: **~3 h** (test writing + iteration on parity
debug if kernel has edge cases with per-peak ZZon).

### Phase 6 — Documentation updates

- Update `docs/killarney_deployment.md` §5 (lightcone mode) and §7
  (optimization roadmap) to note item #1 is complete.
- Update `docs/gpu_session_log.org` with the per-stage profile for
  ievol=1 before and after.
- Update the TL;DR wall-time projections:
  `9 h → ~4.5 h` per octant on lightcone.

Estimated effort: **~30 min**.

## 5. Total effort

~1 day of focused work, broken out as above. Phases 1–4 are the
code change (~4.5 h), Phase 5 is validation (~3 h), Phase 6 is
docs (~0.5 h).

## 6. Risks and unknowns

1. **Register pressure from per-peak array reads**. The
   `_post_process_kernel!` is already near the shmem limit
   (17 × `_MAX_SHELLS_GPU` floats). Adding two device-vector reads
   adds only a couple of register loads per peak, but check the
   compiled PTX for register spill with
   `CUDA.@device_code_warntype` and `CUDA.registers`.

2. **`peak_redshift` distance math**. `peak_redshift` calls
   `chi_to_z` interpolation which is a binary search + linear
   interp. Fine on CPU (called once per peak). Not worth moving to
   GPU — peak counts per tile are ≤ a few hundred.

3. **`z_pk > z_max` boundary peaks**. Small numerical differences
   between CPU and GPU in `peak_redshift` are possible because CPU
   does the lookup in Float64 while GPU currently operates in
   Float32 throughout. Decision: do the `z_pk` computation on **CPU
   in Float64**, upload `Float32(ZZon_pk)` only to the kernel. This
   keeps the filter decision deterministic vs CPU.

4. **Mask side-effects**. In CPU `ievol=1`, peaks with `z > z_max`
   never call `analyse_peak`, so their mask cells are not written.
   Pre-filtering on host preserves this behaviour exactly.

5. **`fcrit_override` interaction**. Some callers pass
   `fcrit_override`. When `ievol=1`, this override should still
   apply (same scalar for all peaks). Keep the existing kwarg
   semantics — when `fcrit_override` is given, broadcast it to the
   per-peak array regardless of ievol.

## 7. Order of operations

Recommended work order to minimise churn and keep tests green at
each step:

1. **Phase 1** (plumbing only, no behaviour change). Verify all
   existing tests pass bit-identically.
2. **Phase 2+3** (host-side assembly + filter). Still no
   behaviour change for `ievol=0`.
3. **Phase 4** (flip the gate). New tests from **Phase 5** catch
   any `ievol=1` regressions.
4. **Phase 5** (validation + bench).
5. **Phase 6** (docs).

At every phase boundary, run:
```bash
CUDA_VISIBLE_DEVICES=0 julia --project=/tmp/pp_gpu_sandbox --threads=32 \
    test/test_multiresolution_gpu.jl
julia --project=/tmp/pp_gpu_sandbox --threads=8 \
    test/test_multigpu_dispatch.jl
```

After Phase 5, also run:
```bash
CUDA_VISIBLE_DEVICES=0 julia --project=/tmp/pp_gpu_sandbox --threads=32 \
    test/test_multiresolution_gpu_ievol1.jl
```

## 8. Expected payoff

- Single-octant lightcone on 4× L40S: **~9 h → ~4.5 h** wall.
- Full 8-octant lightcone mock (8 parallel jobs): **~4.5 h** real
  time (down from ~9 h).
- Full sweep over all optimisations in the Killarney guide roadmap
  then brings lightcone octants to **~3–3.5 h** each.

This is the single highest-leverage change remaining for the
lightcone production workflow.
