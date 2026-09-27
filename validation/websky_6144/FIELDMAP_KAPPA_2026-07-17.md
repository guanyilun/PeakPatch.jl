# Field-matter lightcone κ map: Phase-A exit result (2026-07-17)

Full-matter (halo+field, all lattice cells) Born κ map of octant 000, finecell config
(box 5236 Mpc/h, cellsize 0.85221 Mpc/h), Nside 2048, `run_multitile_fieldmap` on
4× L40S (job 4297504, 2h10m). Compared to Websky `kap.fits` and CAMB-verified Limber
linear C_ℓ with `compare_fieldmap_kappa.jl` (6 flat-sky caps inside the octant).

## The 8× bug and its fix (job 4280770 → 4297504)

The first octant run (job 4280770) painted exactly 2^(n_workers−1) = 8× too much mass:
the worker-task accumulator name `my_acc` was also assigned in the single-worker `else`
branch of `run_multitile_fieldmap`; `if/else` introduces no scope in Julia, so all four
`Threads.@spawn` closures captured ONE shared boxed variable — every worker painted into
the same map and the reduction `total = worker_maps[1]; total .+= worker_maps[wid]`
self-added the aliased array three times. Diagnosed entirely on the login node: the
painter replays at ratio 1.000000 per tile against an exact geometric cell count, and the
job ratio was 7.98882 ≈ 8. Fixed in 668bee0 (`local wacc`, alias-checked copy-reduce,
`cpu_workers` kwarg so the smoke test covers multi-worker dispatch without a GPU).

After the fix (job 4297504):

- painted mass = **6.371121e21 Msun/h = the independent exact geometric count to all
  printed digits** (119,659,142,676 cells × ρ̄·a_latt³)
- mean κ (full sky, unsubtracted) 0.1663 vs analytic (∫W_κ dχ)/8 ≈ 0.170
- sky fraction 0.1247 (octant = 0.125)
- smoke: cpu_workers=4 vs 1 = 1.00000000; gpu_paint vs CPU = 1.00000000, 0 pixels differ

## Phase-A exit C_ℓ (6-cap mean)

| ℓ    | ours/CAMB-linear | ours/kap.fits | kap.fits/linear |
|------|------------------|---------------|-----------------|
| 231  | 0.95             | 0.735         | 1.32            |
| 452  | 0.94             | 0.641         | 1.47            |
| 884  | 0.95             | 0.536         | 1.78            |
| 1730 | 0.94             | 0.389         | 2.42            |
| 3383 | 0.77             | 0.187         | 4.11            |

**Interpretation: the field map tracks CAMB linear theory to ~5% over ℓ ≈ 230–1730** —
exactly what a 2LPT-displaced uniform-mass lattice should produce. The growing deficit
vs kap.fits at high ℓ is the *expected* missing nonlinear/1-halo power: in the Websky
construction that power comes from replacing collapsed regions with NFW-painted resolved
halos (our XGPaint-fork halo painting already matches kap.fits' 1-halo term to ~1%, see
KAPPA_PAINTED_2026-07-16.md). kap.fits/linear rising 1.3 → 4.1 over ℓ 230–3400 matches
the anchor numbers recorded when the Limber code was CAMB-verified.

## Conclusion & next step

Phase A is DONE: deterministic re-paint, exact mass bookkeeping, GPU bit-exact painting,
and linear-level κ power validated. Full kap.fits reproduction = halo κ (XGPaint NFW,
done) + field κ with halos excluded + compensation (Phase C `exclude_halos`) + z>4.5
Gaussian tail. Phase C also unlocks the kSZ field component (v_r kernel).
