# GPU shell early exit (2026-09-27)

**Why.** Lifting the halo mass cap needs nbuff=25, which gives nhunt=24 cells, a
20.4 Mpc/h R_TH cap and about 3e15 Msun/h (`validation/paper/V2_RESULTS_2026-09-27.md`).
But `_shell_gather_full_kernel!` accumulated every shell up to nhunt for every peak. At
nhunt 24 against 15 that is about 4× the cells per peak, and the oct000 test ran at
~20 s per tile against ~3 s in v2 (≈13 h per octant).

**What.** Thread 1 of each peak's block replays `_post_process_kernel!`'s Fbar recurrence
shell by shell, with the same reduced sums and the same Float32 expression. It also
replays the m0 logic: m0 is the first shell with r² > ir2min = (1.75 R_f/a)². At the
first shell s ≥ m0 with Fbar < fcrit·(1 − 1e-4), it sets r_stop = rad_s + 3, and the
gather stops once rad ≥ r_stop.

**Why it is exact.** Post-processing finds the first outward crossing from m0. The
strain, gradient and virialization stencils reach at most 2 cells beyond the crossing
(or beyond m0 for inward crossings). The filter-scale gradient reaches R_f + 2 < m0,
and the formation redshift is taken inside R_TH/2. The −1e-4 margin means that any
rounding difference moves the detected crossing *outward*, never inward. Unvisited
shells keep their zero initialisation, which does not change the first crossing.

- **Where:** production multi-R_f path only, and only when rmax2rs = 0.
- **Switch:** `PEAKPATCH_SHELL_EARLY_EXIT=0`, or `CUDAExt._SHELL_EARLY_EXIT[]`, turns it off.

**Validation.**
- **Example octant** (examples/config_gpu_octant.toml; job 5708923): with the exit off
  and on, the sorted catalogs are identical. That is 0 of 610,713 records differing at
  nbuff 16 and 0 of 610,229 at nbuff 25.
- **nbuff 25 in the same box:** max M rises to 1.9e15 (R_TH 17.4 Mpc/h), from
  7.4e14 capped. The mass cap is gone.
- **New test** `test_multiresolution_gpu.jl` "GPU shell early exit is exact": the
  sorted records are `isequal`. The first version used `==` and falsely failed on NaN
  fields; job 5711281 passes.
- **Speed:** the example is too small to show it (wall time is dominated by
  compilation). It is measured by the oct000 nbuff=25 test, job 5710456.

**Side finding.** `gradrf_x/y/z` (catalog fields 31–33) are NaN for about 0.005% of
production halos: 26 in the first 500k of v2 oct000, and many in small test boxes. Masses
and positions are unaffected. Still to trace (probably a zero-radius mrf shell).
