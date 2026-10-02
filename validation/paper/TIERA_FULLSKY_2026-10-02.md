# Tier-A against the FULL Websky halo catalog (2026-10-02)

**Data and method.**
- **Websky:** the full catalog `websky_ref/halos.pksc`, downloaded 2026-10-01: 862,923,142 halos,
  10 floats each, sorted by mass, minimum 1.29e12 M☉ (the 10-particle cut). It replaces the
  single 10°×10° patch.
- **Ours:** v4fs.
- **Script:** `tierA_fullsky.jl`, job 5846126, 56 min, with 0 GPUs.
- **Results:** `results/tierA_fullsky_v4fs.txt`, log `results/tierA_fullsky_job5846126.log`.

The same **305 caps** (half-angle 5.00°, centres ≥10.5° apart and ≥5.5° from the octant planes;
32–43 per octant) are used in both skies, with the Tier-A shell (1600–2000 Mpc/h), bins and windows
of `tierA_v3fs.jl`.
- **Error:** standard error of the mean over caps, per side, combined in quadrature.
- **Independence:** caps are independent at the shell. The shells around different box corners never
  overlap, since 5236 > 2 × 2000 Mpc/h.
- **Size:** the errors are ~10× smaller than the single-patch Tier-A.

## Results (W/ours, full sky)

| Statistic | W/ours | z | Verdict |
|---|---|---|---|
| N(>M) after AM, M ≥ 3e12, z < 1.5 | 0.997–1.001 | ≤ 1.0 | agree (AM-enforced) |
| N(>M), M ≥ 3e12, z 1.5–3 | 1.003–1.010 | 2–5.5 | ours 0.3–1% low |
| N(>1e14), z < 1.5 | 1.000 / 1.005 / 1.002 | ≤ 0.9 | agree |
| N(>1.2e12), z < 0.5 | 0.935 | −28 | Websky has a 1.29e12 floor (10-particle cut); ours has no cut |
| dN/dz M>3e12, z < 2.5 | 0.996–1.006 | ≤ 2.7 | agree to < 1% |
| dN/dz M>3e12, z 2.75–4 | 1.012–1.030 | 4.7–6.6 | ours 1–3% low at high z |
| dN/dz M>1e13, z 2.25–3 | 1.010–1.017 | 2.7–3.7 | **the patch's 0.90–0.93 was patch variance** |
| σ_vr, 5 mass bins | 0.984–0.994 | −3 to −5.6 | ours 0.6–1.6% hotter |
| **ξ(r), M>5e12, mean 3–15 Mpc/h** | **0.948 ± 0.002** | **−20** | **ours ~5% more clustered** |
| ξ(r) at r = 1.2 / 1.8 / 2.7 Mpc/h | 0.894 / 0.937 / 0.968 | −12 to −21 | |
| ξ(r) at r = 5.8 / 8.6 / 12.7 / 18.8 / 27.8 / 41 Mpc/h | 0.943 / 0.911 / 0.922 / 0.956 / 0.916 / 0.897 | −3 to −27 | the offset persists to the largest r (a linear-bias offset) |
| b(M) ∝ √ξ̄(6–18), M 5e12–1.3e13 / 1.3–3.2e13 / 3.2–7.9e13 / 0.8–2.5e14 | 0.960 / 0.946 / 0.955 / 0.995 | −20 / −21 / −11 / −0.5 | ours 4–5.5% higher below 8e13, equal above |
| v12, mean 5–30 Mpc/h | 0.977 ± 0.005 | +4.8 | ours 2% stronger infall |

**Full sky, no caps, against Tinker08** (`tierA_fullsky_v4fs.txt`, last table):
- **Websky:** 0.992–1.002 for M ≥ 3e12 at all z < 4.5, and 1.000–1.001 for M > 1e13.
- **Ours:** 0.987–1.003 for M ≥ 3e12 at z < 3, and 0.97–0.99 at z > 2 for M ≥ 3e12.
- Both are abundance-matched to the same function, so this is agreement by construction.

## Conclusions

1. **The 0.958 ± 0.026 of the patch-based Tier-A was real.** With 10× smaller errors, ours is
   **~5% more clustered in ξ (≈2.5% in bias)** than Websky at M < 8e13. The offset extends from
   3 to 40 Mpc/h, i.e. it is a linear-bias offset, not only close pairs. The volume reduction
   (v4) removed the larger, mostly small-scale v3fs excess; this remainder is the
   "0.96 ± 0.03" the reduction A/B left unexplained.
   - Absolute calibration: in the z = 0.7 snapshot (`merge_ab.txt`) ours with reduction is
     0.98–1.00 × Tinker10 at low mass. Websky's patch bias was 0.93–0.99 × Tinker10
     (`tierA_diag_v3fs.txt`). Ours is the closer of the two to the N-body-calibrated bias, but
     Tinker10 itself is good to only a few percent.
   - Candidate causes, untested, remaining code-level differences in the merge
     (`CLUSTERING_EXCESS_2026-09-28.md`; as read from the local clone, not validated by the PI):
     - (a) Fortran reduces **equal-radius pairs twice**. About 22% of peaks tie in R_TH, which
       gives stronger reduction.
     - (b) Fortran halos with r³ < ΔV get a **NaN radius**. If they are dropped downstream,
       that is more neighbour removal.
     - (c) the **filter bank**: Websky's production bank was never located.
     - (d) the Fortran hash-window truncation for the smallest halos (this works the opposite
       way).

     (a) and (b) both lower Websky's low-mass bias in the observed direction and can be tested
     with the same-raw-catalog A/B (`merge_ab.jl`) by emulating them.
2. **dN/dz at z 2.25–3: the patch result is withdrawn.** It was patch variance. Full sky, the two
   agree to 1–2%, and both match Tinker to ≤1%. The 2026-10-01 statement "the difference is on
   Websky's side" (`results/dndz_wsky_vs_tinker.txt`) used the 10° patch only and is superseded.
3. **σ_vr: ours is 0.6–1.6% hotter.** This is significant now, and small. It is consistent with
   the slightly higher bias (more halos in denser, faster environments).
