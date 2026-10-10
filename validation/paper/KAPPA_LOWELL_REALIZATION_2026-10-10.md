# κ at low ℓ: the field excess over linear theory is our realization (2026-10-10)

**Question.** In v5 (Huntian v0.1), the field κ map is 3–12% above Limber linear theory per bin at
ℓ ≈ 60–115 (band 1.06), while it is within 1–2% of linear at ℓ = 150–600. Is this the pipeline or the sky?

**Test 1 (mode amplitudes, job 7213141, `kappa_realization_check.jl`).** The realized |δ_k|² of the production
coarse noise (seed 12345, M = 512) is within its mode-count scatter of 1 at every k. Projected with linear Limber
this gives κ factors 1.002 / 0.999 / 1.000 / 1.000 for ℓ 50–150 / 150–300 / 300–600 / 600–1000. Amplitudes alone
do not explain the excess, but this test ignores phases (cross terms between 3D modes on the sphere).

**Test 2 (map level, job 7213274, `kappa_linear_toy.jl`).** A linear Born κ map is ray-traced (Nside 256) from
the production coarse noise with the observer at the box corner, and from 8 independent seeds on the same grid
and with the same code. R_ℓ = C_ℓ(seed 12345) / ⟨C_ℓ(other seeds)⟩ is our realization's linear-sky factor.
Results: `results/kappa_linear_toy.txt`.

| band | R (toy, our seed) | field/lin (production) | corr(toy, field) |
|------|------|------|------|
| ℓ 50–150 | 1.029 | 1.028 | 0.971 |
| ℓ 60–115 | 1.057 | 1.062 | 0.978 |
| ℓ 150–300 | 1.010 | 0.996 | 0.921 |

R tracks field/lin bin by bin, including the fluctuations at ℓ < 50 (e.g. ℓ 27: 1.182 vs 1.185; ℓ 30:
0.826 vs 0.792; ℓ 103: 1.080 vs 1.108). The toy and the production field map are the same sky (r ≈ 0.99 at
ℓ < 50, 0.97 at ℓ ≈ 100, falling with ℓ as the toy's 10 Mpc/h grid and the nonlinear field diverge).

## Conclusion

- **The low-ℓ excess is the sky variance of our one realization, not a pipeline error.** One periodic box seen
  by 8 octants is one sky; seed 12345 happens to be ~3% high at ℓ 50–150 and ~6% high at ℓ 60–115. The
  Gaussian per-band estimate (~1%) understates this variance; the map-level test is the right one.
- **ℓ 150–600:** field κ is 0.99–1.03 × linear theory; the IC P(k) and χ* match the theory table to 0.6%.
  Growth, kernel and P(k) in the field projection are validated at the 1–2% level.
- **ℓ 300–1000 vs Websky (ours/W 1.06–1.11) is a separate, unchanged item:** ours total is 1.01–1.05 ×
  Halofit, Websky `kap_lt4.5` is 0.94–0.95 × Halofit (KAPPA_RESOLUTION_RULE_2026-10-01.md). Our realization factor
  there is ~1.00 (test 1), so it is not the sky. It remains a reference-side construction difference.
- For low-ℓ κ comparisons against Websky or theory, quote the realization factor R_ℓ (or the ratio to the
  toy) rather than the Gaussian error.
