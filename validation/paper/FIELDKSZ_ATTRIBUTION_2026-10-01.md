# Field-kSZ excess v2 → v3: it was the tile buffer, not the Ω_m(a) fix (2026-10-01)

**Question** (report §7 item 1). Field kSZ against linear theory was 1.13–1.21 at ℓ = 70–200 in v2
and 1.01–1.10 in v3. Two changes coincide: the Ω_m(a) 2LPT fix (058207a) and the tile layout
(n 416 → 434, nbuff 16 → 25).

**Test** (`fieldksz_attribution.jl`, job 5845337, 79 min). We regenerated the oct000 field map with
the **v2 layout** (`production_v2/config_v2_oct000.toml`, nbuff 16) and the **current code** (fix in),
then compared it on the octant-000 apodized mask with the v3 oct000 field map (nbuff 25, fix in).
Results: `results/fieldksz_attribution.txt`, log `results/fieldksz_attribution_job5845337.log`.

| ℓ | kSZ, nbuff 16 / nbuff 25 | full-sky v2 / v3 (`theory_v{2,3}_ksz_field.txt`) | κ, nbuff 16 / 25 |
|---|---|---|---|
| 71–104 | 1.07–1.18 | 1.08–1.16 | 1.000–1.001 |
| 114–185 | 1.11–1.14 | 1.10–1.13 | 1.00 |
| 203–297 | 1.06–1.10 | 1.09–1.13 | 1.00–1.01 |
| 434–636 | 1.02–1.05 | 1.03–1.04 | 1.00–1.01 |
| 771–989 | 1.01–1.02 | 1.01–1.02 | 1.00 |

**Conclusion.**
- Changing the layout alone, with the fix present in both maps, reproduces the whole v2 → v3 change.
  **The tile buffer (nbuff 16 → 25) removed the field-kSZ excess; the Ω_m(a) fix did not**, as
  expected from its 0.15–0.34% velocity effect.
- κ is unchanged, so density is not affected, only the velocity / displacement field.
- Mechanism (interpretation): ψ ∝ δ/k² is non-local. The tile-local residual convolution truncates
  it at the buffer edge, which adds spurious velocity power on tile scales (~350 Mpc/h ↔ ℓ ≈ 100–300
  at χ ≈ 2–3 Gpc/h). A larger buffer reduces this.

**Implication / next test.** v3/v4 field kSZ is still 1.03–1.10 × linear theory at ℓ = 65–200.
Part of that may be the same artifact, not yet converged. Deciding test: one oct000 field map at
nbuff ≈ 40 (n ≈ 464) against nbuff 25. This is also an A8 convergence item.

## Follow-up 2026-10-02: nbuff 40 vs 25 — not a monotonic convergence (job 5848321)

The test was one oct000 field map at nbuff 40 (n = 464; `configs/config_nb40_oct000.toml`), on the
same mask, as A/B against the v3 nbuff-25 map (`results/fieldksz_nbuff40_vs_25.txt`, log
`results/fieldksz_nbuff40_job5848321.log`):

| ℓ | kSZ nbuff40/nbuff25 | kSZ nbuff16/nbuff25 (above) | κ nbuff40/25 |
|---|---|---|---|
| 65–104 | 1.04–1.15 | 1.02–1.18 | 1.000 |
| 114–185 | 1.02–1.09 | 1.11–1.14 | 1.00 |
| 203–359 | 0.96–1.01 | 1.06–1.10 | 1.00–1.01 |
| ≥ 395 | 0.99–1.00 | 1.01–1.05 | 1.00 |

**The buffer effect is not a monotonic convergence.** nbuff 16 *and* nbuff 40 both give 5–15% more
field-kSZ power at ℓ = 65–190 than nbuff 25, and κ is identical in all three. Therefore:
- the v2 → v3 change is a tile-layout effect, as concluded above. But it is **not** "a larger buffer
  reduces a truncation error", and nbuff 25 being closest to linear theory looks fortuitous;
- **the low-ℓ field kSZ carries a ~±10% numerical systematic that depends on the tile layout**
  (n, nbuff), at ℓ ≈ 65–200, i.e. k ≈ 0.02–0.1 h/Mpc at χ ≈ 2–3 Gpc/h. At those k the velocity should
  come from the coarse grid, so the layout dependence points to the tile-local residual displacement
  (ψ ∝ δ/k², isolated convolution) leaking power to large scales.
- **Next test:** P_ψ(k) and P_v(k) at k < 0.1 h/Mpc, tiled against exact global FFT, for several
  (n, nbuff), on a small box. `validation/tiling/field_split_test.jl` already measures the δ/ψ split
  accuracy; extend it to the velocity power at low k.
- Until resolved, the field-kSZ comparison at ℓ < 200 should carry a ±10% systematic.
