# Multires split: ψ aliasing error — literature and fix design (2026-10-03)

Context: MATCHED_FORTRAN_2026-10.md §8.

- **Problem.** At the production block (b = 12), the split ψ1 has ~30% error power just above the coarse
  Nyquist k_N, and 14% rms error. δ is fine (≤1.2%).
- **Consequence.** The R_TH, positions and halo clustering change (+3.7% ξ).
- **Status.** Review by a delegated agent, from paper text. Items marked [S] were checked from summaries only.

## What the literature says

**GRAFIC2 (Bertschinger 2001, astro-ph/0103301).** This is the same noise split as ours: block-mean coarse
noise spread over the block, plus fine (ξ − ξ̄) by isolated FFT.

- The images above k_N are real signal of the block-mean noise.
- They are reproduced exactly only if the spread coarse field is convolved with a *field-specific*
  anti-aliasing filter, W = T(k)/T(k₀), with k₀ = k folded into the first Brillouin zone (eqs. 26–28).
- Our D/T compensation is only the diagonal (m = 0) term.
- W_ψ is long-ranged, so with truncation the velocity errors stay at 3–7% rms (§3.3–3.4). Only larger buffers
  help.
- **This confirms our diagnosis, and shows that "fixing the images" in a noise-level ψ split tops out at a few
  percent.**

**Our numbers fit the mechanism.**

- An interpolator tuned for δ makes ψ images too large by k/k₀, because T_ψ = T_δ/k. At k = 1.5 k_N that
  factor is 3.
- The block mean carries sinc²(3π/4) ≈ 9% of the noise there.
- The predicted error power is ≈ (3 − 1)²·9% ≈ 35%. We measured ~30%.

**MUSIC (Hahn & Abel 2011, 1103.6031).**

- The split is by region, not by scale.
- Displacements and velocities are **never** made by convolving noise with a ψ kernel. They come from the
  accurate δ through a Poisson solve (multigrid, with flux-matched coarse–fine boundaries).
- Velocity errors ~1.5e-3 σ_v, about 100× better than GRAFIC2 on the same setup.

**Pen 1997 (astro-ph/9709261).** The kernel is split into a long part and a compact short part (r < R). The
short part acts on the full fine noise, and the buffer equals R.

**Panphasia (Jenkins 2013, 1306.5968) [S].**

- Block-mean bases give large-scale errors ∝ k². Adding gradient moments gives ∝ k⁴.
- Convolved coarse fields are never interpolated.

**TreePM [S].**

- Gaussian long/short split: long part e^{−k² r_s²} on the mesh, short part erfc in real space.
- Typical choice r_s ≈ 1–1.25 mesh cells, short-range cutoff 4.5–5 r_s, force errors ≲ 1% rms, 1–2% max.
- At r_s ≥ Δ the long field at k_N is ~1e-7, so interpolation images are irrelevant.

**monofonIC, Websky / Peak Patch.**

- Both use global FFTs; there is no precedent there for avoiding them.
- monofonIC de-aliases 2LPT products with Orszag 3/2 padding.

## Design (to prototype)

1. **Long/short split by scale, not by noise level.**
   - The coarse level carries only G_L(k) = e^{−k² r_s²} times each field: δ_L, φ_L, ψ_L = ∇φ_L, all from the
     same coarse noise.
   - With r_s ~ Δ_coarse they are negligible at k_N, so the interpolation images vanish.
   - Interpolate φ_L, or all three fields, consistently so that δ_L and ψ_L stay a Poisson pair.
2. **Short part from the tile δ (MUSIC principle).**
   - δ_S = δ_tile − δ_L.
   - ψ_S comes from the isolated (2×-padded) tile solve of ∇²φ_S = δ_S, then ψ_S = −∇φ_S (sign as in our
     LPT convention).
   - The tile δ is already accurate to ~1%, so ψ inherits it.
3. **Gaussian handoff, not sharp.** A sharp-k short kernel has ~1/r² oscillatory tails that a finite buffer
   cannot contain.
4. **Scale vs buffer, the key trade-off.**
   - TreePM practice r_s ≈ Δ_coarse = 12 fine cells needs a cutoff of ~55–60 cells, far more than our
     25-cell buffer.
   - Options:
     - (a) smaller block b = 6–8, with r_s ≈ Δ and cutoff ~30 cells;
     - (b) b = 12 with buffer ~60;
     - (c) r_s ≈ 0.75 Δ (≈ 9 cells; e^{−(0.75π)²} ≈ 4e-3 at k_N), with cutoff ~40.
   - To be settled by measurement.
5. **2LPT.**
   - Build the total φ,ij on the tile and form the products with 3/2 de-aliasing.
   - Solve for ψ2 on the tile. Long modes of the source from short-mode coupling are then missing (ψ2 is
     subdominant; to be measured).

## Acceptance tests (same noise, same field)

- `matched/split_field_error.jl`: ψ1 error power in the band above k_N ≲ 1e-3 (now 0.30); ψ1 rms ≲ 1–2%
  (now 14%); δ no worse than now. Scan r_s, block size and buffer.
- `matched/split_vs_exact.jl` / `matched_split_scan.jl`: split/exact ξ(3–15) and b_E within ~0.5% of 1
  (now 1.037 / 1.017); peaks in the same cell ≳ 99%; Eulerian position error ≪ 0.8 Mpc/h.
- Only then: cost on GPU, and a decision on rerunning production.
