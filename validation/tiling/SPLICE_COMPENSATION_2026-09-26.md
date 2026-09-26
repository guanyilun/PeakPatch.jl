# Multires splice P(k) bump: diagnosis and D/T compensation (2026-09-26)

Follow-up to `TILING_INVARIANCE_2026-09-26.md` §(b). There, the stitched field showed
P_split/P_exact = +11% at 0.5 k_Nc and −7% at 1.2 k_Nc (δ), and ψ was +4–10% at
k = 0.06–0.25 h/Mpc with r ≈ 0.998, i.e. a coherent power boost. This does not depend
on tiling.

## Diagnosis (no free parameters)

The tile field is the Catmull–Rom (tricubic) interpolation of the coarse field, plus K∗r
for the Hoffman–Ribak residual r = n − spread_pc(block mean), where spread_pc is
piecewise-constant (`_generate_extended_residual`). Below the coarse Nyquist, per axis:

    split(k) = K(k)·ñ(k)·[1 + D(k)·(T(k) − D(k))]

- D(θ) = sin(bθ/2)/(b·sin(θ/2)) is the block-average transfer (cell-centred), θ = k·dx_fine.
- T(θ) = ⟨Σ_d w_d(t)·e^{−ibθ(f−i−d)}⟩ is the Catmull–Rom transfer, projected on the
  fine-grid mode and averaged over the b fine sub-positions.
- In 3D, D and T are the products over axes, averaged over directions.

Production block 12 (cell 0.852 Mpc/h, k_Nc = 0.307 h/Mpc):

| k [h/Mpc] | 0.05 | 0.075 | 0.12 | 0.19 | 0.29 | 0.46 |
|---|---|---|---|---|---|---|
| predicted P_split/P_exact | 1.020 | 1.043 | 1.087 | 1.111 | 0.991 | 0.876 |
| measured (δ, ntile 1, `field_N480_b12`) | 1.021 | 1.047 | 1.084 | 1.115 | 0.996 | 0.927 |

The model accounts for the bump below the Nyquist. Above it, what remains is an
aliasing-image mismatch (interpolation images against piecewise-constant images), which is
not diagonal in k.

**Consequences for the frozen campaign.**
- ψ and δ carry +2–11% of power at k = 0.05–0.25 h/Mpc. This likely explains the field-kSZ
  excess over linear theory at ℓ = 73–245 (1.14–1.19): the OV term scales as P_δ·P_v, and
  both are boosted about 8% at k ~ 0.1–0.2.
- It likely also explains part of field κ = 1.02–1.09 × linear Limber
  (`validation/paper_theory/THEORY_ANCHORS_2026-09-26.md`).
- It explains the σ(R) +3.5% bump at R ≈ one coarse cell. The mass function then absorbs it
  through abundance matching.

## Fix (opt-in: `[run] coarse_compensation = true`)

Multiply every coarse kernel that is interpolated into the tiles by C(k) = Π_axes D/T
(`MultiResolution._splice_compensation(M, block)`, passed as `comp=` to
`_periodic_convolve!`). With it, the low-k factor is 1 + D(T·D/T − D) = 1.

- **Where it applies:** δ, ψ1 and ∇²δ in `run_multitile_split`, and δ, ψ and the potential
  in `run_multitile_fieldmap`. The same coarse arrays feed the GPU tiles.
- **Range:** C lies between 0.96 and 1.31 over the band and is 1 at k=0.
- **Default off:** the frozen and held v2 runs are unaffected. It is a physics change, so
  the user decides.
- **Unit test:** `test/test_multiresolution.jl`, "Splice compensation D/T": C(0)=1, C is even,
  and it matches an independent evaluation of D/T.

## Validation

`field_split_test.jl` with `FIELD_COMP=1`, block 12, run as job 5696250, gives
`results/field_N480_b12_comp/`. Only ntile=1 ran: `sbatch --export` splits values on commas.
The effect does not depend on tiling. Same noise and region as the uncompensated
`field_N480_b12`; P_split/P_exact (cross-correlation r) below.

| k [h/Mpc] | 0.05 | 0.075 | 0.12 | 0.15 | 0.19 | 0.23 | 0.29 | 0.37 | 0.46 |
|---|---|---|---|---|---|---|---|---|---|
| δ uncompensated | 1.021 | 1.047 | 1.084 | 1.109 | 1.115 | 1.081 | 0.996 | 0.928 | 0.927 |
| **δ compensated** | **1.002** | **1.003** | **0.999** | **1.001** | **1.000** | **1.000** | **1.000** | 0.979 | 0.943 |
| ψx uncompensated | 1.024 | 1.039 | 1.071 | 1.098 | 1.106 | 1.080 | 1.004 | 0.959 | 0.969 |
| **ψx compensated** | 1.007 | 0.999 | 0.995 | 1.001 | 0.996 | 1.000 | 0.998 | 1.000 | 0.980 |

- **δ, cross-correlation:** r = 1.0000 below the coarse Nyquist, against 0.9990–0.9998 before.
- **δ, rms error:** 4.50% of σ, against 5.10%.
- **ψx, rms error:** 15.2%, unchanged. It is dominated by the near-Nyquist decoherence
  (r = 0.81 at 1.2 k_Nc), the image effect, which a diagonal factor cannot fix.
- **σ(R) (top-hat), stitched/exact:**

  | R [Mpc/h] | 1.4 | 2.1 | 3.3 | 4.9 | 7.5 | 11.4 | 17.4 | 26.5 |
  |---|---|---|---|---|---|---|---|---|
  | compensated | 0.995 | 0.994 | 0.993 | 0.995 | 0.999 | 1.000 | 1.001 | 1.001 |

  Before, it was 0.999–1.000 at R ≤ 2, rising to 1.035 at 11.4. The +3.5% bump is gone;
  a −0.5 to −0.7% deficit remains at R ≲ 5 from the leftover image deficit just above the
  Nyquist. Abundance matching absorbs it.
- **Remaining:** above the coarse Nyquist, δ is 0.94–0.98 at k = 0.37–0.6 h/Mpc (before:
  0.92–0.95). Removing it would need consistent interpolation and residual images, for
  example spreading the residual with the same kernel. That is a larger redesign.

**Verdict:** compensation makes the multires field exact to ≤0.5% in power (r = 1.0000)
below the coarse Nyquist, the range that carries the kSZ and velocity signal. It costs
nothing at run time.

## Side finding (no effect on outputs)

`run_multitile_split` computes a coarse ψ2 (`psi2_coarse`) through `_periodic_convolve!`,
which multiplies the 2LPT *source* by √P (wrong units). The array is never used: ψ2 is
computed per tile from the full tile δ. It is dead code, so nothing changes. Note that the
per-tile periodic 2LPT solve on a non-periodic tile misses ψ2 modes longer than the tile,
an approximation to state in the paper (ψ2 is subdominant on those scales).
