# Abundance-matching top-halo fix — impact on the frozen production catalogs (2026-09-25)

## The bug (fixed in 9f92778, `build_abundance_table`)
The AM table is tabulated at mass-bin edges. Bins above a z-bin's most massive halo have
N_PP(>M) = 0 and were left at the default IDENTITY mapping (M_target = M_TH). The linear
lookup then blended the top halo(s) of every z-bin with their raw mass. Found by
`test/test_finalize_am.jl` (rank-order inversions); fixed by extending the last matched
value upward.

## Impact on the frozen campaign (octant 000; jobs 5639877, 5682203)
`check_am_tophalo_fix.jl` re-matched the raw oct000 catalog with the fixed code (output to
scratch only: `scratch/websky_6144/catalog_prod_oct000_AMfix.pksc`) and compared with the
committed `_AM` catalog:
- 319 of 195,490,586 halos change (|ΔM/M| > 1e-4), all massive (2e13–8e14); median
  ΔM/M = −8.4%, maximum +22%. N(>1e14) 62,540 → 62,538; N(>5e14) 764 → 764.
- In production the raw top masses sit ABOVE the Tinker target, so the identity node
  pulled top halos UP; the fix lowers them (the synthetic test had the opposite sign).

tSZ-weighted moments (`am_fix_tsz_moments.jl`, halos with M > 1e13), fixed / old:

| z | Σ M^{5/3} (∝ mean y) | Σ M^{10/3} (∝ 1-halo tSZ power) |
|---|---|---|
| 0.0–0.5 | 0.9933 | 0.9214 |
| 0.5–1.0 | 0.9994 | 0.9839 |
| 1.0–1.5 | 0.9999 | 0.9929 |
| 1.5–2.0 | 0.9999 | 0.9936 |
| 2.0–3.0 | 0.9999 | 0.9985 |
| 3.0–4.5 | 1.0000 | 0.9985 |
| all | 0.9982 | 0.9543 |

→ mean y −0.2% (negligible); 1-halo tSZ power ~−5% overall, ~−8% at z < 0.5 (crude
proxy: ignores profile and distance weighting, but low-ℓ tSZ power is dominated by exactly
these low-z massive clusters). κ / kSZ halo terms (∝ ~M) ≪ 1%.

## Decision needed (campaign is frozen)
Re-running AM with the fixed code for all 8 octants costs ~16 min each (AM only; raw
catalogs unchanged). Recommended before any tSZ C_ℓ enters the paper; κ, kSZ, CIB results
unaffected at the quoted precision.
