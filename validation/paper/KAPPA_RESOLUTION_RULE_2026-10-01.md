# κ: does Websky's "resolved halo" rule explain our κ offset? (2026-10-01)

**Background.** Our full-sky κ is 1.06–1.13× the released `kap_lt4.5` at ℓ = 300–3000. Ours is
1.00–1.05× Halofit there, while Websky is 0.94–0.99×. Stein+2020 §3.2.3 (`2001.08787.txt:933-938`):
"halos whose virial radii subtend a solid angle less than twice that of a pixel in the map are
considered unresolved and are included in the field component". Our construction B instead paints
*every* halo as a compensated (zero-net) profile on top of the full-matter field.

**Test.** `paint_octant.jl` with `KAPPA_RESOLVED_ONLY=1` and κ only. Halos with πθ(r200m)² < 2 Ω_pix
at Nside 4096 (θ < 0.687′) are not painted; their matter stays in the field map. This is
"v4kres": v4fs catalogs and v3 field maps.
- Jobs 5845278–82 and 5846121–23 (paint), 5846124 (assembly, auto vs `kap_lt4.5`, theory).
- 100.5M of the 112.5M halos per octant are unresolved; 11.95M are painted.
- Results: `results/auto_v4kres_kappa_lt4.5.txt`, `results/theory_v4kres_kappa.txt`, log
  `results/kres_analysis.log`.

| ℓ | ours/W, construction B | ours/W, resolved-only | ours/Halofit B | ours/Halofit resolved-only | W/Halofit |
|---|---|---|---|---|---|
| 104 | 1.192 | 1.191 | — | — | — |
| 297 | 1.073 | 1.072 | 1.022 | 1.021 | 0.953 |
| 771 | 1.098 | 1.092 | 1.027 | 1.021 | 0.936 |
| 1026 | 1.099 | 1.088 | 1.032 | 1.022 | 0.939 |
| 1366 | 1.114 | 1.094 | 1.048 | 1.030 | 0.941 |
| 1999 | 1.111 | 1.070 | — | — | — |
| 2926 / 2894 | 1.134 | **1.049** | 1.134 | **1.050** | 0.999 |
| 5185 | 1.073 | 0.830 | — | — | — |
| 8071 | 0.998 | 0.597 | — | — | — |

The band average over ℓ = 30–1700 is 1.085 ± 0.0014 for resolved-only, against 1.099 for B; 19/43
bands are >4σ in both.

## Conclusion

- **ℓ ≈ 1500–3000:** the resolution rule accounts for most of the high-ℓ excess (ℓ 2900:
  1.13 → 1.05). That part of the difference is construction (resolved vs unresolved halos), not
  catalog.
- **ℓ > 4000:** with the rule, ours falls *below* Websky (0.6–0.8). Websky keeps more small-scale
  power, plausibly from the field sub-volume splitting described in §3.1.3 (lattice sites split into
  up to 5³ points so that they are smaller than a pixel). Our field painter deposits each lattice cell
  as one point.
- **ℓ < 1000:** the rule changes nothing. We stay 7–9% above Websky and within 2–3% of Halofit, while
  Websky is 5–6% *below* Halofit. That large-scale part must come from how Websky builds its field
  component ("matter exterior to resolved halos", with a double-counting correction like the kSZ
  electrons). The released maps do not separate it, so it cannot be tested from them.
- **Reading for the paper:** our κ is the one that agrees with the theory anchor (Halofit) at
  ℓ < 1650. The ℓ < 1000 offset against the released map is a reference-side construction difference
  (neutral framing; not validated by the PI). The ℓ ≳ 1500 part is the resolution rule, which we can
  reproduce on request (`KAPPA_RESOLVED_ONLY=1`).
