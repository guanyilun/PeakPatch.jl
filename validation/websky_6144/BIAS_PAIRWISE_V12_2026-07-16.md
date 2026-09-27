# Halo bias b(M) and pairwise velocity v12(r) vs Websky (2026-07-16)

Tier-A catalog-level tests from the Stein+2020 comparison plan (designed 2026-06-29):
these are the direct drivers of the 2-halo map amplitude (bias) and the pairwise-kSZ
signal (v12), so they predict whether painted-map spectra can match before any painting.

Catalog: `catalog_websky_6144_oct000_finecell_AM.pksc` (oct000, finecell 0.852 Mpc/h,
abundance-matched; Eulerian positions + velocities reconstructed at read time from the
old-convention fields). Reference: Websky public 10°×10° patch. Matched geometry: same
inscribed angular cap (5.01°) for both, thin shell r=1600–2000 Mpc/h (z≈0.73).

## 1. Halo bias b(M) — `compare_bias.jl`

Landy-Szalay ξ_hh per mass bin, window 6–18 Mpc/h; absolute b via linear ξ_mm from the
CAMB P(k) (σ₈ check: 0.8101), D(0.727)=0.689; Tinker 2010 (Δ=200) overlay.

| M_bin [M☉/h] | N_our | N_wsk | b_ours | b_wsky | **b_rel** | b_Tinker | ours/Tinker |
|---|---|---|---|---|---|---|---|
| 7.94e12 | 14715 | 14780 | 1.653 | 1.615 | **1.024** | 1.708 | 0.968 |
| 2.00e13 | 4978 | 5058 | 2.107 | 2.095 | **1.006** | 2.217 | 0.950 |
| 5.01e13 | 1366 | 1419 | 2.726 | 2.782 | **0.980** | 3.031 | 0.899 |
| 1.41e14 | 306 | 274 | 3.999 | 4.202 | **0.952** | 4.614 | 0.867 |

**b(M) matches Websky to ≤5% in every bin** (top bin has only ~300 halos — noise-dominated).
Per-bin counts also agree to ~1–3% (the mass function again, per-bin). Both catalogs sit
~5–13% below analytic Tinker *by the same amount* (b_wsky/b_Tink = 0.95, 0.95, 0.92, 0.91)
— a shared peak-patch-vs-Tinker offset, not a discrepancy between the catalogs.

## 2. Mean pairwise velocity v12(r) — `compare_pairwise_velocity.jl`

v12(r) = ⟨(v_i − v_j)·r̂_ij⟩, M>1e13, same matched cap/shell. n_ours=9026, n_wsky=9077.
Sanity σ_vr: ours 267.5 vs Websky 248.0 km/s (ratio 1.079 in this small cap).

| r [Mpc/h] | v12_ours | v12_wsky | ratio |
|---|---|---|---|
| 3.9 | −224.1 | −231.6 | 0.968 |
| 7.8 | −253.9 | −244.1 | 1.040 |
| 11.7 | −238.9 | −233.1 | 1.025 |
| 15.5 | −214.3 | −206.9 | 1.036 |
| 19.4 | −182.9 | −182.7 | 1.001 |
| 23.3 | −154.6 | −163.6 | 0.945 |
| 27.1 | −136.6 | −145.7 | 0.937 |
| 31.0 | −112.2 | −124.5 | 0.901 |
| 38.7 | −76.2 | −96.8 | 0.788 |
| 46.5 | −47.5 | −80.5 | 0.590 |
| 58.1 | −14.4 | −58.2 | 0.247 |

**Mean ratio over the kSZ-relevant window (5–30 Mpc/h) = 0.997.** The correlated infall
that sources the pairwise-kSZ signal reproduces Websky essentially exactly where it matters.

**Caveat — large-r tail (r ≳ 35 Mpc/h):** our v12 decays faster than Websky's (ratio
0.79 at 39 → 0.25 at 58 Mpc/h). Both catalogs see identical geometry/estimator, so this
is a real difference in large-scale velocity coherence, not an artifact of the comparison.
Plausible suspects: sample variance in a single 5° cap × 400 Mpc/h shell (v12 is small,
tens of km/s, and cap-scale bulk modes don't fully cancel); or genuine attenuation of
long-wavelength velocity modes in our tiled coarse/fine reconstruction. Not investigated
further — it is outside the pairwise-kSZ window and does not affect the Tier-A verdict,
but worth revisiting with more caps/octants if large-scale velocity statistics ever matter.

## Verdict

Both remaining Tier-A statistics pass. Together with the mass function (~1%), dN/dz,
σ_vr (~1–2%), and ξ(r) (bias ratio 1.06), the catalog now matches Websky in every
statistic that drives the painted maps. Next step per the plan: Tier-B painting,
starting with lensing κ C_ℓ (pure mass projection, no gas physics).
