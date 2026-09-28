# Production-v3 campaign (nbuff=25) — submitted 2026-09-27

This is production-v2 (`../production_v2/README.md`) plus:

1. **nbuff = 25, n = 434.** nhunt = min(nbuff−1, 1.75·Rf_max/a) = 24, which raises the R_TH
   cap from 12.8 to 20.4 Mpc/h and the M cap from 7.5e14 to about 3e15.
   - At nbuff=16 no halo above 1e15 was found anywhere, and that is what drove the tSZ C_ℓ
     deficit.
   - The oct000 "v3test" gives N(>1e15, z<1) = 52 against Tinker's 51.5, with no pile-up, and
     tSZ at ℓ = 300–1400 goes from 0.90–0.94 to 1.06–1.21 (`../../paper/V2_RESULTS_2026-09-27.md`).
   - nhunt = 24 is the largest the GPU shell kernel allows.
2. **Code** (as of 058207a):
   - The GPU shell early exit, which is exact and makes nbuff=25 cost about the same as v2.
   - **The Ω_m(a) fix in the 2LPT coefficient** (a³ → a⁻³; ψ2 coefficient was 1–4% high).
     It also affects the field maps, so these are regenerated as well.
3. **Storage:** the raw catalogs go to **scratch**,
   `/home/yguan/scratch/websky_6144/catalogs_v3/`, which is purgeable; copy them before any
   purge. The AM catalogs and all maps go to /project. This fits the quota without deleting
   anything.

All other settings are unchanged: seed 12345, N=6144, box 5236 Mpc/h, cell 0.852213,
ntile 16, periodic_cores, cf 32 plus coarse_compensation, z_max 4.5, 2LPT, ioutshear 1,
and the finecell filter bank.

## Jobs (`bash submit_all.sh all`, submitted 2026-09-27; log `scratch/.../logs/v3_submit.txt`)

| oct | catalog+AM | paint | CIB | fieldmap |
|---|---|---|---|---|
| 000 | 5716328 | 5716329 | 5716330 | 5716331 |
| 001 | 5716332 | 5716333 | 5716334 | 5716335 |
| 010 | 5716336 | 5716337 | 5716338 | 5716339 |
| 011 | 5716340 | 5716341 | 5716342 | 5716343 |
| 100 | 5716344 | 5716345 | 5716346 | 5716347 |
| 101 | 5716348 | 5716349 | 5716350 | 5716351 |
| 110 | 5716352 | 5716353 | 5716354 | 5716355 |
| 111 | 5716356 | 5716357 | 5716358 | 5716359 |

Outputs:
- AM catalogs: `/project/.../websky/catalog_websky_6144_v3_octZYX_AM.pksc`
- maps: `fieldmaps_v3/`, `halomaps_v3/`, `cibmaps_v3/`

**The jobs run from the live working tree.** Do not edit `src/` until every catalog and
fieldmap job has started. Better still, wait until they have all finished, because Julia
loads the code when each job starts.

**Afterwards:**
- `CAMPAIGN=v3` full-sky assembly and spectra, and the tail check (`check_am_tail.jl`) on
  every octant.
- Rerun the Tier-A catalog statistics on v3.
- Compare v3 oct000 with v3test oct000 to size the Ω_m(a) fix.
