#!/usr/bin/env python3
"""Generate the 8 production-v3 octant configs (2026-09-27): production-v2 plus
  * nbuff = 25 (n = 434): nhunt = min(nbuff-1, 1.75 Rf_max/a) = 24 -> R_TH cap 20.4 Mpc/h,
    M cap ~3e15 (was 12.8 Mpc/h / 7.5e14 at nbuff=16, which removed every M>1e15 halo;
    ../../paper/V2_RESULTS_2026-09-27.md). nhunt=24 is the GPU shell-kernel maximum.
  * code: GPU shell early exit (exact) and the Omega_m(a) 2LPT-coefficient fix (058207a).
  * raw catalogs go to scratch (space); AM catalogs to /project (run_catalog.slurm).
Validated on oct000 (v3test): N(>1e15, z<1) = 52 vs Tinker 51.5, no pile-up.

Original v2 description:

Differences from the frozen campaign (../production/, 2026-07-27):
  * periodic_cores = true, n = 416 (nsub = 384 = 6144/16): the tile cores tile the
    whole periodic 6144^3 box and the buffers wrap. The frozen layout (n=414, nsub=382,
    N = nsub*ntile + 2nbuff) left 32 cells per axis that were never core, i.e. an
    unsimulated 13.6 Mpc/h slab next to each octant plane
    (../../paper/FULLSKY_COMPARISON_2026-09-26.md).
  * coarse_compensation = true: coarse kernels x D/T so the multires splice is exact below
    the coarse Nyquist (removes the +2-11% P(k) bump; ../../tiling/SPLICE_COMPENSATION_2026-09-26.md).
  * code fixes: deterministic merge order, data-sized merge hash, GPU/CPU tile-local
    phi_ij zeroes only k=0 (2LPT trace identity).
Unchanged: seed 12345, N=6144, cellsize 5236/6144 = 0.852213 Mpc/h, box 5236 Mpc/h
(= 7700 Mpc), ntile=16, nbuff=16, cf=32 (coarse M = ntile*cf = 512, block 12; nsub = 32
blocks), z_max=4.5, 2LPT,
ioutshear=1, finecell filter bank. Observer bits octZYX: bit=1 -> +2618 on that axis.
The generated configs are COMMITTED — regenerate only to change the set.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
N, NTILE, NBUFF = 6144, 16, 25
NSUB = N // NTILE
NMESH = NSUB + 2 * NBUFF
CELL = 5236.0 / N
OBS = N * CELL / 2          # box corner = core corner (periodic_cores)

TEMPLATE = """# Production-v3 config — nbuff=25 rerun 2026-09-27 (production_v3/generate_configs.py).
# Octant {oct}: observer at ({cenx}, {ceny}, {cenz}) Mpc/h (centered coords) = a box corner.
# Geometry: ntile={ntile}, n={nmesh}, nbuff={nbuff}; nsub={nsub}; periodic_cores -> N = nsub*ntile = {N};
# cellsize = 5236/6144 = {cell:.6f} Mpc/h; box = 5236 Mpc/h = 7700 Mpc. Do not hand-edit.

[cosmology]
Om   = 0.31
OB   = 0.049
OL   = 0.69
h    = 0.68

[grid]
n       = {nmesh}
boxsize = {tile:.5f}
nbuff   = {nbuff}
periodic_cores = true
cenx    = {cenx}
ceny    = {ceny}
cenz    = {cenz}

[run]
z_out     = 0.0
z_max     = 4.5
ievol     = 1
ilpt      = 2
ioutshear = 1
wsmooth   = 1
rmax2rs   = 0.0
NonGauss  = 0
fNL       = 0.0
seed      = 12345
ntile     = {ntile}
coarse_factor = 32
coarse_compensation = true
generate_table = false

[files]
pk         = "validation/websky_6144/data/pk_websky.dat"
filterbank = "validation/websky_6144/data/filters_websky_finecell.dat"
homeltab   = "validation/websky_6144/data/HomelTab_websky.dat"
output     = "catalog_websky_6144_v3_oct{oct}.pksc"

[output]
format = "pksc"
path   = "/home/yguan/scratch/websky_6144/catalogs_v3"
"""

for z in "01":
    for y in "01":
        for x in "01":
            oct_ = z + y + x
            cen = {k: (OBS if b == "1" else -OBS) for k, b in (("cenx", x), ("ceny", y), ("cenz", z))}
            cen = {k: f"{v:.4f}" for k, v in cen.items()}
            path = os.path.join(HERE, f"config_v3_oct{oct_}.toml")
            with open(path, "w") as f:
                f.write(TEMPLATE.format(oct=oct_, ntile=NTILE, nmesh=NMESH, nbuff=NBUFF, nsub=NSUB,
                                        N=N, cell=CELL, tile=NMESH * CELL, **cen))
            print(f"wrote {path}  obs=({cen['cenx']},{cen['ceny']},{cen['cenz']})")
