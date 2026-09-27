#!/usr/bin/env python3
"""
Generate TOML configs for all 8 Websky 6144^3 lightcone octants.

Each octant places the observer at a different corner of the box (5236 Mpc/h = 7700 Mpc).
The naming convention (000..111) is a binary flag: bit set means the observer
sits at the far end of that axis (x=bit0, y=bit1, z=bit2).

IMPORTANT: run_multitile_split uses centered tile coordinates
(tile_center returns (it - (ntile+1)/2) * dcore_box), so the observer
position must be expressed in the same centered coordinate system.
The box extends from -boxsize_full/2 to +boxsize_full/2 in centered coords.

UNITS (see CONVENTIONS.md): the pipeline is Mpc/h. Websky's "(7.7 Gpc)^3" is 7700 *Mpc*
(physical comoving), so the box in Mpc/h is 7700*h = 5236 Mpc/h. The pre-2026-06-15 version
of this file set BOX_FULL = 7700 directly (interpreting Mpc as Mpc/h) -> cells 1/h = 1.47x too
coarse, mass floor ~3e12 instead of Websky's ~1e12. FIXED below: BOX_FULL = BOX_FULL_MPC * h.

Geometry (nbuff=16, factor-h corrected):
  ntile=16, nmesh=414, nbuff=16 ;  nsub = 382 ;  N = 382*16 + 32 = 6144
  cellsize = 5236/6144 = 0.85221 Mpc/h  (= 7700/6144 Mpc, matches Websky)
  per-tile boxsize = 414 * 0.85221 = 352.8164 Mpc/h

Usage:
    python generate_octant_configs.py
"""
import os

# UNITS: Mpc/h pipeline; Websky (7.7 Gpc)^3 = 7700 Mpc -> box_[Mpc/h] = 7700*h. See CONVENTIONS.md.
H = 0.68
BOX_FULL_MPC = 7700.0                  # Websky physical comoving box, Mpc
BOX_FULL = BOX_FULL_MPC * H            # = 5236.0 Mpc/h  <- value the pipeline needs
HALF_BOX = BOX_FULL / 2.0              # = 2618.0 Mpc/h (octant observer offset)
OUTPUT_PATH = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"

# Tile geometry with nbuff=16
NTILE = 16
NBUFF = 16
CELLSIZE = BOX_FULL / 6144.0          # = 0.85221 Mpc/h (= 7700/6144 Mpc, matches Websky)
NMESH = (6144 - 2 * NBUFF) // NTILE + 2 * NBUFF  # = 414
BOXSIZE = round(NMESH * CELLSIZE, 4)  # = 352.8164 Mpc/h (per-tile)

# Octant names: binary flags (bit0=x, bit1=y, bit2=z)
# 000 = origin corner, 111 = far corner, etc.
# In centered coords: bit=0 -> -HALF_BOX, bit=1 -> +HALF_BOX (= -/+2618 Mpc/h, factor-h corrected)
OCTANT_NAMES = ["000", "001", "010", "011", "100", "101", "110", "111"]

for oct_name in OCTANT_NAMES:
    bits = int(oct_name, 2)
    cenx = HALF_BOX if (bits & 1) else -HALF_BOX
    ceny = HALF_BOX if (bits & 2) else -HALF_BOX
    cenz = HALF_BOX if (bits & 4) else -HALF_BOX

    content = f"""# Websky 6144^3 GPU Multi-Resolution — Octant {oct_name}
#
# Killarney Standard Compute (4x L40S 48 GB)
# Observer at centered coord ({cenx}, {ceny}, {cenz}) Mpc/h
# (= box corner in physical coords)
#
# Geometry:
#   ntile={NTILE}, nmesh={NMESH}, nbuff={NBUFF}
#   nsub = {NMESH} - {2*NBUFF} = {NMESH - 2*NBUFF}
#   N = {NMESH - 2*NBUFF}*{NTILE} + {2*NBUFF} = 6144
#   cellsize = {BOX_FULL:.1f}/6144 = {CELLSIZE:.5f} Mpc/h  (= 7700/6144 Mpc, MATCHES Websky)
#   per-tile boxsize = {NMESH} * {CELLSIZE:.5f} = {BOXSIZE} Mpc/h
#   total box = {BOX_FULL:.1f} Mpc/h = {BOX_FULL_MPC:.0f} Mpc  (Websky (7.7 Gpc)^3)
#   UNITS: pipeline is Mpc/h; box = 7700 Mpc * h. See CONVENTIONS.md.

[cosmology]
Om   = 0.31
OB   = 0.049
OL   = 0.69
h    = 0.68

[grid]
n       = {NMESH}
boxsize = {BOXSIZE}
nbuff   = {NBUFF}
cenx    = {cenx}
ceny    = {ceny}
cenz    = {cenz}

[run]
z_out     = 0.0
z_max     = 4.6
ievol     = 1
ilpt      = 2
ioutshear = 1
wsmooth   = 1
rmax2rs   = 0.0
NonGauss  = 0
fNL       = 0.0
seed      = 12345
ntile     = {NTILE}
coarse_factor = 4
generate_table = true

[files]
pk         = "validation/websky_6144/data/pk_websky.dat"
# filter bank MUST match the (corrected) cellsize 0.852 Mpc/h: Rf_min=1.65*0.852=1.406 Mpc/h.
# Regenerate via filter_gen for a different cellsize. (filters_websky.dat is the OLD coarse bank.)
filterbank = "validation/websky_6144/data/filters_websky_finecell.dat"
homeltab   = "validation/websky_6144/data/HomelTab_websky.dat"
output     = "catalog_websky_6144_oct{oct_name}.pksc"

[output]
format = "pksc"
path   = "{OUTPUT_PATH}"
"""

    outdir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(outdir, f"config_websky_6144_oct{oct_name}.toml")
    with open(path, "w") as f:
        f.write(content)
    print(f"Wrote {path}  (observer: ({cenx}, {ceny}, {cenz}))")

print(f"\nDone — {len(OCTANT_NAMES)} octant configs")
