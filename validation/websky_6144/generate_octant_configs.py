#!/usr/bin/env python3
"""
Generate TOML configs for all 8 Websky 6144^3 lightcone octants.

Each octant places the observer at a different corner of the 7700 Mpc/h box.
The naming convention (000..111) is a binary flag: bit set means the observer
sits at the far end of that axis (x=bit0, y=bit1, z=bit2).

IMPORTANT: run_multitile_split uses centered tile coordinates
(tile_center returns (it - (ntile+1)/2) * dcore_box), so the observer
position must be expressed in the same centered coordinate system.
The box extends from -boxsize_full/2 to +boxsize_full/2 in centered coords.

Usage:
    python generate_octant_configs.py
"""
import os

BOX_FULL = 7700.0
HALF_BOX = BOX_FULL / 2.0  # = 3850.0
OUTPUT_PATH = "/home/yguan/projects/aip-aspuru-ab/yguan/websky"

# Octant names: binary flags (bit0=x, bit1=y, bit2=z)
# 000 = origin corner, 111 = far corner, etc.
# In centered coords: bit=0 -> -3850, bit=1 -> +3850
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
#   ntile=16, nmesh=399, nbuff=8
#   nsub = 399 - 16 = 383
#   N = 383*16 + 16 = 6144
#   cellsize = 7700/6144 = 1.25326 Mpc/h
#   per-tile boxsize = 399 * 1.25326 = 500.05 Mpc/h
#   total box = 6144 * 1.25326 = 7700 Mpc/h

[cosmology]
Om   = 0.31
OB   = 0.049
OL   = 0.69
h    = 0.68

[grid]
n       = 399
boxsize = 500.05
nbuff   = 8
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
rmax2rs   = 1.0
NonGauss  = 0
fNL       = 0.0
seed      = 12345
ntile     = 16
coarse_factor = 4
generate_table = true

[files]
pk         = "validation/websky_6144/data/pk_websky.dat"
filterbank = "validation/websky_6144/data/filters_websky.dat"
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
