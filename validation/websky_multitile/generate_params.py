#!/usr/bin/env python3
"""Generate hpkvd_params.bin for the multi-tile validation run (ntile=2, n1=256)."""
import numpy as np

# Geometry: same tile size as 256^3 single-tile validation
nmesh = 256
nbuff = 22
ntile = 2
nsub = nmesh - 2 * nbuff  # = 212
cellsize = 512.0 / nsub   # = 2.41509... Mpc/h (same as 256^3 run)
dcore_box = nsub * cellsize  # = 512.0 Mpc/h (per-tile core)
dL_box = nmesh * cellsize    # = 618.264... Mpc/h (per-tile total)
# Total grid: N = nsub*ntile + 2*nbuff = 212*2 + 44 = 468
# Total box: dcore_box*(ntile-1) + dL_box = 512 + 618.264 = 1130.264 Mpc/h
next_val = nsub * ntile + 2 * nbuff  # = 468

# Cosmology (Planck 2018, same as single-tile runs)
Omx = 0.2645
OmB = 0.0493
Omvac = 0.6862  # 1 - Omx - OmB (flat)
h = 0.6735

# Run parameters
seed = 13579
ireadfield = 0
ioutshear = 0
global_redshift = 0.0
maximum_redshift = 0.0
num_redshifts = 1
ievol = 0
ilpt = 2
wsmooth = 1
rmax2rs = 0.0
ioutfield = 0
NonGauss = 0
fNL = 0.0
A_nG = 0.0
B_nG = 0.0
R_nG = 0.0
iwant_field_part = 0
largerun = 0

# Collapse parameters
ivir_strat = 2
fcoll_3 = 0.171
fcoll_2 = 0.171
fcoll_1 = 0.01
dcrit = 200.0
iforce_strat = 4

# Table interpolation
TabInterpNx = 50
TabInterpNy = 20
TabInterpNz = 20
TabInterpX1 = 1.5
TabInterpX2 = 8.0
TabInterpY1 = 0.0
TabInterpY2 = 0.5
TabInterpZ1 = -1 + 1e-4
TabInterpZ2 = 1 - 1e-4

# String parameters
fielddir = 'fields/'
densfilein = 'websky_mt2'
densfileout = 'websky_mt2'
pkfile = 'tables/planck18_intermittent.dat'
filterfile = 'filter.dat'
fileout = 'output/catalog_raw.pksc'
TabInterpFile = 'HomelTab.dat'

# Build binary
hpkvd_list = [
    np.int32(ireadfield), np.int32(ioutshear),
    np.float32(global_redshift), np.float32(maximum_redshift),
    np.int32(num_redshifts), np.float32(Omx),
    np.float32(OmB), np.float32(Omvac),
    np.float32(h), np.int32(ntile),  # nlx
    np.int32(ntile), np.int32(ntile),  # nly, nlz
    np.float32(dcore_box), np.float32(dL_box),
    np.float32(0.0), np.float32(0.0),  # cenx, ceny
    np.float32(0.0), np.int32(nbuff),  # cenz, nbuff
    np.int32(next_val), np.int32(ievol),
    np.int32(ivir_strat), np.float32(fcoll_3),
    np.float32(fcoll_2), np.float32(fcoll_1),
    np.float32(dcrit), np.int32(iforce_strat),
    np.int32(TabInterpNx), np.int32(TabInterpNy),
    np.int32(TabInterpNz), np.float32(TabInterpX1),
    np.float32(TabInterpX2), np.float32(TabInterpY1),
    np.float32(TabInterpY2), np.float32(TabInterpZ1),
    np.float32(TabInterpZ2), np.int32(wsmooth),
    np.float32(rmax2rs), np.int32(ioutfield),
    np.int32(NonGauss), np.float32(fNL),
    np.float32(A_nG), np.float32(B_nG),
    np.float32(R_nG), np.int32(ilpt),
    np.int32(iwant_field_part), np.int32(largerun),
]

# String params: each is int32(len) + bytes
for s in [fielddir, densfilein, densfileout, pkfile, filterfile, fileout, TabInterpFile]:
    hpkvd_list.append(np.int32(len(s)))
    hpkvd_list.append(bytes(s, 'ascii'))

with open('fortran_run/hpkvd_params.bin', 'wb') as f:
    for item in hpkvd_list:
        f.write(item)

print(f"Generated hpkvd_params.bin")
print(f"  nmesh={nmesh}, nbuff={nbuff}, ntile={ntile}")
print(f"  nsub={nsub}, cellsize={cellsize:.6f}")
print(f"  dcore_box={dcore_box:.5f}, dL_box={dL_box:.5f}")
print(f"  next={next_val}")
print(f"  Total grid: {next_val}^3")
print(f"  Total box: {dcore_box*(ntile-1) + dL_box:.3f} Mpc/h")
