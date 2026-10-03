#!/usr/bin/env python3
"""hpkvd_params.bin + merge_params.txt for the same-field Fortran run (MATCHED_FORTRAN_2026-10.md).

Geometry must equal configs/matched_fortran.toml AND the compiled hpkvd (arrays.f90 n1 = nmesh = 303):
  nsub = nmesh - 2*nbuff = 251, N = nsub*ntile + 2*nbuff = 1056, cell = 0.8522135 Mpc/h.
All lengths are Mpc/h (the Julia convention; the snapshot run has no chi(z), so Fortran is unit-agnostic
as long as box, P(k) and filters agree). ireadfield = 1 reads fields/Fvec_jmatched (written by
matched_julia.jl field). Parameter order follows hpkvd.f90 read_parameters (46 words + 7 strings).
  usage: python3 gen_fortran_inputs.py <rundir>
"""
import struct, sys, os

run = sys.argv[1]
nmesh, nbuff, ntile, cell, z = 303, 26, 4, 258.2207 / 303, 0.7
Om, OB, OL, h = 0.31, 0.049, 0.69, 0.68
nsub = nmesh - 2 * nbuff
N = nsub * ntile + 2 * nbuff
dcore, dL = nsub * cell, nmesh * cell
assert N == 1056 and abs(dcore / nsub - dL / nmesh) < 1e-9

i = lambda x: struct.pack('<i', int(x))
f = lambda x: struct.pack('<f', float(x))
b = i(1) + i(1) + f(z) + f(z) + i(1)                         # ireadfield ioutshear z zmax nz
b += f(Om - OB) + f(OB) + f(OL) + f(h)
b += i(ntile) * 3 + f(dcore) + f(dL) + f(0) + f(0) + f(0)
b += i(nbuff) + i(N) + i(0)                                  # nbuff next ievol
b += i(2) + f(0.171) + f(0.171) + f(0.01) + f(200.0) + i(4)  # ivir fcoll_3 fcoll_2 fcoll_1 dcrit iforce
b += i(50) + i(20) + i(20) + f(1.5) + f(8.0) + f(0.0) + f(0.5) + f(-1 + 1e-4) + f(1 - 1e-4)  # table (file header wins)
b += i(1) + f(0.0) + i(0) + i(0) + f(0) + f(0) + f(0) + f(0)  # wsmooth rmax2rs ioutfield NonGauss fNL A B R
b += i(2) + i(0) + i(0)                                      # ilpt iwant_field_part largerun
assert len(b) == 4 * 46
for s in ['fields/', 'jmatched', 'jmatched', 'pk_websky_nc.dat', 'filters_websky_finecell.dat',
          'output/fortran_raw.pksc', 'HomelTab_websky.dat']:
    b += i(len(s)) + s.encode('ascii')
open(os.path.join(run, 'hpkvd_params.bin'), 'wb').write(b)

# merge_pkvd: boxsize = core extent (dcore_box = boxsize/ntile); iLexc..iFmrg are read but unused;
# maximum_redshift is implicitly INTEGER in merge_pkvd (list read fails on "0.7"), unused at ievol = 0
mp = [nsub * ntile * cell, ntile, 0, 0.0, 0.0, 0.0, 1, Om - OB, OB, OL, h, 2, 3, 0, 0, 0, 1, 0, 0, 0.0, 0]
fmt = "%f\n%d\n%d\n%f\n%f\n%f\n%d\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%f\n%d\n%s\n%s\n"
open(os.path.join(run, 'merge_params.txt'), 'w').write(fmt % tuple(mp + ['output/fortran_raw.pksc',
                                                                          'output/fortran_merge.pksc']))
print(f"N={N} nsub={nsub} cell={cell:.7f} dcore={dcore:.4f} dL={dL:.4f} core extent={nsub*ntile*cell:.4f}")
