#!/usr/bin/env python3
"""Generate a CONSISTENT hpkvd_params.bin for cross-code (Julia-vs-Fortran) tests.

WHY THIS EXISTS — the cellsize trap (2026-06-14):
  hpkvd's grid dimension is a COMPILE-TIME constant: `src/hpkvd/arrays.f90` has
  `integer, parameter :: n1 = N_REPLACE`, substituted at build time (this binary
  was built with n1 = nmesh = 256). So the binary ALWAYS uses
      nsub      = nmesh - 2*nbuff           (= 256 - 44 = 212 for this build)
      n_global  = nsub*ntile + 2*nbuff      (= 468)
  REGARDLESS of the `next` field in the param. The cellsize the run actually uses is
      cellsize  = dcore_box / nsub          (== dL_box / nmesh, must be consistent)
  A param sized for a DIFFERENT nmesh (e.g. 150) silently runs at the wrong cellsize
  (dcore_box/212 instead of dcore_box/106 -> half), which doubled the grid-floor peak
  count and produced a bogus "Fortran finds 2x more peaks". See
  FORTRAN_COMPARISON_2026-06-14.md.

USAGE:
  python3 gen_fortran_params.py --cellsize 1.2533 --nmesh 256 --nbuff 22 --ntile 2 \
      --pk tables/pk_websky_nc.dat --filters filters_websky.dat \
      --out output/catalog.pksc [--ireadfield 0 --densfilein name --ievol 0]

ALWAYS cross-check after launch:
  * the logged "Slab decomposition: n = <N>" equals nsub*ntile + 2*nbuff
  * sigma(R) at a couple of filter scales agrees with Julia to ~1% (physical invariant)
"""
import struct, argparse, sys

def build(cellsize, nmesh, nbuff, ntile, pk, filters, out, ireadfield, densfilein,
          ievol, Om, OB, OL, h, homeltab, fielddir):
    nsub = nmesh - 2*nbuff
    n_global = nsub*ntile + 2*nbuff
    dcore_box = nsub * cellsize
    dL_box    = nmesh * cellsize
    # consistency guard — the exact assertion the bad 2026-06-14 param violated
    c1, c2 = dcore_box/nsub, dL_box/nmesh
    assert abs(c1 - c2) < 1e-6, f"INCONSISTENT cellsize: dcore/nsub={c1} != dL/nmesh={c2}"
    Omx = Om - OB
    i = lambda x: struct.pack('<i', int(x))
    f = lambda x: struct.pack('<f', float(x))
    buf  = i(ireadfield)+i(0)+f(0.0)+f(0.0)+i(1)+f(Omx)+f(OB)+f(OL)+f(h)
    buf += i(ntile)+i(ntile)+i(ntile)+f(dcore_box)+f(dL_box)+f(0.0)+f(0.0)+f(0.0)
    buf += i(nbuff)+i(n_global)+i(ievol)
    buf += i(2)+f(0.171)+f(0.171)+f(0.01)+f(200.0)+i(4)          # collapse params
    buf += i(50)+i(20)+i(20)+f(1.5)+f(8.0)+f(0.0)+f(0.5)+f(-1+1e-4)+f(1-1e-4)
    buf += i(1)+f(0.0)+i(0)+i(0)+f(0.0)+f(0.0)+f(0.0)+f(0.0)+i(2)+i(0)+i(0)  # wsmooth=1,ilpt=2
    for s in [fielddir, densfilein, densfilein, pk, filters, out, homeltab]:
        buf += i(len(s)) + s.encode('ascii')
    open('hpkvd_params.bin','wb').write(buf)
    print(f"hpkvd_params.bin written ({len(buf)} bytes)")
    print(f"  REQUIRES binary compiled with nmesh={nmesh} (check arrays_gen.f90 / log n=={n_global})")
    print(f"  nsub={nsub}  n_global={n_global}  cellsize={c1:.5f} Mpc/h  box={n_global*cellsize:.2f} Mpc/h")
    print(f"  dcore_box={dcore_box:.4f}  dL_box={dL_box:.4f}  ievol={ievol} ireadfield={ireadfield}")

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--cellsize', type=float, required=True)
    p.add_argument('--nmesh', type=int, default=256, help="MUST match the compiled binary (arrays_gen.f90)")
    p.add_argument('--nbuff', type=int, default=22)
    p.add_argument('--ntile', type=int, default=2)
    p.add_argument('--pk', required=True); p.add_argument('--filters', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ireadfield', type=int, default=0); p.add_argument('--densfilein', default='wk')
    p.add_argument('--ievol', type=int, default=0)
    p.add_argument('--Om', type=float, default=0.31); p.add_argument('--OB', type=float, default=0.049)
    p.add_argument('--OL', type=float, default=0.69); p.add_argument('--h', type=float, default=0.68)
    p.add_argument('--homeltab', default='HomelTab.dat'); p.add_argument('--fielddir', default='fields/')
    a = p.parse_args()
    build(a.cellsize, a.nmesh, a.nbuff, a.ntile, a.pk, a.filters, a.out, a.ireadfield,
          a.densfilein, a.ievol, a.Om, a.OB, a.OL, a.h, a.homeltab, a.fielddir)
