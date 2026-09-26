#!/usr/bin/env python3
"""Frozen-campaign FIELD maps (all measured octants) vs the theory anchors.

kappa: field map (full matter, 2LPT lattice, z<4.5, chi*=9656 Mpc/h as painted) vs Limber
       linear (our P0, growth) and Halofit/HMcode, all z<4.5 with the SAME chi*.
kSZ:   field map vs linear Doppler (exact LOS, sharp z=4.5 cut) + linear OV.
Inputs: results/fieldmap_cl_prod_oct???.txt (measure_fieldmap_cl.jl), results/ksz_*.txt.
Writes results/fieldmaps_vs_theory_prod.txt and prints the table.
"""
import os
import numpy as np
import camb
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid
import theory_anchors as T

R = T.OUT
CHISTAR = 14200.0 * T.h       # run_fieldmap_octant.jl: Websky/pks2map hardwired source plane

import glob
files = sorted(glob.glob(os.path.join(R, "fieldmap_cl_prod_oct???.txt")))
octs = [os.path.basename(f)[-7:-4] for f in files]
FM = np.array([np.loadtxt(f) for f in files])          # (noct, nband, 5)
l0, l1, leff = FM[0, :, 0], FM[0, :, 1], FM[0, :, 2]
print("octants:", " ".join(octs))

p, _ = T.camb_pars(False)
zs = np.linspace(0, 5, 101)
curves = {}
for name, nl in (("lin", None), ("halofit", "takahashi"), ("hmcode", "mead2020")):
    q = p.copy()
    if nl is None:
        q.NonLinear = camb.model.NonLinear_none
    else:
        q.NonLinear = camb.model.NonLinear_pk
        q.NonLinearModel.set_params(halofit_version=nl)
    curves[name] = camb.get_matter_power_interpolator(q, zs=zs, kmax=200.0, nonlinear=nl is not None,
                                                      hubble_units=True, k_hunit=True)

chi = np.linspace(1.7, float(T.chi_of_z(4.5)), 6000)       # rmin = 2 cells, as painted
z = T.z_of_chi(chi)
W = 1.5 * T.Om * T.H0_C**2 * (1 + z) * chi * (1 - chi / CHISTAR)

def band(fun, a, b):
    ls = np.arange(int(a), int(b) + 1)
    return np.mean([fun(l) for l in ls[:: max(1, len(ls) // 20)]])

def limber(l, name):
    k = (l + 0.5) / chi
    if name == "linOURS":
        pk = T.P0(k) * T.D_of_z(z) ** 2
    else:
        pk = curves[name].P(z, k, grid=False)
    return np.trapz(W**2 / chi**2 * np.where(k < 200, pk, 0.0), chi)

dop = np.loadtxt(os.path.join(R, "ksz_doppler.txt"))
hdr = open(os.path.join(R, "ksz_doppler.txt")).readlines()[1].split()[1:]
jf = hdr.index("Dl_fullLOS_zmax4.5"); ja = hdr.index("Dl_boundary+bulkLimber_zmax4.5")
ldop = dop[:, 0]
dfull = np.where(np.isfinite(dop[:, jf]), dop[:, jf], dop[:, ja])      # exact <=1000, IBP above
ov = np.loadtxt(os.path.join(R, "ksz_ov.txt"))
def interp_log(x, xs, ys):
    return np.exp(np.interp(np.log(x), np.log(xs), np.log(ys)))

rows = []
print(f"{'l_eff':>7} | {'k/lin':>6} {'sd':>5} {'k/hfit':>6} {'k/hmc':>6} | {'Dksz':>6} {'Dopp':>6} {'OV':>6} {'ksz/th':>6} {'sd':>5} {'min':>5} {'max':>5}")
for ib, (a, b, le) in enumerate(zip(l0, l1, leff)):
    tl = band(lambda l: limber(l, "linOURS"), a, b)
    th = band(lambda l: limber(l, "halofit"), a, b)
    tm = band(lambda l: limber(l, "hmcode"), a, b)
    dd = interp_log(le, ldop, dfull); do = interp_log(le, ov[:, 0], ov[:, 1]) if le >= ov[0, 0] else 0.0
    tot = dd + do
    rk = FM[:, ib, 3] / tl; rz = FM[:, ib, 4] / tot
    ck = FM[:, ib, 3].mean(); dz = FM[:, ib, 4].mean()
    rows.append([le, ck, tl, th, tm, ck / tl, rk.std(ddof=1) if len(rk) > 1 else np.nan, ck / th, ck / tm,
                 dz, dd, do, dz / tot, rz.std(ddof=1) if len(rz) > 1 else np.nan, rz.min(), rz.max()])
    r = rows[-1]
    print(f"{le:7.0f} | {r[5]:6.3f} {r[6]:5.3f} {r[7]:6.3f} {r[8]:6.3f} | {dz:6.3f} {dd:6.3f} {do:6.3f} {r[12]:6.3f} {r[13]:5.3f} {r[14]:5.2f} {r[15]:5.2f}")
np.savetxt(os.path.join(R, "fieldmaps_vs_theory_prod.txt"), np.array(rows),
           header="prod field maps (octants %s; mean over octants, sd = octant-to-octant) vs theory (z<4.5, chi*=%.1f)\n"
                  "l_eff Cl_kap_map Cl_linOURS Cl_halofit Cl_hmcode map/lin sd map/halofit map/hmcode "
                  "Dl_ksz_map Dl_doppler_lin Dl_OV_lin map/(dopp+OV) sd min max" % (",".join(octs), CHISTAR))
