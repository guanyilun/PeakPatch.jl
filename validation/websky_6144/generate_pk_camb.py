#!/usr/bin/env python3
"""
Generate CAMB matter power spectrum for Websky cosmology.

Websky parameters from Stein+ 2020 (arXiv:2001.08787):
  Om=0.31, OB=0.049, OL=0.69, h=0.68, ns=0.965, sigma8=0.81

Requires: pip install camb

Usage:
    python generate_pk_camb.py [output_path]
"""
import sys
import os
import camb
import numpy as np

# --- Websky cosmology (Planck 2018) ---
Om_total = 0.31
OB = 0.049
OL = 0.69
h = 0.68
ns = 0.965
sigma8 = 0.81

Omh2 = Om_total * h**2
Obh2 = OB * h**2

# --- Output path ---
default_out = os.path.join(os.path.dirname(__file__), "data", "pk_websky.dat")
pk_path = sys.argv[1] if len(sys.argv) > 1 else default_out
os.makedirs(os.path.dirname(pk_path), exist_ok=True)

# --- Run CAMB ---
pars = camb.CAMBparams()
pars.set_cosmology(H0=100 * h, ombh2=Obh2, omch2=Omh2 - Obh2,
                   omk=1.0 - Om_total - OL)
pars.InitPower.set_params(ns=ns, As=2e-9)
pars.set_matter_power(redshifts=[0.0], kmax=100.0)
pars.NonLinear = camb.model.NonLinear_none
pars.WantLensing = False

# Iteratively scale As to match target sigma8
results = camb.get_results(pars)
s8 = results.get_sigma8_0()
As_new = 2e-9 * (sigma8 / s8) ** 2

pars.InitPower.set_params(ns=ns, As=As_new)
results = camb.get_results(pars)
s8_final = results.get_sigma8_0()

# --- Get P(k) ---
kh, z, pk = results.get_matter_power_spectrum(minkh=1e-5, maxkh=100.0, npoints=500)

# --- Write ---
with open(pk_path, "w") as f:
    f.write(f"# CAMB power spectrum for Websky cosmology\n")
    f.write(f"# Om={Om_total}, OB={OB}, OL={OL}, h={h}, ns={ns}, sigma8={s8_final:.4f}\n")
    f.write(f"# k [h/Mpc]    P(k) [(Mpc/h)^3]\n")
    for i in range(len(kh)):
        f.write(f"{kh[i]:.10e}  {pk[0][i]:.10e}\n")

print(f"Wrote {len(kh)} k-points to {pk_path}")
print(f"sigma8 = {s8_final:.4f}")
print(f"k range: {kh[0]:.2e} to {kh[-1]:.2e} h/Mpc")
print(f"P(k) at k=0.1 h/Mpc: {pk[0][np.argmin(np.abs(kh - 0.1))]:.2e}")
