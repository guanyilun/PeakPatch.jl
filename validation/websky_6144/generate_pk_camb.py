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

# --- Peak Patch normalization convention ---
# The field-generation convolution uses amp = sqrt(P * dk^3 * n^3) = sqrt(P/dx^3)*(2pi)^1.5,
# so the P(k) FILE must be PRE-DIVIDED by (2*pi)^3 for the generated field to have the
# correct power spectrum. (The official peakpatch/tools/powerspectrum_create.py divides
# by (2*pi*h)^3 because its CAMB output is in physical Mpc; CAMB's get_matter_power_spectrum
# here returns (Mpc/h)^3 with k in h/Mpc, so the factor is (2*pi)^3, NO h. Verified
# empirically 2026-06-14 via probe_pk_norm.jl: ÷(2pi)^3 gives field sigma(R)/theory ~1.00.)
pp_norm = (2.0 * np.pi) ** 3

# --- Write ---
with open(pk_path, "w") as f:
    f.write(f"# CAMB power spectrum for Websky cosmology\n")
    f.write(f"# Om={Om_total}, OB={OB}, OL={OL}, h={h}, ns={ns}, sigma8={s8_final:.4f}\n")
    f.write(f"# P(k) PRE-DIVIDED by (2*pi)^3 for Peak Patch convolution convention\n")
    f.write(f"# k [h/Mpc]    P_pp(k) = P(k)/(2pi)^3 [(Mpc/h)^3]\n")
    for i in range(len(kh)):
        f.write(f"{kh[i]:.10e}  {pk[0][i] / pp_norm:.10e}\n")

print(f"Wrote {len(kh)} k-points to {pk_path}")
print(f"sigma8 = {s8_final:.4f}")
print(f"k range: {kh[0]:.2e} to {kh[-1]:.2e} h/Mpc")
print(f"P(k) at k=0.1 h/Mpc: {pk[0][np.argmin(np.abs(kh - 0.1))]:.2e}")
