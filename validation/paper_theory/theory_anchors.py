#!/usr/bin/env python3
"""Theory anchor curves for the methods paper (docs/paper_comparison_plan_2026-09.md, A4).

Outputs (validation/paper_theory/results/):
  kappa_limber.txt  C_ell^kk (Limber, Websky eq 3.21) for linear / Halofit(Takahashi) /
                    HMcode2020 P(k,z); ranges 0->z*, 0->4.5, 4.5->z*
  ksz_doppler.txt   linear Doppler kSZ: Websky eq 4.3 boundary term (z_max 4.5, 4.6) and the
                    exact linear line-of-sight integral (sharp cutoff at z_max)
  ksz_ov.txt        linear Ostriker-Vishniac (Websky eq 4.4), Limber, z < 4.5
  pk_check.txt      CAMB linear P(k,0) vs the pk_websky.dat used for the ICs

Conventions match src/FieldMap.jl so the curves overlay our field maps directly:
  lengths Mpc/h, flat LCDM (Om=0.31, OL=0.69, NO radiation) for chi(z), D(z), f(z);
  P0(k) = pk_websky.dat * (2pi)^3 (the exact spectrum the ICs were drawn from);
  dtau/dchi = sigT_ne0 * x_e(z) * (1+z)^2, sigT_ne0 = sigma_T * (rho_c/m_p) * f_e * Ob h^2 [per Mpc/h],
  f_e = 0.9, Y_He = 0.245, He doubly ionized at z<3 and singly at z>=3;
  v = a H f D delta0 k_hat/k [km/s], Delta T/T = -int dchi (dtau/dchi) v_r / c.
Nonlinear P(k,z) comes from CAMB with the generate_pk_camb.py setup (As scaled to sigma8=0.81).

Run: module load python/3.11 scipy-stack; source /home/yguan/scratch/theory_env/bin/activate
     python theory_anchors.py            (~10 min, <2 GB)
"""
import os, sys, time
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import interp1d
from scipy.special import spherical_jn, jv
import camb
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
PKFILE = os.path.join(HERE, "..", "websky_6144", "data", "pk_websky.dat")

Om, OB, OL, h, ns, sigma8 = 0.31, 0.049, 0.69, 0.68, 0.965, 0.81
C_KMS = 299792.458
TCMB_UK = 2.7255e6
F_E, Y_HE = 0.9, 0.245
SIGT_NE0 = 6.65246e-29 * 11.2299 * F_E * OB * h**2 * 3.0857e22 / h   # per Mpc/h (FieldMap.jl)
H0_C = 1.0 / 2997.92458                                              # (Mpc/h)^-1

def x_e(z):
    return np.where(z < 3, 1 - Y_HE / 2, 1 - 3 * Y_HE / 4)

# ---------------- background: flat LCDM, no radiation (= PeakPatch Cosmology) -------------
E = lambda z: np.sqrt(Om * (1 + z) ** 3 + OL)
zg = np.concatenate([np.linspace(0, 10, 20001), np.linspace(10.01, 1200, 20000)])
chig = np.concatenate([[0.0], np.cumsum(0.5 * (1 / E(zg[1:]) + 1 / E(zg[:-1])) * np.diff(zg))]) * 2997.92458
chi_of_z = interp1d(zg, chig)
z_of_chi = interp1d(chig, zg)

def growth():
    # linear growth D(a), normalized D(1)=1, and f = dlnD/dlna, from the standard ODE
    def rhs(lna, y):
        a = np.exp(lna); Oma = Om / (Om + OL * a**3)
        D, dD = y
        return [dD, -(2 - 1.5 * Oma) * dD + 1.5 * Oma * D]
    lna = np.linspace(np.log(1e-3), 0, 4000)
    s = solve_ivp(rhs, (lna[0], 0), [1e-3, 1e-3], t_eval=lna, rtol=1e-10, atol=1e-13)
    D = s.y[0] / s.y[0][-1]; f = s.y[1] / s.y[0]
    z = np.exp(-lna) - 1
    return interp1d(z[::-1], D[::-1]), interp1d(z[::-1], f[::-1])
D_of_z, f_of_z = growth()

# ---------------- P0(k): the IC spectrum ----------------
pkdat = np.loadtxt(PKFILE)
k_file, P_file = pkdat[:, 0], pkdat[:, 1] * (2 * np.pi) ** 3
_lp = interp1d(np.log(k_file), np.log(P_file), bounds_error=False, fill_value="extrapolate")
def P0(k):
    return np.exp(_lp(np.log(k)))

# ---------------- CAMB (generate_pk_camb.py setup) ----------------
def camb_pars(nonlin):
    Omh2 = Om * h**2; Obh2 = OB * h**2
    p = camb.CAMBparams()
    p.set_cosmology(H0=100 * h, ombh2=Obh2, omch2=Omh2 - Obh2, omk=1.0 - Om - OL)
    p.InitPower.set_params(ns=ns, As=2e-9)
    p.set_matter_power(redshifts=[0.0], kmax=100.0)
    p.NonLinear = camb.model.NonLinear_none
    r = camb.get_results(p)
    As = 2e-9 * (sigma8 / r.get_sigma8_0()) ** 2
    p.InitPower.set_params(ns=ns, As=As)
    return p, As

# ---------------- kSZ Doppler helpers (module level for multiprocessing) ----------------
def G(chi):   # dtau/dchi * a H f D / c    [(Mpc/h)^-1]
    z = z_of_chi(chi); a = 1 / (1 + z)
    return SIGT_NE0 * x_e(z) * (1 + z) ** 2 * a * 100 * E(z) * f_of_z(z) * D_of_z(z) / C_KMS

def sph_j(l, x):
    """spherical j_l via jv (5x faster than spherical_jn at high l); j_l(x) set to 0 for x < l/2."""
    out = np.zeros_like(x)
    m = x > 0.5 * l
    out[m] = jv(l + 0.5, x[m]) * np.sqrt(np.pi / (2 * x[m]))
    return out

def doppler_boundary(l, zm):
    # Websky eq 4.3: C_l = (2/pi) G(chi_m)^2 int dk P0(k)/k^2 j_l(k chi_m)^2
    chim = float(chi_of_z(zm))
    kmax = max(0.05, (3.0 * l + 300.0) / chim)
    k = np.arange(1e-5, kmax, min(5e-5, 0.05 / chim))
    return 2 / np.pi * G(chim) ** 2 * np.trapz(P0(k) / k**2 * sph_j(l, k * chim) ** 2, k)

CHIM45 = float(chi_of_z(4.5))
_CHI = np.arange(1.0, CHIM45, 1.0)
_GC = G(_CHI)
_GP = np.gradient(_GC, _CHI)           # includes the He-reionization step at z=3 as a 1-cell spike

def doppler_full(l):
    # exact linear LOS: Delta_l(k) = int G j_l'(k chi)/k dchi, done by parts:
    #   = G(chi_m) j_l(k chi_m)/k^2 - (1/k^2) int G'(chi) j_l(k chi) dchi
    kmax = max(0.05, (3.0 * l + 150.0) / CHIM45)
    kk = np.arange(2e-5, kmax, 5e-5)
    acc = np.empty(len(kk))
    for i0 in range(0, len(kk), 100):
        ks = kk[i0:i0 + 100]
        J = sph_j(l, np.outer(ks, _CHI))
        acc[i0:i0 + 100] = (_GC[-1] * sph_j(l, ks * CHIM45) - np.trapz(J * _GP[None, :], _CHI, axis=1)) / ks**2
    return 2 / np.pi * np.trapz(kk**2 * P0(kk) * acc**2, kk)

def doppler_bulk_limber(l):
    nu = l + 0.5
    return np.trapz(_CHI**2 * _GP**2 * P0(nu / _CHI) / nu**4, _CHI)

def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)

def main():
    t0 = time.time()
    p, As = camb_pars(False)
    res = camb.get_results(p)
    log(f"CAMB {camb.__version__}  As={As:.4e}  sigma8={res.get_sigma8_0():.4f}")
    kh, _, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=10.0, npoints=300)
    rat = pk[0] / P0(kh)
    np.savetxt(os.path.join(OUT, "pk_check.txt"), np.c_[kh, pk[0], P0(kh), rat],
               header="k[h/Mpc]  P_camb_lin(z=0)  P0_icfile*(2pi)^3  ratio   [CAMB %s vs 1.6.6 file]" % camb.__version__)
    log(f"P_camb/P_file: median {np.median(rat):.4f}, range {rat.min():.4f}-{rat.max():.4f}")

    if not os.environ.get("SKIP_KAPPA"):
        # ---------------- kappa Limber ----------------
        chistar_nr = float(chi_of_z(1089.0))                       # FieldMap default (no radiation)
        chistar_camb = res.comoving_radial_distance(1089.0) * h    # Mpc -> Mpc/h, with radiation
        log(f"chi*(1089): no-radiation {chistar_nr:.1f}, CAMB {chistar_camb:.1f} Mpc/h")
        zs = np.concatenate([np.linspace(0, 5, 101)[:-1], np.geomspace(5, 1089, 60)])
        interps = {}
        for name, nl in (("lin", camb.model.NonLinear_none), ("halofit", "takahashi"), ("hmcode", "mead2020")):
            q = p.copy()
            if name == "lin":
                q.NonLinear = camb.model.NonLinear_none
            else:
                q.NonLinear = camb.model.NonLinear_pk
                q.NonLinearModel.set_params(halofit_version=nl)
            interps[name] = camb.get_matter_power_interpolator(q, zs=zs, kmax=200.0, nonlinear=(name != "lin"),
                                                               hubble_units=True, k_hunit=True)
            log(f"interpolator {name} ready")
        ells = np.unique(np.round(np.geomspace(2, 10000, 400)).astype(int))

        def limber(PK, zlo, zhi, chistar, lin_ours=False):
            chi = np.linspace(float(chi_of_z(zlo)) + 1e-3, float(chi_of_z(zhi)), 6000)
            z = z_of_chi(chi)
            W = 1.5 * Om * H0_C**2 * (1 + z) * chi * (1 - chi / chistar)
            cl = np.empty(len(ells))
            for i, l in enumerate(ells):
                k = (l + 0.5) / chi
                if lin_ours:
                    pkv = P0(k) * D_of_z(np.minimum(z, 999.0)) ** 2
                else:
                    pkv = PK.P(z, k, grid=False)
                pkv = np.where(k < 200.0, pkv, 0.0)
                cl[i] = np.trapz(W**2 / chi**2 * pkv, chi)
            return cl

        cols = [ells]; hdr = ["ell"]
        zstar = 1089.0
        for name in ("lin", "halofit", "hmcode"):
            for (zlo, zhi, tag) in ((0.0, zstar, "0-zs"), (0.0, 4.5, "0-4.5"), (4.5, zstar, "4.5-zs")):
                cols.append(limber(interps[name], zlo, zhi, chistar_camb)); hdr.append(f"{name}_{tag}")
        # the linear z<4.5 curve with OUR P0, growth and chi* (field-map-matched)
        cols.append(limber(None, 0.0, 4.5, chistar_nr, lin_ours=True)); hdr.append("linOURS_0-4.5")
        cols.append(limber(interps["halofit"], 0.0, 4.5, chistar_nr)); hdr.append("halofit_0-4.5_chistarNR")
        np.savetxt(os.path.join(OUT, "kappa_limber.txt"), np.column_stack(cols),
                   header="C_ell^kk, Limber k=(l+1/2)/chi, Om=0.31 kernel. chi* = CAMB chi(1089)=%.1f Mpc/h "
                          "unless noted; NR = no-radiation chi*=%.1f (FieldMap default)\n" % (chistar_camb, chistar_nr)
                          + "  ".join(hdr))
        log("kappa done")

    # ---------------- kSZ Doppler ----------------
    lD = np.unique(np.round(np.geomspace(2, 3000, 60)).astype(int))
    cols = [lD]; hdr = ["ell"]
    for zm in (4.5, 4.6):
        clb = np.array([doppler_boundary(l, zm) for l in lD])
        cols.append(clb * lD * (lD + 1) / (2 * np.pi) * TCMB_UK**2); hdr.append(f"Dl_eq4.3_zmax{zm}")
    log("Doppler boundary terms done")
    # exact linear LOS integral (sharp cutoff at z_max=4.5 = the field map geometry), l<=1000
    lfull = lD[lD <= 1000]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(int(os.environ.get("NPROC", "4"))) as pool:
        clf = pool.map(doppler_full, [int(l) for l in lfull])
    clf = np.array(list(clf) + [np.nan] * (len(lD) - len(lfull)))
    cols.append(clf * lD * (lD + 1) / (2 * np.pi) * TCMB_UK**2); hdr.append("Dl_fullLOS_zmax4.5")
    # bulk term in Limber (integration by parts): boundary + int dchi chi^2 G'^2 P0(nu/chi)/nu^4
    cll = np.array([doppler_boundary(l, 4.5) + doppler_bulk_limber(l) for l in lD])
    cols.append(cll * lD * (lD + 1) / (2 * np.pi) * TCMB_UK**2); hdr.append("Dl_boundary+bulkLimber_zmax4.5")
    np.savetxt(os.path.join(OUT, "ksz_doppler.txt"), np.column_stack(cols),
               header="linear Doppler kSZ D_ell [uK^2]; f_e=0.9 inside n_e0 (enters squared); "
                      "eq4.3 = Websky boundary term; fullLOS = exact linear integral 0<chi<chi(4.5) (l<=1000); "
                      "boundary+bulkLimber = IBP approximation (no cross term)\n"
                      + "  ".join(hdr))
    for l, a, b in zip(lD, clf, cll):
        if l in (2, 10, 30, 100, 300) or (np.isfinite(a) and l == lfull[-1]):
            log(f"  l={l}: full/approx = {a/b:.3f}")
    log("Doppler done")

    # ---------------- Ostriker-Vishniac (linear, Limber) ----------------
    kq = np.geomspace(1e-3, 30.0, 160)
    kp = np.geomspace(1e-5, 100.0, 3000)
    mu, wmu = np.polynomial.legendre.leggauss(256)
    I = np.empty(len(kq))
    for i, k in enumerate(kq):
        K, M = np.meshgrid(kp, mu, indexing="ij")
        q2 = k**2 + K**2 - 2 * k * K * M
        q2 = np.maximum(q2, 1e-30)
        integ = P0(np.sqrt(q2)) * P0(K) * (1 - M**2) * k * (k - 2 * K * M) / (K**2 * q2)
        inner = integ @ wmu
        I[i] = np.trapz(inner * kp**2, kp) / (4 * np.pi**2)
    Iint = interp1d(np.log(kq), np.log(np.maximum(I, 1e-300)), bounds_error=False, fill_value=-700)
    lO = np.unique(np.round(np.geomspace(10, 10000, 80)).astype(int))
    chi = np.linspace(5.0, float(chi_of_z(4.5)), 5000); z = z_of_chi(chi); a = 1 / (1 + z)
    g = SIGT_NE0 * x_e(z) * (1 + z) ** 2
    amp = (a * 100 * E(z) * f_of_z(z) / C_KMS) ** 2 * D_of_z(z) ** 4
    clo = np.array([0.5 * np.trapz(g**2 / chi**2 * amp * np.exp(Iint(np.log((l + 0.5) / chi))), chi) for l in lO])
    np.savetxt(os.path.join(OUT, "ksz_ov.txt"), np.c_[lO, clo * lO * (lO + 1) / (2 * np.pi) * TCMB_UK**2],
               header="linear Ostriker-Vishniac D_ell [uK^2], z<4.5, Limber, P_qperp from linear P0 x D^4\nell  Dl_OV_lin")
    log(f"OV done; D_3000={np.interp(3000, lO, clo*lO*(lO+1)/2/np.pi*TCMB_UK**2):.3f} uK^2")
    log(f"total {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
