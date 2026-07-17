#!/usr/bin/env python3
"""
15_neff_vs_z.py -- Galaxy surface density vs redshift from the UV luminosity function.

Computes dN/dz per arcmin^2 for NIR (Y/J) limiting magnitudes, by integrating a
Schechter UV luminosity function with observed redshift evolution.  This is the
Malmquist-free way to get the number density a survey of a given depth WOULD see
-- it cannot be obtained by cutting a real catalog fainter than its own detection
limit (e.g. COSMOS-Web saturates at ~27.5).

Method
------
  M_lim(z) = m_lim - DM(z) + 2.5*log10(1+z)        # flat-f_nu k-correction
  n(<M_lim, z) = \int_{-inf}^{M_lim} Phi(M) dM     # Schechter, integrated numerically
  dN/dz/arcmin^2 = n(z) * (dV/dz/dOmega)           # comoving volume element

Schechter form (Harikane+2022 Eq. 35; standard):
  Phi(M) = (ln10/2.5) phi* 10^{-0.4(M-M*)(alpha+1)} exp(-10^{-0.4(M-M*)})

Luminosity function references
-----------------------------
  z = 2.1 - 8.9 : Bouwens et al. 2021, AJ 162, 47 (arXiv:2102.07775), Table 5
                  rest-frame UV LF Schechter parameters.
  z = 9, 12, 16 : Harikane et al. 2023, ApJS 265, 5 (arXiv:2208.01612), Table 8
                  Schechter fits (alpha fixed to -2.35).
  Both papers define M_UV assuming a flat rest-frame UV continuum, which is the
  same k-correction convention adopted here.

CAVEATS (read before using in a proposal)
-----------------------------------------
  * This is a rest-UV-selected number density throughout.  At z < ~4 an observed
    Y/J band does NOT probe rest-frame 1500-1600A, so the flat-f_nu k-correction
    is what makes the low-z end self-consistent -- it is exact only for a
    flat-f_nu SED.  A J-selected sample at z~1 (rest-frame optical) would differ.
  * All z<2 nodes are deep HST (UVCANDELS, UVUDF) rather than shallow
    ground-based data.  Mehta z~2.2 (M*=-19.92, phi*=5.0e-3, a=-1.32) is
    consistent with Bouwens z~2.1 (-20.28, 4.0e-3, -1.52) at the join.
  * Small cosmology differences between sources (Mehta uses Planck15; others
    Om=0.3,h=0.7) are not renormalized -- a few-percent level effect.
  * No LF node below z~0.7, so the curve starts there.
  * The result at m_lim = 29.3 is dominated by the faint-end slope alpha, which
    is an extrapolation well below the magnitudes where the LF is constrained.

Usage:
    python3 15_neff_vs_z.py        # writes outputs/neff_vs_z.{pdf,png}
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

# ---------------------------------------------------------------- Times fonts
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"

COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)
SR_PER_ARCMIN2 = (180.0 / np.pi * 60.0) ** 2          # arcmin^2 per steradian

# ---------------------------------------------------------------- LF table
# z, M*_UV, phi* [Mpc^-3], alpha
LF_NODES = np.array([
    # --- Bouwens et al. 2021 (AJ 162, 47), Table 5 ---
    #     rest-frame UV LF Schechter parameters (alpha is already monotonic
    #     across these rows, so plain interpolation is faithful here).
    [2.1, -20.28, 4.0e-3,   -1.52],
    [2.9, -20.87, 2.1e-3,   -1.61],
    [3.8, -20.93, 1.69e-3,  -1.69],
    [4.9, -21.10, 0.79e-3,  -1.74],
    [5.9, -20.93, 0.51e-3,  -1.93],
    [6.8, -21.15, 0.19e-3,  -2.06],
    [7.9, -20.93, 0.09e-3,  -2.23],
    [8.9, -21.15, 0.021e-3, -2.33],
    # --- Harikane et al. 2023 (ApJS 265, 5), Table 8 (Schechter, alpha fixed) ---
    [9.0,  -21.24, 10**-4.83, -2.35],
    [12.0, -20.47, 10**-5.06, -2.35],
    [16.0, -20.80, 10**-5.84, -2.35],
])
Z_MIN, Z_MAX = LF_NODES[0, 0], LF_NODES[-1, 0]


def lf_params(z):
    """Interpolate (M*, phi*, alpha) at redshift z; log-interpolate phi*."""
    zn = LF_NODES[:, 0]
    mstar = np.interp(z, zn, LF_NODES[:, 1])
    logphi = np.interp(z, zn, np.log10(LF_NODES[:, 2]))
    alpha = np.interp(z, zn, LF_NODES[:, 3])
    return mstar, 10.0**logphi, alpha


def schechter_ndens(m_lim_abs, mstar, phistar, alpha, m_bright=-30.0, n=3000):
    """Number density [Mpc^-3] of galaxies brighter than absolute mag m_lim_abs."""
    if m_lim_abs <= m_bright:
        return 0.0
    M = np.linspace(m_bright, m_lim_abs, n)
    x = 10.0 ** (-0.4 * (M - mstar))
    phi = (np.log(10) / 2.5) * phistar * x ** (alpha + 1.0) * np.exp(-x)
    return float(np.trapezoid(phi, M))


def dNdz_per_arcmin2(z, m_lim):
    """dN/dz per arcmin^2 for an observed-frame (Y/J) limiting magnitude."""
    dm = COSMO.distmod(z).value
    M_lim = m_lim - dm + 2.5 * np.log10(1.0 + z)      # flat-f_nu k-correction
    mstar, phistar, alpha = lf_params(z)
    n = schechter_ndens(M_lim, mstar, phistar, alpha)
    dV = COSMO.differential_comoving_volume(z).value  # Mpc^3 / sr / dz
    return n * dV / SR_PER_ARCMIN2


# ---------------------------------------------------------------- compute
DEPTHS = [(27.5, "COSMOS-Web  (27.5)",        "#1f6fd6"),
          (28.2, "HLWAS UD  (28.2)",          "#e08a1e"),
          (29.3, "HLTDS Deep  (29.3)",        "#e12b25")]

zgrid = np.linspace(Z_MIN, 10.0, 300)

fig, ax = plt.subplots(figsize=(8.6, 5.8))
for m_lim, label, col in DEPTHS:
    y = np.array([dNdz_per_arcmin2(z, m_lim) for z in zgrid])
    pk = zgrid[np.argmax(y)]
    tot = np.trapezoid(y, zgrid)
    ax.plot(zgrid, y, lw=2.4, color=col,
            label=f"{label}   peak $z$ = {pk:.2f}")
    print(f"m_lim={m_lim}:  peak z={pk:.2f}   max dN/dz={y.max():.1f} /arcmin^2   "
          f"integrated N({Z_MIN:.1f}<z<10)={tot:.1f} /arcmin^2")

ax.set_ylim(0, 250)
ax.set_xlim(2, 5.2)
ax.set_xlabel("Redshift", fontsize=17)
ax.set_ylabel(r"d$N$/d$z$  (arcmin$^{-2}$)", fontsize=17)
ax.tick_params(axis="both", which="both", direction="in", labelsize=14)
ax.grid(True, which="major", color="0.88", lw=0.7)
ax.legend(fontsize=13, framealpha=0.95)
ax.set_title("Galaxy surface density from the rest-UV LF\n"
             "(Bouwens+2021; Harikane+2023)", fontsize=12)
fig.tight_layout()

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
os.makedirs(os.path.join(OUTDIR, "png"), exist_ok=True)
pdf = os.path.join(OUTDIR, "neff_vs_z.pdf")
png = os.path.join(OUTDIR, "png", "neff_vs_z.png")
fig.savefig(pdf); fig.savefig(png, dpi=140)
print("wrote", pdf); print("wrote", png)
