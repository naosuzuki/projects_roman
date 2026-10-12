#!/usr/bin/env python3
"""
02_aat_success.py -- AAT 2dF/AAOmega redshift-success "ETC".

The AAT counterpart of the PFS ETC (../programs_pfs/01,02). Instead of a
continuum S/N, it predicts the probability of a secure redshift (QOP = 3, 4)
as a function of magnitude and exposure time, calibrated on the HSC-AAT host
campaign (~/github/papers/24_hscsnhost_aat, Fig. 2; data
~/github/projects_hsc/specz/csvfiles_aat/{g,r,i,z,y}_magfraction.csv):

  * Template: a logistic fit to the measured HSC-AAT success rate per HSC
    band at the run's effective exposure T_REF = 2 h per pointing,
        S(m, T_REF) = 1 / (1 + exp((m - m50) / w)).
    HSC-z is used for Roman Z087 (~877 nm), the magnitude tabulated by the
    HOURGLASS2 simulation for both SNe and hosts.
  * Exposure scaling (background limited, as in the HSC-AAT paper, which
    shows it reconciles the 2 h HSC-AAT curve with the ~5-10x longer OzDES
    one): at fixed success the limiting magnitude shifts by
        dm = 1.25 log10(t / T_REF),  i.e.  S(m, t) = S(m - dm, T_REF).

    python 02_aat_success.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
HSCAAT_DIR = os.path.expanduser("~/github/projects_hsc/specz/csvfiles_aat")

T_REF = 2.0        # h, effective exposure per pointing of the HSC-AAT run
FIT_RANGE = (18.0, 24.5)   # fainter bins are small-number / mis-ID dominated
# logistic fit to HSC-z (Z087 proxy), frozen here so other scripts can import;
# main() refits and warns if these drift
M50_Z, W_Z = 21.30, 0.97


def logistic(m, m50, w):
    return 1.0 / (1.0 + np.exp((np.asarray(m, dtype=float) - m50) / w))


def success(mag, hours, m50=M50_Z, w=W_Z):
    """Probability of a secure (QOP>=3) redshift for magnitude `mag` after `hours`."""
    hours = np.maximum(np.asarray(hours, dtype=float), 1e-6)
    return logistic(np.asarray(mag) - 1.25 * np.log10(hours / T_REF), m50, w)


def hours_for(mag, frac, m50=M50_Z, w=W_Z):
    """Exposure (h) to reach success probability `frac` at magnitude `mag`."""
    m_frac = m50 + w * np.log(1.0 / frac - 1.0)      # S(m_frac, T_REF) = frac
    return T_REF * 10 ** ((np.asarray(mag) - m_frac) / 1.25)


def fit_logistic(mag, frac):
    sel = (mag >= FIT_RANGE[0]) & (mag <= FIT_RANGE[1]) & np.isfinite(frac)
    best = (np.inf, None, None)
    for m50 in np.arange(19.0, 25.0, 0.01):
        for w in np.arange(0.2, 2.0, 0.01):
            chi = np.sum((logistic(mag[sel], m50, w) - frac[sel]) ** 2)
            if chi < best[0]:
                best = (chi, m50, w)
    return best[1], best[2]


def main():
    import csv
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    meas = {}
    for b in "grizy":
        fn = os.path.join(HSCAAT_DIR, f"{b}_magfraction.csv")
        if os.path.exists(fn):
            d = np.genfromtxt(fn, delimiter=",", names=True)
            meas[b] = (d["mag"], d[f"{b}_fraction"])
    fits = {b: fit_logistic(*meas[b]) for b in meas}
    print(f"logistic fits to HSC-AAT success at T_REF = {T_REF} h:")
    for b, (m50, w) in fits.items():
        print(f"  HSC-{b}: m50 = {m50:.2f}  w = {w:.2f}")
    if "z" in fits and (abs(fits["z"][0] - M50_Z) > 0.02 or abs(fits["z"][1] - W_Z) > 0.02):
        print(f"  NOTE: refit z = {fits['z']} differs from frozen ({M50_Z}, {W_Z}); update")

    Zs = np.array([20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0])
    csvf = os.path.join(CSV_DIR, "02_aat_success.csv")
    with open(csvf, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["Z", "S_1h", "S_2h", "S_4h", "S_8h", "S_16h", "hours_50pct", "hours_80pct"])
        print("\n   Z    S(1h) S(2h) S(4h) S(8h) S(16h)   t50[h]   t80[h]")
        for z in Zs:
            s = [float(success(z, t)) for t in (1, 2, 4, 8, 16)]
            t50, t80 = float(hours_for(z, 0.5)), float(hours_for(z, 0.8))
            w.writerow([z] + [f"{v:.3f}" for v in s] + [f"{t50:.2f}", f"{t80:.2f}"])
            print(f"  {z:4.1f}  " + " ".join(f"{v:5.2f}" for v in s) + f"   {t50:7.1f}  {t80:7.1f}")
    print("table ->", csvf)

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, ax = plt.subplots(figsize=(9, 6.5), constrained_layout=True)
    hours = [1, 2, 4, 8, 16, 32]
    cmap = plt.get_cmap("rainbow_r")
    m = np.linspace(18, 27, 901)
    for i, t in enumerate(hours):
        ax.plot(m, success(m, t), color=cmap(i / (len(hours) - 1)), lw=2.4,
                label=f"{t} h")
    if "z" in meas:
        mz, fz = meas["z"]
        ax.plot(mz, fz, "ks", ms=7, mfc="none", mew=1.5,
                label="HSC-AAT, HSC-$z$ (2 h)")
    ax.axhline(0.5, color="0.5", ls=":", lw=1)
    ax.axhline(0.8, color="0.5", ls=":", lw=1)
    ax.set_xlim(18, 27)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Roman Z087 Magnitude (AB)", fontsize=18)
    ax.set_ylabel("Redshift Success Rate (QOP $\\geq$ 3)", fontsize=18)
    ax.set_title("AAT 2dF/AAOmega Redshift Success vs Exposure Time", fontsize=18)
    ax.tick_params(labelsize=14, direction="in", which="both")
    ax.legend(fontsize=13, title="Exposure", title_fontsize=13, loc="lower left")
    ax.grid(True, alpha=0.3)
    png = os.path.join(PNG_DIR, "02_aat_success_zband.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
