#!/usr/bin/env python3
"""
08_hsc_etc.py -- Subaru/HSC imaging exposure-time calculator.

Point-source S/N (native 1x1 pixels) as a function of total exposure time,
for HSC g / r2 / i2, at source AB magnitudes 24-28. One figure per band, with
60/90/120-min guides; PFS-aligned plot style.

Empirical sky-background-limited model anchored to the observed depth: with the
5-sigma 1-hour limiting magnitude m5_1h,
    S/N(mag, t) = 5 * 10^(-0.4 (mag - m5_1h)) * sqrt(t / 1 h).
m5_1h is set to ~26 (S/N=5 for a 26th-mag point source in 1 h), matching HSC
observing experience; edit it per band at the top.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")

# --- telescope / detector ---
# Empirical sky-background-limited imaging model. The point-source S/N is
# anchored to the 5-sigma limiting AB magnitude reached in 1 hour (m5_1h) and
# scales as  S/N(mag, t) = 5 * 10^(-0.4 (mag - m5_1h)) * sqrt(t / 1 h),
# i.e. S/N is proportional to source flux and to sqrt(exposure time). m5_1h is
# anchored to ~26 mag at S/N = 5 in 1 h (HSC observing experience).
HSC = {
    "g":  dict(lam=477.0, m5_1h=26.4, color="#d62728"),
    "r2": dict(lam=623.0, m5_1h=26.2, color="#1f77b4"),
    "i2": dict(lam=771.0, m5_1h=26.0, color="#2ca02c"),
}
MAGS = [24.0, 25.0, 26.0, 27.0, 28.0]
VLINES = [60, 90, 120]       # minutes


def sn(mag, t_min, b):
    """Point-source imaging S/N (sky-background limited)."""
    return 5.0 * 10 ** (-0.4 * (mag - b["m5_1h"])) * np.sqrt(np.asarray(t_min) / 60.0)


def main():
    os.makedirs(PNG_DIR, exist_ok=True); os.makedirs(CSV_DIR, exist_ok=True)
    print("HSC imaging ETC (empirical, sky-limited; "
          "S/N = 5*10^(-0.4(m-m5))*sqrt(t/1h))")
    csv = os.path.join(CSV_DIR, "08_hsc_etc.csv")
    with open(csv, "w") as fo:
        fo.write("band,mag,t_min,SN\n")
        for bn, b in HSC.items():
            print(f"  {bn}: 5sigma(1h) depth = {b['m5_1h']:.1f};  "
                  f"mag 26 in 60 min -> S/N {sn(26.0, 60.0, b):.1f}")
            for mag in MAGS:
                for tm in (10, 20, 30, 60, 90, 120, 240, 360):
                    fo.write(f"{bn},{mag},{tm},{sn(mag, tm, b):.3f}\n")
    print("table ->", csv)

    # --- plots: one per band ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    from matplotlib.ticker import FixedLocator, FuncFormatter
    cmap = plt.get_cmap("rainbow_r")                          # same as PFS plots
    cols = [cmap(x) for x in np.linspace(0, 1, len(MAGS))]    # 24 red -> 28 violet
    tmk = np.array([5, 7.5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 240, 300, 360])
    tref = np.logspace(np.log10(tmk.min()), np.log10(tmk.max()), 100)
    yticks = [0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500]

    for bn, b in HSC.items():
        fig, ax = plt.subplots(figsize=(8, 6))
        for mag, c in zip(MAGS, cols):
            ax.plot(tmk, sn(mag, tmk, b), "o-", color=c, lw=1.8,
                    ms=5, label=f"{bn} = {mag:.0f} AB")
        # sqrt(t) reference anchored to the middle magnitude
        mref = MAGS[len(MAGS) // 2]
        s0 = sn(mref, tmk[0], b)
        ax.plot(tref, s0 * np.sqrt(tref / tmk[0]), "k--", lw=1.2, alpha=0.7,
                label=r"$\propto\sqrt{t}$ (bkg-limited)")
        for tv in VLINES:
            ax.axvline(tv, color="0.35", ls=":", lw=1.3)
            ax.text(tv, 0.985, f"{tv} min", transform=ax.get_xaxis_transform(),
                    rotation=90, va="top", ha="right", fontsize=13, color="0.3")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.yaxis.set_major_locator(FixedLocator(yticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_xlabel("Total Exposure Time  [min]", fontsize=17)
        ax.set_ylabel("Imaging S/N (Point Source)", fontsize=17)
        ax.set_title(f"Subaru/HSC Imaging ETC — {bn}-band S/N vs Exposure Time (log–log)\n"
                     f"5$\\sigma$ Depth $\\approx$ {b['m5_1h']:.1f} mag in 1 h, "
                     f"$\\lambda\\approx{b['lam']:.0f}$ nm", fontsize=15)
        ax.tick_params(labelsize=13)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(frameon=True, fontsize=12, loc="upper left")
        png = os.path.join(PNG_DIR, f"08_hsc_etc_{bn}band.png")
        fig.tight_layout(); fig.savefig(png, dpi=140)
        print("plot ->", png)


if __name__ == "__main__":
    main()
