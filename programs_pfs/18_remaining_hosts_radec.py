#!/usr/bin/env python3
"""
18_remaining_hosts_radec.py -- sky distribution (RA, Dec) of the remaining
(incomplete) SN-host galaxies at the end of the ELAIS-N1 survey, color-coded by
host Roman Z magnitude (bright = blue, faint = red).

Remaining = Phase III host target (Z < 25.5) not integrated to S/N = 5 by the
survey end. Reads the program-SN catalog written by 07_fiber_budget.py.

    python 18_remaining_hosts_radec.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "07_program_sne_ELAIS-N1.csv")


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    remaining = (htar == 1) & (comp == 0)
    ra = np.asarray(d["host_ra"]).astype(float)[remaining]
    dec = np.asarray(d["host_dec"]).astype(float)[remaining]
    zmag = np.asarray(d["host_Z"]).astype(float)[remaining]
    print(f"remaining hosts: {remaining.sum()}  "
          f"(Z {zmag.min():.2f}-{zmag.max():.2f}, median {np.median(zmag):.2f})")

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

    cosd = np.cos(np.radians(np.median(dec)))
    fig, ax = plt.subplots(figsize=(9, 8))
    # bright (low Z) -> blue, faint (high Z) -> red
    sc = ax.scatter(ra, dec, c=zmag, cmap="coolwarm", vmin=20.5, vmax=25.5,
                    s=11, alpha=0.7, linewidths=0)
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Host Roman $Z$ (AB Mag)", fontsize=16)
    cb.ax.tick_params(labelsize=12)

    ax.set_aspect(1.0 / cosd)
    ax.invert_xaxis()
    ax.set_xlabel("RA (deg)", fontsize=18)
    ax.set_ylabel("Dec (deg)", fontsize=18)
    ax.set_title(f"ELAIS-N1 Remaining (Incomplete) SN-Host Galaxies  "
                 f"($N={remaining.sum()}$)", fontsize=15)
    ax.tick_params(labelsize=13)
    ax.grid(True, alpha=0.25)

    png = os.path.join(PNG_DIR, "18_remaining_hosts_radec_ELAIS-N1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
