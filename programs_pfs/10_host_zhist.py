#!/usr/bin/env python3
"""
10_host_zhist.py -- redshift histogram of ELAIS-N1 SN Ia host galaxies, each
bar split three ways by host Roman Z-band magnitude:
  blue  = Z < 24,  green = 24 <= Z < 25.5,  red = Z >= 25.5.

    python 10_host_zhist.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "03_snia_host_radec_ELAIS-N1.csv")


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True)
    z = d["z_cmb"]
    zm = d["peakmag_Z"]                       # host HOSTGAL_MAG_Z
    v = (zm > 0) & (zm < 90)
    blue = v & (zm < 24.0)
    green = v & (zm >= 24.0) & (zm < 25.5)
    red = ~blue & ~green                     # Z >= 25.5 (incl. any undefined)
    print(f"ELAIS-N1 SN Ia hosts: {len(z)}  "
          f"(Z<24: {blue.sum()}, 24-25.5: {green.sum()}, >=25.5: {red.sum()})")

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

    bins = np.arange(0.0, 3.0 + 1e-9, 0.1)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.hist([z[blue], z[green], z[red]], bins=bins, stacked=True,
            color=["blue", "limegreen", "red"], edgecolor="white", linewidth=0.3,
            label=[f"Host $Z<24$  ($N={blue.sum():,}$)",
                   f"Host $24\\leq Z<25.5$  ($N={green.sum():,}$)",
                   f"Host $Z\\geq25.5$  ($N={red.sum():,}$)"])
    ax.set_xlabel("Redshift", fontsize=18)
    ax.set_ylabel("Number of Host Galaxies", fontsize=18)
    ax.set_title(f"ELAIS-N1 SN Ia Host Galaxies vs Redshift  ($N={len(z):,}$)", fontsize=16)
    ax.tick_params(labelsize=13)
    ax.set_xlim(0, 3.0)
    ax.legend(fontsize=13, loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    png = os.path.join(PNG_DIR, "10_host_zhist_ELAIS-N1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
