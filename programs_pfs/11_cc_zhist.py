#!/usr/bin/env python3
"""
11_cc_zhist.py -- redshift histogram of ELAIS-N1 core-collapse (non-Ia) SNe,
each bar split by peak Roman Z-band magnitude: blue = Z < 24 (PFS-observable),
red = Z >= 24 (too faint at peak; undefined/high-z sentinels included).

    python 11_cc_zhist.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "03_cc_radec_ELAIS-N1.csv")


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True)
    z = d["z_cmb"]
    zm = d["peakmag_Z"]
    blue = (zm > 0) & (zm < 24.0)        # valid peak Z < 24
    red = ~blue                          # Z >= 24 or undefined (-9 sentinel)
    print(f"ELAIS-N1 CC SNe: {len(z)}  (Z<24: {blue.sum()}, Z>=24/undef: {red.sum()})")

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
    ax.hist([z[blue], z[red]], bins=bins, stacked=True,
            color=["blue", "red"], edgecolor="white", linewidth=0.3,
            label=[f"Peak $Z<24$  ($N={blue.sum():,}$)",
                   f"Peak $Z\\geq24$  ($N={red.sum():,}$)"])
    ax.set_xlabel("Redshift", fontsize=18)
    ax.set_ylabel("Number of CC SNe", fontsize=18)
    ax.set_title(f"ELAIS-N1 Core-Collapse SNe vs Redshift  ($N={len(z):,}$)", fontsize=16)
    ax.tick_params(labelsize=13)
    ax.set_xlim(0, 3.0)
    ax.legend(fontsize=13, loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    png = os.path.join(PNG_DIR, "11_cc_zhist_ELAIS-N1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
