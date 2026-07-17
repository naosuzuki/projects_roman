#!/usr/bin/env python3
"""
13_snia_specz_zhist_v2.py -- as 13_snia_specz_zhist.py (redshift histogram of
SN Ia with a successful live PFS spec-z, for the 100/80/70% weather cases), but
with the red and blue colors SWAPPED: 100% = red, 70% = blue (80% stays green).

    python 13_snia_specz_zhist_v2.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "07_program_sne_ELAIS-N1.csv")

CASES = [(1.0, "100% (no weather loss)", "red"),
         (0.8, "80% weather",            "g"),
         (0.7, "70% weather",            "blue")]


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    typ = np.asarray(d["type"]).astype(str)
    obs = np.asarray(d["observed"]).astype(int)
    z = np.asarray(d["z"]).astype(float)
    m = (typ == "Ia") & (obs == 1)
    zo = z[m]
    print(f"successful SN Ia spec-z: {len(zo)}  (z {zo.min():.3f}-{zo.max():.3f})")

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

    bins = np.arange(0.0, np.ceil(zo.max() * 10) / 10 + 0.051, 0.05)
    counts, edges = np.histogram(zo, bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = (edges[1] - edges[0]) * 0.93         # slightly narrower than the bin -> small gap
    fig, ax = plt.subplots(figsize=(9, 6))
    # nested bars: draw 100% first, then 80%, then 70% on top, so each bin shows
    # blue (0-70%), a green band (70-80%) and a red band (80-100%)
    for i, (w, lab, col) in enumerate(CASES):
        n = w * len(zo)
        ax.bar(centers, counts * w, width=bw, color=col, edgecolor="white",
               linewidth=0.4, zorder=i + 1, label=f"{lab}  ($N={n:.0f}$)")

    ax.set_xlabel("Redshift", fontsize=18)
    ax.set_ylabel("Number of SN Ia (Successful Spec-$z$)", fontsize=18)
    ax.set_title("ELAIS-N1 SN Ia with Successful PFS Spec-$z$ vs Redshift",
                 fontsize=16)
    ax.tick_params(labelsize=13)
    ax.set_xlim(0, bins[-1])
    ax.legend(fontsize=13, loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)

    png = os.path.join(PNG_DIR, "13_snia_specz_zhist_v2_ELAIS-N1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
