#!/usr/bin/env python3
"""
26_snia_specz_sn_host_zhist.py -- combined redshift histogram: the live SN Ia
PFS spec-z distribution for the 100/80/70% weather cases (as Fig 14 /
13_snia_specz_zhist.py), OVERLAID with the SN Ia host-galaxy spec-z distribution
for the 80% weather case (from Fig 15 / 14_snia_host_specz_zhist.py) in light
green with alpha=0.5.

Reads the program-SN catalog written by 07_fiber_budget.py:
  live spec-z : `observed`=1  (caught while Z<24)
  host spec-z : `host_completed`=1 (integrated to S/N=5 at 10-pix)

    python 26_snia_specz_sn_host_zhist.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "07_program_sne_ELAIS-N1.csv")

CASES = [(1.0, "SNIa 100% (no weather loss)", "blue"),
         (0.8, "SNIa 80% weather",            "g"),
         (0.7, "SNIa 70% weather",            "red")]
W_HOST = 0.8                                  # host case: 80% weather


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    typ = np.asarray(d["type"]).astype(str)
    obs = np.asarray(d["observed"]).astype(int)
    hcomp = np.asarray(d["host_completed"]).astype(int)
    z = np.asarray(d["z"]).astype(float)
    ia = typ == "Ia"
    zo = z[ia & (obs == 1)]                    # live SN spec-z
    zh = z[ia & (hcomp == 1)]                  # host spec-z
    print(f"live SN Ia spec-z: {len(zo)};  host spec-z: {len(zh)} "
          f"(x{W_HOST} -> {W_HOST*len(zh):.0f})")

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

    zmax = max(zo.max(), zh.max())
    bins = np.arange(0.0, np.ceil(zmax * 10) / 10 + 0.051, 0.05)
    counts, edges = np.histogram(zo, bins)
    hcounts, _ = np.histogram(zh, bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = (edges[1] - edges[0]) * 0.93         # slightly narrower than the bin -> small gap
    fig, ax = plt.subplots(figsize=(9, 6))
    # host 80% case in the BACKGROUND (lowest zorder), light green alpha=0.6;
    # the opaque SN bars are drawn on top, so the SN colors stay pure.
    nh = W_HOST * len(zh)
    ax.bar(centers, hcounts * W_HOST, width=bw, color="lightgreen", alpha=0.6,
           edgecolor="white", linewidth=0.4, zorder=1,
           label=f"Host 80% weather  ($N={nh:.0f}$)")
    # nested SN bars: 100% first, then 80%, then 70% on top (as Fig 14)
    for i, (w, lab, col) in enumerate(CASES):
        n = w * len(zo)
        ax.bar(centers, counts * w, width=bw, color=col, edgecolor="white",
               linewidth=0.4, zorder=i + 2, label=f"{lab}  ($N={n:.0f}$)")

    ax.set_xlabel("Redshift", fontsize=21)
    ax.set_ylabel("Expected Number of Successful Spec-$z$", fontsize=21)
    ax.set_title("ELAIS-N1 SN Ia and Host-Galaxy PFS Spec-$z$ vs Redshift",
                 fontsize=18)
    ax.tick_params(labelsize=16)
    ax.set_xlim(0, bins[-1])
    ax.legend(fontsize=14, loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)

    png = os.path.join(PNG_DIR, "26_snia_specz_sn_host_zhist_ELAIS-N1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
