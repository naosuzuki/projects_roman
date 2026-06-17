#!/usr/bin/env python3
"""
23_central_fill.py -- with Config E+F at r=1.5 deg (24 pointings) leaving a
central hexagonal hole (radius ~0.93 deg), find the MINIMUM number of PFS
pointings to fill that hole and observe the central remaining hosts.

Results (central-hole remaining hosts, r<0.93): n=1 center -> 55%; n=3 triangle
(d=0.43) -> 95%; n=4 square (d=0.60) -> 99%. Minimum useful fill is 3 pointings;
4 gives near-complete coverage.

Reads 07_program_sne_ELAIS-N1.csv. Produces a 3-example comparison figure.

    python 23_central_fill.py
"""
import os
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "07_program_sne_ELAIS-N1.csv")
RA0, DEC0 = 242.498, 54.497
FOV = 1.3
R = FOV / 2.0
M = 12
R_RING = 1.5


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def ef_centers():
    return [(R_RING * np.cos(np.radians(p + 360 * k / M)),
             R_RING * np.sin(np.radians(p + 360 * k / M)))
            for p in (0.0, 15.0) for k in range(M)]


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    cosd = np.cos(np.radians(DEC0))
    xi = (hra - RA0) * cosd; eta = hdec - DEC0
    r = np.hypot(xi, eta)
    rem = (htar == 1) & (comp == 0)
    pts = np.column_stack([xi, eta])
    EF = ef_centers()
    covEF = np.zeros(len(pts), bool)
    for xc, yc in EF:
        covEF |= Path(hexagon(xc, yc)).contains_points(pts)
    hole = rem & ~covEF & (r < 1.0)
    hp = pts[hole]; nh = int(hole.sum())

    # candidate fills
    def tri(d): return [(d * np.cos(np.radians(90 + 120 * k)), d * np.sin(np.radians(90 + 120 * k))) for k in range(3)]
    def sq(d): return [(d * np.cos(np.radians(45 + 90 * k)), d * np.sin(np.radians(45 + 90 * k))) for k in range(4)]
    EXAMPLES = [("1 pointing (centre)", [(0.0, 0.0)]),
                ("3 pointings (triangle, $d=0.43^\\circ$)", tri(0.425)),
                ("Config G: 4 pointings (square, $d=0.60^\\circ$)", sq(0.60))]

    def covmask(cens):
        c = np.zeros(nh, bool)
        for xc, yc in cens:
            c |= Path(hexagon(xc, yc)).contains_points(hp)
        return c
    print(f"central-hole remaining hosts (r<1.0, not in E+F): {nh}")
    for lab, cens in EXAMPLES:
        print(f"  {lab}: {100*covmask(cens).mean():.0f}%")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib import font_manager as fm
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

    def t2sky(x, y):
        return RA0 + x / cosd, DEC0 + y

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    hra_, hdec_ = t2sky(hp[:, 0], hp[:, 1])
    for ax, (lab, cens) in zip(axes, EXAMPLES):
        cm = covmask(cens)
        # E+F ring (light context)
        for xc, yc in EF:
            sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
            ax.add_patch(Polygon(sky, closed=True, fill=False, edgecolor="0.8", lw=0.8))
        ax.scatter(hra_[cm], hdec_[cm], s=7, c="#2ca02c", alpha=0.6, linewidths=0, label="Covered")
        ax.scatter(hra_[~cm], hdec_[~cm], s=9, c="#d62728", alpha=0.8, linewidths=0, label="Missed")
        for xc, yc in cens:
            sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
            ax.add_patch(Polygon(sky, closed=True, fill=False, edgecolor="#1f6fe0", lw=2.0))
        ax.set_aspect(1.0 / cosd); ax.invert_xaxis()
        ax.set_xlim(RA0 + 1.65 / cosd, RA0 - 1.65 / cosd)
        ax.set_ylim(DEC0 - 1.65, DEC0 + 1.65)
        ax.set_xlabel("RA (deg)", fontsize=14)
        ax.set_title(f"{lab}\n{100*cm.mean():.0f}\\% of {nh} central hosts", fontsize=13)
        ax.tick_params(labelsize=11); ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("Dec (deg)", fontsize=14)
    axes[0].legend(fontsize=11, loc="upper right")

    png = os.path.join(PNG_DIR, "23_central_fill_ELAIS-N1.png")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
