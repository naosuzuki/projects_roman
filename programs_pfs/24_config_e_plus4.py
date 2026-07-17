#!/usr/bin/env python3
"""
24_config_e_plus4.py -- Configuration E (12-pointing ring, r=1.5 deg) plus the
4 central pointings (square, d=0.6 deg) = 16 fields, as a complete tiling for the
remaining SN-host galaxies of ELAIS-N1. Shows the coverage map and the recovery
(hosts integrated to S/N=5) vs total 1-hour exposures.

Model: each round observes all 16 fields (16 exposures); a host covered by n of
them gains n hours/round (overlaps integrate faster) and completes when the
accumulated time reaches its remaining requirement.

Reads 07_program_sne_ELAIS-N1.csv.   python 24_config_e_plus4.py
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
R_FOOT = 2.09


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def config_centers():
    ring = [(R_RING * np.cos(np.radians(360 * k / M)),
             R_RING * np.sin(np.radians(360 * k / M))) for k in range(M)]
    central = [(0.60 * np.cos(np.radians(45 + 90 * k)),
                0.60 * np.sin(np.radians(45 + 90 * k))) for k in range(4)]
    return ring, central


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    vneed = np.asarray(d["host_visits_needed"]).astype(float)
    vdone = np.asarray(d["host_visits_done"]).astype(float)
    cosd = np.cos(np.radians(DEC0))
    xi = (hra - RA0) * cosd; eta = hdec - DEC0
    rneed = np.maximum(1, np.ceil(vneed - vdone)).astype(int)
    rem = (htar == 1) & (comp == 0)
    pts = np.column_stack([xi, eta])[rem]
    rn = rneed[rem]
    n = len(pts)

    ring, central = config_centers()
    allc = ring + central
    ncov = np.zeros(n, int)
    for xc, yc in allc:
        ncov += Path(hexagon(xc, yc)).contains_points(pts).astype(int)
    covered = ncov > 0
    print(f"remaining hosts: {n};  Config E(12)+4 central = {len(allc)} fields")
    print(f"  coverage: {covered.sum()} ({100*covered.sum()/n:.0f}%)")

    # recovery vs exposures: per round +ncov h; complete when rounds*ncov >= rneed
    rounds = np.arange(0, 26)
    rounds_needed = np.where(covered, np.ceil(rn / np.maximum(ncov, 1)), 1e9)
    rec = np.array([int((covered & (rounds_needed <= Rn)).sum()) for Rn in rounds])
    exp = len(allc) * rounds
    for Rn in (2, 4, 6, 8, 12):
        print(f"  {len(allc)*Rn:4d} exp ({Rn} rounds): {rec[Rn]:5d} ({100*rec[Rn]/n:.0f}%)")

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
    ra, dec = t2sky(pts[:, 0], pts[:, 1])

    fig, (axm, axr) = plt.subplots(1, 2, figsize=(15, 6.5))
    axm.scatter(ra[covered], dec[covered], s=5, c="#2ca02c", alpha=0.5, linewidths=0, label="Covered")
    axm.scatter(ra[~covered], dec[~covered], s=6, c="#d62728", alpha=0.6, linewidths=0, label="Missed")
    def draw_map(ax, lfs=16, tfs=12, tts=14, lgs=11):
        ax.scatter(ra[covered], dec[covered], s=5, c="#2ca02c", alpha=0.5,
                   linewidths=0, label="Covered")
        ax.scatter(ra[~covered], dec[~covered], s=6, c="#d62728", alpha=0.6,
                   linewidths=0, label="Missed")
        for grp, col, lab in ((ring, "#1f6fe0", "Config E ring (12)"),
                              (central, "#9467bd", "Config G central (4)")):
            first = True
            for xc, yc in grp:
                sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
                ax.add_patch(Polygon(sky, closed=True, fill=False, edgecolor=col,
                                     lw=1.7, alpha=0.9, label=(lab if first else None)))
                first = False
        th = np.linspace(0, 2 * np.pi, 240)
        fra, fdec = t2sky(R_FOOT * np.cos(th), R_FOOT * np.sin(th))
        ax.plot(fra, fdec, color="0.4", ls="--", lw=1.1, label="Footprint edge")
        ax.set_aspect(1.0 / cosd); ax.invert_xaxis()
        ax.set_xlabel("RA (deg)", fontsize=lfs); ax.set_ylabel("Dec (deg)", fontsize=lfs)
        ax.set_title(f"Config E (12) + Config G (4 Central) = 16 Fields  "
                     f"(coverage {100*covered.sum()/n:.0f}%)", fontsize=tts)
        ax.tick_params(labelsize=tfs); ax.legend(fontsize=lgs, loc="upper right")
        ax.grid(True, alpha=0.25)

    # (points are drawn inside draw_map; clear the pre-drawn ones)
    axm.clear()
    draw_map(axm)

    axr.plot(exp, 100 * rec / n, "-o", color="#1f6fe0", ms=4)
    axr.axhline(100 * covered.sum() / n, color="0.6", ls=":", lw=1,
                label=f"Coverage ceiling {100*covered.sum()/n:.0f}\\%")
    axr.set_xlabel("Total 1-Hour Exposures", fontsize=16)
    axr.set_ylabel("Remaining Hosts Recovered (%)", fontsize=16)
    axr.set_title("Time Evolution: Hosts Recovered vs Exposure", fontsize=15)
    axr.set_xlim(0, 100); axr.set_ylim(0, 100)
    axr.tick_params(labelsize=12); axr.legend(fontsize=12, loc="lower right"); axr.grid(True, alpha=0.3)

    png = os.path.join(PNG_DIR, "24_config_e_plus4_ELAIS-N1.png")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)

    # standalone version of the coverage map: "Configuration alpha", painted
    # blue (Fig 13 style: light fill + solid outline), bigger host points
    fig2, axs = plt.subplots(figsize=(9, 9))
    axs.scatter(ra[covered], dec[covered], s=11, c="#2ca02c", alpha=0.6,
                linewidths=0, label="Covered")
    axs.scatter(ra[~covered], dec[~covered], s=13, c="#d62728", alpha=0.75,
                linewidths=0, label="Missed")
    first = True
    for xc, yc in ring + central:
        sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
        axs.add_patch(Polygon(sky, closed=True, facecolor="#1f6fe0",
                              edgecolor="none", alpha=0.10))
        axs.add_patch(Polygon(sky, closed=True, fill=False, edgecolor="#1f6fe0",
                              lw=1.8, alpha=0.9,
                              label=("Configuration $\\alpha$ (16)" if first else None)))
        first = False
    th2 = np.linspace(0, 2 * np.pi, 240)
    fra2, fdec2 = t2sky(R_FOOT * np.cos(th2), R_FOOT * np.sin(th2))
    axs.plot(fra2, fdec2, color="0.4", ls="--", lw=1.1, label="Footprint edge")
    axs.set_aspect(1.0 / cosd); axs.invert_xaxis()
    axs.set_xlabel("RA (deg)", fontsize=20); axs.set_ylabel("Dec (deg)", fontsize=20)
    axs.set_title(f"Configuration $\\alpha$ = E Ring (12) + G Central (4)  "
                  f"(Coverage {100*covered.sum()/n:.0f}%)", fontsize=17)
    axs.tick_params(labelsize=15); axs.legend(fontsize=14, loc="upper right")
    axs.grid(True, alpha=0.25)
    png2 = os.path.join(PNG_DIR, "24_config_e_plus4_map_ELAIS-N1.png")
    fig2.tight_layout(); fig2.savefig(png2, dpi=140)
    print("plot ->", png2)


if __name__ == "__main__":
    main()
