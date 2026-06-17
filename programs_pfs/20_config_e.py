#!/usr/bin/env python3
"""
20_config_e.py -- Configuration E: an efficient overlapping ring of pointings to
recover the OUTSKIRT remaining SN-host galaxies of ELAIS-N1 (those beyond the
A/B core flower), without re-pointing the centre.

Config C's outer ring sat at r ~ 2.25 deg, half outside the footprint (~50%
filled). The outskirt hosts actually live in a thin annulus r = 1.6-2.2 deg
(footprint edge 2.09 deg). Config E places 12 pointings on a ring at r = 1.6 deg
so each FoV is ~95% inside the footprint, overlapping so the annulus is fully
tiled and hosts in the overlaps integrate faster.

Reads 07_program_sne_ELAIS-N1.csv. Produces the layout + recovery figure.

    python 20_config_e.py
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
R_RING = 1.6        # config E ring radius [deg]
M = 12              # config E pointings
R_FOOT = 2.09       # footprint radius (99th pct of host targets) [deg]


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def core_centers():            # 7-pointing A/B flower
    sp = np.sqrt(3) * R
    c = [(0.0, 0.0)]
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang); c.append((sp * np.cos(a), sp * np.sin(a)))
    return c


def e_centers():               # config E ring
    return [(R_RING * np.cos(2 * np.pi * k / M), R_RING * np.sin(2 * np.pi * k / M))
            for k in range(M)]


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    vneed = np.asarray(d["host_visits_needed"]).astype(float)
    vdone = np.asarray(d["host_visits_done"]).astype(float)
    pAh = np.asarray(d["pA_host"]).astype(int)
    pBh = np.asarray(d["pB_host"]).astype(int)
    cosd = np.cos(np.radians(DEC0))
    xi = (hra - RA0) * cosd
    eta = hdec - DEC0
    rneed = np.maximum(1, np.ceil(vneed - vdone)).astype(int)
    rem = (htar == 1) & (comp == 0)
    outer = rem & ~((pAh >= 0) | (pBh >= 0))      # outskirt (uncovered) remaining hosts
    pts = np.column_stack([xi, eta])

    ec = e_centers()
    # number of E-pointings covering each outskirt host (overlap -> faster integ.)
    ncov = np.zeros(len(pts), int)
    for (xc, yc) in ec:
        ncov += Path(hexagon(xc, yc)).contains_points(pts).astype(int)
    covered = outer & (ncov > 0)
    print(f"outskirt remaining hosts: {outer.sum()};  Config E covers "
          f"{covered.sum()} ({100*covered.sum()/outer.sum():.0f}%)")

    # fill fraction of each E pointing (FoV grid inside footprint disk)
    gx, gy = np.meshgrid(np.linspace(-R, R, 41), np.linspace(-R, R, 41))
    gx, gy = gx.ravel(), gy.ravel()
    inhex = Path(hexagon(0, 0)).contains_points(np.column_stack([gx, gy]))
    gx, gy = gx[inhex], gy[inhex]
    fills = [np.mean(np.hypot(gx + xc, gy + yc) <= R_FOOT) for (xc, yc) in ec]
    print(f"  mean FoV fill inside footprint: {100*np.mean(fills):.0f}%  "
          f"(vs ~50% for Config C's outer ring at r=2.25)")

    # recovery vs exposures: per round each E pointing adds 1 h; an overlap host
    # gets ncov hours/round -> completes after ceil(rneed/ncov) rounds. 12 exp/round.
    rounds_needed = np.where(covered, np.ceil(rneed / np.maximum(ncov, 1)), 1e9)
    Rmax = 40
    rounds = np.arange(0, Rmax + 1)
    rec = np.array([int((covered & (rounds_needed <= r)).sum()) for r in rounds])
    exp = M * rounds
    for r in (3, 5, 8, 12):
        print(f"  E depth {r:2d} rounds = {M*r:3d} exp -> recover {rec[r]:4d} "
              f"({100*rec[r]/outer.sum():.0f}% of outskirts)")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Circle
    from matplotlib import font_manager as fm
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, (axm, axr) = plt.subplots(1, 2, figsize=(15, 6.5))

    def t2sky(x, y):
        return RA0 + x / cosd, DEC0 + y
    ra, dec = t2sky(xi, eta)
    axm.scatter(ra[rem & ~outer], dec[rem & ~outer], s=5, c="0.78", alpha=0.5,
                linewidths=0, label="Remaining (core)")
    axm.scatter(ra[outer], dec[outer], s=6, c="#d62728", alpha=0.55,
                linewidths=0, label="Remaining (outskirts)")
    for (xc, yc) in core_centers():
        sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
        axm.add_patch(Polygon(sky, closed=True, fill=False, edgecolor="#1f6fe0", lw=1.4, alpha=0.7))
    first = True
    for (xc, yc) in ec:
        sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
        axm.add_patch(Polygon(sky, closed=True, fill=False, edgecolor="#2ca02c",
                              lw=1.8, alpha=0.95, label=("Config E ring (12)" if first else None)))
        first = False
    fra, fdec = t2sky(R_FOOT, 0)
    axm.add_patch(Circle((RA0, DEC0), R_FOOT / cosd, fill=False, ec="0.4", ls="--", lw=1))
    axm.set_aspect(1.0 / cosd); axm.invert_xaxis()
    axm.set_xlabel("RA (deg)", fontsize=16); axm.set_ylabel("Dec (deg)", fontsize=16)
    axm.set_title(f"Configuration E: Overlapping Outskirt Ring (r={R_RING}$^\\circ$, {M} pointings)",
                  fontsize=14)
    axm.tick_params(labelsize=12); axm.legend(fontsize=11, loc="upper right"); axm.grid(True, alpha=0.25)

    axr.plot(exp, 100 * rec / outer.sum(), "-o", color="#2ca02c", ms=4)
    axr.axhline(100 * covered.sum() / outer.sum(), color="0.6", ls=":", lw=1,
                label=f"Coverage ceiling {100*covered.sum()/outer.sum():.0f}%")
    for r, lab in ((5, "E: 60 exp"),):
        axr.axvline(M * r, color="0.6", ls=":", lw=1)
        axr.text(M * r, 3, f"{M*r} exp", rotation=90, va="bottom", ha="right", fontsize=11, color="0.3")
    axr.set_xlabel("Total 1-Hour Exposures (Config E)", fontsize=16)
    axr.set_ylabel("Outskirt Hosts Recovered (%)", fontsize=16)
    axr.set_title("Config E Recovery vs Exposure", fontsize=15)
    axr.set_xlim(0, M * Rmax); axr.set_ylim(0, 100)
    axr.tick_params(labelsize=12); axr.legend(fontsize=12, loc="lower right"); axr.grid(True, alpha=0.3)

    png = os.path.join(PNG_DIR, "20_config_e_ELAIS-N1.png")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
