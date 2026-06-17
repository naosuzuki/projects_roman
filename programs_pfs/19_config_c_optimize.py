#!/usr/bin/env python3
"""
19_config_c_optimize.py -- design a follow-up pointing configuration ("C", with
optional deep tier "D") to recover the remaining (incomplete) SN-host galaxies
of the ELAIS-N1 survey, maximizing hosts recovered per 1-hour exposure.

Finding (see 07): of 2270 remaining hosts, ~1265 sit in the OUTER ANNULUS beyond
the single 7-pointing flower (A/B), so rotating the flower does not help. Adding
a second ring of 12 pointings -> a 19-point hexagonal mosaic ("config C") covers
100% of them. Integration is steeply front-loaded (median 2 h more), so a
shallow/moderate pass recovers most; the faint Z>24.5 tail (config D, deep) is
not cost-effective.

Reads 07_program_sne_ELAIS-N1.csv. Produces the mosaic + recovery-curve figure.

    python 19_config_c_optimize.py
"""
import os
import heapq
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
CSV = os.path.join(CSV_DIR, "07_program_sne_ELAIS-N1.csv")
RA0, DEC0 = 242.498, 54.497
FOV = 1.3
R = FOV / 2.0
SP = np.sqrt(3) * R


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def mosaic_centers():
    """19-point flower: ring0 (center) + ring1 (6) = A/B core; ring2 (12) = config C."""
    c = [(0.0, 0.0)]
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang); c.append((SP * np.cos(a), SP * np.sin(a)))
    outer = []
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang); outer.append((2 * SP * np.cos(a), 2 * SP * np.sin(a)))
    for ang in (60, 120, 180, 240, 300, 0):
        a = np.radians(ang); outer.append((np.sqrt(3) * SP * np.cos(a), np.sqrt(3) * SP * np.sin(a)))
    return c, outer            # 7 core, 12 outer


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    vneed = np.asarray(d["host_visits_needed"]).astype(float)
    vdone = np.asarray(d["host_visits_done"]).astype(float)
    rem = (htar == 1) & (comp == 0)
    cosd = np.cos(np.radians(DEC0))
    xi = ((hra - RA0) * cosd)[rem]
    eta = (hdec - DEC0)[rem]
    rneed = np.maximum(1, np.ceil(vneed - vdone)).astype(int)[rem]    # >=1 more hour
    pts = np.column_stack([xi, eta])
    n = len(pts)

    core, outer = mosaic_centers()
    allc = core + outer
    assign = np.full(n, -1, int)
    for i, (xc, yc) in enumerate(allc):
        inside = Path(hexagon(xc, yc)).contains_points(pts)
        assign[(assign < 0) & inside] = i
    P = len(allc)
    print(f"remaining hosts: {n};  19-pt mosaic covers {(assign>=0).sum()} ({100*(assign>=0).sum()/n:.0f}%)")

    # ---- recovery curves ----
    Ts = np.arange(0, 61)
    uni_rec = np.array([((rneed <= T) & (assign >= 0)).sum() for T in Ts])
    uni_exp = P * Ts
    # per-pointing greedy-optimal depth
    maxd = int(rneed.max())
    cumc = []
    for p in range(P):
        nb = np.sort(rneed[assign == p])
        cumc.append(np.searchsorted(nb, np.arange(0, maxd + 2), side="right"))
    depth = np.zeros(P, int)
    def gain(p):
        t = depth[p]
        return cumc[p][t + 1] - cumc[p][t] if t + 1 <= maxd else 0
    heap = [(-gain(p), p) for p in range(P)]; heapq.heapify(heap)
    opt_exp, opt_rec, rec = [0], [0], 0
    while len(opt_exp) < 700:
        g, p = heapq.heappop(heap); g = -g
        if g != gain(p):
            heapq.heappush(heap, (-gain(p), p)); continue
        if g == 0:
            break
        depth[p] += 1; rec += g
        opt_exp.append(opt_exp[-1] + 1); opt_rec.append(rec)
        heapq.heappush(heap, (-gain(p), p))
    opt_exp, opt_rec = np.array(opt_exp), np.array(opt_rec)

    for T in (5, 8, 12, 20):
        i = T
        print(f"  uniform T={T:2d}h: {uni_exp[i]:4d} exp -> {uni_rec[i]:4d} ({100*uni_rec[i]/n:.0f}%)")
    print(f"  optimal knee ~{opt_rec.max()} hosts ({100*opt_rec.max()/n:.0f}%) by ~{opt_exp[np.argmax(opt_rec>=opt_rec.max()-1)]} exp")

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

    fig, (axm, axr) = plt.subplots(1, 2, figsize=(15, 6.5))

    # left: mosaic over remaining hosts
    def t2sky(x, y):
        return RA0 + x / cosd, DEC0 + y
    ra = RA0 + xi / cosd; dec = DEC0 + eta
    axm.scatter(ra, dec, s=5, c="0.55", alpha=0.5, linewidths=0)
    for grp, col, lab in ((core, "#1f6fe0", "A/B core (7)"),
                          (outer, "#d62728", "Config C extension (12)")):
        first = True
        for (xc, yc) in grp:
            sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
            axm.add_patch(Polygon(sky, closed=True, fill=False, edgecolor=col,
                                  lw=1.8, alpha=0.9, label=(lab if first else None)))
            first = False
    axm.set_aspect(1.0 / cosd); axm.invert_xaxis()
    axm.set_xlabel("RA (deg)", fontsize=16); axm.set_ylabel("Dec (deg)", fontsize=16)
    axm.set_title(f"19-Point Mosaic Over Remaining Hosts ($N={n}$)", fontsize=15)
    axm.tick_params(labelsize=12); axm.legend(fontsize=12, loc="upper right"); axm.grid(True, alpha=0.25)

    # right: recovery vs exposures
    axr.plot(uni_exp, 100 * uni_rec / n, "-o", color="#1f6fe0", ms=4,
             label="Uniform depth (all 19)")
    axr.plot(opt_exp, 100 * opt_rec / n, "-", color="#d62728", lw=2.2,
             label="Per-pointing optimized")
    for T, lab in ((8, "C: ~8 h"), (20, "deeper")):
        axr.axvline(P * T, color="0.6", ls=":", lw=1)
        axr.text(P * T, 3, f"{P*T} exp", rotation=90, va="bottom", ha="right",
                 fontsize=11, color="0.3")
    axr.set_xlabel("Total 1-Hour Exposures", fontsize=16)
    axr.set_ylabel("Remaining Hosts Recovered (%)", fontsize=16)
    axr.set_title("Recovery vs Exposure Time", fontsize=15)
    axr.set_xlim(0, 600); axr.set_ylim(0, 100)
    axr.tick_params(labelsize=12); axr.legend(fontsize=12, loc="lower right"); axr.grid(True, alpha=0.3)

    png = os.path.join(PNG_DIR, "19_config_c_optimize_ELAIS-N1.png")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
