#!/usr/bin/env python3
"""
32_config_beta.py -- Configuration beta = Configuration alpha rotated by 15 deg
(the half-period of the 12-pointing ring's 30-deg symmetry, hence the optimal
offset; a 30-deg rotation maps the ring onto itself and gains nothing).

alpha = E ring (12 pointings, r=1.5 deg) + G central square (4, d=0.6 deg).
beta  = the same 16 fields rigidly rotated by 15 deg.
Coverage of the 2270 remaining hosts: alpha 92%, alpha U beta 97%.
Alternating visits (alpha, beta, alpha, ...) recover 934 / 1140 / 1239 hosts
after 1 / 2 / 3 visits (vs 934 / 1072 / 1160 for alpha-only).

Reads 07_program_sne_ELAIS-N1.csv. Produces layout + recovery figure.

    python 32_config_beta.py
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
ROT_B = 15.0


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def config(rot):
    """Configuration alpha rigidly rotated by `rot` deg (16 fields)."""
    c = [(R_RING * np.cos(np.radians(rot + 360 * k / M)),
          R_RING * np.sin(np.radians(rot + 360 * k / M))) for k in range(M)]
    c += [(0.60 * np.cos(np.radians(rot + 45 + 90 * k)),
           0.60 * np.sin(np.radians(rot + 45 + 90 * k))) for k in range(4)]
    return c


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int) == 1
    comp = np.asarray(d["host_completed"]).astype(int) == 1
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    vneed = np.asarray(d["host_visits_needed"]).astype(float)
    vdone = np.asarray(d["host_visits_done"]).astype(float)
    cosd = np.cos(np.radians(DEC0))
    pts = np.column_stack([(hra - RA0) * cosd, hdec - DEC0])
    rneed = np.maximum(1, np.ceil(vneed - vdone))
    rem = htar & ~comp
    n = int(rem.sum())

    def ncov(cens):
        m = np.zeros(len(pts), int)
        for xc, yc in cens:
            m += Path(hexagon(xc, yc)).contains_points(pts).astype(int)
        return m
    A, B = config(0.0), config(ROT_B)
    nA, nB = ncov(A), ncov(B)
    covA = rem & (nA > 0)
    covU = rem & ((nA + nB) > 0)
    print(f"remaining hosts: {n}")
    print(f"  alpha: {covA.sum()} ({100*covA.sum()/n:.1f}%)   "
          f"alpha U beta: {covU.sum()} ({100*covU.sum()/n:.1f}%)")

    # alternating alpha/beta recovery
    Rmax = 12
    rec_ab, rec_a = [], []
    for v in range(0, Rmax + 1):
        ka = (v + 1) // 2; kb = v // 2
        rec_ab.append(int((rem & (ka * nA + kb * nB >= rneed)).sum()))
        rec_a.append(int((rem & (v * nA >= rneed)).sum()))
    rec_ab, rec_a = np.array(rec_ab), np.array(rec_a)
    exp = 16 * np.arange(0, Rmax + 1)
    for v in (1, 2, 3):
        print(f"  visits={v}: alt {rec_ab[v]} ({100*rec_ab[v]/n:.0f}%)  "
              f"alpha-only {rec_a[v]} ({100*rec_a[v]/n:.0f}%)")

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

    def draw_map(ax, lfs=16, tfs=12, tts=15, lgs=11):
        ax.scatter(ra[rem], dec[rem], s=5, c="0.55", alpha=0.5, linewidths=0,
                   label="Remaining hosts")
        for cens, col, lab in ((A, "#1f6fe0", "Configuration $\\alpha$ (16)"),
                               (B, "#d62728", "Configuration $\\beta$ ($+15^\\circ$)")):
            first = True
            for xc, yc in cens:
                sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
                ax.add_patch(Polygon(sky, closed=True, fill=False, edgecolor=col,
                                     lw=1.5, alpha=0.85, label=(lab if first else None)))
                first = False
        th = np.linspace(0, 2 * np.pi, 240)
        fra, fdec = t2sky(R_FOOT * np.cos(th), R_FOOT * np.sin(th))
        ax.plot(fra, fdec, color="0.4", ls="--", lw=1.1, label="Footprint edge")
        ax.set_aspect(1.0 / cosd); ax.invert_xaxis()
        ax.set_xlabel("RA (deg)", fontsize=lfs); ax.set_ylabel("Dec (deg)", fontsize=lfs)
        ax.set_title(f"Configurations $\\alpha$ + $\\beta$  "
                     f"(coverage {100*covA.sum()/n:.0f}% $\\to$ {100*covU.sum()/n:.0f}%)",
                     fontsize=tts)
        ax.tick_params(labelsize=tfs); ax.legend(fontsize=lgs, loc="upper right")
        ax.grid(True, alpha=0.25)

    fig, (axm, axr) = plt.subplots(1, 2, figsize=(15, 6.5))
    draw_map(axm)

    axr.plot(exp, 100 * rec_ab / n, "-o", color="#1f6fe0", ms=4,
             label="$\\alpha/\\beta$ alternating")
    axr.plot(exp, 100 * rec_a / n, "-s", color="0.5", ms=4, label="$\\alpha$ only")
    axr.axhline(100 * covU.sum() / n, color="0.6", ls=":", lw=1,
                label=f"$\\alpha\\cup\\beta$ coverage ceiling {100*covU.sum()/n:.0f}%")
    axr.set_xlabel("Total 1-Hour Exposures", fontsize=16)
    axr.set_ylabel("Remaining Hosts Recovered (%)", fontsize=16)
    axr.set_title("Recovery: $\\alpha/\\beta$ Alternating vs $\\alpha$ Only", fontsize=15)
    axr.set_xlim(0, 100); axr.set_ylim(0, 100)
    axr.tick_params(labelsize=12); axr.legend(fontsize=12, loc="lower right")
    axr.grid(True, alpha=0.3)

    png = os.path.join(PNG_DIR, "32_config_beta_ELAIS-N1.png")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)

    # standalone version of the alpha+beta FoV map (left panel)
    fig2, axs = plt.subplots(figsize=(9, 9))
    draw_map(axs, lfs=20, tfs=15, tts=17, lgs=14)
    png2 = os.path.join(PNG_DIR, "32_config_beta_map_ELAIS-N1.png")
    fig2.tight_layout(); fig2.savefig(png2, dpi=140)
    print("plot ->", png2)

    # standalone Configuration-beta-only map (companion to the alpha map of 24)
    covB = rem & (nB > 0)
    fig3, axb = plt.subplots(figsize=(9, 9))
    rb = covB[rem]                                  # covered mask within remaining
    axb.scatter(ra[rem][rb], dec[rem][rb], s=5, c="#2ca02c", alpha=0.5,
                linewidths=0, label="Covered")
    axb.scatter(ra[rem][~rb], dec[rem][~rb], s=6, c="#d62728", alpha=0.6,
                linewidths=0, label="Missed")
    for grp, col, lab in ((B[:M], "#1f6fe0", "E ring $+15^\\circ$ (12)"),
                          (B[M:], "#9467bd", "G central $+15^\\circ$ (4)")):
        first = True
        for xc, yc in grp:
            sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
            axb.add_patch(Polygon(sky, closed=True, fill=False, edgecolor=col,
                                  lw=1.7, alpha=0.9, label=(lab if first else None)))
            first = False
    th = np.linspace(0, 2 * np.pi, 240)
    fra, fdec = t2sky(R_FOOT * np.cos(th), R_FOOT * np.sin(th))
    axb.plot(fra, fdec, color="0.4", ls="--", lw=1.1, label="Footprint edge")
    axb.set_aspect(1.0 / cosd); axb.invert_xaxis()
    axb.set_xlabel("RA (deg)", fontsize=20); axb.set_ylabel("Dec (deg)", fontsize=20)
    axb.set_title(f"Configuration $\\beta$ = $\\alpha$ Rotated $15^\\circ$  "
                  f"(Coverage {100*covB.sum()/n:.0f}%)", fontsize=17)
    axb.tick_params(labelsize=15); axb.legend(fontsize=14, loc="upper right")
    axb.grid(True, alpha=0.25)
    png3 = os.path.join(PNG_DIR, "32_config_beta_only_map_ELAIS-N1.png")
    fig3.tight_layout(); fig3.savefig(png3, dpi=140)
    print("plot ->", png3)


if __name__ == "__main__":
    main()
