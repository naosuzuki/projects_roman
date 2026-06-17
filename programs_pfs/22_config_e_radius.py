#!/usr/bin/env python3
"""
22_config_e_radius.py -- empirically optimize the radius of the Config E ring
(12 pointings) for recovering the outskirt SN-host galaxies of ELAIS-N1.

Trade-off: a larger ring reaches the outermost hosts (higher coverage) but pushes
each field of view past the footprint edge (lower fill = more empty sky); a
smaller ring is fully inside the footprint (fill 100%) but misses the outer
hosts. Coverage and fill therefore pull in opposite directions; their product
(a balanced figure of merit) peaks near r = 1.6 deg.

Reads 07_program_sne_ELAIS-N1.csv. Produces coverage/fill vs radius figure.

    python 22_config_e_radius.py
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
R_FOOT = 2.09
M = 12


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    pAh = np.asarray(d["pA_host"]).astype(int)
    pBh = np.asarray(d["pB_host"]).astype(int)
    cosd = np.cos(np.radians(DEC0))
    xi = (hra - RA0) * cosd; eta = hdec - DEC0
    rem = (htar == 1) & (comp == 0)
    outer = rem & ~((pAh >= 0) | (pBh >= 0))
    pts = np.column_stack([xi, eta])[outer]
    nout = int(outer.sum())

    gx, gy = np.meshgrid(np.linspace(-R, R, 41), np.linspace(-R, R, 41))
    gx, gy = gx.ravel(), gy.ravel()
    inhex = Path(hexagon(0, 0)).contains_points(np.column_stack([gx, gy]))
    gx, gy = gx[inhex], gy[inhex]

    radii = np.round(np.arange(1.25, 2.001, 0.025), 3)
    cov, fill = [], []
    for rr in radii:
        cen = [(rr * np.cos(np.radians(p + 360 * k / M)), rr * np.sin(np.radians(p + 360 * k / M)))
               for p in (0.0, 15.0) for k in range(M)]      # Config E + Config F = 24 pointings
        ncov = np.zeros(len(pts), int); fr = []
        for xc, yc in cen:
            ncov += Path(hexagon(xc, yc)).contains_points(pts).astype(int)
            fr.append(np.mean(np.hypot(gx + xc, gy + yc) <= R_FOOT))
        cov.append(100 * (ncov > 0).sum() / nout)
        fill.append(100 * np.mean(fr))
    cov, fill = np.array(cov), np.array(fill)
    prod = cov * fill / 100.0
    r_nospill = radii[fill >= 99.5].max()                # largest ring with no empty space
    r_bal = radii[np.argmax(prod)]                       # balanced (coverage*fill) optimum
    print(f"outskirt hosts: {nout}  (Config E+F = 24 pointings)")
    print(f"  no empty space (fill=100%), max coverage: r = {r_nospill:.2f} deg  "
          f"(cov {cov[radii==r_nospill][0]:.0f}%, fill 100%)")
    print(f"  balanced optimum (max coverage*fill):     r = {r_bal:.2f} deg  "
          f"(cov {cov[radii==r_bal][0]:.0f}%, fill {fill[radii==r_bal][0]:.0f}%)")

    # ---- figure ----
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

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(radii, cov, "-o", color="#1f6fe0", ms=3.5, label="Coverage of outskirt hosts")
    ax.plot(radii, fill, "-s", color="#d62728", ms=3.5, label="Mean FoV fill (inside footprint)")
    ax.plot(radii, prod, "-", color="#2ca02c", lw=2.2, label="Coverage $\\times$ fill (figure of merit)")
    ax.axvline(r_nospill, color="0.4", ls="--", lw=1.3)
    ax.text(r_nospill, 4, f"  no empty space, $r={r_nospill:.2f}^\\circ$ (cov {cov[radii==r_nospill][0]:.0f}\\%)",
            rotation=90, va="bottom", ha="left", fontsize=12, color="0.3")
    ax.set_xlabel("Config E+F Ring Radius (deg)", fontsize=17)
    ax.set_ylabel("Percent", fontsize=17)
    ax.set_title("Optimizing the Config E+F Ring Radius (24 Pointings)", fontsize=15)
    ax.set_ylim(0, 105); ax.set_xlim(radii.min(), radii.max())
    ax.tick_params(labelsize=13); ax.legend(fontsize=12, loc="lower center"); ax.grid(True, alpha=0.3)

    png = os.path.join(PNG_DIR, "22_config_e_radius_ELAIS-N1.png")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
