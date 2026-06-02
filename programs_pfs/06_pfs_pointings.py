#!/usr/bin/env python3
"""
06_pfs_pointings.py -- lay out 7 PFS pointings over ELAIS-N1 in a hexagonal
("flower") configuration: one central pointing surrounded by six, packed
edge-to-edge (honeycomb).

The PFS field of view is a hexagon ~1.3 deg in diameter (circumradius
R = 0.65 deg). Edge-sharing neighbors sit at spacing sqrt(3)*R = 1.126 deg.
We draw the flower over the simulated SN Ia distribution, mark which SNe are
covered, and report per-pointing counts (the seed of the fiber budget).

    python 06_pfs_pointings.py
    python 06_pfs_pointings.py --fov 1.3
"""
import os
import argparse
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
SNE_CSV = os.path.join(CSV_DIR, "03_snia_radec_ELAIS-N1.csv")


def hexagon(xc, yc, R, phase_deg=0.0):
    """Regular hexagon vertices (tangent-plane deg) about (xc,yc); phase rotates it."""
    a = np.radians(np.arange(0, 360, 60) + phase_deg)
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def flower_centers(R, rotate_deg=0.0):
    """7 edge-sharing pointing centers (center + 6 ring), rigidly rotated."""
    spacing = np.sqrt(3) * R
    centers = [(0.0, 0.0)]
    for ang in (90, 150, 210, 270, 330, 30):     # across the 6 flat edges
        a = np.radians(ang + rotate_deg)
        centers.append((spacing * np.cos(a), spacing * np.sin(a)))
    return centers


def coverage(centers, R, phase_deg, pts):
    """Return (covered_mask, per_pointing_counts) for SNe points `pts`."""
    covered = np.zeros(len(pts), bool)
    counts = []
    for (xc, yc) in centers:
        inside = Path(hexagon(xc, yc, R, phase_deg)).contains_points(pts)
        counts.append(int(inside.sum()))
        covered |= inside
    return covered, counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fov", type=float, default=1.3, help="PFS FoV diameter [deg]")
    ap.add_argument("--rotate", type=float, default=0.0,
                    help="rigid rotation of the flower [deg] (A=0, B=30)")
    ap.add_argument("--label", default=None, help="config label for title/filename")
    ap.add_argument("--sne-csv", default=SNE_CSV)
    ap.add_argument("--host-csv",
                    default=os.path.join(CSV_DIR, "03_snia_host_radec_ELAIS-N1.csv"))
    args = ap.parse_args()
    if args.label is None:
        args.label = {0.0: "A", 30.0: "B"}.get(args.rotate, f"{args.rotate:g}deg")
    os.makedirs(PNG_DIR, exist_ok=True)

    d = np.genfromtxt(args.sne_csv, delimiter=",", names=True)
    ra, dec = d["RA_deg"], d["DEC_deg"]
    ra0, dec0 = float(np.median(ra)), float(np.median(dec))
    cosd = np.cos(np.radians(dec0))
    xi = (ra - ra0) * cosd          # tangent-plane offsets [deg]
    eta = dec - dec0

    R = args.fov / 2.0              # hex circumradius
    spacing = np.sqrt(3) * R        # edge-sharing (flat-to-flat) neighbor distance

    pts = np.column_stack([xi, eta])
    centers = flower_centers(R, args.rotate)
    covered, counts = coverage(centers, R, args.rotate, pts)

    # A (0 deg) and B (30 deg) for the combined two-visit coverage
    covA, _ = coverage(flower_centers(R, 0.0), R, 0.0, pts)
    covB, _ = coverage(flower_centers(R, 30.0), R, 30.0, pts)
    comb = covA | covB
    n = len(xi)
    ncov = int(covered.sum())
    print(f"ELAIS-N1 center RA0={ra0:.3f} Dec0={dec0:.3f}")
    print(f"PFS FoV {args.fov} deg, spacing {spacing:.3f} deg (edge-sharing), "
          f"config {args.label} (rotate {args.rotate:g} deg)")
    print(f"per-pointing SN counts (P0=center): {counts}")
    print(f"this config ({args.label}): {ncov}/{n} SNe = {100*ncov/n:.1f}%")
    print(f"A alone: {int(covA.sum())} ({100*covA.sum()/n:.1f}%)  "
          f"B alone: {int(covB.sum())} ({100*covB.sum()/n:.1f}%)  "
          f"A∪B: {int(comb.sum())} ({100*comb.sum()/n:.1f}%)")

    # ---- coverage for several target samples (same 7-pointing geometry) ----
    def cov_fracs(rr, dd):
        x = (rr - ra0) * cosd
        y = dd - dec0
        p = np.column_stack([x, y])
        a = coverage(flower_centers(R, 0.0), R, 0.0, p)[0]
        b = coverage(flower_centers(R, 30.0), R, 30.0, p)[0]
        return len(x), int(a.sum()), int(b.sum()), int((a | b).sum())

    C = CSV_DIR
    SAMPLES = [
        ("SNe Ia (all)",       args.sne_csv),
        ("SNe Ia (Z<24)",      os.path.join(C, "03_snia_radec_ELAIS-N1_Zlt24.csv")),
        ("CC SNe (all)",       os.path.join(C, "03_cc_radec_ELAIS-N1.csv")),
        ("CC SNe (Z<24)",      os.path.join(C, "03_cc_radec_ELAIS-N1_Zlt24.csv")),
        ("TDE (all)",          os.path.join(C, "03_tde_radec_ELAIS-N1.csv")),
        ("TDE (Z<24)",         os.path.join(C, "03_tde_radec_ELAIS-N1_Zlt24.csv")),
        ("Ia host (all)",      args.host_csv),
        ("Ia host (Z<24.5)",   os.path.join(C, "03_snia_host_radec_ELAIS-N1_Zlt24.5.csv")),
        ("CC host (all)",      os.path.join(C, "03_cc_host_radec_ELAIS-N1.csv")),
        ("CC host (Z<24.5)",   os.path.join(C, "03_cc_host_radec_ELAIS-N1_Zlt24.5.csv")),
    ]
    cov_csv = os.path.join(C, "06_coverage_ELAIS-N1.csv")
    print("\n  sample                total     A            B          A∪B")
    with open(cov_csv, "w") as fo:
        fo.write("sample,total,A,B,AunionB,A_pct,B_pct,AunionB_pct\n")
        for nm, path in SAMPLES:
            if not os.path.exists(path):
                print(f"  {nm:18s} (missing: {os.path.basename(path)})")
                continue
            cat = np.genfromtxt(path, delimiter=",", names=True)
            if cat.size == 0:
                print(f"  {nm:18s} (empty)")
                continue
            rr = np.atleast_1d(cat["RA_deg"])
            dd = np.atleast_1d(cat["DEC_deg"])
            tot, na, nb, nab = cov_fracs(rr, dd)
            fo.write(f"{nm},{tot},{na},{nb},{nab},"
                     f"{100*na/tot:.1f},{100*nb/tot:.1f},{100*nab/tot:.1f}\n")
            print(f"  {nm:18s} {tot:6d}  {na:5d}({100*na/tot:4.1f}%) "
                  f"{nb:5d}({100*nb/tot:4.1f}%) {nab:5d}({100*nab/tot:4.1f}%)")
    print("coverage table ->", cov_csv)

    # tangent -> RA/Dec for plotting
    def t2sky(x, y):
        return ra0 + x / cosd, dec0 + y

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Circle
    from matplotlib import font_manager as fm

    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    LABEL_FS, TICK_FS, TITLE_FS = 16, 13, 16

    fig, ax = plt.subplots(figsize=(9, 9))

    # SNe: lost (gray) and covered (blue)
    raL, decL = ra[~covered], dec[~covered]
    raC, decC = ra[covered], dec[covered]
    ax.scatter(raL, decL, s=4, c="0.7", alpha=0.5, linewidths=0, label="not covered")
    ax.scatter(raC, decC, s=4, c="#2c6fbb", alpha=0.5, linewidths=0, label="covered")

    # 7 hexagonal pointings
    for i, (xc, yc) in enumerate(centers):
        sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc, R, args.rotate)])
        ax.add_patch(Polygon(sky, closed=True, fill=True, facecolor="#d62728",
                             edgecolor="#d62728", alpha=0.10, lw=2.0))
        rc, dc = t2sky(xc, yc)
        ax.plot(rc, dc, "+", color="#d62728", ms=10, mew=2)
        ax.text(rc, dc, f"  P{i}\n  ({counts[i]})", color="#7a1417",
                fontsize=11, ha="left", va="center")

    ax.set_aspect(1.0 / cosd)
    ax.invert_xaxis()
    ax.set_xlabel("RA [deg]", fontsize=LABEL_FS)
    ax.set_ylabel("Dec [deg]", fontsize=LABEL_FS)
    ax.tick_params(labelsize=TICK_FS)
    ax.set_title(f"PFS 7-pointing hexagonal config {args.label} "
                 f"(rotate {args.rotate:g}$^\\circ$) over ELAIS-N1\n"
                 f"FoV {args.fov}$^\\circ$, {ncov}/{len(xi)} SNe covered "
                 f"({100*ncov/len(xi):.0f}\\%);  A$\\cup$B = "
                 f"{100*comb.sum()/len(xi):.0f}\\%", fontsize=TITLE_FS)
    ax.legend(fontsize=11, loc="upper right", markerscale=2)
    ax.grid(True, alpha=0.25)
    ax.text(0.02, 0.02, "P$i$ (n) = pointing index (SNe in field)",
            transform=ax.transAxes, fontsize=10, color="0.3",
            bbox=dict(fc="white", ec="0.7", alpha=0.8))

    png = os.path.join(PNG_DIR, f"06_pfs_pointings_hex7{args.label}_ELAIS-N1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
