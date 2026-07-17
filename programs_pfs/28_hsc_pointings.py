#!/usr/bin/env python3
"""
28_hsc_pointings.py -- HSC version of the PFS 7-pointing figure (Fig 12):
lay out 7 Subaru/HSC pointings (circular FoV, 1.77 deg^2 -> radius 0.75 deg)
over the ELAIS-N1 SN Ia field, one central + 6 on a ring, with overlap, to
maximize the SN coverage. Configurations A (ring at 0 deg phase) and B
(rotated 30 deg) are drawn as two separate figures, (a) and (b).

The ring radius is optimized empirically against the simulated SN Ia
distribution (a scan around the geometric 7-circle covering optimum
d = sqrt(3) R = 1.30 deg, which fully covers r < 2R = 1.5 deg).

    python 28_hsc_pointings.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "03_snia_radec_ELAIS-N1.csv")

FOV_AREA = 1.77                          # HSC FoV [deg^2]
R = np.sqrt(FOV_AREA / np.pi)            # -> circular FoV radius 0.7506 deg


def ring_centers(d, rot):
    """Center + 6 ring pointings at ring radius d, rigid rotation rot [deg]."""
    c = [(0.0, 0.0)]
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang + rot)
        c.append((d * np.cos(a), d * np.sin(a)))
    return c


def order_ccw(centers):
    """P0 center; ring numbered P1..P6 from the westmost (screen-right), CCW."""
    ring = centers[1:]
    sang = [np.degrees(np.arctan2(e, -x)) % 360 for (x, e) in ring]
    return [centers[0]] + [ring[k] for k in np.argsort(sang)]


def coverage(centers, pts):
    """covered mask + per-pointing counts (overlaps count in each pointing)."""
    covered = np.zeros(len(pts), bool)
    counts = []
    for (xc, yc) in centers:
        inside = np.hypot(pts[:, 0] - xc, pts[:, 1] - yc) <= R
        counts.append(int(inside.sum()))
        covered |= inside
    return covered, counts


def main():
    d0 = np.genfromtxt(CSV, delimiter=",", names=True)
    ra, dec = d0["RA_deg"], d0["DEC_deg"]
    ra0, dec0 = float(np.median(ra)), float(np.median(dec))
    cosd = np.cos(np.radians(dec0))
    pts = np.column_stack([(ra - ra0) * cosd, dec - dec0])
    n = len(pts)
    print(f"HSC FoV {FOV_AREA} deg^2 -> circular radius R = {R:.3f} deg")

    # ---- optimize the ring radius against the SN distribution (config A) ----
    best = (0, None)
    for d in np.round(np.arange(1.00, 1.61, 0.025), 3):
        cov, _ = coverage(ring_centers(d, 0.0), pts)
        if cov.sum() > best[0]:
            best = (int(cov.sum()), d)
    D = best[1]
    print(f"optimal ring radius (max SN coverage, config A): d = {D:.3f} deg "
          f"(geometric covering optimum sqrt(3)R = {np.sqrt(3)*R:.3f})")

    cenA = order_ccw(ring_centers(D, 0.0))
    cenB = order_ccw(ring_centers(D, 30.0))
    covA, cntA = coverage(cenA, pts)
    covB, cntB = coverage(cenB, pts)
    comb = covA | covB
    print(f"A: {covA.sum()}/{n} ({100*covA.sum()/n:.1f}%)  "
          f"B: {covB.sum()}/{n} ({100*covB.sum()/n:.1f}%)  "
          f"A∪B: {comb.sum()}/{n} ({100*comb.sum()/n:.1f}%)")

    # ---- plots ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

    def t2sky(x, y):
        return ra0 + x / cosd, dec0 + y

    for label, cen, covered, counts, ncov, fc, lc in (
            ("A", cenA, covA, cntA, int(covA.sum()), "#1f6fe0", "red"),
            ("B", cenB, covB, cntB, int(covB.sum()), "red", "blue")):
        fig, ax = plt.subplots(figsize=(9, 9))
        raL, decL = ra[~covered], dec[~covered]
        raC, decC = ra[covered], dec[covered]
        ax.scatter(raL, decL, s=4, c="0.7", alpha=0.5, linewidths=0, label="Not covered")
        ax.scatter(raC, decC, s=4, c="#2c6fbb", alpha=0.5, linewidths=0, label="Covered")
        th = np.linspace(0, 2 * np.pi, 181)
        for i, (xc, yc) in enumerate(cen):
            cra, cdec = t2sky(xc + R * np.cos(th), yc + R * np.sin(th))
            ax.fill(cra, cdec, facecolor=fc, alpha=0.10)
            ax.plot(cra, cdec, color=fc, lw=2.0, alpha=0.9)
            rc, dc = t2sky(xc, yc)
            ax.plot(rc, dc, "+", color=lc, ms=12, mew=2.2)
            ax.text(rc, dc, f"  P{i}", color=lc, fontsize=19, fontweight="bold",
                    ha="left", va="bottom")
            ax.text(rc, dc, f"  ({counts[i]})", color=lc, fontsize=16,
                    ha="left", va="top")
        ax.set_aspect(1.0 / cosd)
        ax.invert_xaxis()
        ax.set_xlabel("RA [deg]", fontsize=20)
        ax.set_ylabel("Dec [deg]", fontsize=20)
        ax.tick_params(labelsize=15)
        ax.set_title(f"HSC 7-pointing config {label} (ring $d={D:.2f}^\\circ$) "
                     f"over ELAIS-N1\nFoV {FOV_AREA} deg$^2$, {ncov}/{n} SNe covered "
                     f"({100*ncov/n:.0f}%);  A$\\cup$B = {100*comb.sum()/n:.0f}%",
                     fontsize=18)
        ax.legend(fontsize=17, loc="upper right", markerscale=2)
        ax.grid(True, alpha=0.25)
        ax.text(0.02, 0.02, "P$i$ (n) = pointing index (SNe in field)",
                transform=ax.transAxes, fontsize=12, color="0.3",
                bbox=dict(fc="white", ec="0.7", alpha=0.8))
        png = os.path.join(PNG_DIR, f"28_hsc_pointings_{label}_ELAIS-N1.png")
        fig.tight_layout()
        fig.savefig(png, dpi=140)
        print("plot ->", png)


if __name__ == "__main__":
    main()
