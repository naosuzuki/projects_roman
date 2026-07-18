#!/usr/bin/env python3
"""
06_desi_pointings.py -- DESI counterpart of programs_pfs/06_pfs_pointings.py.

Lay out DESI pointings over ELAIS-N1 with the circular 3.2 deg-diameter DESI
field of view (radius R = 1.6 deg), in two configurations analogous to the PFS
A/B flowers:

  * Config A -- a 3-pointing equilateral TRIANGLE.
  * Config B -- a 6-pointing hexagonal RING.

The target field is the live-SN + host program in ELAIS-N1 (Phase II live
supernovae AND Phase III host galaxies): all program SNe positions plus their
valid host-galaxy positions. The configurations are deliberately COMPACT --
small ring radii that pile the circular fields toward the center, fully covering
the high-S/N core (r < 1.5 deg, ~60% of the targets) while letting the sparse,
low-S/N outer rim fall outside. Unlike the PFS flower, which gains area through
field OVERLAP, the DESI case gains through the wide field: one 3.2 deg pointing's
5000 fibers observe the whole local population at once.

    python 06_desi_pointings.py                      # compact default radii
    python 06_desi_pointings.py --d3 0.7 --d6 0.6    # override ring radii

Field center RA0 = 242.498, Dec0 = +54.497 (as in programs_pfs/25_case_summary).
"""
import os
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
HOST_CSV = os.path.join(HERE, "..", "programs_pfs", "outputs", "csv",
                        "07_program_sne_ELAIS-N1.csv")

RA0, DEC0 = 242.498, 54.497      # ELAIS-N1 field center (deg)
FOV = 3.2                        # DESI field of view diameter (deg)
R = FOV / 2.0                    # DESI field radius (deg)
RCORE = 1.5                      # high-S/N core radius (deg)
RFIELD = 2.09                    # ELAIS-N1 field boundary (r99 of live-SN+host field)

# Compact, rim-INSIDE ring radii: every DESI field of view stays entirely within
# the ELAIS-N1 boundary, requiring d + R <= RFIELD, i.e. d <= 0.49 deg. We go
# smaller still to concentrate the fields on the high-S/N core (r<1.5 deg, which
# stays 100% covered), accepting loss of the sparse low-S/N rim.
DEFAULT_D3 = 0.40                # 3-pointing triangle  (rim at 2.00 deg)
DEFAULT_D6 = 0.40                # 6-pointing hexagon   (rim at 2.00 deg)


def triangle_centers(d, phase_deg=90.0):
    """3 pointing centers on a ring of radius d, 120 deg apart (one vertex up)."""
    a = np.radians(phase_deg + 120.0 * np.arange(3))
    return np.column_stack([d * np.cos(a), d * np.sin(a)])


def hexagon_centers(d, phase_deg=0.0):
    """6 pointing centers on a ring of radius d, 60 deg apart."""
    a = np.radians(phase_deg + 60.0 * np.arange(6))
    return np.column_stack([d * np.cos(a), d * np.sin(a)])


def coverage(centers, pts):
    """(covered_mask, per_pointing_counts) for tangent-plane points `pts`."""
    covered = np.zeros(len(pts), bool)
    counts = []
    for xc, yc in centers:
        inside = np.hypot(pts[:, 0] - xc, pts[:, 1] - yc) <= R
        counts.append(int(inside.sum()))
        covered |= inside
    return covered, counts


def optimize_radius(make_centers, pts, lo=0.0, hi=1.8, step=0.025):
    """Ring radius d (in [lo,hi]) maximizing covered fraction of `pts`."""
    best_d, best_c = lo, -1.0
    for d in np.arange(lo, hi + 1e-9, step):
        c = coverage(make_centers(d), pts)[0].mean()
        if c > best_c + 1e-9:
            best_d, best_c = d, c
    return float(best_d), float(best_c)


def load_samples():
    """Tangent-plane points for the live-SN + host samples; field center fixed.

    Returns (samples, field, cosd) where `field` is the combined live-SN + host
    target set (the population the DESI tiling must cover) and `samples` are the
    individual sub-samples for the coverage table. Sentinel/invalid host
    coordinates are filtered out.
    """
    d = np.genfromtxt(HOST_CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    cosd = np.cos(np.radians(DEC0))
    htar = d["host_target"].astype(int)
    comp = d["host_completed"].astype(int)
    hra, hdec = d["host_ra"].astype(float), d["host_dec"].astype(float)
    ra, dec = d["ra"].astype(float), d["dec"].astype(float)

    def tp(rr, dd):
        return np.column_stack([(rr - RA0) * cosd, dd - DEC0])

    # valid host position: characterized host, finite, within 5 deg of center
    okh = ((htar == 1) & np.isfinite(hra) & np.isfinite(hdec)
           & (np.abs(hra - RA0) < 5) & (np.abs(hdec - DEC0) < 5))
    rem = okh & (comp == 0)

    sne = tp(ra, dec)                       # live SNe (Phase II): all program SNe
    hosts = tp(hra[okh], hdec[okh])         # host galaxies (Phase III)
    remaining = tp(hra[rem], hdec[rem])     # remaining (incomplete) hosts
    field = np.vstack([sne, hosts])         # the live-SN + host program field

    samples = {
        "Live SNe": sne,
        "Host galaxies": hosts,
        "Remaining hosts": remaining,
    }
    return samples, field, cosd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--d3", type=float, default=DEFAULT_D3,
                    help=f"triangle ring radius [deg] (default {DEFAULT_D3})")
    ap.add_argument("--d6", type=float, default=DEFAULT_D6,
                    help=f"hexagon ring radius [deg] (default {DEFAULT_D6})")
    ap.add_argument("--tri-phase", type=float, default=90.0)
    ap.add_argument("--hex-phase", type=float, default=0.0)
    args = ap.parse_args()

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    samples, field, cosd = load_samples()
    d3, d6 = args.d3, args.d6
    cen3 = triangle_centers(d3, args.tri_phase)
    cen6 = hexagon_centers(d6, args.hex_phase)

    dmax = RFIELD - R     # largest ring radius with the FoV rim inside ELAIS-N1
    print(f"DESI FoV {FOV} deg (R={R} deg) over ELAIS-N1 (RA0={RA0}, Dec0={DEC0})")
    print(f"ELAIS-N1 boundary RFIELD={RFIELD} deg; rim-inside requires d <= {dmax:.2f} deg")
    print(f"Config A triangle: 3 pointings, d = {d3:.3f} deg, rim = {d3+R:.2f} deg "
          f"({'INSIDE' if d3+R <= RFIELD else 'OUTSIDE'} field)")
    print(f"Config B hexagon : 6 pointings, d = {d6:.3f} deg, rim = {d6+R:.2f} deg "
          f"({'INSIDE' if d6+R <= RFIELD else 'OUTSIDE'} field)")
    print(f"live-SN + host field: {len(field)} targets; "
          f"core (r<{RCORE} deg): {int((np.hypot(field[:,0],field[:,1])<RCORE).sum())}")

    # ---- coverage table over the samples (whole-field and high-S/N core) ----
    cov_csv = os.path.join(CSV_DIR, "06_desi_coverage_ELAIS-N1.csv")
    print("\n  sample                 total   A(3-tri)        B(6-hex)")
    with open(cov_csv, "w") as fo:
        fo.write("config,ring_radius_deg,sample,total,covered,covered_pct\n")
        rows = list(samples.items()) + [
            ("Field (SNe+hosts)", field),
            (f"  core r<{RCORE:g}", field[np.hypot(field[:, 0], field[:, 1]) < RCORE]),
        ]
        for nm, pts in rows:
            tot = len(pts)
            na = int(coverage(cen3, pts)[0].sum())
            nb = int(coverage(cen6, pts)[0].sum())
            fo.write(f"A_triangle,{d3:.3f},{nm.strip()},{tot},{na},{100*na/tot:.1f}\n")
            fo.write(f"B_hexagon,{d6:.3f},{nm.strip()},{tot},{nb},{100*nb/tot:.1f}\n")
            print(f"  {nm:20s} {tot:6d}  {na:5d}({100*na/tot:5.1f}%)  {nb:5d}({100*nb/tot:5.1f}%)")
    print("coverage table ->", cov_csv)

    # per-pointing target pool (cumulative live SNe + hosts falling in each FoV;
    # instantaneous fiber demand is far lower -- that's the fiber-budget sim).
    _, cnt3 = coverage(cen3, field)
    _, cnt6 = coverage(cen6, field)
    print(f"\n  per-pointing live-SN+host target pool (cumulative, heavy overlap):")
    print(f"    Config A (3): {cnt3}   max {max(cnt3)}")
    print(f"    Config B (6): {cnt6}   max {max(cnt6)}")

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib import font_manager as fm

    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    LABEL_FS, TICK_FS, TITLE_FS = 20, 15, 18

    # combined live-SN + host field sky coords for plotting
    fra = RA0 + field[:, 0] / cosd
    fdec = DEC0 + field[:, 1]

    def t2sky(x, y):
        return RA0 + x / cosd, DEC0 + y

    def circle_sky(xc, yc, rad, n=240):
        a = np.linspace(0, 2 * np.pi, n)
        return np.array([t2sky(xc + rad * np.cos(t), yc + rad * np.sin(t)) for t in a])

    def draw_config(ax, centers, counts, color, label):
        # field reference boundary and high-S/N core
        fb = circle_sky(0, 0, RFIELD)
        ax.plot(fb[:, 0], fb[:, 1], color="0.55", lw=1.2, ls="--",
                label="ELAIS-N1 ($r_{99}=2.09^\\circ$)")
        cc = circle_sky(0, 0, RCORE)
        ax.plot(cc[:, 0], cc[:, 1], color="green", lw=1.3, ls=":",
                label=f"High-S/N core ($r<{RCORE:g}^\\circ$)")
        covered, _ = coverage(centers, field)
        ax.scatter(fra[~covered], fdec[~covered], s=4, c="0.7", alpha=0.5,
                   linewidths=0, label="Not covered (rim)")
        ax.scatter(fra[covered], fdec[covered], s=4, c="#2c6fbb", alpha=0.45,
                   linewidths=0, label="Covered (SNe + hosts)")
        for i, (xc, yc) in enumerate(centers):
            sky = circle_sky(xc, yc, R)
            ax.add_patch(Polygon(sky, closed=True, facecolor=color,
                                 edgecolor="none", alpha=0.10))
            ax.add_patch(Polygon(sky, closed=True, fill=False,
                                 edgecolor=color, lw=2.0, alpha=0.9))
            rc, dc = t2sky(xc, yc)
            ax.plot(rc, dc, "+", color="red", ms=13, mew=2.3)
            ax.text(rc, dc, f"  P{i+1}", color="red", fontsize=18,
                    fontweight="bold", ha="left", va="bottom")
            ax.text(rc, dc, f"  ({counts[i]})", color="red", fontsize=14,
                    ha="left", va="top")
        ax.set_aspect(1.0 / cosd)
        ax.invert_xaxis()
        ax.set_xlabel("RA [deg]", fontsize=LABEL_FS)
        ax.set_ylabel("Dec [deg]", fontsize=LABEL_FS)
        ax.tick_params(labelsize=TICK_FS)
        ax.set_title(label, fontsize=TITLE_FS)
        ax.legend(fontsize=11, loc="upper right", markerscale=2.5, framealpha=0.9)
        ax.grid(True, alpha=0.25)

    n = len(field)
    ncovA = int(coverage(cen3, field)[0].sum())
    ncovB = int(coverage(cen6, field)[0].sum())

    # --- Config A: 3-pointing triangle ---
    figA, axA = plt.subplots(figsize=(8.6, 8.6))
    draw_config(axA, cen3, cnt3, "#1f6fe0",
                f"DESI Config A: 3-pointing triangle (FoV $3.2^\\circ$, $d={d3:.2f}^\\circ$)\n"
                f"live SNe + hosts: {100*ncovA/n:.0f}\\% covered, core $r<1.5^\\circ$ 100\\%")
    pngA = os.path.join(PNG_DIR, "06_desi_pointings_tri3_ELAIS-N1.png")
    figA.tight_layout(); figA.savefig(pngA, dpi=140)
    print("plot ->", pngA)

    # --- Config B: 6-pointing hexagon ---
    figB, axB = plt.subplots(figsize=(8.6, 8.6))
    draw_config(axB, cen6, cnt6, "#d1410c",
                f"DESI Config B: 6-pointing hexagon (FoV $3.2^\\circ$, $d={d6:.2f}^\\circ$)\n"
                f"live SNe + hosts: {100*ncovB/n:.0f}\\% covered, core $r<1.5^\\circ$ 100\\%")
    pngB = os.path.join(PNG_DIR, "06_desi_pointings_hex6_ELAIS-N1.png")
    figB.tight_layout(); figB.savefig(pngB, dpi=140)
    print("plot ->", pngB)


if __name__ == "__main__":
    main()
