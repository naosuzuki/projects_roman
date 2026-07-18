#!/usr/bin/env python3
"""
31_lofar_radec.py -- sky positions (RA, Dec) of LOFAR LoTSS DR2 radio sources
in the ELAIS-N1 field, colored by radio source strength (144 MHz total flux).

Companion to 03_snia_radec.py (same style: Times New Roman, true angular
aspect, RA increasing left, rainbow colormap). Reads the full-sky LoTSS DR2
source catalog (radio.astron.nl / lofar-surveys.org, ~4.4M sources) and keeps
sources within --radius of the field center. Flux spans ~4 dex, so color is
logarithmic (mJy).

    python 31_lofar_radec.py
    python 31_lofar_radec.py --catalog ~/Downloads/LoTSS_DR2_v110_masked.srl.fits \
        --center 242.5 54.5 --radius 2.1 --fmin 0.3 --fmax 100
"""
import os
import argparse
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")

DEFAULT_CATALOG = os.path.expanduser("~/Downloads/LoTSS_DR2_v110_masked.srl.fits")
# ELAIS-N1 center/radius matched to the Roman HLTDS field of 03_snia_radec.py
DEFAULT_CENTER = (242.5, 54.5)
DEFAULT_RADIUS = 2.1


def load_field(catalog, ra0, dec0, radius):
    """Sources within `radius` deg of (ra0, dec0): name, ra, dec, total & peak
    flux (mJy). A flat-sky cut prefilters the 4.4M-row table before the exact
    angular-distance selection."""
    with fits.open(catalog, memmap=True) as h:
        d = h[1].data
        ra, dec = np.asarray(d["RA"], float), np.asarray(d["DEC"], float)
        cosd = np.cos(np.deg2rad(dec0))
        near = (np.abs(dec - dec0) < radius) & (np.abs(ra - ra0) * cosd < 1.2 * radius)
        d = d[near]
        ra, dec = ra[near], dec[near]
        dra = np.deg2rad(ra - ra0)
        sep = np.rad2deg(np.arccos(np.clip(
            np.sin(np.deg2rad(dec0)) * np.sin(np.deg2rad(dec)) +
            np.cos(np.deg2rad(dec0)) * np.cos(np.deg2rad(dec)) * np.cos(dra), -1, 1)))
        keep = sep < radius
        return (np.asarray(d["Source_Name"])[keep].astype(str),
                ra[keep], dec[keep],
                np.asarray(d["Total_flux"], float)[keep],
                np.asarray(d["Peak_flux"], float)[keep])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", default=DEFAULT_CATALOG,
                    help=f"LoTSS DR2 source list (default {DEFAULT_CATALOG})")
    ap.add_argument("--center", nargs=2, type=float, default=DEFAULT_CENTER,
                    metavar=("RA", "DEC"), help="field center in deg")
    ap.add_argument("--radius", type=float, default=DEFAULT_RADIUS,
                    help=f"selection radius in deg (default {DEFAULT_RADIUS})")
    ap.add_argument("--cmap", default="rainbow",
                    help="flux colormap (default rainbow: faint blue -> bright red)")
    ap.add_argument("--flux", choices=["total", "peak"], default="total",
                    help="color by Total_flux or Peak_flux (default total)")
    ap.add_argument("--fmin", type=float, default=0.3,
                    help="colorbar min flux in mJy (default 0.3)")
    ap.add_argument("--fmax", type=float, default=100.0,
                    help="colorbar max flux in mJy (default 100)")
    ap.add_argument("--msize", type=float, default=4.0,
                    help="scatter marker size (default 4)")
    ap.add_argument("--alpha", type=float, default=0.7,
                    help="scatter marker opacity (default 0.7)")
    args = ap.parse_args()
    ra0, dec0 = args.center

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    name, ra, dec, ftot, fpeak = load_field(args.catalog, ra0, dec0, args.radius)
    flux = ftot if args.flux == "total" else fpeak
    print(f"ELAIS-N1  : N={len(ra):6d}  "
          f"RA [{ra.min():.2f},{ra.max():.2f}]  Dec [{dec.min():.2f},{dec.max():.2f}]  "
          f"flux [{flux.min():.2f},{flux.max():.1f}] mJy")

    base = "31_lofar_radec_ELAIS-N1"
    fsuffix = "_byPeak" if args.flux == "peak" else ""

    # ---- CSV catalog ----
    csv = os.path.join(CSV_DIR, f"{base}{fsuffix}.csv")
    with open(csv, "w") as fo:
        fo.write("Source_Name,RA_deg,DEC_deg,Total_flux_mJy,Peak_flux_mJy\n")
        for s, r, dd, ft, fp in zip(name, ra, dec, ftot, fpeak):
            fo.write(f"{s},{r:.5f},{dd:.5f},{ft:.4f},{fp:.4f}\n")
    print("catalog ->", csv)

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib import font_manager as fm

    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf",
               "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    LABEL_FS, TICK_FS, TITLE_FS = 16, 12, 15

    # brightest sources drawn last so they sit on top
    order = np.argsort(flux)
    ra, dec, cval = ra[order], dec[order], np.clip(flux[order], args.fmin, args.fmax)

    fig, ax = plt.subplots(figsize=(8.5, 7.5), constrained_layout=True)
    sc = ax.scatter(ra, dec, c=cval, cmap=args.cmap,
                    norm=LogNorm(vmin=args.fmin, vmax=args.fmax),
                    s=args.msize, alpha=args.alpha, linewidths=0)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(dec0)))   # true angular scale
    ax.invert_xaxis()                               # RA increases left
    ax.set_xlabel("RA  [deg]", fontsize=LABEL_FS)
    ax.set_ylabel("Dec  [deg]", fontsize=LABEL_FS)
    ax.tick_params(labelsize=TICK_FS)
    ax.set_title(f"ELAIS-N1   (N = {len(ra):,})\n"
                 f"RA≈{ra0:.1f}°, Dec≈{dec0:+.1f}°", fontsize=TITLE_FS)
    ax.grid(True, alpha=0.25)

    fluxname = "Total" if args.flux == "total" else "Peak"
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f"144 MHz {fluxname} Flux  [mJy]", fontsize=LABEL_FS)
    cbar.ax.tick_params(labelsize=TICK_FS)
    fig.suptitle(f"LOFAR LoTSS DR2 Radio Sources — Sky Positions  "
                 f"(r < {args.radius:g}°, N = {len(ra):,})", fontsize=13)

    png = os.path.join(PNG_DIR, f"{base}{fsuffix}.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
