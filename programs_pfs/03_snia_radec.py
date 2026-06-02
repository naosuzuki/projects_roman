#!/usr/bin/env python3
"""
03_snia_radec.py -- sky positions (RA, Dec) of all Type Ia SNe in the
Roman HLTDS HOURGLASS2 simulation, colored by redshift.

The Ia sample (SNANA model SNIaMODEL00) spans two widely-separated deep
fields, so we draw one panel per field (ELAIS-N1 north, EDF-S south) with
the proper RA/cos(Dec) aspect and astronomical (RA increasing left) axes.

Reads the split field directories under DATADIR; handles both .FITS and
.FITS.gz (astropy detects gzip), skipping macOS ._ AppleDouble files.

    python 03_snia_radec.py
    python 03_snia_radec.py --datadir /path/to/ltcv_sim
"""
import os
import glob
import argparse
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")

DEFAULT_DATADIR = "/Volumes/exdisk1/data/Roman/ltcv_sim"
FIELDS = ["ELAIS-N1", "EDF-S"]          # north, south

# What to plot: the SN itself or its host galaxy. Each maps to HEAD columns
# (position, redshift for color, and the Roman Z-band mag used for the cut).
TARGETS = {
    "sn":   dict(ra="RA", dec="DEC", z="SIM_REDSHIFT_CMB", mag="SIM_PEAKMAG_Z",
                 label="Type Ia SNe", magname="peak Z"),
    "host": dict(ra="HOSTGAL_RA", dec="HOSTGAL_DEC", z="SIM_REDSHIFT_CMB",
                 mag="HOSTGAL_MAG_Z", label="Type Ia host galaxies", magname="host Z"),
}


def field_head_files(datadir, field, model="SNIaMODEL00"):
    """All HEAD files for `model` in a field dir (.FITS or .FITS.gz, no ._)."""
    files = []
    for ext in ("FITS", "FITS.gz"):
        files += glob.glob(os.path.join(datadir, field, f"*{model}-*_HEAD.{ext}"))
    files = [f for f in files if not os.path.basename(f).startswith("._")]
    return sorted(files)


def load_field(datadir, field, cols, model="SNIaMODEL00"):
    """Concatenate SNID/RA/Dec/z/zmag over all HEAD files, using the column
    names in `cols` (the SN or its host galaxy)."""
    snid, ra, dec, z, zmag = [], [], [], [], []
    for f in field_head_files(datadir, field, model):
        with fits.open(f) as h:
            d = h[1].data
            if d is None or len(d) == 0:
                continue
            snid.append(np.asarray(d["SNID"]).astype(str))
            ra.append(np.asarray(d[cols["ra"]], float))
            dec.append(np.asarray(d[cols["dec"]], float))
            z.append(np.asarray(d[cols["z"]], float))
            zmag.append(np.asarray(d[cols["mag"]], float))
    if not ra:
        return (np.array([]),) * 5
    return (np.concatenate(snid), np.concatenate(ra), np.concatenate(dec),
            np.concatenate(z), np.concatenate(zmag))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datadir", default=DEFAULT_DATADIR,
                    help=f"ltcv_sim dir with field subdirs (default {DEFAULT_DATADIR})")
    ap.add_argument("--cmap", default="rainbow",
                    help="redshift colormap (default rainbow: low-z blue -> high-z red)")
    ap.add_argument("--fields", nargs="+", default=FIELDS, choices=FIELDS,
                    help="which field(s) to plot (default: both)")
    ap.add_argument("--zmag-cut", type=float, default=None, dest="zmag_cut",
                    help="keep only SNe with peak Roman Z-band AB mag brighter "
                         "than this (e.g. 24.5)")
    ap.add_argument("--msize", type=float, default=4.0,
                    help="scatter marker size (default 4)")
    ap.add_argument("--alpha", type=float, default=0.7,
                    help="scatter marker opacity (default 0.7)")
    ap.add_argument("--target", choices=list(TARGETS), default="sn",
                    help="plot the SN ('sn') or its host galaxy ('host')")
    ap.add_argument("--zmin", type=float, default=0.0, help="colorbar min redshift")
    ap.add_argument("--zmax", type=float, default=None,
                    help="colorbar max redshift (fix it to compare samples; "
                         "default: auto from data)")
    ap.add_argument("--color-by", choices=["redshift", "zmag"], default="redshift",
                    dest="color_by", help="quantity mapped to color (default redshift)")
    ap.add_argument("--cmin", type=float, default=None, help="color-axis min override")
    ap.add_argument("--cmax", type=float, default=None, help="color-axis max override")
    ap.add_argument("--model", default="SNIaMODEL00",
                    help="SNANA model: SNIaMODEL00 (Ia), NONIaMODEL06 (CC), "
                         "NONIaMODEL02 (TDE), ...")
    ap.add_argument("--tag", default="snia",
                    help="short type label for output filenames (e.g. snia, cc, tde)")
    args = ap.parse_args()
    fields = args.fields
    cols = TARGETS[args.target]

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    data = {}
    for fld in fields:
        snid, ra, dec, z, zmag = load_field(args.datadir, fld, cols, args.model)
        if len(ra):   # drop sentinel positions (e.g. unmatched hosts at -999)
            ok = np.isfinite(ra) & np.isfinite(dec) & (ra >= 0) & (ra <= 360) & (np.abs(dec) <= 90)
            snid, ra, dec, z, zmag = (a[ok] for a in (snid, ra, dec, z, zmag))
        if args.zmag_cut is not None and len(ra):
            keep = np.isfinite(zmag) & (zmag > 0) & (zmag < 90) & (zmag < args.zmag_cut)
            snid, ra, dec, z, zmag = (a[keep] for a in (snid, ra, dec, z, zmag))
        data[fld] = (snid, ra, dec, z, zmag)
        if len(ra):
            print(f"{fld:10s}: N={len(ra):6d}  "
                  f"RA [{ra.min():.2f},{ra.max():.2f}]  "
                  f"Dec [{dec.min():.2f},{dec.max():.2f}]  "
                  f"z [{z.min():.2f},{z.max():.2f}]")
    ntot = sum(len(d[1]) for d in data.values())
    base = f"03_{args.tag}_host_radec" if args.target == "host" else f"03_{args.tag}_radec"
    fsuffix = "" if fields == FIELDS else "_" + "_".join(fields)
    if args.zmag_cut is not None:
        fsuffix += f"_Zlt{args.zmag_cut:g}"
    if args.color_by == "zmag":
        fsuffix += "_byZmag"

    # ---- CSV catalog ----
    csv = os.path.join(CSV_DIR, f"{base}{fsuffix}.csv")
    with open(csv, "w") as fo:
        fo.write("field,SNID,RA_deg,DEC_deg,z_cmb,peakmag_Z\n")
        for fld in fields:
            snid, ra, dec, z, zmag = data[fld]
            for s, r, dd, zz, zm in zip(snid, ra, dec, z, zmag):
                fo.write(f"{fld},{s},{r:.5f},{dd:.5f},{zz:.4f},{zm:.3f}\n")
    print("catalog ->", csv)

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
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

    # what gets mapped to color: redshift (index 3) or Z-band mag (index 4)
    if args.color_by == "zmag":
        cidx, clabel = 4, f"{cols['magname']}  (AB mag)"
        cmin = args.cmin if args.cmin is not None else 20.0
        cmax = args.cmax if args.cmax is not None else 24.5
    else:
        cidx, clabel = 3, "redshift  $z_{\\rm CMB}$"
        zall = np.concatenate([data[f][3] for f in fields if len(data[f][3])])
        cmin = args.cmin if args.cmin is not None else args.zmin
        cmax = (args.cmax if args.cmax is not None else
                (args.zmax if args.zmax is not None else float(np.nanmax(zall))))

    n = len(fields)
    figw = 8.5 if n == 1 else 6.5 * n
    fig, axes = plt.subplots(1, n, figsize=(figw, 7.5), constrained_layout=True,
                             squeeze=False)
    axes = axes.ravel()
    sc = None
    for ax, fld in zip(axes, fields):
        snid, ra, dec, z, zmag = data[fld]
        cval = zmag if args.color_by == "zmag" else z
        sc = ax.scatter(ra, dec, c=cval, cmap=args.cmap, vmin=cmin, vmax=cmax,
                        s=args.msize, alpha=args.alpha, linewidths=0)
        dec0 = np.median(dec)
        ax.set_aspect(1.0 / np.cos(np.deg2rad(dec0)))   # true angular scale
        ax.invert_xaxis()                               # RA increases left
        ax.set_xlabel("RA  [deg]", fontsize=LABEL_FS)
        ax.set_ylabel("Dec  [deg]", fontsize=LABEL_FS)
        ax.tick_params(labelsize=TICK_FS)
        ax.set_title(f"{fld}   (N = {len(ra):,})\n"
                     f"RA≈{np.median(ra):.1f}°, Dec≈{dec0:+.1f}°", fontsize=TITLE_FS)
        ax.grid(True, alpha=0.25)

    cbar = fig.colorbar(sc, ax=axes, fraction=0.046, pad=0.04)
    cbar.set_label(clabel, fontsize=LABEL_FS)
    cbar.ax.tick_params(labelsize=TICK_FS)
    cutnote = f", {cols['magname']} < {args.zmag_cut:g}" if args.zmag_cut is not None else ""
    fig.suptitle(f"Roman HLTDS {cols['label']} — sky positions  "
                 f"(HOURGLASS2 sim{cutnote}, N = {ntot:,})", fontsize=13)

    png = os.path.join(PNG_DIR, f"{base}{fsuffix}.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
