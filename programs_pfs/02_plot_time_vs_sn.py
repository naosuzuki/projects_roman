#!/usr/bin/env python3
"""
plot_time_vs_sn.py -- PFS continuum S/N vs total exposure time, per band.

Runs the real PFS ETC over a grid of (magnitude, exposure time) and plots the
median continuum S/N per pixel in a chosen broadband filter (default: i).

Curves are color-coded as a red->blue rainbow from BRIGHT (red) to FAINT (blue).

Examples
--------
    python plot_time_vs_sn.py                       # i-band, linear, default grid
    python plot_time_vs_sn.py --loglog              # log-log version
    python plot_time_vs_sn.py --loglog --from-csv   # restyle from cached grid (no ETC re-run)
    python plot_time_vs_sn.py --band z --mags 23 24 25
    python plot_time_vs_sn.py --seeing 0.6 --moon 0.25

Outputs land in outputs/png/ and outputs/csv/ next to this script.
Conditions not overridden use the ETC defaults
(seeing 0.8", zenith 45 deg, new moon, field edge 0.675 deg, REFF 0.3").
"""
import os
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# The ETC wrapper lives in 01_pfs_etc.py; a leading digit isn't a valid module
# name for `import`, so load it by file path.
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("pfs_etc", os.path.join(HERE, "01_pfs_etc.py"))
pfs_etc = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(pfs_etc)
run_etc, band_sn, band_resolution, BANDS = (
    pfs_etc.run_etc, pfs_etc.band_sn, pfs_etc.band_resolution, pfs_etc.BANDS)

PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")


def compute_grid(args, overrides):
    """Run the ETC over the (mag, exp_num) grid -> {mag: (times_min, sn)}."""
    results = {}
    for mag in args.mags:
        times, sns = [], []
        for n in args.exp_nums:
            res = run_etc(mag=mag, exp_time=args.single_exp, exp_num=n, **overrides)
            sn = band_sn(res, args.band)
            times.append(args.single_exp * n / 60.0)
            sns.append(sn)
            print(f"  {args.band}={mag:5.1f}  N={n:3d}  "
                  f"t={args.single_exp*n/60:6.1f} min  S/N/pix={sn:7.2f}")
        results[mag] = (np.array(times), np.array(sns))
    return results


def write_csv(path, results, mags, exp_nums, band, cond):
    with open(path, "w") as f:
        f.write(f"# PFS ETC continuum S/N per pixel in {band}-band; {cond}\n")
        f.write("mag_AB,exp_num,total_time_min,SN_per_pixel\n")
        for mag in mags:
            t, s = results[mag]
            for n, tt, ss in zip(exp_nums, t, s):
                f.write(f"{mag},{n},{tt:.1f},{ss:.4f}\n")


def load_csv(path):
    """Read a cached grid CSV -> ({mag: (times_min, sn)}, condition_string).

    Columns: mag_AB, exp_num, total_time_min, SN_per_pixel.
    """
    with open(path) as f:
        lines = f.readlines()
    cond = ""
    if lines and lines[0].startswith("#"):
        cond = lines[0].lstrip("# ").split(";", 1)[-1].strip()
    hdr = next(i for i, l in enumerate(lines) if l.startswith("mag_AB"))
    d = np.genfromtxt(path, delimiter=",", skip_header=hdr + 1)
    mag, total, sn = d[:, 0], d[:, 2], d[:, 3]
    results = {}
    for m in np.unique(mag):
        sel = mag == m
        order = np.argsort(total[sel])
        results[float(m)] = (total[sel][order], sn[sel][order])
    return results, cond


def make_plot(results, mags, band, cond, loglog, cmap_name, single_exp, outpath,
              vlines=None, npix=1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Times-Roman serif font throughout (incl. mathtext for sqrt/R_eff).
    # Explicitly register the system Times faces so matplotlib can't fall back.
    from matplotlib import font_manager as fm
    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf",
               "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    LABEL_FS = 24   # axis label font size
    TICK_FS = 24    # tick label font size
    TITLE_FS = 24   # title font size
    LEGEND_FS = 24  # legend font size
    VLINE_FS = 24   # vertical-line annotation font size

    # bright (lowest mag) -> red, faint (highest mag) -> blue
    mags = sorted(mags)
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0.0, 1.0, len(mags)))

    fig, ax = plt.subplots(figsize=(11, 6.5))
    for c, mag in zip(colors, mags):
        t, s = results[mag]
        ax.plot(t, s, "o-", color=c, lw=1.8, ms=5, label=f"{band} = {mag:.1f} AB")

    # sqrt(t) reference anchored to the middle curve
    ref_mag = mags[len(mags) // 2]
    tr, sr = results[ref_mag]
    tref = np.linspace(tr.min(), tr.max(), 100)
    ax.plot(tref, sr[0] * np.sqrt(tref / tr[0]), "k--", lw=1.2, alpha=0.7,
            label=r"$\propto\sqrt{t}$ (bkg-limited)")

    lo, hi, _ = BANDS[band]
    binunit = "pixel" if npix == 1 else f"{npix}-pix bin"
    ax.set_xlabel(f"Total Exposure Time [min]  ({single_exp:.0f} s frames)",
                  fontsize=LABEL_FS)
    ax.set_ylabel(f"Continuum S/N per {binunit}\n"
                  f"(Median, {lo:.0f}–{hi:.0f} nm)", fontsize=LABEL_FS)
    ax.tick_params(labelsize=TICK_FS)
    scale = "log–log" if loglog else "linear"
    ax.set_title(f"Subaru PFS ETC — {band}-band ({scale})",
                 fontsize=TITLE_FS, pad=28)
    ax.text(0.5, 1.012, cond, transform=ax.transAxes, fontsize=13.5,
            ha="center", va="bottom", color="0.25")
    ax.legend(frameon=True, fontsize=LEGEND_FS, loc="center left",
              bbox_to_anchor=(1.01, 0.5), borderpad=0.3, labelspacing=0.35,
              handlelength=1.4, handletextpad=0.5)

    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
    else:
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        # keep the small log-log orientation inset on the linear version
        axins = ax.inset_axes([0.58, 0.12, 0.38, 0.38])
        for c, mag in zip(colors, mags):
            t, s = results[mag]
            axins.loglog(t, s, "o-", color=c, ms=3, lw=1.2)
        axins.loglog(tref, sr[0] * np.sqrt(tref / tr[0]), "k--", lw=1, alpha=0.6)
        axins.set_title("log–log", fontsize=8)
        axins.tick_params(labelsize=7)
        axins.grid(True, which="both", alpha=0.3)

    # vertical reference lines at requested exposure times
    if vlines:
        trans = ax.get_xaxis_transform()
        for tv in vlines:
            ax.axvline(tv, color="0.35", ls=":", lw=1.3, alpha=0.85, zorder=0)
            ax.text(tv, 0.985, f"{tv:g} min", transform=trans, rotation=90,
                    va="top", ha="right", fontsize=VLINE_FS, color="0.3")

    # explicit y-axis tick labels at useful S/N values (plain numbers)
    if loglog:
        from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
        yticks = [0.5, 1, 2, 3, 4, 5, 10, 20, 50, 100]  # only those in range are drawn
        ax.yaxis.set_major_locator(FixedLocator(yticks))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))

    fig.tight_layout()
    fig.savefig(outpath, dpi=140)
    print("plot ->", outpath)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band", default="i", choices=list(BANDS),
                    help="broadband filter for the S/N summary (default: i)")
    ap.add_argument("--mags", type=float, nargs="+",
                    default=[22.0, 22.5, 23.0, 23.5, 24.0], help="AB magnitudes")
    ap.add_argument("--exp-nums", type=int, nargs="+",
                    default=[1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48],
                    help="number of frames (total time = single_exp * N)")
    ap.add_argument("--single-exp", type=float, default=450.0,
                    help="single-frame exposure time [s] (default: 450)")
    ap.add_argument("--seeing", type=float, default=None, help="seeing FWHM [arcsec]")
    ap.add_argument("--zenith", type=float, default=None, help="zenith angle [deg]")
    ap.add_argument("--moon", type=float, default=None, help="moon phase [0-0.5]")
    ap.add_argument("--reff", type=float, default=None, help="effective radius [arcsec]")
    ap.add_argument("--loglog", action="store_true",
                    help="log-log axes (and drop the inset); writes a *_loglog.png")
    ap.add_argument("--cmap", default="rainbow_r",
                    help="matplotlib colormap, bright->faint (default rainbow_r = red->blue)")
    ap.add_argument("--from-csv", action="store_true",
                    help="replot from the cached CSV for this band (skip ETC runs)")
    ap.add_argument("--vlines", type=float, nargs="*", default=[60, 90, 120],
                    help="vertical reference lines at these total times [min] "
                         "(default 60 90 120; pass with no values to disable)")
    ap.add_argument("--bin", type=int, default=1, dest="bin",
                    help="pixel binning M: S/N scaled by sqrt(M), R=lambda/(M*disp) "
                         "(default 1)")
    args = ap.parse_args()
    args.mags = sorted(args.mags)

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)
    tag = f"{args.band}band"
    csv = os.path.join(CSV_DIR, f"time_vs_sn_{tag}.csv")
    suffix = "_loglog" if args.loglog else ""
    if args.bin != 1:
        suffix += f"_bin{args.bin}"
    png = os.path.join(PNG_DIR, f"time_vs_sn_{tag}{suffix}.png")

    if args.from_csv:
        if not os.path.exists(csv):
            ap.error(f"no cached CSV at {csv} -- run once without --from-csv first")
        results, cond = load_csv(csv)
        mags = sorted(results)
        # recover single-frame exposure from the cache for the axis label
        any_t, _ = results[mags[0]]
        single_exp = args.single_exp
        print(f"loaded {len(mags)} magnitudes x {len(any_t)} times from {csv}")
    else:
        overrides = {}
        if args.seeing is not None: overrides["SEEING"] = args.seeing
        if args.zenith is not None: overrides["ZENITH_ANG"] = args.zenith
        if args.moon is not None:   overrides["MOON_PHASE"] = args.moon
        if args.reff is not None:   overrides["REFF"] = args.reff
        cond = (f"seeing {overrides.get('SEEING', 0.8)}″, "
                f"zenith {overrides.get('ZENITH_ANG', 45)}°, "
                f"moon {overrides.get('MOON_PHASE', 0)}, "
                f"R$_{{eff}}$={overrides.get('REFF', 0.3)}″")
        results = compute_grid(args, overrides)
        mags = args.mags
        single_exp = args.single_exp
        write_csv(csv, results, mags, args.exp_nums, args.band, cond)
        print("table ->", csv)

    # pixel binning: scale continuum S/N by sqrt(M) and record effective R
    if args.bin != 1:
        f = np.sqrt(args.bin)
        results = {m: (t, s * f) for m, (t, s) in results.items()}
    lam, R = band_resolution(args.band, args.bin)
    binnote = f"{args.bin}-pix bin, R≈{R:,.0f} @ {lam:.0f} nm"
    cond = f"{cond}  |  {binnote}"

    make_plot(results, mags, args.band, cond, args.loglog, args.cmap, single_exp, png,
              vlines=args.vlines, npix=args.bin)


if __name__ == "__main__":
    main()
