#!/usr/bin/env python3
"""
04_obs_calendar.py -- visualize the Roman HLTDS observation cadence for a field
as calendar dates: when each band was observed over the survey.

The visit epochs are the MJDs in the PHOT tables (shared across all SNe in a
field via the SIMLIB). We take the union of unique epochs per band, convert
MJD -> calendar date, and draw:
  * top:    histogram of all observations per ~month (observing seasons)
  * bottom: an event raster, one row per band (the cadence)

    python 04_obs_calendar.py                       # ELAIS-N1
    python 04_obs_calendar.py --field EDF-S
"""
import os
import glob
import argparse
import warnings
import numpy as np
from astropy.io import fits
from astropy.time import Time

warnings.simplefilter("ignore")          # silence ERFA dubious-year warnings

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
DEFAULT_DATADIR = "/Volumes/exdisk1/data/Roman/ltcv_sim"

# Roman bands present in the sim, ordered by effective wavelength (nm)
BAND_WL = {"R062-R": 620, "Z087-Z": 870, "Y106-Y": 1060,
           "J129-J": 1290, "H158-H": 1580, "F184-F": 1840}


def collect_epochs(datadir, field, model):
    """Union of (rounded) observation MJDs per band across a field's PHOT files."""
    phot = [f for f in glob.glob(os.path.join(datadir, field, f"*{model}-*_PHOT.FITS*"))
            if not os.path.basename(f).startswith("._")]
    per_band = {}
    for f in phot:
        with fits.open(f) as h:
            d = h[1].data
            if d is None or len(d) == 0:
                continue
            m = np.asarray(d["MJD"], float)
            b = np.asarray(d["BAND"]).astype(str)
        good = m > 0                      # drop -777 separator rows
        m, b = m[good], b[good]
        for bn in np.unique(b):
            per_band.setdefault(bn, set()).update(np.round(m[b == bn], 3).tolist())
    return {bn: np.array(sorted(v)) for bn, v in per_band.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datadir", default=DEFAULT_DATADIR)
    ap.add_argument("--field", default="ELAIS-N1")
    ap.add_argument("--model", default="SNIaMODEL00",
                    help="SNANA model whose PHOT defines the cadence (default SNIaMODEL00)")
    args = ap.parse_args()

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    epochs = collect_epochs(args.datadir, args.field, args.model)
    bands = [b for b in BAND_WL if b in epochs]            # wavelength order
    allmjd = np.concatenate([epochs[b] for b in bands])
    t0, t1 = allmjd.min(), allmjd.max()
    print(f"{args.field}: {len(bands)} bands, "
          f"{Time(t0,format='mjd').iso[:10]} -> {Time(t1,format='mjd').iso[:10]}")
    for b in bands:
        print(f"  {b:10s}: {len(epochs[b]):4d} epochs")

    # ---- CSV (unique nights per band) ----
    csv = os.path.join(CSV_DIR, f"04_obs_calendar_{args.field}.csv")
    with open(csv, "w") as fo:
        fo.write("band,MJD,date\n")
        for b in bands:
            for m in epochs[b]:
                fo.write(f"{b},{m:.3f},{Time(m,format='mjd').iso[:10]}\n")
    print("epochs ->", csv)

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib import font_manager as fm

    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    LABEL_FS, TICK_FS, TITLE_FS = 16, 12, 16

    # color by wavelength (blue=short -> red=long)
    wl = np.array([BAND_WL[b] for b in bands], float)
    cmap = plt.get_cmap("rainbow")
    norm = (wl - wl.min()) / (wl.max() - wl.min())
    colors = [cmap(x) for x in norm]

    dnum = {b: mdates.date2num(Time(epochs[b], format="mjd").to_datetime()) for b in bands}
    alln = mdates.date2num(Time(allmjd, format="mjd").to_datetime())

    fig, (axtop, axev) = plt.subplots(
        2, 1, figsize=(13, 7), height_ratios=[1, 3], sharex=True,
        constrained_layout=True)

    # top: observations per ~month
    nb = max(1, int((alln.max() - alln.min()) / 30))
    axtop.hist(alln, bins=nb, color="0.4")
    axtop.set_ylabel("obs / month", fontsize=LABEL_FS)
    axtop.tick_params(labelsize=TICK_FS)
    axtop.set_title(f"Roman HLTDS observation cadence — {args.field}  "
                    f"({Time(t0,format='mjd').iso[:7]} to {Time(t1,format='mjd').iso[:7]})",
                    fontsize=TITLE_FS)

    # bottom: event raster per band
    axev.eventplot([dnum[b] for b in bands],
                   lineoffsets=np.arange(len(bands)), linelengths=0.8,
                   colors=colors)
    axev.set_yticks(np.arange(len(bands)))
    axev.set_yticklabels(bands, fontsize=TICK_FS)
    axev.set_ylabel("band", fontsize=LABEL_FS)
    axev.set_xlabel("calendar date", fontsize=LABEL_FS)
    axev.tick_params(labelsize=TICK_FS)
    axev.grid(True, axis="x", alpha=0.3)

    axev.xaxis.set_major_locator(mdates.YearLocator())
    axev.xaxis.set_minor_locator(mdates.MonthLocator((1, 4, 7, 10)))
    axev.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    png = os.path.join(PNG_DIR, f"04_obs_calendar_{args.field}.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
