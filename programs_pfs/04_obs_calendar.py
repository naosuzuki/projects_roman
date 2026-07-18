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
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun

warnings.simplefilter("ignore")          # silence ERFA dubious-year warnings

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
DEFAULT_DATADIR = "/Volumes/exdisk1/data/Roman/ltcv_sim"
ELAIS_N1_RADEC = (242.5, 54.5)           # field center (deg)


def subaru_site():
    try:
        return EarthLocation.of_site("Subaru")
    except Exception:
        return EarthLocation(lat=19.8255 * u.deg, lon=-155.4761 * u.deg, height=4139 * u.m)


def nightly_visible_hours(date_mjds, target, site, elev=30.0, twilight=-18.0, dt=0.5):
    """Visible hours per night (target elevation > elev during Sun < twilight)
    for each evening-date MJD."""
    utc_h = np.arange(4.0, 16.0 + 1e-6, dt)            # UTC 04-16 = HST 18:00-06:00
    base = Time(np.asarray(date_mjds, float), format="mjd")
    times = (base[:, None] + 1 * u.day) + utc_h[None, :] * u.hour
    flat = times.ravel()
    aa = AltAz(obstime=flat, location=site)
    talt = target.transform_to(aa).alt.deg.reshape(times.shape)
    salt = get_sun(flat).transform_to(aa).alt.deg.reshape(times.shape)
    vis = (salt < twilight) & (talt > elev)
    return vis.sum(axis=1) * dt

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
    ap.add_argument("--shift-to-date", default=None, dest="shift_to",
                    help="shift the cadence so --main-start maps to this date (e.g. 2028-06-01)")
    ap.add_argument("--main-start", default="2029-01-01",
                    help="current main-survey start date in the sim (default 2029-01-01)")
    ap.add_argument("--shade-visible", action="store_true",
                    help="shade months ELAIS-N1 is visible > --vis-hours/night from Subaru")
    ap.add_argument("--vis-hours", type=float, default=2.0,
                    help="visibility threshold in hours/night (default 2)")
    ap.add_argument("--label", default="",
                    help="filename suffix to keep variants separate (e.g. _shift2028)")
    ap.add_argument("--pilot-start", default=None, dest="pilot_start",
                    help="re-shift the early pilot group to start on this date "
                         "(e.g. 2027-04-15), leaving the main survey unchanged")
    ap.add_argument("--pilot-gap", type=float, default=90.0,
                    help="min gap (days) separating the pilot group from the main survey")
    ap.add_argument("--phase-labels", action="store_true", dest="phase_labels",
                    help="open a gap between the panels and mark the Subaru phase "
                         "starts (Phase-I HSC, Phase-II PFS, Phase-III PFS)")
    args = ap.parse_args()

    shift_days = 0.0
    if args.shift_to is not None:
        shift_days = Time(args.shift_to).mjd - Time(args.main_start).mjd

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    epochs = collect_epochs(args.datadir, args.field, args.model)
    bands = [b for b in BAND_WL if b in epochs]            # wavelength order
    if shift_days:
        print(f"shifting cadence by {shift_days:.0f} days "
              f"({args.main_start} -> {args.shift_to})")

    # optional: re-shift the early "pilot" group to start on --pilot-start,
    # leaving the main survey untouched.
    pilot_extra, pilot_end = 0.0, -1e18
    if args.pilot_start is not None:
        nights = np.unique(np.round(
            np.concatenate([epochs[b] for b in bands]) + shift_days, 0))
        nights.sort()
        g = np.where(np.diff(nights) > args.pilot_gap)[0]
        pilot_end = nights[g[0]] if len(g) else nights[-1]
        pilot_extra = Time(args.pilot_start).mjd - nights[0]
        print(f"  pilot group {Time(nights[0],format='mjd').iso[:10]}.."
              f"{Time(pilot_end,format='mjd').iso[:10]} -> +{pilot_extra:.0f} d "
              f"(starts {args.pilot_start})")

    def shmjd(m):
        """Final shifted MJD: global shift + extra pilot shift for pilot epochs."""
        s = np.asarray(m, float) + shift_days
        return s + pilot_extra * (s <= pilot_end + 0.5)

    sh = {b: shmjd(epochs[b]) for b in bands}
    sh_all = np.concatenate([sh[b] for b in bands])
    t0, t1 = sh_all.min(), sh_all.max()
    print(f"{args.field}: {len(bands)} bands, "
          f"{Time(t0,format='mjd').iso[:10]} -> {Time(t1,format='mjd').iso[:10]}")
    for b in bands:
        print(f"  {b:10s}: {len(epochs[b]):4d} epochs")

    # ---- CSV (unique nights per band) ----
    csv = os.path.join(CSV_DIR, f"04_obs_calendar_{args.field}{args.label}.csv")
    with open(csv, "w") as fo:
        fo.write("band,MJD,date\n")
        for b in bands:
            for m, s in zip(epochs[b], sh[b]):
                fo.write(f"{b},{m:.3f},{Time(s,format='mjd').iso[:10]}\n")
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
    LABEL_FS, TICK_FS, TITLE_FS = 24, 19, 17

    # color by wavelength (blue=short -> red=long)
    wl = np.array([BAND_WL[b] for b in bands], float)
    cmap = plt.get_cmap("rainbow")
    norm = (wl - wl.min()) / (wl.max() - wl.min())
    colors = [cmap(x) for x in norm]

    dnum = {b: mdates.date2num(Time(sh[b], format="mjd").to_datetime()) for b in bands}
    alln = mdates.date2num(Time(sh_all, format="mjd").to_datetime())

    fig, (axtop, axev) = plt.subplots(
        2, 1, figsize=(13, 7), height_ratios=[1, 3], sharex=True,
        constrained_layout=True)

    # top: observations per ~month
    nb = max(1, int((alln.max() - alln.min()) / 30))
    axtop.hist(alln, bins=nb, color="0.4")
    axtop.set_ylabel("Obs/Month", fontsize=LABEL_FS)
    axtop.tick_params(labelsize=TICK_FS)
    subtitle = (f"Main Survey from {args.shift_to}" if shift_days
                else f"{Time(t0,format='mjd').iso[:7]} to {Time(t1,format='mjd').iso[:7]}")
    axtop.set_title(f"Roman HLTDS Observation Cadence — {args.field}  ({subtitle})",
                    fontsize=LABEL_FS)

    # bottom: event raster per band
    axev.eventplot([dnum[b] for b in bands],
                   lineoffsets=np.arange(len(bands)), linelengths=0.8,
                   colors=colors)
    axev.set_yticks(np.arange(len(bands)))
    axev.set_yticklabels(bands, fontsize=TICK_FS)
    axev.set_ylabel("Band", fontsize=LABEL_FS)
    axev.set_xlabel("Calendar Date", fontsize=LABEL_FS)
    axev.tick_params(labelsize=TICK_FS)
    axev.grid(True, axis="x", alpha=0.3)

    axev.xaxis.set_major_locator(mdates.YearLocator())
    axev.xaxis.set_minor_locator(mdates.MonthLocator((1, 4, 7, 10)))
    axev.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    # shade the months ELAIS-N1 is visible > vis_hours/night from Subaru
    if args.shade_visible:
        site = subaru_site()
        target = SkyCoord(ELAIS_N1_RADEC[0] * u.deg, ELAIS_N1_RADEC[1] * u.deg)
        d0 = int(np.floor(sh_all.min()))
        d1 = int(np.ceil(sh_all.max()))
        days = np.arange(d0, d1 + 1)
        vh = nightly_visible_hours(days, target, site, elev=30.0, twilight=-18.0)
        mask = vh > args.vis_hours
        ddn = mdates.date2num(Time(days, format="mjd").to_datetime())
        first, i = True, 0
        while i < len(mask):
            if mask[i]:
                j = i
                while j + 1 < len(mask) and mask[j + 1]:
                    j += 1
                lbl = ("Visibility from Subaru" if first else None)
                for ax in (axtop, axev):
                    ax.axvspan(ddn[i], ddn[j], color="#2ca02c", alpha=0.13, lw=0,
                               zorder=0, label=(lbl if ax is axev else None))
                first, i = False, j + 1
            else:
                i += 1
        axev.legend(loc="upper right", fontsize=18, framealpha=0.9)

    # mark the (shifted) survey start (dotted line only, no text)
    if shift_days:
        x0 = mdates.date2num(Time(args.shift_to).to_datetime())
        for ax in (axtop, axev):
            ax.axvline(x0, color="k", ls=":", lw=1.3, zorder=3)

    # mark the Subaru phase starts in a gap between the panels
    if args.phase_labels:
        fig.get_layout_engine().set(hspace=0.14)
        PHASES = (("2027-04-15", "Phase-I HSC"),
                  ("2028-06-01", "Phase-II PFS"),
                  ("2030-07-01", "Phase-III PFS"))
        tr = axev.get_xaxis_transform()          # x = data (date), y = axes fraction
        for dstr, lab in PHASES:
            xd = mdates.date2num(Time(dstr).to_datetime())
            axev.plot([xd, xd], [1.0, 1.05], transform=tr, color="k", lw=1.5,
                      clip_on=False)
            axev.text(xd, 1.06, " " + lab, transform=tr, fontsize=19,
                      fontweight="bold", ha="left", va="bottom", clip_on=False)

    png = os.path.join(PNG_DIR, f"04_obs_calendar_{args.field}{args.label}.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
