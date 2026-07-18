#!/usr/bin/env python3
"""
04_desi_obs_calendar.py -- DESI counterpart of programs_pfs/04_obs_calendar.py.

Reproduces the Roman HLTDS observation-cadence calendar for ELAIS-N1 (the
analog of Figure 10 / fig:schedule_shift in the Roman-Subaru paper), but shades
the **observing-time window** computed for DESI / the Mayall 4-m at Kitt Peak
(KPNO) instead of Subaru, in a light-orange tone.

The Roman cadence (per-band visit MJDs) is identical to the Subaru study -- it
is a property of the survey, not the follow-up telescope -- so we read the
already-committed epoch list

    ../programs_pfs/outputs/csv/04_obs_calendar_ELAIS-N1.csv   (band, MJD, date)

and recompute only the visibility shading from KPNO. No external disk / FITS
needed. The cadence shift and pilot-group relocation match the Subaru figure so
the two panels are directly comparable:

    python 04_desi_obs_calendar.py            # shift to 2028-06-01, pilot 2027-04-15

Output: outputs/png/04_desi_obs_calendar_ELAIS-N1_shift2028.png
"""
import os
import argparse
import warnings
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun

warnings.simplefilter("ignore")          # silence ERFA dubious-year warnings

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
# Roman cadence (per-band visit MJDs) committed by the PFS study:
EPOCHS_CSV = os.path.join(HERE, "..", "programs_pfs", "outputs", "csv",
                          "04_obs_calendar_ELAIS-N1.csv")
ELAIS_N1_RADEC = (242.5, 54.5)           # field center (deg)

# Roman bands present in the sim, ordered by effective wavelength (nm)
BAND_WL = {"R062-R": 620, "Z087-Z": 870, "Y106-Y": 1060,
           "J129-J": 1290, "H158-H": 1580, "F184-F": 1840}

# Light-orange shade for the DESI observing-time window
DESI_SHADE = "#FFB347"


def kpno_site():
    """Mayall / Kitt Peak National Observatory."""
    try:
        return EarthLocation.of_site("kpno")
    except Exception:
        return EarthLocation(lat=31.9633 * u.deg, lon=-111.5997 * u.deg,
                             height=2120 * u.m)


def nightly_visible_hours(date_mjds, target, site, elev=30.0, twilight=-18.0, dt=0.5):
    """Visible hours per night (target elevation > elev during Sun < twilight)
    for each evening-date MJD. KPNO is UTC-7 (MST), so the local night spans
    roughly UTC 02-15; we sample a generous UTC 01-16 window."""
    utc_h = np.arange(1.0, 16.0 + 1e-6, dt)
    base = Time(np.asarray(date_mjds, float), format="mjd")
    times = (base[:, None] + 1 * u.day) + utc_h[None, :] * u.hour
    flat = times.ravel()
    aa = AltAz(obstime=flat, location=site)
    talt = target.transform_to(aa).alt.deg.reshape(times.shape)
    salt = get_sun(flat).transform_to(aa).alt.deg.reshape(times.shape)
    vis = (salt < twilight) & (talt > elev)
    return vis.sum(axis=1) * dt


def collect_epochs(csv):
    """Union of (rounded) observation MJDs per band from the committed CSV."""
    d = np.genfromtxt(csv, delimiter=",", names=True, dtype=None, encoding="utf-8")
    band = np.asarray(d["band"]).astype(str)
    mjd = np.asarray(d["MJD"]).astype(float)
    per_band = {}
    for bn in np.unique(band):
        per_band[bn] = np.array(sorted(set(np.round(mjd[band == bn], 3).tolist())))
    return per_band


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--epochs-csv", default=EPOCHS_CSV)
    ap.add_argument("--field", default="ELAIS-N1")
    ap.add_argument("--shift-to-date", default="2028-06-01", dest="shift_to",
                    help="shift the cadence so --main-start maps to this date")
    ap.add_argument("--main-start", default="2029-01-01",
                    help="current main-survey start date in the sim")
    ap.add_argument("--vis-hours", type=float, default=2.0,
                    help="visibility threshold in hours/night (default 2)")
    ap.add_argument("--no-shade", action="store_true",
                    help="omit the DESI visibility shading")
    ap.add_argument("--label", default="_shift2028",
                    help="filename suffix (default _shift2028)")
    ap.add_argument("--pilot-start", default="2027-04-15", dest="pilot_start",
                    help="re-shift the early pilot group to start on this date, "
                         "leaving the main survey unchanged (set '' to disable)")
    ap.add_argument("--pilot-gap", type=float, default=90.0,
                    help="min gap (days) separating the pilot group from main survey")
    args = ap.parse_args()

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    shift_days = 0.0
    if args.shift_to:
        shift_days = Time(args.shift_to).mjd - Time(args.main_start).mjd

    epochs = collect_epochs(args.epochs_csv)
    bands = [b for b in BAND_WL if b in epochs]            # wavelength order
    if shift_days:
        print(f"shifting cadence by {shift_days:.0f} days "
              f"({args.main_start} -> {args.shift_to})")

    # optional: re-shift the early "pilot" group to start on --pilot-start,
    # leaving the main survey untouched (matches the Subaru figure).
    pilot_extra, pilot_end = 0.0, -1e18
    if args.pilot_start:
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

    dnum = {b: mdates.date2num(Time(sh[b], format="mjd").to_datetime()) for b in bands}
    alln = mdates.date2num(Time(sh_all, format="mjd").to_datetime())

    fig, (axtop, axev) = plt.subplots(
        2, 1, figsize=(13, 7), height_ratios=[1, 3], sharex=True,
        constrained_layout=True)

    # top: observations per ~month
    nb = max(1, int((alln.max() - alln.min()) / 30))
    axtop.hist(alln, bins=nb, color="0.4")
    axtop.set_ylabel("obs / month", fontsize=LABEL_FS)
    axtop.tick_params(labelsize=TICK_FS)
    subtitle = (f"main survey from {args.shift_to}" if shift_days
                else f"{Time(t0,format='mjd').iso[:7]} to {Time(t1,format='mjd').iso[:7]}")
    axtop.set_title(f"Roman HLTDS observation cadence — {args.field}  ({subtitle})",
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

    # shade the months ELAIS-N1 is visible > vis_hours/night from KPNO (DESI)
    if not args.no_shade:
        site = kpno_site()
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
                lbl = (f"DESI/KPNO: visible $>{args.vis_hours:g}$ h/night" if first else None)
                for ax in (axtop, axev):
                    ax.axvspan(ddn[i], ddn[j], color=DESI_SHADE, alpha=0.15, lw=0,
                               zorder=0, label=(lbl if ax is axev else None))
                first, i = False, j + 1
            else:
                i += 1
        axev.legend(loc="upper right", fontsize=10, framealpha=0.9)

    # mark the (shifted) survey start
    if shift_days:
        x0 = mdates.date2num(Time(args.shift_to).to_datetime())
        for ax in (axtop, axev):
            ax.axvline(x0, color="k", ls=":", lw=1.3, zorder=3)
        axev.text(x0, len(bands) - 0.4, f" survey start\n {args.shift_to}",
                  fontsize=9, va="top", ha="left")

    png = os.path.join(PNG_DIR, f"04_desi_obs_calendar_{args.field}{args.label}.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
