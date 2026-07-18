#!/usr/bin/env python3
"""
05_desi_visibility.py -- DESI counterpart of programs_pfs/05_visibility.py.

Visibility of a field (default ELAIS-N1) from the DESI / Mayall 4-m site at
Kitt Peak (KPNO, lat +31.96 deg) over a year -- the analog of Figure 13
(fig:vis) in the Roman-Subaru paper, recomputed for KPNO instead of Maunakea.

Because ELAIS-N1 (Dec +54.5 deg) lies much closer to KPNO's latitude than to
Subaru's (+19.83 deg), it transits far higher from Kitt Peak (~67.5 deg,
airmass ~1.08, vs 55.3 deg / 1.22 from Subaru) and stays observable a larger
fraction of the year -- a key part of the DESI host-follow-up case.

Computes target altitude/airmass and Sun altitude on a (date, time-of-night)
grid and shows:
  * top:    observable dark hours per night (airmass < limit, Sun < twilight)
  * bottom: airmass map over the year vs MST clock time, with twilight edges
and a dedicated stacked-bar "visible hours" figure (the Fig 13 analog).

    python 05_desi_visibility.py                       # ELAIS-N1, 2029
    python 05_desi_visibility.py --ra 242.5 --dec 54.5 --year 2030
"""
import os
import argparse
import warnings
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun

warnings.simplefilter("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")

# ELAIS-N1 field center (Roman HLTDS north), deg
ELAIS_N1 = (242.5, 54.5)
MST_OFFSET = 7.0            # MST = UTC - 7 (Arizona, no DST)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ra", type=float, default=ELAIS_N1[0], help="RA [deg]")
    ap.add_argument("--dec", type=float, default=ELAIS_N1[1], help="Dec [deg]")
    ap.add_argument("--name", default="ELAIS-N1")
    ap.add_argument("--year", type=int, default=2029)
    ap.add_argument("--elev-limit", type=float, default=30.0, dest="elevlim",
                    help="minimum elevation [deg] defining visibility (default 30)")
    ap.add_argument("--twilight", type=float, default=-18.0,
                    help="Sun altitude defining dark time [deg] (default -18)")
    args = ap.parse_args()

    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)

    try:
        site = EarthLocation.of_site("kpno")
    except Exception:
        site = EarthLocation(lat=31.9633 * u.deg, lon=-111.5997 * u.deg, height=2120 * u.m)
    target = SkyCoord(args.ra * u.deg, args.dec * u.deg)
    alt_lim = args.elevlim                                      # visibility: elevation > alt_lim
    am_at_lim = 1.0 / np.cos(np.radians(90.0 - alt_lim))       # = 2.0 for 30 deg

    # grid: evening dates over the year x time-of-night (UTC hours covering the night)
    base = Time(f"{args.year}-01-01 00:00:00")
    day_off = np.arange(0, 365, 1)
    utc_h = np.arange(1.0, 15.001, 0.25)              # UTC 01-15 = MST 18:00-08:00
    dates = base + day_off * u.day                    # evening date (UTC midnight)
    # night belongs to UTC (date+1) 01:00..15:00
    times = (dates[:, None] + 1 * u.day) + utc_h[None, :] * u.hour
    flat = times.ravel()
    aa = AltAz(obstime=flat, location=site)
    talt = target.transform_to(aa).alt.deg.reshape(times.shape)
    salt = get_sun(flat).transform_to(aa).alt.deg.reshape(times.shape)

    with np.errstate(invalid="ignore"):
        airmass = 1.0 / np.cos(np.radians(90.0 - talt))
    airmass[talt <= 0] = np.nan

    dark = salt < args.twilight
    observable = dark & (talt > alt_lim)
    dt = 0.25
    hrs_obs = (observable * dt).sum(axis=1)            # per night, elevation > alt_lim
    hrs_dark = (dark * dt).sum(axis=1)

    # monthly aggregates
    months = np.array([d.month for d in dates.datetime])
    MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    mon_mean = np.array([hrs_obs[months == m].mean() if (months == m).any() else 0
                         for m in range(1, 13)])
    mon_tot = np.array([hrs_obs[months == m].sum() for m in range(1, 13)])

    # elevation-band bins during dark: red <45, green 45-50, blue >50 deg
    # (same thresholds as the Subaru figure for a direct comparison; from KPNO
    # the field reaches ~67.5 deg, so the >50 band dominates -- the visual point.)
    E1, E2 = 45.0, 50.0
    red_h = ((dark & (talt > alt_lim) & (talt < E1)) * dt).sum(axis=1)   # 30-45
    grn_h = ((dark & (talt >= E1) & (talt < E2)) * dt).sum(axis=1)       # 45-50
    blu_h = ((dark & (talt >= E2)) * dt).sum(axis=1)                     # >50

    def _monthly(arr):
        return np.array([arr[months == m].mean() if (months == m).any() else 0
                         for m in range(1, 13)])
    mon_red, mon_grn, mon_blu = _monthly(red_h), _monthly(grn_h), _monthly(blu_h)

    # MST clock for y-axis (monotonic 18->32)
    mst = utc_h - MST_OFFSET
    mst = np.where(mst < 12, mst + 24, mst)

    # stats
    bestidx = int(np.argmax(hrs_obs))
    transit_alt = 90 - abs(args.dec - site.lat.deg)
    print(f"{args.name}  (RA {args.ra}, Dec {args.dec})  from KPNO (Kitt Peak), {args.year}")
    print(f"  transit altitude ~ {transit_alt:.1f} deg "
          f"(min airmass ~ {1/np.cos(np.radians(abs(args.dec-site.lat.deg))):.2f})")
    print(f"  visibility = elevation > {alt_lim:.0f} deg (airmass < {am_at_lim:.2f}) during astro. night")
    print(f"  nights with >0 h visible: {int((hrs_obs>0).sum())}/365")
    print(f"  nights with >4 h visible: {int((hrs_obs>4).sum())}")
    print(f"  best night: {Time(dates[bestidx]).iso[:10]}  ->  {hrs_obs[bestidx]:.1f} h visible")
    print(f"  total visible hours / yr: {hrs_obs.sum():.0f}")
    print("  month   mean h/night   total h/month")
    for i, mn in enumerate(MON):
        print(f"    {mn}       {mon_mean[i]:5.1f}          {mon_tot[i]:6.0f}")

    # CSV
    csv = os.path.join(CSV_DIR, f"05_desi_visibility_{args.name}_{args.year}.csv")
    with open(csv, "w") as fo:
        fo.write("evening_date,dark_hours,observable_hours\n")
        for d, hd, ho in zip(dates, hrs_dark, hrs_obs):
            fo.write(f"{Time(d).iso[:10]},{hd:.2f},{ho:.2f}\n")
    print("table ->", csv)

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

    dnum = mdates.date2num(Time(dates).to_datetime())

    fig, (axtop, axmap) = plt.subplots(
        2, 1, figsize=(13, 8), height_ratios=[1, 2.2], sharex=True,
        constrained_layout=True)

    axtop.fill_between(dnum, hrs_obs, color="#e8820c", alpha=0.85,
                       label=f"elevation > {alt_lim:.0f}$^\\circ$")
    axtop.plot(dnum, hrs_dark, color="0.4", lw=1, ls="--", label="astro. dark hours")
    axtop.set_ylabel("visible\nhours / night", fontsize=LABEL_FS)
    axtop.tick_params(labelsize=TICK_FS)
    axtop.legend(fontsize=11, loc="upper right")
    axtop.set_title(f"Visibility of {args.name} from DESI/KPNO (Kitt Peak), {args.year}  "
                    f"(elevation $>{alt_lim:.0f}^\\circ$)", fontsize=TITLE_FS)
    axtop.grid(True, alpha=0.3)

    # bottom: airmass map (only dark, target up); twilight edges via Sun alt
    am_plot = np.ma.masked_where(~(dark & (talt > 0)), airmass)
    X, Y = np.meshgrid(dnum, mst, indexing="ij")
    pm = axmap.pcolormesh(X, Y, np.clip(am_plot, 1, 3), cmap="viridis_r",
                          vmin=1, vmax=3, shading="nearest")
    axmap.contour(X, Y, np.where(dark, talt, np.nan), levels=[alt_lim],
                  colors="white", linewidths=1.3)
    axmap.contour(X, Y, salt, levels=[0.0], colors="orange", linewidths=1.2)
    axmap.contour(X, Y, salt, levels=[args.twilight], colors="red", linewidths=1.0)

    axmap.set_ylabel("MST clock time", fontsize=LABEL_FS)
    axmap.set_xlabel("date", fontsize=LABEL_FS)
    axmap.set_yticks([18, 21, 24, 27, 30])
    axmap.set_yticklabels(["18:00", "21:00", "00:00", "03:00", "06:00"], fontsize=TICK_FS)
    axmap.tick_params(labelsize=TICK_FS)
    axmap.xaxis.set_major_locator(mdates.MonthLocator())
    axmap.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    cbar = fig.colorbar(pm, ax=axmap, pad=0.01)
    cbar.set_label("airmass (sec z)", fontsize=LABEL_FS)
    cbar.ax.tick_params(labelsize=TICK_FS)
    axmap.text(0.01, 0.02, f"white: elevation = {alt_lim:.0f}$^\\circ$ (airmass {am_at_lim:.1f});  "
               "orange: sunset/rise;  red: astro. twilight",
               transform=axmap.transAxes, fontsize=9, color="0.2",
               bbox=dict(fc="white", ec="0.7", alpha=0.8))

    png = os.path.join(PNG_DIR, f"05_desi_visibility_{args.name}_{args.year}.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)

    # ---- dedicated "observable hours" figure (the Fig 13 analog) ----
    fig2, (axn, axm) = plt.subplots(2, 1, figsize=(13, 8), constrained_layout=True)

    s1 = red_h
    s2 = red_h + grn_h
    s3 = red_h + grn_h + blu_h
    axn.fill_between(dnum, 0, s1, color="#d62728", alpha=0.85,
                     label="Elevation $<45^\\circ$")
    axn.fill_between(dnum, s1, s2, color="#2ca02c", alpha=0.85,
                     label="Elevation 45--50$^\\circ$")
    axn.fill_between(dnum, s2, s3, color="#3b7dd8", alpha=0.85,
                     label="Elevation $>50^\\circ$")
    axn.plot(dnum, hrs_dark, color="0.35", lw=1.2, ls="--", label="Astro. dark hours")
    axn.set_ylabel("Visible hours / night", fontsize=18)
    axn.set_xlabel("Date", fontsize=18)
    axn.tick_params(labelsize=TICK_FS)
    axn.legend(fontsize=12, loc="upper right")
    axn.grid(True, alpha=0.3)
    axn.set_title(f"Visible hours for {args.name} from DESI/KPNO (Kitt Peak), {args.year}  "
                  f"(elevation $>{alt_lim:.0f}^\\circ$)", fontsize=TITLE_FS)
    axn.xaxis.set_major_locator(mdates.MonthLocator())
    axn.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    axn.set_xlim(dnum.min(), dnum.max())

    x = np.arange(12)
    axm.bar(x, mon_red, width=0.62, color="#d62728",
            label="Elevation $<45^\\circ$")
    axm.bar(x, mon_grn, width=0.62, bottom=mon_red, color="#2ca02c",
            label="Elevation 45--50$^\\circ$")
    axm.bar(x, mon_blu, width=0.62, bottom=mon_red + mon_grn, color="#3b7dd8",
            label="Elevation $>50^\\circ$")
    for i in range(12):
        if mon_mean[i] > 0:
            axm.text(i, mon_mean[i] + 0.1, f"{mon_tot[i]:.0f}h", ha="center",
                     va="bottom", fontsize=13, color="0.3")
    axm.set_xticks(x)
    axm.set_xticklabels(MON, fontsize=TICK_FS)
    axm.set_ylabel("Mean visible hours / night", fontsize=18)
    axm.set_xlabel("Month  (number above bar = total visible hours)", fontsize=18)
    axm.tick_params(labelsize=TICK_FS)
    axm.legend(fontsize=11, loc="upper right")
    axm.grid(True, axis="y", alpha=0.3)

    png2 = os.path.join(PNG_DIR, f"05_desi_visibility_hours_{args.name}_{args.year}.png")
    fig2.savefig(png2, dpi=140)
    print("plot ->", png2)


if __name__ == "__main__":
    main()
