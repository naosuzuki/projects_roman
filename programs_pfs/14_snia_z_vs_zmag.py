#!/usr/bin/env python3
"""
14_snia_z_vs_zmag.py -- redshift vs peak Roman Z-band (Z087) magnitude for
Type Ia SNe in the HOURGLASS2 WIDE+DEEP simulation
(../data/OUT_ROMAN_WIDE+DEEP.TEXT), split by Subaru (Hawaii) observability.

Selects true SNe Ia (SIM_GENTYPE = 10) with a valid PEAKMAG_Z (> 0; the
-888/-9 sentinels mean the band was not observed or the SN failed cuts).

The sim uses a dummy sky position and an arbitrary MJD zero-point, so per
project convention the sample is treated as ELAIS-N1 on the adopted survey
timeline: SIM_PEAKMJD is shifted so the survey start maps to 2028-06-01
(2-year survey through 2030-05). A SN counts as observable from Hawaii if
ELAIS-N1 is visible from Subaru for > VIS_HOURS h (elevation > 30 deg during
astronomical dark) on its peak night -- the Feb-Aug observing season.
Observable SNe are drawn solid; invisible ones smaller with alpha = 0.3.

    python 14_snia_z_vs_zmag.py
"""
import os
import warnings
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun

warnings.simplefilter("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
DATA = os.path.join(HERE, "..", "data", "OUT_ROMAN_WIDE+DEEP.TEXT")

ELAIS_N1 = (242.5, 54.5)          # field center (deg)
SURVEY_START = "2028-06-01"       # adopted (shifted) main-survey start
VIS_HOURS = 2.0                   # observable if > this many visible h/night
ELEV_LIM = 30.0                   # deg
TWILIGHT = -18.0                  # Sun altitude for astronomical dark


def subaru_site():
    try:
        return EarthLocation.of_site("Subaru")
    except Exception:
        return EarthLocation(lat=19.8255 * u.deg, lon=-155.4761 * u.deg,
                             height=4139 * u.m)


def nightly_visible_hours(night_mjds, target, site, elev=ELEV_LIM,
                          twilight=TWILIGHT, dt=0.5):
    """Visible hours per night (target elevation > elev during Sun < twilight)
    for each evening-date MJD.  Same recipe as 04_obs_calendar.py."""
    utc_h = np.arange(4.0, 16.0 + 1e-6, dt)            # UTC 04-16 = HST 18:00-06:00
    base = Time(np.asarray(night_mjds, float), format="mjd")
    times = (base[:, None] + 1 * u.day) + utc_h[None, :] * u.hour
    flat = times.ravel()
    aa = AltAz(obstime=flat, location=site)
    talt = target.transform_to(aa).alt.deg.reshape(times.shape)
    salt = get_sun(flat).transform_to(aa).alt.deg.reshape(times.shape)
    vis = (salt < twilight) & (talt > elev)
    return vis.sum(axis=1) * dt


def main():
    df = pd.read_csv(DATA, comment="#", sep=r"\s+")
    df = df.drop(columns=["VARNAMES:"], errors="ignore")

    ia = df[(df["SIM_GENTYPE"] == 10) & (df["PEAKMAG_Z"] > 0)].copy()
    print(f"SNe Ia with valid peak Z mag: {len(ia)}")

    # map the sim clock onto the adopted timeline: survey start -> 2028-06-01
    shift = Time(SURVEY_START).mjd - ia["SIM_PEAKMJD"].min()
    ia["PEAKMJD_CAL"] = ia["SIM_PEAKMJD"] + shift
    t0, t1 = ia["PEAKMJD_CAL"].min(), ia["PEAKMJD_CAL"].max()
    print(f"peak epochs shifted by {shift:+.0f} d -> "
          f"{Time(t0, format='mjd').iso[:10]} .. {Time(t1, format='mjd').iso[:10]}")

    # Subaru visibility of ELAIS-N1 on each peak night
    nights = np.arange(int(np.floor(t0)), int(np.ceil(t1)) + 1)
    site = subaru_site()
    target = SkyCoord(ELAIS_N1[0] * u.deg, ELAIS_N1[1] * u.deg)
    vh = nightly_visible_hours(nights, target, site)
    vis_of_night = dict(zip(nights, vh))
    ia["VIS_HOURS"] = [vis_of_night[int(np.floor(m))] for m in ia["PEAKMJD_CAL"]]
    obs = ia["VIS_HOURS"] > VIS_HOURS

    print(f"observable from Hawaii (> {VIS_HOURS:g} h/night at peak): "
          f"{obs.sum()}  |  invisible: {(~obs).sum()}")
    for f in ("WIDE", "DEEP"):
        s = ia["FIELD"] == f
        print(f"  {f}: {int((obs & s).sum())} observable / {int(s.sum())}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, ax = plt.subplots(figsize=(9, 6))
    for f, c in (("WIDE", "blue"), ("DEEP", "red")):
        s = ia["FIELD"] == f
        vis, inv = ia[s & obs], ia[s & ~obs]
        ax.scatter(vis["zCMB"], vis["PEAKMAG_Z"], s=8, c=c, alpha=1.0,
                   label=f"{f} Observable  ($N={len(vis):,}$)", rasterized=True)
        ax.scatter(inv["zCMB"], inv["PEAKMAG_Z"], s=3, c=c, alpha=0.3,
                   label=f"{f} Invisible  ($N={len(inv):,}$)", rasterized=True)
    ax.axhline(24.0, color="k", linestyle=":", linewidth=1.2)
    ax.text(0.03, 24.08, "$Z=24$ (PFS Observable)", fontsize=12, va="bottom")

    ax.set_xlabel("Redshift", fontsize=18)
    ax.set_ylabel("Peak Z-Band Magnitude (Z087, AB)", fontsize=18)
    ax.set_title(f"Roman HLTDS Type Ia SNe: Redshift vs Peak Z-Band Magnitude  "
                 f"($N={len(ia):,}$)", fontsize=15)
    ax.tick_params(labelsize=13)
    ax.set_xlim(0, 2.0)
    ax.legend(fontsize=12, loc="lower right")
    ax.grid(True, alpha=0.3)

    png = os.path.join(PNG_DIR, "14_snia_z_vs_zmag.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
