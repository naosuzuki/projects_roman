#!/usr/bin/env python3
"""
07_fiber_budget.py -- PFS fiber-budget simulation for the ELAIS-N1 SN program.

Survey model (see ms.tex):
  * 7-pointing hexagonal flower; A (0deg) at start of each dark run, B (30deg)
    at the end. 1 h per pointing, 7 h per configuration.
  * Twice per month, ONLY in months ELAIS-N1 is observable from Subaru
    (>2 h/night, Feb-Aug). Survey: 2028-06-01 .. 2030-05-31.
  * Targets = SNe Ia + CC with Z-band (Roman Z087) light-curve mag < 24 (the
    transient is a fiber target on every visit while Z<24).
  * When a SN fades past Z=24, observe its host if host Z<25, integrating until
    S/N=5 at 10-pixel binning (S/N accumulates as sqrt(hours); hours-to-5 from
    the PFS ETC, with host i-band continuum mag approximated by its Z mag).

Outputs: per-pointing fibers needed (SN/host) and spare fibers over the survey,
and the remaining (incomplete) Z<25 hosts at survey end -- as plots + tables.
"""
import os
import glob
import importlib.util as ilu
import warnings
import numpy as np
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun, get_body
from matplotlib.path import Path

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
DATADIR = "/Volumes/exdisk1/data/Roman/ltcv_sim"
FIELD = "ELAIS-N1"
MODELS = {"Ia": "SNIaMODEL00", "CC": "NONIaMODEL06"}

FOV = 1.3
R = FOV / 2.0
FIBERS_AVAIL = 1850          # science fibers/pointing (2394 - ~400 sky - ~150 fluxstd)
SN_ZCUT = 24.0
HOST_ZCUT = 25.0
SN_TARGET = 5.0              # target host S/N at 10-pix binning
NPIX_BIN = 10
SHIFT_DAYS = Time("2028-06-01").mjd - Time("2029-01-01").mjd   # -214
SURVEY_START = Time("2028-06-01").mjd
SURVEY_END = Time("2030-05-31").mjd
OBS_MONTHS = {2, 3, 4, 5, 6, 7, 8}     # ELAIS-N1 visible >2h/night from Subaru
ELAIS = (242.5, 54.5)

# ---------------------------------------------------------------- ETC helper
_pe_spec = ilu.spec_from_file_location("pfs_etc", os.path.join(HERE, "01_pfs_etc.py"))
pe = ilu.module_from_spec(_pe_spec); _pe_spec.loader.exec_module(pe)


def build_hours_to_target():
    """Interpolator: AB mag -> hours (1h visits) to reach S/N=SN_TARGET at 10-pix."""
    mags = np.arange(21.0, 26.01, 0.25)
    sn10 = []
    for m in mags:
        res = pe.run_etc(mag=float(m), exp_time=450, exp_num=8)   # 1 hour
        sn10.append(pe.band_sn(res, "i") * np.sqrt(NPIX_BIN))
    sn10 = np.array(sn10)
    hrs = (SN_TARGET / sn10) ** 2
    return lambda mag: np.interp(mag, mags, hrs)


# ---------------------------------------------------------------- geometry
def hexagon(xc, yc, rr, phase):
    a = np.radians(np.arange(0, 360, 60) + phase)
    return np.column_stack([xc + rr * np.cos(a), yc + rr * np.sin(a)])


def flower_centers(rr, rot):
    sp = np.sqrt(3) * rr
    c = [(0.0, 0.0)]
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang + rot)
        c.append((sp * np.cos(a), sp * np.sin(a)))
    return c


def assign_pointing(xi, eta, rot):
    """Return pointing index 0..6 for each (xi,eta), or -1 if uncovered."""
    out = np.full(len(xi), -1, int)
    pts = np.column_stack([xi, eta])
    for i, (xc, yc) in enumerate(flower_centers(R, rot)):
        inside = Path(hexagon(xc, yc, R, rot)).contains_points(pts)
        out[(out < 0) & inside] = i
    return out


# ---------------------------------------------------------------- data
def field_files(model, kind):
    fs = []
    for ext in ("FITS", "FITS.gz"):
        fs += glob.glob(f"{DATADIR}/{FIELD}/*{model}-*_{kind}.{ext}")
    return sorted(f for f in fs if not os.path.basename(f).startswith("._"))


def load_program_sne():
    """Return dict arrays for SNe whose Z087 light curve reaches <24 (program SNe):
    type, ra,dec (SN), hra,hdec (host), hostZ, t1,t2 (active window, shifted MJD)."""
    rec = dict(typ=[], ra=[], dec=[], hra=[], hdec=[], hostz=[], t1=[], t2=[])
    for typ, model in MODELS.items():
        for hf in field_files(model, "HEAD"):
            pf = hf.replace("HEAD", "PHOT")
            with fits.open(hf) as fh, fits.open(pf) as fp:
                H = fh[1].data; P = fp[1].data
                if H is None or len(H) == 0:
                    continue
                band = np.asarray(P["BAND"]).astype(str)
                mjd = np.asarray(P["MJD"], float)
                smag = np.asarray(P["SIM_MAGOBS"], float)
                pmin = np.asarray(H["PTROBS_MIN"], int)
                pmax = np.asarray(H["PTROBS_MAX"], int)
                for k in range(len(H)):
                    a, b = pmin[k], pmax[k]
                    sl = slice(a - 1, b)
                    zb = band[sl] == "Z087-Z"
                    zm = smag[sl][zb]; zt = mjd[sl][zb]
                    good = (zm > 0) & (zm < 90) & (zm < SN_ZCUT)
                    if not good.any():
                        continue
                    t = zt[good]
                    rec["typ"].append(typ)
                    rec["ra"].append(float(H["RA"][k]))
                    rec["dec"].append(float(H["DEC"][k]))
                    rec["hra"].append(float(H["HOSTGAL_RA"][k]))
                    rec["hdec"].append(float(H["HOSTGAL_DEC"][k]))
                    rec["hostz"].append(float(H["HOSTGAL_MAG_Z"][k]))
                    rec["t1"].append(t.min() + SHIFT_DAYS)
                    rec["t2"].append(t.max() + SHIFT_DAYS)
    for k in rec:
        rec[k] = np.array(rec[k]) if k != "typ" else np.array(rec[k], dtype=object)
    return rec


# ---------------------------------------------------------------- schedule
def new_moons(mjd0, mjd1):
    days = np.arange(int(mjd0), int(mjd1) + 1, 1.0)
    t = Time(days, format="mjd")
    sep = get_sun(t).separation(get_body("moon", t)).deg
    illum = (1 - np.cos(np.radians(sep))) / 2
    nm = days[1:-1][(illum[1:-1] < illum[:-2]) & (illum[1:-1] < illum[2:]) & (illum[1:-1] < 0.05)]
    return nm


def visit_schedule():
    """List of (mjd, config) for A/B visits in observable months, 2028-06..2030-05."""
    nm = new_moons(SURVEY_START - 40, SURVEY_END + 40)
    visits = []
    for m in nm:
        dt = Time(m, format="mjd").datetime
        if dt.month not in OBS_MONTHS:
            continue
        a, b = m - 6.0, m + 6.0          # start / end of the ~2-week dark run
        for mjd, cfg in ((a, "A"), (b, "B")):
            if SURVEY_START <= mjd <= SURVEY_END:
                visits.append((mjd, cfg))
    return sorted(visits)


# ---------------------------------------------------------------- simulation
def main():
    os.makedirs(PNG_DIR, exist_ok=True); os.makedirs(CSV_DIR, exist_ok=True)
    print("Building ETC host-integration interpolator ...")
    hours_to_target = build_hours_to_target()

    print("Loading program SNe (Z087 light curve < 24) ...")
    sn = load_program_sne()
    n = len(sn["ra"])
    print(f"  program SNe: {n}  (Ia={np.sum(sn['typ']=='Ia')}, CC={np.sum(sn['typ']=='CC')})")

    ra0, dec0 = 242.498, 54.497
    cosd = np.cos(np.radians(dec0))
    # SN and host pointing assignment for configs A(0) and B(30)
    sxi, seta = (sn["ra"] - ra0) * cosd, sn["dec"] - dec0
    hxi, heta = (sn["hra"] - ra0) * cosd, sn["hdec"] - dec0
    pA_sn = assign_pointing(sxi, seta, 0.0)
    pB_sn = assign_pointing(sxi, seta, 30.0)
    pA_ho = assign_pointing(hxi, heta, 0.0)
    pB_ho = assign_pointing(hxi, heta, 30.0)

    host_ok = (sn["hostz"] > 0) & (sn["hostz"] < HOST_ZCUT)
    visits_needed = np.ceil(hours_to_target(np.clip(sn["hostz"], 21, 26))).astype(int)
    visits_needed[~host_ok] = 10**9          # non-targets never "needed"

    visits = visit_schedule()
    print(f"  visits: {len(visits)}  "
          f"({Time(visits[0][0],format='mjd').iso[:10]} .. {Time(visits[-1][0],format='mjd').iso[:10]})")

    host_obs = np.zeros(n, int)      # 1-h visits accumulated on each host
    host_done = np.zeros(n, bool)
    # per-visit per-pointing counts
    nv = len(visits)
    sn_fib = np.zeros((nv, 7), int)
    ho_fib = np.zeros((nv, 7), int)
    vmjd = np.array([v[0] for v in visits])
    vcfg = [v[1] for v in visits]

    for vi, (mjd, cfg) in enumerate(visits):
        psn = pA_sn if cfg == "A" else pB_sn
        pho = pA_ho if cfg == "A" else pB_ho
        active_sn = (mjd >= sn["t1"]) & (mjd <= sn["t2"]) & (psn >= 0)
        # host active: SN already faded (past t2), host qualifies, not done, covered
        faded = mjd > sn["t2"]
        active_ho = faded & host_ok & (~host_done) & (pho >= 0)
        for p in range(7):
            sn_fib[vi, p] = np.sum(active_sn & (psn == p))
            ho_fib[vi, p] = np.sum(active_ho & (pho == p))
        # accumulate host integration for those observed this visit
        obs_now = np.where(active_ho)[0]
        host_obs[obs_now] += 1
        host_done[host_obs >= visits_needed] = True

    tot_fib = sn_fib + ho_fib
    spare = FIBERS_AVAIL - tot_fib

    # ---- remaining hosts (Z<25) not completed by survey end ----
    started = host_obs > 0
    completed = host_done
    # a host is "remaining" if it qualifies (Z<25) and not completed by survey end
    remaining = host_ok & (~completed)
    print(f"\nHosts (Z<25) of program SNe: {int(host_ok.sum())}")
    print(f"  completed (reached S/N=5): {int(completed.sum())}")
    print(f"  remaining (incomplete) at 2030-05-31: {int(remaining.sum())}")
    print(f"    of which started but unfinished: {int((remaining & started).sum())}")
    print(f"    not yet started (SN still bright / faded late): {int((remaining & ~started).sum())}")

    # ---- per-pointing table ----
    csv = os.path.join(CSV_DIR, "07_fiber_budget_ELAIS-N1.csv")
    with open(csv, "w") as fo:
        fo.write("pointing,peak_SN,peak_host,peak_total,mean_total,min_spare,mean_spare\n")
        for p in range(7):
            fo.write(f"P{p},{sn_fib[:,p].max()},{ho_fib[:,p].max()},{tot_fib[:,p].max()},"
                     f"{tot_fib[:,p].mean():.1f},{spare[:,p].min()},{spare[:,p].mean():.1f}\n")
    print("\n  pointing  peakSN  peakHost  peakTot  meanTot  minSpare  meanSpare")
    for p in range(7):
        print(f"  P{p}       {sn_fib[:,p].max():5d}   {ho_fib[:,p].max():6d}   "
              f"{tot_fib[:,p].max():6d}   {tot_fib[:,p].mean():6.1f}  {spare[:,p].min():7d}  {spare[:,p].mean():8.1f}")
    print("table ->", csv)

    # =================== plots ===================
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib import font_manager as fm
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    LFS, TFS, TT = 15, 12, 15
    dn = mdates.date2num(Time(vmjd, format="mjd").to_datetime())

    def gapped(x, y, gap=60.0):          # break lines across the dead seasons
        xs, ys = [x[0]], [y[0]]
        for i in range(1, len(x)):
            if vmjd[i] - vmjd[i - 1] > gap:
                xs.append(np.nan); ys.append(np.nan)
            xs.append(x[i]); ys.append(y[i])
        return np.array(xs), np.array(ys)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9), sharex=True,
                                   constrained_layout=True)
    # panel 1: central pointing P0 SN+host stacked, and ring mean
    ax1.bar(dn, sn_fib[:, 0], width=5, color="#1f77b4", label="P0 SNe")
    ax1.bar(dn, ho_fib[:, 0], width=5, bottom=sn_fib[:, 0], color="#ff7f0e", label="P0 hosts")
    ring = tot_fib[:, 1:].mean(axis=1)
    ax1.plot(*gapped(dn, ring), "k.-", lw=1, ms=4, label="ring mean (P1-P6) total")
    ax1.set_ylabel("fibers needed / pointing", fontsize=LFS)
    ax1.tick_params(labelsize=TFS)
    ax1.legend(fontsize=10, loc="upper left")
    ax1.set_title("PFS fiber demand in ELAIS-N1 (central pointing P0)  "
                  f"[{FIBERS_AVAIL} science fibers available]", fontsize=TT)
    ax1.grid(True, alpha=0.3)

    # panel 2: spare fibers (min over pointings = worst case, and P0)
    ax2.plot(*gapped(dn, spare.min(axis=1)), "-o", color="#d62728", ms=4,
             label="min spare (worst pointing)")
    ax2.plot(*gapped(dn, spare[:, 0].astype(float)), "-s", color="#2ca02c", ms=3, label="P0 spare")
    ax2.axhline(FIBERS_AVAIL, color="0.6", ls=":", lw=1)
    ax2.set_ylabel("spare fibers / pointing", fontsize=LFS)
    ax2.set_xlabel("date", fontsize=LFS)
    ax2.tick_params(labelsize=TFS)
    ax2.legend(fontsize=10, loc="lower left")
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_locator(mdates.MonthLocator((1, 4, 7, 10)))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    for lab in ax2.get_xticklabels():
        lab.set_rotation(30); lab.set_ha("right")

    png = os.path.join(PNG_DIR, "07_fiber_budget_ELAIS-N1.png")
    fig.savefig(png, dpi=140); print("plot ->", png)

    # remaining hosts histogram by Z mag
    fig2, axh = plt.subplots(figsize=(9, 6))
    hz = sn["hostz"][host_ok]
    axh.hist(hz, bins=np.arange(20, 25.01, 0.25), color="0.7", label="all Z<25 hosts")
    axh.hist(sn["hostz"][remaining], bins=np.arange(20, 25.01, 0.25),
             color="#d62728", alpha=0.8, label="remaining (incomplete) at 2030-05-31")
    axh.set_xlabel("host Roman $Z$ (AB mag)", fontsize=LFS)
    axh.set_ylabel("number of hosts", fontsize=LFS)
    axh.set_title("ELAIS-N1 SN-host galaxies: completed vs remaining", fontsize=TT)
    axh.tick_params(labelsize=TFS); axh.legend(fontsize=11); axh.grid(True, alpha=0.3)
    png2 = os.path.join(PNG_DIR, "07_remaining_hosts_ELAIS-N1.png")
    fig2.savefig(png2, dpi=140, bbox_inches="tight"); print("plot ->", png2)


if __name__ == "__main__":
    main()
