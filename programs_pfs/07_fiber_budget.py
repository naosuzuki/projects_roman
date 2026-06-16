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
MODELS = {"Ia": "SNIaMODEL00", "CC": "NONIaMODEL06", "TDE": "NONIaMODEL02"}

FOV = 1.3
R = FOV / 2.0
FIBERS_TOTAL = 2394          # PFS science fibers per pointing
SKY = 400                    # reserved sky fibers (observatory-provided)
FLUXSTD = 200                # reserved flux-standard fibers
FIBERS_AVAIL = FIBERS_TOTAL - SKY - FLUXSTD   # = 1794 for program targets
SN_ZCUT = 24.0
HOST_ZCUT = 25.5
SN_TARGET = 5.0              # target host S/N at 10-pix binning
NPIX_BIN = 10
WEATHER = 0.70               # fraction of usable (clear) nights at Maunakea
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
    ntotal = {t: 0 for t in MODELS}        # all simulated Roman transients per class
    for typ, model in MODELS.items():
        for hf in field_files(model, "HEAD"):
            pf = hf.replace("HEAD", "PHOT")
            with fits.open(hf) as fh, fits.open(pf) as fp:
                H = fh[1].data; P = fp[1].data
                if H is None or len(H) == 0:
                    continue
                ntotal[typ] += len(H)
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
    return rec, ntotal


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
    sn, ntotal = load_program_sne()
    n = len(sn["ra"])
    CLASSES = ["Ia", "CC", "TDE"]
    print("  program SNe: %d  (%s)" % (n, ", ".join(
        f"{c}={int(np.sum(sn['typ']==c))}" for c in CLASSES)))

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

    typ = sn["typ"]
    sn_c = {c: np.zeros((nv, 7), int) for c in CLASSES}   # live-SN fibers per class
    ho_c = {c: np.zeros((nv, 7), int) for c in CLASSES}   # host fibers per class
    sn_obs_A = np.zeros(n, bool)     # observed live (Z<24, covered) on >=1 A visit
    sn_obs_B = np.zeros(n, bool)     # observed live on >=1 B visit
    for vi, (mjd, cfg) in enumerate(visits):
        psn = pA_sn if cfg == "A" else pB_sn
        pho = pA_ho if cfg == "A" else pB_ho
        active_sn = (mjd >= sn["t1"]) & (mjd <= sn["t2"]) & (psn >= 0)
        if cfg == "A":
            sn_obs_A |= active_sn
        else:
            sn_obs_B |= active_sn
        # host active: SN already faded (past t2), host qualifies, not done, covered
        faded = mjd > sn["t2"]
        active_ho = faded & host_ok & (~host_done) & (pho >= 0)
        for p in range(7):
            asn = active_sn & (psn == p)
            aho = active_ho & (pho == p)
            for c in CLASSES:
                cm = (typ == c)
                sn_c[c][vi, p] = np.sum(asn & cm)
                ho_c[c][vi, p] = np.sum(aho & cm)
        # accumulate host integration for those observed this visit
        obs_now = np.where(active_ho)[0]
        host_obs[obs_now] += 1
        host_done[host_obs >= visits_needed] = True

    sn_fib = sum(sn_c[c] for c in CLASSES)
    ho_fib = sum(ho_c[c] for c in CLASSES)
    tot_fib = sn_fib + ho_fib
    spare = FIBERS_AVAIL - tot_fib

    # ---- remaining hosts (Z<25) not completed by survey end ----
    started = host_obs > 0
    completed = host_done
    # a host is "remaining" if it qualifies (Z<25) and not completed by survey end
    remaining = host_ok & (~completed)
    print(f"\nfiber budget/pointing: {FIBERS_TOTAL} total = {SKY} sky + {FLUXSTD} FluxStd "
          f"+ demand + spare  ({FIBERS_AVAIL} available for program)")
    print(f"Hosts (Z<{HOST_ZCUT}) of program SNe: {int(host_ok.sum())}")
    print(f"  completed (reached S/N=5): {int(completed.sum())}")
    print(f"  remaining (incomplete) at 2030-05-31: {int(remaining.sum())}")
    print(f"    of which started but unfinished: {int((remaining & started).sum())}")
    print(f"    not yet started (SN still bright / faded late): {int((remaining & ~started).sum())}")

    # ---- successful PFS spec-z of live transients (observed >=1x while Z<24) ----
    sn_observed = sn_obs_A | sn_obs_B
    # "visible from Subaru" = the Z<24 window overlaps an observable (Feb-Aug)
    # night within the survey -> the correct denominator for the success rate.
    surv_days = np.arange(int(SURVEY_START), int(SURVEY_END) + 1)
    surv_mon = np.array([Time(int(d), format="mjd").datetime.month for d in surv_days])
    obs_mjds = surv_days[np.isin(surv_mon, list(OBS_MONTHS))]
    visible = (np.searchsorted(obs_mjds, sn["t2"], "right")
               > np.searchsorted(obs_mjds, sn["t1"], "left"))

    ycsv = os.path.join(CSV_DIR, "07_specz_yield_ELAIS-N1.csv")
    print("\nPFS spec-z yield (live transient observed >=1 visit while Z<24, covered):")
    print(f"  class  Roman  reachZ<24  visible  observed  success  net(x{WEATHER:.2f} wx)")
    with open(ycsv, "w") as fo:
        fo.write("class,N_roman,N_program,N_visible,observed_A,observed_B,observed,"
                 "success_rate,net_success\n")
        for c in CLASSES:
            cm = (typ == c)
            nrom = int(ntotal[c])
            npr, nvis = int(cm.sum()), int((visible & cm).sum())
            oa, ob, ou = (int((sn_obs_A & cm).sum()), int((sn_obs_B & cm).sum()),
                          int((sn_observed & cm).sum()))
            sr = ou / nvis if nvis else 0.0
            net = sr * WEATHER
            fo.write(f"{c},{nrom},{npr},{nvis},{oa},{ob},{ou},{sr:.4f},{net:.4f}\n")
            print(f"  {c:3s}  {nrom:6d}  {npr:7d}  {nvis:6d}  {ou:6d}   "
                  f"{100*sr:5.1f}%   {100*net:5.1f}%")
    print("specz yield ->", ycsv)

    # ---- reproducible CSVs ----
    # (1) program transient catalog (one row per program SN)
    pcsv = os.path.join(CSV_DIR, "07_program_sne_ELAIS-N1.csv")
    with open(pcsv, "w") as fo:
        fo.write("id,type,ra,dec,host_ra,host_dec,host_Z,t1_mjd,t2_mjd,"
                 "pA_sn,pB_sn,pA_host,pB_host,host_target,host_visits_needed,"
                 "observed_A,observed_B,observed\n")
        for k in range(n):
            vn = int(visits_needed[k]) if host_ok[k] else -1
            fo.write(f"{k},{typ[k]},{sn['ra'][k]:.5f},{sn['dec'][k]:.5f},"
                     f"{sn['hra'][k]:.5f},{sn['hdec'][k]:.5f},{sn['hostz'][k]:.3f},"
                     f"{sn['t1'][k]:.2f},{sn['t2'][k]:.2f},"
                     f"{pA_sn[k]},{pB_sn[k]},{pA_ho[k]},{pB_ho[k]},"
                     f"{int(host_ok[k])},{vn},"
                     f"{int(sn_obs_A[k])},{int(sn_obs_B[k])},{int(sn_observed[k])}\n")
    print("program SNe ->", pcsv)

    # (2) per-visit per-pointing demand (one row per visit x pointing)
    vcsv = os.path.join(CSV_DIR, "07_visit_demand_ELAIS-N1.csv")
    with open(vcsv, "w") as fo:
        fo.write("visit,mjd,date,config,pointing,sn_Ia,sn_CC,sn_TDE,"
                 "host_Ia,host_CC,host_TDE,sn_total,host_total,total,spare\n")
        for vi in range(nv):
            for p in range(7):
                st, ht = int(sn_fib[vi, p]), int(ho_fib[vi, p])
                fo.write(f"{vi},{vmjd[vi]:.2f},{Time(vmjd[vi],format='mjd').iso[:10]},"
                         f"{vcfg[vi]},P{p},"
                         f"{sn_c['Ia'][vi,p]},{sn_c['CC'][vi,p]},{sn_c['TDE'][vi,p]},"
                         f"{ho_c['Ia'][vi,p]},{ho_c['CC'][vi,p]},{ho_c['TDE'][vi,p]},"
                         f"{st},{ht},{st+ht},{FIBERS_AVAIL-st-ht}\n")
    print("visit demand ->", vcsv)

    # (3) host completion ledger (one row per qualifying host)
    hcsv = os.path.join(CSV_DIR, "07_host_completion_ELAIS-N1.csv")
    with open(hcsv, "w") as fo:
        fo.write("id,type,host_Z,visits_needed,visits_done,started,completed\n")
        for k in np.where(host_ok)[0]:
            fo.write(f"{k},{typ[k]},{sn['hostz'][k]:.3f},{int(visits_needed[k])},"
                     f"{int(host_obs[k])},{int(host_obs[k]>0)},{int(host_done[k])}\n")
    print("host completion ->", hcsv)

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
    # panel 1: target-fiber demand composition for the central pointing P0
    # (sky/FluxStd are subtracted from the 1794 available, not shown)
    w = 9.0
    b = np.zeros(nv, float)
    for arr, col, lab in ((sn_c["Ia"][:, 0],  "#d62728", "SN Ia"),
                          (sn_c["CC"][:, 0],  "#1f77b4", "CC SN"),
                          (sn_c["TDE"][:, 0], "#6a3d9a", "TDE"),
                          (ho_c["Ia"][:, 0],  "#ff7f0e", "Ia host"),
                          (ho_c["CC"][:, 0],  "#9ecae1", "CC host"),
                          (ho_c["TDE"][:, 0], "#cab2d6", "TDE host")):
        ax1.bar(dn, arr, width=w, bottom=b, color=col, label=lab)
        b = b + arr
    ax1.set_ylabel("Fibers / pointing", fontsize=18)
    ax1.tick_params(labelsize=TFS)
    ax1.legend(fontsize=12, loc="upper right", ncol=2)
    ax1.set_title("PFS target-fiber demand in ELAIS-N1 (central pointing P0)  "
                  f"[{FIBERS_AVAIL} available $= 2394-400$ sky$-200$ FluxStd]", fontsize=TT)
    ax1.grid(True, axis="y", alpha=0.3)

    # panel 2: spare fibers (min over pointings = worst case, and P0)
    ax2.plot(*gapped(dn, spare.max(axis=1).astype(float)), "-o", color="#d62728", ms=4,
             label="Largest spare (of 7 pointings)")
    ax2.plot(*gapped(dn, spare.min(axis=1).astype(float)), "-o", color="#1f77b4", ms=4,
             label="Smallest spare (of 7 pointings)")
    ax2.plot(*gapped(dn, spare[:, 0].astype(float)), "-s", color="#2ca02c", ms=3,
             label="P0 spare")
    ax2.axhline(FIBERS_AVAIL, color="0.6", ls=":", lw=1, label=f"Available = {FIBERS_AVAIL}")
    ax2.set_ylabel("Spare fibers / pointing", fontsize=18)
    ax2.set_xlabel("Date", fontsize=18)
    ax2.tick_params(labelsize=TFS)
    ax2.set_ylim(spare.min() - 25, FIBERS_AVAIL + 55)     # headroom for upper-right legend
    ax2.legend(fontsize=12, loc="upper right", ncol=2)
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_locator(mdates.MonthLocator((1, 4, 7, 10)))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    for lab in ax2.get_xticklabels():
        lab.set_rotation(30); lab.set_ha("right")

    png = os.path.join(PNG_DIR, "07_fiber_budget_ELAIS-N1.png")
    fig.savefig(png, dpi=140); print("plot ->", png)

    # remaining hosts histogram by Z mag
    fig2, axh = plt.subplots(figsize=(9, 6))
    hbins = np.arange(20, 25.51, 0.25)
    hz = sn["hostz"][host_ok]
    axh.hist(hz, bins=hbins, color="0.7", label=f"all $Z<{HOST_ZCUT}$ hosts")
    axh.hist(sn["hostz"][remaining], bins=hbins,
             color="#d62728", alpha=0.8, label="remaining (incomplete) at 2030-05-31")
    axh.set_xlabel("host Roman $Z$ (AB mag)", fontsize=LFS)
    axh.set_ylabel("number of hosts", fontsize=LFS)
    axh.set_title("ELAIS-N1 SN-host galaxies: completed vs remaining", fontsize=TT)
    axh.tick_params(labelsize=TFS); axh.legend(fontsize=11); axh.grid(True, alpha=0.3)
    png2 = os.path.join(PNG_DIR, "07_remaining_hosts_ELAIS-N1.png")
    fig2.savefig(png2, dpi=140, bbox_inches="tight"); print("plot ->", png2)


if __name__ == "__main__":
    main()
