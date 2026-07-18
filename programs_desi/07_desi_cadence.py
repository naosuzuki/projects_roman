#!/usr/bin/env python3
"""
07_desi_cadence.py -- DESI cadence / fiber-budget simulation for the ELAIS-N1
live-SN + host program, using Configuration B (6-pointing hexagon, d=0.40 deg).

This is the DESI counterpart of programs_pfs/07_fiber_budget.py, adapted to a
DESI observing model in which ONE pointing is observed per visit night, cycling
through the six Configuration-B pointings.  Two scenarios:

  * Scenario 1 (--cadence 5):  one pointing every  5 days -> the 6-pointing
    circle is completed every 30 days (each pointing revisited every 30 d).
  * Scenario 2 (--cadence 10): one pointing every 10 days -> circle every 60 d
    (each pointing revisited every 60 d).

Model / assumptions
  * Each pointing-visit integrates a FIXED 2.0 h, delivered in parallel by DESI's
    5000 fibers to every target in that field (live SNe get their epoch; hosts
    accumulate integration).
  * A live transient gets a spec-z if >=1 visit to its pointing falls within its
    Z087 < 24 active window [t1, t2] (an observable KPNO night).
  * A host completes (spec-z) when its accumulated integration reaches the DESI
    requirement for S/N=5 at 10-pix = APERTURE x (PFS hours), APERTURE=4.2 for
    the 4-m vs 8.2-m aperture (host i-band continuum ~ host Z, as in PFS).
  * Visits happen only on KPNO-observable nights (ELAIS-N1 visible > 2 h/night).
  * Weather: net yields tabulated at KPNO clear-night fractions 0.75 and 0.65.

All inputs come from the committed catalog
../programs_pfs/outputs/csv/07_program_sne_ELAIS-N1.csv (t1/t2, z, host_Z,
host_visits_needed = PFS hours to S/N=5), so no external disk is needed.

    python 07_desi_cadence.py --cadence 5
    python 07_desi_cadence.py --cadence 10

Outputs (suffix _cad{N}): a DESI program catalog CSV (observed / host_completed),
spec-z and host yield CSVs, the fiber-budget figure (Fig 18 analog), and the four
redshift histograms (Figs 14-17 analogs).
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
PROG_CSV = os.path.join(HERE, "..", "programs_pfs", "outputs", "csv",
                        "07_program_sne_ELAIS-N1.csv")

RA0, DEC0 = 242.498, 54.497
FOV = 3.2
R = FOV / 2.0
D_HEX = 0.40                       # Config B ring radius (deg)
FIBERS = 5000                      # DESI fibers per pointing
APERTURE = 4.2                     # DESI(4m) / PFS(8.2m) integration-time penalty
T_VISIT = 2.0                      # fixed hours of integration per pointing-visit
HOST_ZCUT = 25.5
SURVEY_START = Time("2028-06-01").mjd
SURVEY_END = Time("2030-05-31").mjd
WEATHERS = (0.75, 0.65)            # KPNO clear-night fractions
# Full Roman-simulated ELAIS-N1 totals per class (field property; from the
# HOURGLASS2 sim, same as the Roman-Subaru study Table 3).
ROMAN_ALL = {"Ia": 13094, "CC": 18774, "TDE": 21}
CLASSES = ["Ia", "CC", "TDE"]


def kpno_site():
    try:
        return EarthLocation.of_site("kpno")
    except Exception:
        return EarthLocation(lat=31.9633 * u.deg, lon=-111.5997 * u.deg, height=2120 * u.m)


def hexagon_centers(d=D_HEX, phase_deg=0.0):
    a = np.radians(phase_deg + 60.0 * np.arange(6))
    return np.column_stack([d * np.cos(a), d * np.sin(a)])


def assign_pointing(xi, eta, centers):
    """Nearest containing-circle pointing index (0..5) for each point, else -1."""
    out = np.full(len(xi), -1, int)
    best = np.full(len(xi), 1e9)
    for i, (xc, yc) in enumerate(centers):
        dd = np.hypot(xi - xc, eta - yc)
        inside = dd <= R
        take = inside & (dd < best)
        out[take] = i
        best[take] = dd[take]
    return out


def nightly_visible_hours(days, target, site, elev=30.0, twilight=-18.0, dt=0.5):
    utc_h = np.arange(1.0, 16.0 + 1e-6, dt)
    base = Time(np.asarray(days, float), format="mjd")
    times = (base[:, None] + 1 * u.day) + utc_h[None, :] * u.hour
    flat = times.ravel()
    aa = AltAz(obstime=flat, location=site)
    talt = target.transform_to(aa).alt.deg.reshape(times.shape)
    salt = get_sun(flat).transform_to(aa).alt.deg.reshape(times.shape)
    return ((salt < twilight) & (talt > elev)).sum(axis=1) * dt


def setup_fonts():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager as fm
    import matplotlib.pyplot as plt
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    return plt


def weather_hist(plt, z_obs, title, ylabel, png, cad):
    """Nested weather-band redshift histogram (mirror PFS 13/14/16/17)."""
    cases = [(1.0, "100% (no weather loss)", "blue"),
             (WEATHERS[0], f"{int(100*WEATHERS[0])}% weather", "g"),
             (WEATHERS[1], f"{int(100*WEATHERS[1])}% weather", "red")]
    if len(z_obs) == 0:
        zmax = 1.0
    else:
        zmax = np.ceil(z_obs.max() * 10) / 10
    bins = np.arange(0.0, zmax + 0.051, 0.05)
    counts, edges = np.histogram(z_obs, bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = (edges[1] - edges[0]) * 0.93
    fig, ax = plt.subplots(figsize=(9, 6))
    for i, (w, lab, col) in enumerate(cases):
        ax.bar(centers, counts * w, width=bw, color=col, edgecolor="white",
               linewidth=0.4, zorder=i + 1, label=f"{lab}  ($N={w*len(z_obs):.0f}$)")
    ax.set_xlabel("Redshift", fontsize=18)
    ax.set_ylabel(ylabel, fontsize=17)
    ax.set_title(title, fontsize=14)
    ax.tick_params(labelsize=13)
    ax.set_xlim(0, bins[-1] if len(bins) > 1 else 1.0)
    ax.legend(fontsize=12, loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cadence", type=int, default=5, help="days between pointing visits (5 or 10)")
    ap.add_argument("--t-visit", type=float, default=T_VISIT, dest="tvisit",
                    help="fixed integration hours per pointing-visit (default 2.0)")
    ap.add_argument("--aperture", type=float, default=APERTURE,
                    help="DESI/PFS integration-time penalty (default 4.2)")
    args = ap.parse_args()
    cad = args.cadence
    tag = f"cad{cad}"
    os.makedirs(PNG_DIR, exist_ok=True); os.makedirs(CSV_DIR, exist_ok=True)

    # ---- load committed program-SN catalog ----
    d = np.genfromtxt(PROG_CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    typ = np.asarray(d["type"]).astype(str)
    z = d["z"].astype(float)
    ra, dec = d["ra"].astype(float), d["dec"].astype(float)
    hra, hdec = d["host_ra"].astype(float), d["host_dec"].astype(float)
    hostz = d["host_Z"].astype(float)
    t1, t2 = d["t1_mjd"].astype(float), d["t2_mjd"].astype(float)
    host_target = d["host_target"].astype(int)
    vneed_pfs = d["host_visits_needed"].astype(float)   # PFS 1-h visits to S/N=5
    n = len(z)
    cosd = np.cos(np.radians(DEC0))

    # ---- Config B pointing assignment (nearest containing circle) ----
    cen = hexagon_centers()
    psn = assign_pointing((ra - RA0) * cosd, dec - DEC0, cen)        # SN -> pointing
    pho = assign_pointing((hra - RA0) * cosd, hdec - DEC0, cen)      # host -> pointing

    host_ok = (host_target == 1) & np.isfinite(hostz)
    needed_desi = args.aperture * np.where(host_ok, vneed_pfs, 1e18)  # DESI hours to S/N=5

    # ---- KPNO observable nights over the survey ----
    site, target = kpno_site(), SkyCoord(RA0 * u.deg, DEC0 * u.deg)
    days = np.arange(int(SURVEY_START), int(SURVEY_END) + 1)
    vis_h = nightly_visible_hours(days, target, site)
    observable = vis_h > 2.0
    obs_mjds = days[observable]

    # ---- build the visit schedule: 1 pointing per cadence-step, cycling 6 ----
    visits = []          # (mjd, pointing)
    cyc = 0
    for day in range(int(SURVEY_START), int(SURVEY_END) + 1, cad):
        if observable[day - int(SURVEY_START)]:
            visits.append((day, cyc % 6))
            cyc += 1
    nv = len(visits)
    print(f"=== DESI Config B, cadence {cad} d  (visit every {cad} d, circle every {6*cad} d) ===")
    print(f"  observable nights: {int(observable.sum())}/{len(days)};  "
          f"scheduled visits: {nv}  ({Time(visits[0][0],format='mjd').iso[:10]} .. "
          f"{Time(visits[-1][0],format='mjd').iso[:10]})")
    pv = np.bincount([p for _, p in visits], minlength=6)
    print(f"  visits per pointing P1..P6: {list(pv)}  (revisit ~ every {6*cad} d in season)")

    # ---- simulate ----
    observed = np.zeros(n, bool)          # live transient caught (Z<24, in pointing)
    host_acc = np.zeros(n, float)         # accumulated DESI integration hours
    host_done = np.zeros(n, bool)
    dem = []                               # per-visit (mjd, p, n_live, n_host)
    for (mjd, p) in visits:
        live = (psn == p) & (t1 <= mjd) & (mjd <= t2)
        observed |= live
        hint = (pho == p) & host_ok & (mjd > t2) & (~host_done)
        host_acc[hint] += args.tvisit
        host_done[host_ok & (host_acc >= needed_desi)] = True
        dem.append((mjd, p, int(live.sum()), int(hint.sum())))
    dem = np.array(dem, float)

    # ---- "visible" denominators (KPNO observable within active window) ----
    visible = (np.searchsorted(obs_mjds, t2, "right")
               > np.searchsorted(obs_mjds, t1, "left")) & (psn >= 0)
    host_visible = host_ok & (np.searchsorted(obs_mjds, t2, "right") < len(obs_mjds)) & (pho >= 0)

    # ---- DESI program catalog (mirror PFS columns used by histograms) ----
    cat = os.path.join(CSV_DIR, f"07_desi_program_sne_{tag}.csv")
    with open(cat, "w") as fo:
        fo.write("id,type,z,host_Z,p_sn,p_host,host_target,host_visits_needed_pfs,"
                 "host_hours_needed_desi,host_hours_done,observed,host_completed\n")
        for k in range(n):
            fo.write(f"{k},{typ[k]},{z[k]:.5f},{hostz[k]:.3f},{psn[k]},{pho[k]},"
                     f"{host_target[k]},{vneed_pfs[k]:.0f},"
                     f"{needed_desi[k] if host_ok[k] else -1:.1f},{host_acc[k]:.1f},"
                     f"{int(observed[k])},{int(host_done[k])}\n")
    print("catalog ->", cat)

    # ---- spec-z yield (Table 3 data) ----
    def yield_rows(num_mask, den_mask, observed_mask):
        rows = []
        for c in CLASSES:
            cm = (typ == c)
            nrom = ROMAN_ALL[c]
            npr = int((num_mask & cm).sum())
            nvis = int((den_mask & cm).sum())
            nob = int((observed_mask & cm).sum())
            sr = nob / nvis if nvis else 0.0
            nets = [(int(round(nob * w)), sr * w) for w in WEATHERS]
            rows.append((c, nrom, npr, nvis, nob, sr, nets))
        return rows

    program = (psn >= 0) | True   # all program SNe reach Z<24 by construction
    srows = yield_rows(np.ones(n, bool), visible, observed)
    ycsv = os.path.join(CSV_DIR, f"07_desi_specz_yield_{tag}.csv")
    print("\nLive-transient spec-z yield:")
    print("  class  Roman  program  visible  observed  success  " +
          "  ".join(f"net(x{w:.2f})" for w in WEATHERS))
    with open(ycsv, "w") as fo:
        wcols = "".join(f",net_obs_w{int(100*w)},net_rate_w{int(100*w)}" for w in WEATHERS)
        fo.write("class,N_roman,N_program,N_visible,N_observed,success_rate" + wcols + "\n")
        for c, nrom, npr, nvis, nob, sr, nets in srows:
            fo.write(f"{c},{nrom},{npr},{nvis},{nob},{sr:.4f}" +
                     "".join(f",{no},{nr:.4f}" for no, nr in nets) + "\n")
            print(f"  {c:3s}  {nrom:6d}  {npr:6d}  {nvis:6d}  {nob:6d}   {100*sr:5.1f}%   " +
                  "  ".join(f"{no} ({100*nr:4.1f}%)" for no, nr in nets))
    print("specz yield ->", ycsv)

    # ---- host-galaxy yield ----
    hrows = yield_rows(host_ok, host_visible, host_done)
    hycsv = os.path.join(CSV_DIR, f"07_desi_host_yield_{tag}.csv")
    print("\nHost-galaxy spec-z yield (completed to S/N=5 at 10-pix):")
    print(f"  class  Roman  hostZ<{HOST_ZCUT}  visible  completed  success  " +
          "  ".join(f"net(x{w:.2f})" for w in WEATHERS))
    with open(hycsv, "w") as fo:
        wcols = "".join(f",net_obs_w{int(100*w)},net_rate_w{int(100*w)}" for w in WEATHERS)
        fo.write("class,N_roman,N_host,N_visible,N_completed,success_rate" + wcols + "\n")
        for c, nrom, nho, nvis, ncomp, sr, nets in hrows:
            fo.write(f"{c},{nrom},{nho},{nvis},{ncomp},{sr:.4f}" +
                     "".join(f",{no},{nr:.4f}" for no, nr in nets) + "\n")
            print(f"  {c:3s}  {nrom:6d}  {nho:7d}  {nvis:6d}  {ncomp:6d}   {100*sr:5.1f}%   " +
                  "  ".join(f"{no} ({100*nr:4.1f}%)" for no, nr in nets))
    print("host yield ->", hycsv)

    nhost = int(host_ok.sum())
    print(f"\nHosts (Z<{HOST_ZCUT}): {nhost}; completed {int(host_done.sum())} "
          f"({100*host_done.sum()/nhost:.0f}%); per-visit 2h, aperture x{args.aperture}")

    # ======================= plots =======================
    plt = setup_fonts()
    import matplotlib.dates as mdates

    # --- Fig 18 analog: fiber budget (per-visit demand vs 5000 fibers) ---
    vmjd = dem[:, 0]; nlive = dem[:, 2]; nhost = dem[:, 3]; ntot = nlive + nhost
    dn = mdates.date2num(Time(vmjd, format="mjd").to_datetime())
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(13, 9), sharex=True, constrained_layout=True)
    a1.bar(dn, nlive, width=cad * 0.8, color="#d62728", label="Live SNe (Phase II)")
    a1.bar(dn, nhost, width=cad * 0.8, bottom=nlive, color="#ff7f0e",
           label="Hosts integrating (Phase III)")
    a1.set_ylabel("Targets in visited pointing", fontsize=16)
    a1.legend(fontsize=12, loc="upper right")
    a1.set_title(f"DESI Config B fiber demand per visit, ELAIS-N1 — "
                 f"cadence {cad} d (1 pointing/{cad} d, circle/{6*cad} d)", fontsize=15)
    a1.grid(True, axis="y", alpha=0.3)
    a2.plot(dn, FIBERS - ntot, "-o", color="#1f77b4", ms=3, label="Spare fibers")
    a2.axhline(FIBERS, color="0.6", ls=":", lw=1, label=f"DESI fibers = {FIBERS}")
    a2.set_ylabel("Spare fibers", fontsize=16)
    a2.set_xlabel("Date", fontsize=16)
    a2.set_ylim(0, FIBERS + 200)
    a2.legend(fontsize=12, loc="lower right")
    a2.grid(True, alpha=0.3)
    a2.xaxis.set_major_locator(mdates.MonthLocator((1, 4, 7, 10)))
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    for lab in a2.get_xticklabels():
        lab.set_rotation(30); lab.set_ha("right")
    png = os.path.join(PNG_DIR, f"07_desi_fiber_budget_{tag}.png")
    fig.savefig(png, dpi=140); print("plot ->", png); plt.close(fig)

    # --- Figs 14-17 analogs: redshift histograms ---
    cl = f"(cadence {cad} d)"
    weather_hist(plt, z[(typ == "Ia") & observed],
                 f"ELAIS-N1 SN Ia with Successful DESI Spec-$z$ vs Redshift  {cl}",
                 "Number of SN Ia (Successful Spec-$z$)",
                 os.path.join(PNG_DIR, f"08_desi_snia_specz_zhist_{tag}.png"), cad)
    weather_hist(plt, z[(typ == "Ia") & host_done],
                 f"ELAIS-N1 SN Ia Host Galaxies with Successful DESI Spec-$z$ vs Redshift  {cl}",
                 "Number of SN Ia Hosts (Successful Spec-$z$)",
                 os.path.join(PNG_DIR, f"08_desi_snia_host_specz_zhist_{tag}.png"), cad)
    weather_hist(plt, z[(typ == "CC") & observed],
                 f"ELAIS-N1 Core-Collapse SNe with Successful DESI Spec-$z$ vs Redshift  {cl}",
                 "Number of CC SNe (Successful Spec-$z$)",
                 os.path.join(PNG_DIR, f"08_desi_cc_specz_zhist_{tag}.png"), cad)
    weather_hist(plt, z[(typ == "CC") & host_done],
                 f"ELAIS-N1 CC SN Host Galaxies with Successful DESI Spec-$z$ vs Redshift  {cl}",
                 "Number of CC SN Hosts (Successful Spec-$z$)",
                 os.path.join(PNG_DIR, f"08_desi_cc_host_specz_zhist_{tag}.png"), cad)


if __name__ == "__main__":
    main()
