#!/usr/bin/env python3
"""
09_desi_phase3.py -- DESI Phase III host-galaxy completeness for ELAIS-N1,
Configuration B (6-pointing hexagon).  The analog of the Roman-Subaru study's
Figures 19 (remaining hosts vs host Z) and 20 (sky map of remaining hosts), plus
a Phase-III planning curve: how many integration hours are needed to complete the
hosts left after the cadence survey.

For each cadence scenario (5 d, 10 d) we read the DESI cadence catalog written by
07_desi_cadence.py (host_completed) and the committed program catalog (host_ra,
host_dec, host_Z, host_visits_needed = PFS hours to S/N=5).  A host's DESI
integration requirement is APERTURE x (PFS hours) for the 4-m vs 8.2-m aperture.

Figures
  * 09_desi_remaining_hosts_cad{N}.png      -- Fig 19: host Z, completed vs remaining
  * 09_desi_remaining_hosts_radec_cad{N}.png-- Fig 20: sky map of remaining hosts
  * 09_desi_hours_vs_remaining.png          -- hosts completed/remaining vs hours
                                               per pointing (Config B), both cadences

    python 09_desi_phase3.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
PROG_CSV = os.path.join(HERE, "..", "programs_pfs", "outputs", "csv",
                        "07_program_sne_ELAIS-N1.csv")

RA0, DEC0 = 242.498, 54.497
FOV = 3.2
R = FOV / 2.0
D_HEX = 0.40
APERTURE = 4.2
HOST_ZCUT = 25.5
CADENCES = (5, 10)


def hexagon_centers(d=D_HEX, phase_deg=0.0):
    a = np.radians(phase_deg + 60.0 * np.arange(6))
    return np.column_stack([d * np.cos(a), d * np.sin(a)])


def assign_pointing(xi, eta, centers):
    out = np.full(len(xi), -1, int)
    best = np.full(len(xi), 1e9)
    for i, (xc, yc) in enumerate(centers):
        dd = np.hypot(xi - xc, eta - yc)
        take = (dd <= R) & (dd < best)
        out[take] = i
        best[take] = dd[take]
    return out


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


def main():
    os.makedirs(PNG_DIR, exist_ok=True)

    # ---- committed program catalog: host positions, Z, PFS hours ----
    d = np.genfromtxt(PROG_CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    hra, hdec = d["host_ra"].astype(float), d["host_dec"].astype(float)
    hostz = d["host_Z"].astype(float)
    host_target = d["host_target"].astype(int)
    vneed_pfs = d["host_visits_needed"].astype(float)
    n = len(hostz)
    cosd = np.cos(np.radians(DEC0))

    cen = hexagon_centers()
    pho = assign_pointing((hra - RA0) * cosd, hdec - DEC0, cen)
    host_ok = (host_target == 1) & np.isfinite(hostz)
    covered = host_ok & (pho >= 0)                  # Config-B Phase III pool
    need_desi = APERTURE * vneed_pfs                # DESI hours to S/N=5 per host

    plt = setup_fonts()

    # =================================================================
    # Fig 19 & 20 per scenario: completed vs remaining after the survey
    # =================================================================
    achieved = {}
    for cad in CADENCES:
        cat = os.path.join(CSV_DIR, f"07_desi_program_sne_cad{cad}.csv")
        dc = np.genfromtxt(cat, delimiter=",", names=True, dtype=None, encoding="utf-8")
        completed = dc["host_completed"].astype(int) == 1
        comp = covered & completed
        rem = covered & (~completed)
        achieved[cad] = int(comp.sum())
        print(f"cadence {cad} d: Config-B hosts {int(covered.sum())}; "
              f"completed {int(comp.sum())} ({100*comp.mean()/covered.mean():.0f}% of covered); "
              f"remaining {int(rem.sum())}")

        # ---- Fig 19: host Z histogram, completed (blue) vs remaining (red) ----
        fig, ax = plt.subplots(figsize=(9, 6))
        hbins = np.arange(20, HOST_ZCUT + 0.01, 0.25)
        ax.hist(hostz[comp], bins=hbins, color="lightskyblue",
                label=f"Completed (S/N$=5$, $N={int(comp.sum())}$)")
        ax.hist(hostz[rem], bins=hbins, color="red", alpha=0.65,
                label=f"Remaining (incomplete, $N={int(rem.sum())}$)")
        ax.set_xlabel("Host Roman $Z$ (AB Mag)", fontsize=16)
        ax.set_ylabel("Number of Hosts", fontsize=16)
        ax.set_title(f"ELAIS-N1 SN-Host Galaxies: Completed vs Remaining  "
                     f"(DESI Config B, cadence {cad} d)", fontsize=14)
        ax.tick_params(labelsize=12); ax.legend(fontsize=12); ax.grid(True, alpha=0.3)
        png = os.path.join(PNG_DIR, f"09_desi_remaining_hosts_cad{cad}.png")
        fig.tight_layout(); fig.savefig(png, dpi=140); print("plot ->", png); plt.close(fig)

        # ---- Fig 20: sky map of remaining hosts colored by host Z ----
        fig, ax = plt.subplots(figsize=(8.4, 8.0))
        sc = ax.scatter(hra[rem], hdec[rem], c=hostz[rem], cmap="rainbow_r",
                        s=7, vmin=22, vmax=25.5, linewidths=0)
        # Config-B field outlines for context
        from matplotlib.patches import Polygon

        def circle_sky(xc, yc, rad, n=200):
            a = np.linspace(0, 2 * np.pi, n)
            return np.column_stack([RA0 + (xc + rad * np.cos(a)) / cosd,
                                    DEC0 + (yc + rad * np.sin(a))])
        for xc, yc in cen:
            ax.add_patch(Polygon(circle_sky(xc, yc, R), closed=True, fill=False,
                                 edgecolor="0.5", lw=1.0, alpha=0.7))
        ax.set_aspect(1.0 / cosd); ax.invert_xaxis()
        ax.set_xlabel("RA [deg]", fontsize=16); ax.set_ylabel("Dec [deg]", fontsize=16)
        ax.set_title(f"ELAIS-N1 Remaining Hosts ($N={int(rem.sum())}$)  "
                     f"(DESI Config B, cadence {cad} d)", fontsize=13)
        ax.tick_params(labelsize=12); ax.grid(True, alpha=0.25)
        cb = fig.colorbar(sc, ax=ax, pad=0.01); cb.set_label("Host Roman $Z$ (AB Mag)", fontsize=14)
        png = os.path.join(PNG_DIR, f"09_desi_remaining_hosts_radec_cad{cad}.png")
        fig.tight_layout(); fig.savefig(png, dpi=140); print("plot ->", png); plt.close(fig)

    # =================================================================
    # Hours-vs-remaining: Phase III recovery curve (Config B, scenario-indep.)
    # If every pointing is integrated to a uniform depth T (hours), a covered
    # host completes when T >= its DESI requirement.  Per-pointing depth (not
    # summed) because the 6 fields integrate in parallel on different nights.
    # =================================================================
    need_cov = need_desi[covered]
    ntot = len(need_cov)
    T = np.linspace(0, 200, 1000)                   # hours per pointing
    completed_vs_T = np.array([(need_cov <= t).sum() for t in T])
    remaining_vs_T = ntot - completed_vs_T

    fig, ax = plt.subplots(figsize=(9.5, 6.3))
    ax.plot(T, completed_vs_T, "-", color="#1f77b4", lw=2.2, label="Hosts completed")
    ax.plot(T, remaining_vs_T, "-", color="#d62728", lw=2.2, label="Hosts remaining")
    # mark a few completeness depths
    for frac, txt in ((0.5, "50%"), (0.8, "80%"), (0.9, "90%")):
        th = np.interp(frac * ntot, completed_vs_T, T)
        ax.axvline(th, color="0.6", ls=":", lw=1)
        ax.text(th + 1.5, ntot * 0.10, f"{txt}\n{th:.0f} h", fontsize=10, color="0.3")
    # mark each cadence survey's achieved completion
    for cad, col in zip(CADENCES, ("#2ca02c", "#9467bd")):
        ax.axhline(achieved[cad], color=col, ls="--", lw=1.4,
                   label=f"Cadence {cad} d survey: {achieved[cad]} completed")
    ax.set_xlabel("Integration Hours per Pointing (Config B, 6 fields)", fontsize=15)
    ax.set_ylabel("Number of Hosts", fontsize=15)
    ax.set_title("DESI Phase III: Host Spec-$z$ Recovery vs Integration Time  "
                 f"({ntot} Config-B hosts, $4.2\\times$ aperture penalty)", fontsize=13)
    secax = ax.secondary_xaxis("top", functions=(lambda h: 6 * h, lambda p: p / 6))
    secax.set_xlabel("Total Exposure (pointing-hours, all 6 fields)", fontsize=13)
    ax.set_xlim(0, 150); ax.set_ylim(0, ntot * 1.02)
    ax.tick_params(labelsize=12); secax.tick_params(labelsize=11)
    ax.legend(fontsize=11, loc="center right"); ax.grid(True, alpha=0.3)
    png = os.path.join(PNG_DIR, "09_desi_hours_vs_remaining.png")
    fig.tight_layout(); fig.savefig(png, dpi=140); print("plot ->", png); plt.close(fig)

    # quick text summary of depths
    print("\nPhase III depth to complete the Config-B host pool:")
    for frac in (0.5, 0.8, 0.9, 0.95):
        th = np.interp(frac * ntot, completed_vs_T, T)
        print(f"  {int(100*frac)}% ({int(frac*ntot)} hosts): {th:.0f} h/pointing  "
              f"({6*th:.0f} pointing-hours total)")


if __name__ == "__main__":
    main()
