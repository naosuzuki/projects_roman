#!/usr/bin/env python3
"""
21_config_ef.py -- Configuration F = Configuration E rotated by 15 deg, and an
alternating E/F observing scenario for the outskirt SN-host galaxies of ELAIS-N1.

Config E: 12 pointings on a ring at r=1.6 deg (angles 0,30,...,330).
Config F: same ring rotated 15 deg (angles 15,45,...,345) -> fills E's gaps.
Scenario: observe E for 1 h (all 12 pointings), then F for 1 h, and repeat the
pair. A host covered by nE E-pointings gains nE h per E round; nF h per F round
(overlaps integrate faster). Interleaving F fills the gaps between E pointings.

Reads 07_program_sne_ELAIS-N1.csv. Produces the layout + recovery figure.

    python 21_config_ef.py
"""
import os
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "07_program_sne_ELAIS-N1.csv")
RA0, DEC0 = 242.498, 54.497
FOV = 1.3
R = FOV / 2.0
R_RING = 1.5        # optimized for E+F: 100% FoV fill (no empty space), 95% coverage
M = 12
R_FOOT = 2.09


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def ring_centers(phase_deg):
    return [(R_RING * np.cos(np.radians(phase_deg + 360 * k / M)),
             R_RING * np.sin(np.radians(phase_deg + 360 * k / M))) for k in range(M)]


def core_centers():            # 7-pointing A/B core flower (covers the inner hosts)
    sp = np.sqrt(3) * R
    c = [(0.0, 0.0)]
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang); c.append((sp * np.cos(a), sp * np.sin(a)))
    return c


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    vneed = np.asarray(d["host_visits_needed"]).astype(float)
    vdone = np.asarray(d["host_visits_done"]).astype(float)
    pAh = np.asarray(d["pA_host"]).astype(int)
    pBh = np.asarray(d["pB_host"]).astype(int)
    cosd = np.cos(np.radians(DEC0))
    xi = (hra - RA0) * cosd; eta = hdec - DEC0
    rneed = np.maximum(1, np.ceil(vneed - vdone)).astype(int)
    rem = (htar == 1) & (comp == 0)
    inA = (pAh >= 0); inB = (pBh >= 0)
    inner = rem & (inA | inB)               # remaining hosts in the A/B core
    outer = rem & ~(inA | inB)              # remaining hosts in the outskirts
    pts = np.column_stack([xi, eta])
    nout = int(outer.sum()); ninn = int(inner.sum())

    EC, FC = ring_centers(0.0), ring_centers(15.0)
    nE = np.zeros(len(pts), int); nF = np.zeros(len(pts), int)
    for (xc, yc) in EC:
        nE += Path(hexagon(xc, yc)).contains_points(pts).astype(int)
    for (xc, yc) in FC:
        nF += Path(hexagon(xc, yc)).contains_points(pts).astype(int)
    covE = outer & (nE > 0)
    covEF = outer & ((nE + nF) > 0)
    print(f"outskirt hosts: {nout}")
    print(f"  Config E covers {covE.sum()} ({100*covE.sum()/nout:.0f}%); "
          f"E+F covers {covEF.sum()} ({100*covEF.sum()/nout:.0f}%)")

    # alternating E,F,E,F... : after R rounds, host hours = nEr*nE + nFr*nF
    Rmax = 24
    rounds = np.arange(0, Rmax + 1)
    def recovered(seq):                 # seq: 'EF' alternating or 'E' only
        rec = []
        for Rn in rounds:
            if seq == "EF":
                nEr = (Rn + 1) // 2; nFr = Rn // 2
            else:
                nEr = Rn; nFr = 0
            hrs = nEr * nE + nFr * nF
            rec.append(int((outer & (hrs >= rneed)).sum()))
        return np.array(rec)
    rec_ef = recovered("EF")
    exp = M * rounds
    # inner hosts recovered by re-observing the A/B core (alternate A,B; 7 pts/round)
    nA, nB = inA.astype(int), inB.astype(int)
    rec_in = np.array([int((inner & (((Rn + 1) // 2) * nA + (Rn // 2) * nB >= rneed)).sum())
                       for Rn in rounds])
    exp_in = 7 * rounds
    print(f"  inner hosts {ninn}; outskirts {nout}")
    print("  outskirts via E/F        inner via A/B core")
    for Rn in (2, 4, 6, 8):
        print(f"   {M*Rn:3d}exp -> {rec_ef[Rn]:4d} ({100*rec_ef[Rn]/nout:2.0f}%)   "
              f"{7*Rn:3d}exp -> {rec_in[Rn]:4d} ({100*rec_in[Rn]/ninn:2.0f}%)")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib import font_manager as fm
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, (axm, axr) = plt.subplots(1, 2, figsize=(15, 6.5))

    def t2sky(x, y):
        return RA0 + x / cosd, DEC0 + y
    ra, dec = t2sky(xi, eta)
    axm.scatter(ra[inner], dec[inner], s=6, c="#1f6fe0", alpha=0.5, linewidths=0,
                label="Remaining (inner)")
    axm.scatter(ra[outer], dec[outer], s=6, c="0.5", alpha=0.5, linewidths=0,
                label="Remaining (outskirts)")
    for (xc, yc) in core_centers():
        sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
        axm.add_patch(Polygon(sky, closed=True, fill=False, edgecolor="#9ecae1", lw=1.2, alpha=0.85))
    for grp, col, lab in ((EC, "#2ca02c", "Config E ring ($0^\\circ$)"),
                          (FC, "#d62728", "Config F ring ($15^\\circ$)")):
        first = True
        for (xc, yc) in grp:
            sky = np.array([t2sky(x, y) for x, y in hexagon(xc, yc)])
            axm.add_patch(Polygon(sky, closed=True, fill=False, edgecolor=col,
                                  lw=1.6, alpha=0.9, label=(lab if first else None)))
            first = False
    th = np.linspace(0, 2 * np.pi, 240)
    fra, fdec = t2sky(R_FOOT * np.cos(th), R_FOOT * np.sin(th))
    axm.plot(fra, fdec, color="0.4", ls="--", lw=1.2, label=f"Footprint edge ($r={R_FOOT}^\\circ$)")
    axm.set_aspect(1.0 / cosd); axm.invert_xaxis()
    axm.set_xlabel("RA (deg)", fontsize=16); axm.set_ylabel("Dec (deg)", fontsize=16)
    axm.set_title("Config E ($0^\\circ$) + Config F ($15^\\circ$) Interleaved Rings", fontsize=14)
    axm.tick_params(labelsize=12); axm.legend(fontsize=11, loc="upper right"); axm.grid(True, alpha=0.25)

    axr.plot(exp, 100 * rec_ef / nout, "-o", color="#d62728", ms=4,
             label=f"Outskirts via E/F ($N={nout}$)")
    axr.plot(exp_in, 100 * rec_in / ninn, "-s", color="#1f6fe0", ms=4,
             label=f"Inner via A/B core ($N={ninn}$)")
    axr.set_xlabel("Total 1-Hour Exposures", fontsize=16)
    axr.set_ylabel("Remaining Hosts Recovered (%)", fontsize=16)
    axr.set_title("Recovery vs Exposure: Inner and Outskirts", fontsize=15)
    axr.set_xlim(0, 100); axr.set_ylim(0, 100)
    axr.tick_params(labelsize=12); axr.legend(fontsize=12, loc="lower right"); axr.grid(True, alpha=0.3)

    png = os.path.join(PNG_DIR, "21_config_ef_ELAIS-N1.png")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
