#!/usr/bin/env python3
"""
03_edfs_targets.py -- EDF-S SN Ia targets for AAT follow-up.

AAT counterpart of ../programs_pfs/03 + 09 for the southern HLTDS field. Uses
the committed HOURGLASS2 catalog ../programs_pfs/outputs/csv/03_snia_radec.csv
(field, SNID, RA, Dec, z_cmb, peakmag_Z). peakmag_Z = -9 is a sentinel (no Z087
peak), treated as fainter than any cut.

The EDF-S footprint is two discs: a dense DEEP disc (NE) and a WIDE disc (SW).

Outputs:
  png/03_edfs_snia_radec.png    sky map colored by redshift, peak Z<23 marked
  png/03_edfs_snia_zhist.png    redshift histogram split by peak Z
  csv/03_edfs_snia.csv          per-SN table with disc + expected AAT success
"""
import os
import importlib.util as ilu
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
CAT = os.path.join(HERE, "..", "programs_pfs", "outputs", "csv", "03_snia_radec.csv")

# disc centres (deg; 2-means, median) and radii (deg; r99.9) fitted to the
# HOURGLASS2 EDF-S footprint
DISCS = {"DEEP": (63.596, -47.599), "WIDE": (58.915, -49.317)}
DISC_R = {"DEEP": 1.42, "WIDE": 1.81}
EDF_S = (62.6, -48.2)          # plotting tangent point
ZCUTS = (22.0, 23.0, 24.0)

_s = ilu.spec_from_file_location("aat_success", os.path.join(HERE, "02_aat_success.py"))
aat = ilu.module_from_spec(_s); _s.loader.exec_module(aat)


def angsep(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def gnomonic(ra, dec, ra0, dec0):
    """Tangent-plane (xi, eta) in deg about (ra0, dec0); xi increases to the East."""
    r, d, r0, d0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
    xi = np.cos(d) * np.sin(r - r0) / cosc
    eta = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0)) / cosc
    return np.degrees(xi), np.degrees(eta)


def load():
    d = np.genfromtxt(CAT, delimiter=",", names=True, dtype=None, encoding="utf-8")
    d = d[d["field"] == "EDF-S"]
    ra, dec = d["RA_deg"].astype(float), d["DEC_deg"].astype(float)
    pz = d["peakmag_Z"].astype(float)
    pz = np.where(pz < 0, 99.0, pz)                       # sentinel -> faint
    sep = {k: angsep(ra, dec, *c) for k, c in DISCS.items()}
    disc = np.where(sep["DEEP"] / DISC_R["DEEP"] < sep["WIDE"] / DISC_R["WIDE"], "DEEP", "WIDE")
    return dict(snid=d["SNID"], ra=ra, dec=dec, z=d["z_cmb"].astype(float),
                pz=pz, disc=disc, sep=sep)


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    for fn in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fn}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    return plt


def main():
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)
    T = load()
    n = len(T["z"])
    print(f"EDF-S SNe Ia: {n}   (no Z087 peak: {(T['pz'] > 90).sum()})")
    print("disc    N     r_max  z_med   " + "  ".join(f"Z<{c:.0f}" for c in ZCUTS))
    for k in DISCS:
        m = T["disc"] == k
        print(f"{k:5s} {m.sum():6d}  {T['sep'][k][m].max():5.2f}  {np.median(T['z'][m]):5.2f}   "
              + "  ".join(f"{(m & (T['pz'] < c)).sum():5d}" for c in ZCUTS))
    print("all   " + f"{n:6d}" + " " * 17 + "  ".join(f"{(T['pz'] < c).sum():5d}" for c in ZCUTS))

    # expected live-SN redshift success for a 2 h AAT visit at peak
    s2 = aat.success(T["pz"], 2.0)
    print(f"expected secure live-SN redshifts if every SN got 2 h at peak: {s2.sum():.0f}")

    csvf = os.path.join(CSV_DIR, "03_edfs_snia.csv")
    with open(csvf, "w") as fo:
        fo.write("SNID,disc,ra,dec,z,peakmag_Z,sep_deep,sep_wide,S_peak_2h\n")
        for i in range(n):
            pz = "" if T["pz"][i] > 90 else f"{T['pz'][i]:.3f}"
            fo.write(f"{T['snid'][i]},{T['disc'][i]},{T['ra'][i]:.5f},{T['dec'][i]:.5f},"
                     f"{T['z'][i]:.4f},{pz},{T['sep']['DEEP'][i]:.4f},"
                     f"{T['sep']['WIDE'][i]:.4f},{s2[i]:.4f}\n")
    print("table ->", csvf)

    plt = style()
    # ---- sky map ----
    xi, eta = gnomonic(T["ra"], T["dec"], *EDF_S)
    fig, ax = plt.subplots(figsize=(10, 8.5), constrained_layout=True)
    sc = ax.scatter(xi, eta, c=T["z"], s=0.6, cmap="rainbow", vmin=0, vmax=2.5,
                    alpha=0.5, rasterized=True)
    br = T["pz"] < 23
    ax.scatter(xi[br], eta[br], c=T["z"][br], s=14, cmap="rainbow", vmin=0, vmax=2.5,
               edgecolors="k", linewidths=0.4, label="Peak $Z<23$")
    for k, (r0, d0) in DISCS.items():
        x0, y0 = gnomonic(r0, d0, *EDF_S)
        ax.text(x0 + (-1.6 if k == "DEEP" else 1.9), y0 + (1.0 if k == "DEEP" else -1.0), k, ha="center",
                fontsize=20, weight="bold")
    ax.invert_xaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("$\\Delta$RA cos(Dec) [deg]", fontsize=18)
    ax.set_ylabel("$\\Delta$Dec [deg]", fontsize=18)
    ax.set_title(f"EDF-S SNe Ia (HOURGLASS2, N = {n:,})", fontsize=20)
    ax.tick_params(labelsize=14, direction="in")
    ax.legend(fontsize=15, loc="lower left", markerscale=1.5)
    cb = fig.colorbar(sc, ax=ax, pad=0.01, shrink=0.85)
    cb.solids.set_alpha(1)
    cb.set_label("Redshift", fontsize=18)
    cb.ax.tick_params(labelsize=14)
    png = os.path.join(PNG_DIR, "03_edfs_snia_radec.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)

    # ---- redshift histogram split by peak Z ----
    fig, ax = plt.subplots(figsize=(9, 6.5), constrained_layout=True)
    bins = np.arange(0, 3.01, 0.1)
    groups = [(T["pz"] < 22, "Peak $Z<22$", "#d62728"),
              ((T["pz"] >= 22) & (T["pz"] < 23), "$22\\leq Z<23$", "#ff7f0e"),
              ((T["pz"] >= 23) & (T["pz"] < 24), "$23\\leq Z<24$", "#2ca02c"),
              (T["pz"] >= 24, "$Z\\geq24$ or no peak", "#1f77b4")]
    ax.hist([T["z"][g] for g, _, _ in groups], bins=bins, stacked=True,
            color=[c for _, _, c in groups],
            label=[f"{l} ({g.sum():,})" for g, l, _ in groups])
    ax.set_xlabel("Redshift", fontsize=18)
    ax.set_ylabel("Number of SNe Ia", fontsize=18)
    ax.set_title("EDF-S SN Ia Redshift Distribution", fontsize=20)
    ax.tick_params(labelsize=14, direction="in")
    ax.legend(fontsize=14)
    ax.grid(True, alpha=0.3)
    png = os.path.join(PNG_DIR, "03_edfs_snia_zhist.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
