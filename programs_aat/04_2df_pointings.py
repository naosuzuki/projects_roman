#!/usr/bin/env python3
"""
04_2df_pointings.py -- AAT 2dF pointing layout over EDF-S.

AAT counterpart of ../programs_pfs/06_pfs_pointings.py. 2dF has a circular
2.1 deg diameter field (R = 1.05 deg) with 392 science fibers. EDF-S is two
discs (DEEP r~1.42 deg, WIDE r~1.81 deg; see 03_edfs_targets.py), so each disc
is tiled separately with N 2dF fields: either N fields on a ring of radius d,
or one central field plus N-1 on the ring. For each N we search d and the ring
rotation to maximize the fraction of SNe Ia covered.

Outputs:
  csv/04_2df_coverage.csv     coverage vs N per disc (all SNe, peak Z<23, Z<24)
  csv/04_2df_pointings.csv    adopted field centres (RA, Dec)
  png/04_2df_pointings.png    adopted layout over the SN map
"""
import os
import importlib.util as ilu
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
R2DF = 1.05                     # deg, 2dF field radius (2.1 deg diameter)
FIBERS = 392                    # science fibers per 2dF configuration
NMAX = 7
COVER_GOAL = 0.99               # adopt the smallest N covering >= 99% of SNe

_s = ilu.spec_from_file_location("edfs", os.path.join(HERE, "03_edfs_targets.py"))
ed = ilu.module_from_spec(_s); _s.loader.exec_module(ed)


def layout(n, center, d, rot):
    """Field centres (tangent-plane deg about the disc centre)."""
    k = n - 1 if center else n
    th = np.radians(rot) + 2 * np.pi * np.arange(k) / max(k, 1)
    pts = [(d * np.cos(t), d * np.sin(t)) for t in th]
    return ([(0.0, 0.0)] if center else []) + pts


def covered(x, y, fields):
    m = np.zeros(len(x), bool)
    for fx, fy in fields:
        m |= np.hypot(x - fx, y - fy) < R2DF
    return m


def best_layout(x, y, n):
    if n == 1:
        f = [(0.0, 0.0)]
        return covered(x, y, f).mean(), f
    best = (-1, None)
    for center in (False, True):
        if center and n < 3:
            continue
        for d in np.arange(0.3, 1.81, 0.02):
            for rot in np.arange(0, 360.0 / (n - center), 3.0):
                f = layout(n, center, d, rot)
                c = covered(x, y, f).mean()
                if c > best[0] + 1e-9:
                    best = (c, f)
    return best


def tangent_to_radec(x, y, ra0, dec0):
    X, Y = np.radians(x), np.radians(y)
    a0, d0 = np.radians(ra0), np.radians(dec0)
    rho = np.hypot(X, Y)
    if rho == 0:
        return ra0, dec0
    c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(d0) + Y * np.sin(c) * np.cos(d0) / rho)
    ra = a0 + np.arctan2(X * np.sin(c), rho * np.cos(d0) * np.cos(c) - Y * np.sin(d0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def main():
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)
    T = ed.load()
    rows, adopted = [], {}
    for disc, (ra0, dec0) in ed.DISCS.items():
        m = T["disc"] == disc
        x, y = ed.gnomonic(T["ra"][m], T["dec"][m], ra0, dec0)
        pz = T["pz"][m]
        print(f"{disc}: N_SN = {m.sum()}")
        print("   N   all    Z<23   Z<24   SNe/field(Z<23)")
        for n in range(1, NMAX + 1):
            c, f = best_layout(x, y, n)
            cov = covered(x, y, f)
            c23 = cov[pz < 23].mean()
            c24 = cov[pz < 24].mean()
            per = (cov & (pz < 23)).sum() / n
            rows.append((disc, n, c, c23, c24))
            print(f"  {n:2d}  {c:5.3f}  {c23:5.3f}  {c24:5.3f}   {per:6.0f}")
            if disc not in adopted and c >= COVER_GOAL:
                adopted[disc] = (n, f)
        if disc not in adopted:
            adopted[disc] = (n, f)

    with open(os.path.join(CSV_DIR, "04_2df_coverage.csv"), "w") as fo:
        fo.write("disc,n_fields,cover_all,cover_Zlt23,cover_Zlt24\n")
        for r in rows:
            fo.write(f"{r[0]},{r[1]},{r[2]:.4f},{r[3]:.4f},{r[4]:.4f}\n")

    fields = []
    with open(os.path.join(CSV_DIR, "04_2df_pointings.csv"), "w") as fo:
        fo.write("field,disc,ra,dec\n")
        for disc, (n, f) in adopted.items():
            for i, (fx, fy) in enumerate(f):
                ra, dec = tangent_to_radec(fx, fy, *ed.DISCS[disc])
                name = f"{disc[0]}{i}"
                fields.append((name, disc, ra, dec))
                fo.write(f"{name},{disc},{ra:.4f},{dec:.4f}\n")
    ntot = len(fields)
    print(f"adopted: " + ", ".join(f"{d} {n} fields" for d, (n, _) in adopted.items())
          + f"  -> {ntot} 2dF fields")

    # per-field target counts over the whole survey
    sep = np.array([ed.angsep(T["ra"], T["dec"], ra, dec) for _, _, ra, dec in fields])
    inf = sep < R2DF
    ncov = inf.sum(axis=0)
    print(f"covered: all {np.mean(ncov > 0):.3f}, Z<23 {np.mean(ncov[T['pz'] < 23] > 0):.3f};"
          f" mean fields per covered SN {ncov[ncov > 0].mean():.2f}")
    for (name, disc, _, _), row in zip(fields, inf):
        print(f"   {name}: {row.sum():5d} SNe, {(row & (T['pz'] < 23)).sum():4d} with Z<23")

    # ---- plot ----
    plt = ed.style()
    xi, eta = ed.gnomonic(T["ra"], T["dec"], *ed.EDF_S)
    cov = ncov > 0
    fig, ax = plt.subplots(figsize=(10, 8.5), constrained_layout=True)
    ax.scatter(xi[~cov], eta[~cov], s=0.5, c="0.6", rasterized=True)
    ax.scatter(xi[cov], eta[cov], s=0.5, c="#3b7dd8", alpha=0.5, rasterized=True)
    t = np.linspace(0, 2 * np.pi, 300)
    for name, disc, ra, dec in fields:
        fx, fy = ed.gnomonic(ra, dec, *ed.EDF_S)
        # draw the true 1.05-deg circle on the sky, projected
        cr = ra + R2DF * np.cos(t) / np.cos(np.radians(dec))
        cd = dec + R2DF * np.sin(t)
        px, py = ed.gnomonic(cr, cd, *ed.EDF_S)
        ax.plot(px, py, color="#d62728" if disc == "DEEP" else "#2ca02c", lw=2)
        n_in = inf[[f[0] for f in fields].index(name)].sum()
        ax.text(fx, fy, f"{name}\n{n_in}", ha="center", va="center", fontsize=13,
                weight="bold")
    ax.invert_xaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("$\\Delta$RA cos(Dec) [deg]", fontsize=18)
    ax.set_ylabel("$\\Delta$Dec [deg]", fontsize=18)
    ax.set_title(f"AAT 2dF Pointings over EDF-S ({ntot} Fields, "
                 f"{100 * cov.mean():.0f}% of SNe Ia Covered)", fontsize=19)
    ax.tick_params(labelsize=14, direction="in")
    png = os.path.join(PNG_DIR, "04_2df_pointings.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
