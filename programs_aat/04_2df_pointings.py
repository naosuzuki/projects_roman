#!/usr/bin/env python3
"""
04_2df_pointings.py -- AAT 2dF pointing layout over EDF-S.

AAT counterpart of ../programs_pfs/06_pfs_pointings.py. 2dF has a circular
2.1 deg diameter field (R = 1.05 deg) with 392 science fibers. EDF-S is two
discs (DEEP r~1.42 deg, WIDE r~1.81 deg; see 03_edfs_targets.py), so each disc
is tiled separately with N 2dF fields: either N fields on a ring of radius d,
or one central field plus N-1 on the ring. For each N we search d and the ring
rotation. We adopt the smallest N that can cover >= COVER_MIN of the disc's
SNe and, among those layouts, the TIGHTEST one: maximal mean number of fields
per SN (central overlap = extra hours per host), which is equivalent to
minimal field area falling on empty sky outside the disc. Chasing the last
~1% of edge SNe would otherwise push the fields outward and waste ~10% of the
2dF area on blank sky.

    python 04_2df_pointings.py                  # COVER_MIN = 0.99
    python 04_2df_pointings.py --cover-min 1.0  # maximal-coverage layout

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
COVER_MIN = 0.99                # required SN coverage per disc (default)

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


def candidates(n):
    if n == 1:
        yield [(0.0, 0.0)]
        return
    for center in (False, True):
        if center and n < 3:
            continue
        for d in np.arange(0.0, 1.81, 0.02):
            for rot in np.arange(0, 360.0 / (n - center), 3.0):
                yield layout(n, center, d, rot)


def best_layout(x, y, n, cover_min):
    """(max coverage, its layout), (tightest layout with coverage >= cover_min or None)."""
    best_cov = (-1, None)
    tight = (-1, None, None)                         # (mean overlap, coverage, layout)
    for f in candidates(n):
        nc = sum((np.hypot(x - fx, y - fy) < R2DF).astype(int) for fx, fy in f)
        c = (nc > 0).mean()
        if c > best_cov[0] + 1e-9:
            best_cov = (c, f)
        if c >= cover_min - 1e-9 and nc.mean() > tight[0] + 1e-9:
            tight = (nc.mean(), c, f)
    return best_cov, (tight if tight[2] is not None else None)


def waste(f, rdisc):
    """Fraction of the summed 2dF field area lying outside the disc (blank sky)."""
    g = np.linspace(-3.5, 3.5, 351)
    gx, gy = np.meshgrid(g, g)
    ind = np.hypot(gx, gy) < rdisc
    tot = inside = 0
    for fx, fy in f:
        inf = np.hypot(gx - fx, gy - fy) < R2DF
        tot += inf.sum(); inside += (inf & ind).sum()
    return 1 - inside / tot


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
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cover-min", type=float, default=COVER_MIN, dest="cover_min")
    args = ap.parse_args()
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)
    T = ed.load()
    rows, adopted = [], {}
    for disc, (ra0, dec0) in ed.DISCS.items():
        m = T["disc"] == disc
        x, y = ed.gnomonic(T["ra"][m], T["dec"][m], ra0, dec0)
        pz = T["pz"][m]
        print(f"{disc}: N_SN = {m.sum()}")
        print("        max-coverage layout          |  tightest layout with coverage >= "
              f"{args.cover_min:.2f}")
        print("   N   cover  overlap  blank-sky   |  cover  Z<23   overlap  blank-sky")
        for n in range(1, NMAX + 1):
            (c, f), tight = best_layout(x, y, n, args.cover_min)
            nc = sum((np.hypot(x - fx, y - fy) < R2DF).astype(int) for fx, fy in f)
            line = f"  {n:2d}  {c:5.3f}  {nc.mean():5.2f}    {waste(f, ed.DISC_R[disc]):5.3f}     |"
            row = [disc, n, c, nc.mean(), waste(f, ed.DISC_R[disc])]
            if tight:
                ov, ct, ft = tight
                cov = covered(x, y, ft)
                wt = waste(ft, ed.DISC_R[disc])
                line += f"  {ct:5.3f}  {cov[pz < 23].mean():5.3f}  {ov:5.2f}     {wt:5.3f}"
                row += [ct, cov[pz < 23].mean(), ov, wt]
                if disc not in adopted:
                    adopted[disc] = (n, ft)
            else:
                row += [np.nan] * 4
            rows.append(row)
            print(line)
        if disc not in adopted:
            adopted[disc] = (n, f)

    with open(os.path.join(CSV_DIR, "04_2df_coverage.csv"), "w") as fo:
        fo.write("disc,n_fields,maxcov_cover,maxcov_overlap,maxcov_blank,"
                 "tight_cover,tight_cover_Zlt23,tight_overlap,tight_blank\n")
        for r in rows:
            fo.write(",".join([r[0], str(r[1])] + [f"{v:.4f}" for v in r[2:]]) + "\n")

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
    blank = {d: waste(f, ed.DISC_R[d]) for d, (n, f) in adopted.items()}
    print("blank-sky fraction of 2dF area: " + ", ".join(f"{d} {w:.3f}" for d, w in blank.items()))

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
    for disc, (r0, d0) in ed.DISCS.items():
        cr = r0 + ed.DISC_R[disc] * np.cos(t) / np.cos(np.radians(d0))
        cd = d0 + ed.DISC_R[disc] * np.sin(t)
        px, py = ed.gnomonic(cr, cd, *ed.EDF_S)
        ax.plot(px, py, color="k", lw=1.2, ls="--")
    for name, disc, ra, dec in fields:
        fx, fy = ed.gnomonic(ra, dec, *ed.EDF_S)
        # draw the true 1.05-deg circle on the sky, projected
        cr = ra + R2DF * np.cos(t) / np.cos(np.radians(dec))
        cd = dec + R2DF * np.sin(t)
        px, py = ed.gnomonic(cr, cd, *ed.EDF_S)
        ax.plot(px, py, color="#d62728" if disc == "DEEP" else "#2ca02c", lw=2)
        n_in = inf[[f[0] for f in fields].index(name)].sum()
        # label pushed outward from the disc centre so overlapping fields stay legible
        cx, cy = ed.gnomonic(*ed.DISCS[disc], *ed.EDF_S)
        vx, vy = fx - cx, fy - cy
        rr = np.hypot(vx, vy)
        lx, ly = (cx + vx / rr * (rr + 0.6), cy + vy / rr * (rr + 0.6)) if rr > 1e-3 else (fx, fy)
        ax.text(lx, ly, f"{name}\n{n_in}", ha="center", va="center", fontsize=13,
                weight="bold", bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
    ax.invert_xaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("$\\Delta$RA cos(Dec) [deg]", fontsize=18)
    ax.set_ylabel("$\\Delta$Dec [deg]", fontsize=18)
    ax.set_title(f"AAT 2dF Pointings over EDF-S ({ntot} Fields, "
                 f"{100 * cov.mean():.0f}% of SNe Ia Covered, "
                 f"{ncov[cov].mean():.1f} Fields per SN)", fontsize=18)
    ax.tick_params(labelsize=14, direction="in")
    png = os.path.join(PNG_DIR, "04_2df_pointings.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
