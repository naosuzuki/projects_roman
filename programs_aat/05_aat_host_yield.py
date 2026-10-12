#!/usr/bin/env python3
"""
05_aat_host_yield.py -- AAT 2dF host-galaxy redshift campaign over EDF-S.

AAT counterpart of the PFS host simulation (../programs_pfs/07, 24, 25).

Inputs
  * EDF-S SNe Ia (03_edfs_targets.py) and the 9-field 2dF layout
    (04_2df_pointings.py).
  * Host Z087 magnitudes: the committed local EDF-S table has none, so each
    EDF-S SN is given the host of a random ELAIS-N1 SN Ia of the same redshift
    (|dz| < 0.05) from the same HOURGLASS2 simulation and SN Ia model
    (../programs_pfs/outputs/csv/03_snia_host_radec_ELAIS-N1.csv). ELAIS-N1
    SNe without a host carry over as hostless.
  * Redshift success S(Z, t) from 02_aat_success.py (HSC-AAT calibrated).

Campaign model
  * A "round" observes each of the 9 fields once for T_VIS hours with
    FIB_SCI target fibers (392 science fibers minus sky fibers).
  * In a field, fibers go to covered, unsecured, unretired hosts with
    Z < HOST_ZCUT, brightest first (OzDES-like: keep a target until it yields
    a redshift). Overlapping fields let a host accumulate hours from each.
  * Each host has a fixed random threshold u ~ U(0,1); it is secured once
    S(Z, hours) >= u. A host is retired after HMAX hours.

Outputs
  csv/05_aat_host_yield.csv      secured redshifts vs total AAT hours
  csv/05_aat_host_final.csv      per-host final state
  png/05_aat_host_yield.png      yield vs hours, and success by host Z / redshift
"""
import os
import argparse
import importlib.util as ilu
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")
PFS_CSV = os.path.join(HERE, "..", "programs_pfs", "outputs", "csv")


def _load(name, fn):
    s = ilu.spec_from_file_location(name, os.path.join(HERE, fn))
    m = ilu.module_from_spec(s); s.loader.exec_module(m)
    return m


aat = _load("aat_success", "02_aat_success.py")
ed = _load("edfs", "03_edfs_targets.py")
R2DF = 1.05


def assign_hosts(z_edfs, rng, dz=0.05):
    """Host Z087 for each EDF-S SN, drawn from ELAIS-N1 SNe Ia at the same z."""
    sn = np.genfromtxt(os.path.join(PFS_CSV, "03_snia_radec.csv"), delimiter=",",
                       names=True, dtype=None, encoding="utf-8")
    sn = sn[sn["field"] == "ELAIS-N1"]
    h = np.genfromtxt(os.path.join(PFS_CSV, "03_snia_host_radec_ELAIS-N1.csv"),
                      delimiter=",", names=True, dtype=None, encoding="utf-8")
    hmag = dict(zip(h["SNID"], h["peakmag_Z"].astype(float)))
    zn = sn["z_cmb"].astype(float)
    hn = np.array([hmag.get(i, np.nan) for i in sn["SNID"]])     # nan = hostless
    order = np.argsort(zn)
    zn, hn = zn[order], hn[order]
    lo = np.searchsorted(zn, z_edfs - dz)
    hi = np.searchsorted(zn, z_edfs + dz)
    out = np.empty(len(z_edfs))
    for i, (a, b) in enumerate(zip(lo, hi)):
        if b <= a:                                   # no neighbour: nearest
            a = min(max(a - 1, 0), len(zn) - 1); b = a + 1
        out[i] = hn[rng.integers(a, b)]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rounds", type=int, default=60, help="rounds of all 9 fields")
    ap.add_argument("--tvis", type=float, default=1.0, help="hours per field visit")
    ap.add_argument("--fibers", type=int, default=360, help="target fibers per field")
    ap.add_argument("--zcut", type=float, default=24.0, help="host Z087 limit")
    ap.add_argument("--hmax", type=float, default=20.0, help="retire a host after HMAX h")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(PNG_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    T = ed.load()
    hz = assign_hosts(T["z"], rng)
    hostless = ~np.isfinite(hz)
    target = (~hostless) & (hz < args.zcut)
    fields = np.genfromtxt(os.path.join(CSV_DIR, "04_2df_pointings.csv"), delimiter=",",
                           names=True, dtype=None, encoding="utf-8")
    nf = len(fields)
    inf = np.array([ed.angsep(T["ra"], T["dec"], r, d) < R2DF
                    for r, d in zip(fields["ra"], fields["dec"])])
    n = len(hz)
    print(f"EDF-S SNe Ia {n}: hostless {hostless.sum()}, host Z<{args.zcut}: {target.sum()}"
          f"  (Z<22 {(target & (hz < 22)).sum()}, Z<23 {(target & (hz < 23)).sum()})")

    u = rng.random(n)
    hours = np.zeros(n)
    secured = np.zeros(n, bool)
    retired = np.zeros(n, bool)
    hzs = np.where(target, hz, 99.0)
    log = []                                  # (total_hours, n_secured, fibers_used_mean)
    tot = 0.0
    for r in range(args.rounds):
        used = []
        for f in range(nf):
            cand = np.where(inf[f] & target & ~secured & ~retired)[0]
            cand = cand[np.argsort(hzs[cand])][:args.fibers]
            used.append(len(cand))
            hours[cand] += args.tvis
            secured[cand] |= aat.success(hz[cand], hours[cand]) >= u[cand]
            retired[cand] |= (~secured[cand]) & (hours[cand] >= args.hmax)
            tot += args.tvis
        log.append((tot, secured.sum(), np.mean(used)))

    print("\n  hours  secured  frac_of_targets  fibers_used/field")
    for t, s, fu in log:
        if t in (9, 18, 45, 90, 180, 270, 360, 450, 540) or t == log[-1][0]:
            print(f"  {t:5.0f}  {s:7d}   {s / target.sum():6.3f}        {fu:6.0f}")

    zb = [(0, 0.5), (0.5, 0.8), (0.8, 1.1), (1.1, 1.5), (1.5, 3.0)]
    print("\n  SN redshift   N_SN  N_target  secured  frac_of_all_SNe")
    for a, b in zb:
        m = (T["z"] >= a) & (T["z"] < b)
        print(f"  {a:.1f}-{b:.1f}     {m.sum():6d}  {(m & target).sum():6d}  "
              f"{(m & secured).sum():6d}   {(m & secured).sum() / m.sum():.3f}")
    print("\n  host Z    N_target  secured  frac")
    for a in np.arange(18, args.zcut, 1.0):
        m = target & (hz >= a) & (hz < a + 1)
        if m.sum():
            print(f"  {a:.0f}-{a + 1:.0f}    {m.sum():6d}  {(m & secured).sum():6d}   "
                  f"{(m & secured).sum() / m.sum():.3f}")

    with open(os.path.join(CSV_DIR, "05_aat_host_yield.csv"), "w") as fo:
        fo.write("total_hours,n_secured,frac_targets,fibers_used_per_field\n")
        for t, s, fu in log:
            fo.write(f"{t:.1f},{s},{s / target.sum():.4f},{fu:.1f}\n")
    with open(os.path.join(CSV_DIR, "05_aat_host_final.csv"), "w") as fo:
        fo.write("SNID,disc,z,host_Z,target,hours,secured,retired\n")
        for i in range(n):
            hzv = "" if hostless[i] else f"{hz[i]:.3f}"
            fo.write(f"{T['snid'][i]},{T['disc'][i]},{T['z'][i]:.4f},{hzv},"
                     f"{int(target[i])},{hours[i]:.1f},{int(secured[i])},{int(retired[i])}\n")
    print("tables ->", os.path.join(CSV_DIR, "05_aat_host_{yield,final}.csv"))

    # ---- plot ----
    plt = ed.style()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 6.2), constrained_layout=True)
    t, s, _ = np.array(log).T
    a1.plot(t, s, color="#3b7dd8", lw=3)
    a1.axhline(target.sum(), color="0.4", ls="--", lw=1.2)
    a1.text(t[-1], target.sum(), f"Hosts with $Z<{args.zcut:.0f}$: {target.sum():,}",
            ha="right", va="bottom", fontsize=14)
    a1.set_xlabel("Total AAT Exposure [Hours]", fontsize=18)
    a1.set_ylabel("Secure Host Redshifts", fontsize=18)
    a1.set_title(f"EDF-S Host Redshifts with 2dF ({nf} Fields)", fontsize=19)
    a1.set_xlim(0, t[-1]); a1.set_ylim(0, target.sum() * 1.1)
    a1.tick_params(labelsize=14, direction="in"); a1.grid(True, alpha=0.3)

    bins = np.arange(0, 3.01, 0.1)
    a2.hist([T["z"][secured], T["z"][target & ~secured], T["z"][~target]], bins=bins,
            stacked=True, color=["#d62728", "#ff7f0e", "#bbbbbb"],
            label=[f"Secured ({secured.sum():,})",
                   f"Targeted, no redshift ({(target & ~secured).sum():,})",
                   f"Host $Z\\geq{args.zcut:.0f}$ or hostless ({(~target).sum():,})"])
    a2.set_xlabel("Redshift", fontsize=18)
    a2.set_ylabel("Number of SNe Ia", fontsize=18)
    a2.set_title(f"After {t[-1]:.0f} h ({t[-1] / nf:.0f} h per Field)", fontsize=19)
    a2.tick_params(labelsize=14, direction="in"); a2.legend(fontsize=13)
    a2.grid(True, alpha=0.3)
    png = os.path.join(PNG_DIR, "05_aat_host_yield.png")
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
