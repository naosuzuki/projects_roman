#!/usr/bin/env python3
"""
15_snia_z_reach.py -- redshift reach of the ELAIS-N1 SN Ia program: per Delta z
= 0.1 bin, how many SNe Ia reach the live-SN trigger (Z087 < 24 at peak), how
many of their hosts are Phase III targets (Z < 25.5), how many hosts complete to
S/N = 5, and the median host Z. Shows that the survey's z ~ 1.2 ceiling is set
by the SN peak-brightness cut (the program-SN counts collapse beyond z ~ 1.1),
not by host faintness (the high-z hosts that exist are bright, ~22.5-23.7 mag).

Reads the program-SN catalog written by 07_fiber_budget.py.

    python 15_snia_z_reach.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(HERE, "outputs", "csv")
CSV = os.path.join(CSV_DIR, "07_program_sne_ELAIS-N1.csv")


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    typ = np.asarray(d["type"]).astype(str)
    z = np.asarray(d["z"]).astype(float)
    hz = np.asarray(d["host_Z"]).astype(float)
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    ia = typ == "Ia"

    out = os.path.join(CSV_DIR, "15_snia_z_reach_ELAIS-N1.csv")
    edges = np.round(np.arange(0.0, 1.5 + 1e-9, 0.1), 1)
    print("  z-bin    progSNe  hostZ<25.5  completed  medianHostZ")
    with open(out, "w") as fo:
        fo.write("z_lo,z_hi,prog_SNe,host_zlt255,completed,median_host_Z\n")
        for lo in edges[:-1]:
            hi = round(lo + 0.1, 1)
            m = ia & (z >= lo) & (z < hi)
            mt = m & (htar == 1)
            mc = m & (comp == 1)
            med = float(np.median(hz[mt])) if mt.sum() else float("nan")
            fo.write(f"{lo},{hi},{int(m.sum())},{int(mt.sum())},{int(mc.sum())},"
                     f"{'' if mt.sum()==0 else f'{med:.2f}'}\n")
            if m.sum():
                print(f"  {lo:.1f}-{hi:.1f}   {m.sum():6d}   {mt.sum():7d}    "
                      f"{mc.sum():6d}      {med:.2f}")
    print("table ->", out)


if __name__ == "__main__":
    main()
