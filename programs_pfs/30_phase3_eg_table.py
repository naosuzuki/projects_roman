#!/usr/bin/env python3
"""
30_phase3_eg_table.py -- Phase III ledger for the E+G configuration (Case 6:
12-pointing ring at r=1.5 deg + 4 central fields = 16 fields): per transient
class, the Table-5 remaining-host decomposition PLUS the number of host spec-z
completed (reaching S/N=5) after Phase III visits 1, 2, and 3, where one visit
observes every E+G field once for 1 h and a host in n overlapping fields gains
n hours per visit.

Reads 07_program_sne_ELAIS-N1.csv.   python 30_phase3_eg_table.py
"""
import os
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(HERE, "outputs", "csv")
CSV = os.path.join(CSV_DIR, "07_program_sne_ELAIS-N1.csv")
RA0, DEC0 = 242.498, 54.497
FOV = 1.3
R = FOV / 2.0
M = 12
R_RING = 1.5
CLASSES = ["Ia", "CC", "TDE"]


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def eg_centers():
    ring = [(R_RING * np.cos(np.radians(360 * k / M)),
             R_RING * np.sin(np.radians(360 * k / M))) for k in range(M)]
    central = [(0.60 * np.cos(np.radians(45 + 90 * k)),
                0.60 * np.sin(np.radians(45 + 90 * k))) for k in range(4)]
    return ring + central


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    typ = np.asarray(d["type"]).astype(str)
    htar = np.asarray(d["host_target"]).astype(int) == 1
    comp = np.asarray(d["host_completed"]).astype(int) == 1
    start = np.asarray(d["host_started"]).astype(int) == 1
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    vneed = np.asarray(d["host_visits_needed"]).astype(float)
    vdone = np.asarray(d["host_visits_done"]).astype(float)
    cosd = np.cos(np.radians(DEC0))
    xi = (hra - RA0) * cosd
    eta = hdec - DEC0
    rneed = np.maximum(1, np.ceil(vneed - vdone))
    rem = htar & ~comp
    pts = np.column_stack([xi, eta])

    ncov = np.zeros(len(pts), int)
    for xc, yc in eg_centers():
        ncov += Path(hexagon(xc, yc)).contains_points(pts).astype(int)

    out = os.path.join(CSV_DIR, "30_phase3_eg_ELAIS-N1.csv")
    print("Class  hostZ<25.5  completed  unobs  SN<5  remaining  P3visit1  P3visit2  P3visit3")
    with open(out, "w") as fo:
        fo.write("class,N_host_zcut,N_completed,N_unobserved,N_started_below,"
                 "N_remaining,EG_visit1,EG_visit2,EG_visit3\n")
        tots = np.zeros(8, int)
        for c in CLASSES:
            m = htar & (typ == c)
            n, nc = int(m.sum()), int((m & comp).sum())
            r = m & ~comp
            nun, nsb = int((r & ~start).sum()), int((r & start).sum())
            v1 = int((r & (1 * ncov >= rneed)).sum())
            v2 = int((r & (2 * ncov >= rneed)).sum())
            v3 = int((r & (3 * ncov >= rneed)).sum())
            tots += (n, nc, nun, nsb, nun + nsb, v1, v2, v3)
            fo.write(f"{c},{n},{nc},{nun},{nsb},{nun+nsb},{v1},{v2},{v3}\n")
            print(f"{c:5s}  {n:9d}  {nc:9d}  {nun:5d}  {nsb:4d}  {nun+nsb:9d}"
                  f"  {v1:8d}  {v2:8d}  {v3:8d}")
        fo.write("total," + ",".join(str(t) for t in tots) + "\n")
        print("total  " + "  ".join(f"{t}" for t in tots))
    print("table ->", out)


if __name__ == "__main__":
    main()
