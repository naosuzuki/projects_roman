#!/usr/bin/env python3
"""
31_phase3_eg_w80.py -- 80%-weather version of the Phase III E+G ledger (Table 9):
the main survey completes only 80% of the hosts (Ia 1329, CC 842, TDE 2,
matching the Net-80% column of the host-yield table), so the Phase III target
list is larger; the E+G configuration (Case 6) is then applied to it.

Model (deterministic, reproducible):
  * Completed(80%) = round(0.8 x Completed) per class. The demoted hosts are the
    weather-fragile ones -- those with the LARGEST integration requirement
    (visits_needed), which depend on the most scheduled hours and so lose the
    most to canceled visits. They re-enter the pool as observed-but-below-S/N=5
    hosts with remaining need = ceil(0.2 x visits_needed) (they received 80% of
    the required time).
  * Unobserved hosts and originally under-integrated hosts are unchanged
    (the latter with remaining need ceil(visits_needed - 0.8 x visits_done)).
  * Phase III visits are counted per executed visit of all 16 E+G fields
    (weather stretches the calendar, not the per-visit yield).

Reads 07_program_sne_ELAIS-N1.csv.   python 31_phase3_eg_w80.py
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
W = 0.8
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
    pts = np.column_stack([(hra - RA0) * cosd, hdec - DEC0])
    n = len(pts)

    ncov = np.zeros(n, int)
    for xc, yc in eg_centers():
        ncov += Path(hexagon(xc, yc)).contains_points(pts).astype(int)

    # demote the weather-fragile completed hosts (largest vneed) to hit 0.8x
    comp80 = comp.copy()
    demoted = np.zeros(n, bool)
    for c in CLASSES:
        mc = np.where(comp & (typ == c))[0]
        n80 = int(round(W * len(mc)))
        ndem = len(mc) - n80
        order = mc[np.argsort(-vneed[mc], kind="stable")]   # largest need first
        dem = order[:ndem]
        comp80[dem] = False
        demoted[dem] = True

    # remaining need per host in the 80%-weather Phase III pool
    rneed = np.full(n, np.nan)
    rneed[htar & ~start] = vneed[htar & ~start]                            # unobserved
    m_sb = htar & start & ~comp                                            # originally below S/N=5
    rneed[m_sb] = np.maximum(1, np.ceil(vneed[m_sb] - W * vdone[m_sb]))
    rneed[demoted] = np.maximum(1, np.ceil((1 - W) * vneed[demoted]))      # weather-lost

    rem = htar & ~comp80
    out = os.path.join(CSV_DIR, "31_phase3_eg_w80_ELAIS-N1.csv")
    print("Class  hostZ<25.5  compl80  unobs  SN<5  remaining  P3visit1  P3visit2")
    with open(out, "w") as fo:
        fo.write("class,N_host_zcut,N_completed80,N_unobserved,N_started_below,"
                 "N_remaining,EG_visit1,EG_visit2\n")
        tots = np.zeros(7, int)
        for c in CLASSES:
            m = htar & (typ == c)
            npool, nc = int(m.sum()), int((m & comp80).sum())
            r = m & ~comp80
            nun = int((r & ~start).sum())
            nsb = int((r & start).sum())          # incl. demoted (they were started)
            v1 = int((r & (1 * ncov >= rneed)).sum())
            v2 = int((r & (2 * ncov >= rneed)).sum())
            tots += (npool, nc, nun, nsb, nun + nsb, v1, v2)
            fo.write(f"{c},{npool},{nc},{nun},{nsb},{nun+nsb},{v1},{v2}\n")
            print(f"{c:5s}  {npool:9d}  {nc:7d}  {nun:5d}  {nsb:4d}  "
                  f"{nun+nsb:9d}  {v1:8d}  {v2:8d}")
        fo.write("total," + ",".join(str(t) for t in tots) + "\n")
        print("total  " + "  ".join(str(t) for t in tots))
    print("table ->", out)


if __name__ == "__main__":
    main()
