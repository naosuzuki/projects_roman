#!/usr/bin/env python3
"""
25_case_summary.py -- summary of the Phase III host-recovery cases: fields,
coverage, and completeness fraction (hosts integrated to S/N=5) after 1, 2, and
3 visits. One visit = every field observed once for 1 h; a host lying in n
overlapping fields receives n hours per visit, so it completes when
(visits x n_cov) >= its remaining-hour requirement. Fractions are of all 2270
remaining hosts.

Reads 07_program_sne_ELAIS-N1.csv.   python 25_case_summary.py
"""
import os
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "outputs", "csv", "07_program_sne_ELAIS-N1.csv")
RA0, DEC0 = 242.498, 54.497
FOV = 1.3
R = FOV / 2.0
SP = np.sqrt(3) * R
M = 12


def hexagon(xc, yc):
    a = np.radians(np.arange(0, 360, 60))
    return np.column_stack([xc + R * np.cos(a), yc + R * np.sin(a)])


def flower(rot=0.0):
    c = [(0.0, 0.0)]
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang + rot); c.append((SP * np.cos(a), SP * np.sin(a)))
    return c


def mosaic19():
    c = flower(0.0)
    for ang in (90, 150, 210, 270, 330, 30):
        a = np.radians(ang); c.append((2 * SP * np.cos(a), 2 * SP * np.sin(a)))
    for ang in (60, 120, 180, 240, 300, 0):
        a = np.radians(ang); c.append((np.sqrt(3) * SP * np.cos(a), np.sqrt(3) * SP * np.sin(a)))
    return c


def ring(rr, phase=0.0):
    return [(rr * np.cos(np.radians(phase + 360 * k / M)),
             rr * np.sin(np.radians(phase + 360 * k / M))) for k in range(M)]


def central4():
    return [(0.60 * np.cos(np.radians(45 + 90 * k)), 0.60 * np.sin(np.radians(45 + 90 * k))) for k in range(4)]


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    htar = np.asarray(d["host_target"]).astype(int)
    comp = np.asarray(d["host_completed"]).astype(int)
    hra = np.asarray(d["host_ra"]).astype(float)
    hdec = np.asarray(d["host_dec"]).astype(float)
    vneed = np.asarray(d["host_visits_needed"]).astype(float)
    vdone = np.asarray(d["host_visits_done"]).astype(float)
    cosd = np.cos(np.radians(DEC0))
    xi = (hra - RA0) * cosd; eta = hdec - DEC0
    rneed = np.maximum(1, np.ceil(vneed - vdone)).astype(int)
    rem = (htar == 1) & (comp == 0)
    pts = np.column_stack([xi, eta])[rem]
    rn = rneed[rem]
    n = len(pts)

    CASES = [
        ("1  C    (mosaic)",     mosaic19()),
        ("3  E    (ring 1.6)",   ring(1.6)),
        ("4  E+F  (ring 1.5)",   ring(1.5, 0) + ring(1.5, 15)),
        ("5  E+F+G (28 fields)", ring(1.5, 0) + ring(1.5, 15) + central4()),
        ("6  E+G  (16 fields)",  ring(1.5, 0) + central4()),
    ]
    print(f"remaining hosts: {n}\n")
    print("Case Config             Fields  Cover   1visit  2visit  3visit")
    for lab, cens in CASES:
        ncov = np.zeros(n, int)
        for xc, yc in cens:
            ncov += Path(hexagon(xc, yc)).contains_points(pts).astype(int)
        cover = (ncov > 0)
        c1 = (1 * ncov >= rn) & cover
        c2 = (2 * ncov >= rn) & cover
        c3 = (3 * ncov >= rn) & cover
        print(f"{lab:22s} {len(cens):4d}   {100*cover.mean():4.0f}%   "
              f"{100*c1.mean():4.0f}%   {100*c2.mean():4.0f}%   {100*c3.mean():4.0f}%")


if __name__ == "__main__":
    main()
