#!/usr/bin/env python3
"""
29_host_remaining_table.py -- per-class ledger of the Phase III target list:
host galaxies (Z < 25.5) of program SNe that, at the end of the main 2028-2030
survey (the Phase II campaign), are still UNOBSERVED (never received a fiber
hour) or OBSERVED BUT BELOW S/N = 5 (started, unfinished).

Reads the program-SN catalog written by 07_fiber_budget.py.

    python 29_host_remaining_table.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(HERE, "outputs", "csv")
CSV = os.path.join(CSV_DIR, "07_program_sne_ELAIS-N1.csv")
CLASSES = ["Ia", "CC", "TDE"]


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    typ = np.asarray(d["type"]).astype(str)
    htar = np.asarray(d["host_target"]).astype(int) == 1
    comp = np.asarray(d["host_completed"]).astype(int) == 1
    start = np.asarray(d["host_started"]).astype(int) == 1

    out = os.path.join(CSV_DIR, "29_host_remaining_ELAIS-N1.csv")
    print("Class  hostZ<25.5  completed  remaining  unobserved  started_SN<5")
    with open(out, "w") as fo:
        fo.write("class,N_host_zcut,N_completed,N_remaining,N_unobserved,N_started_below\n")
        tots = np.zeros(5, int)
        for c in CLASSES:
            m = htar & (typ == c)
            n = int(m.sum())
            nc = int((m & comp).sum())
            rem = m & ~comp
            nun = int((rem & ~start).sum())          # never received a fiber hour
            nsb = int((rem & start).sum())           # started but S/N < 5
            nr = nun + nsb
            tots += (n, nc, nr, nun, nsb)
            fo.write(f"{c},{n},{nc},{nr},{nun},{nsb}\n")
            print(f"{c:5s}  {n:9d}  {nc:9d}  {nr:9d}  {nun:10d}  {nsb:12d}")
        fo.write(f"total,{tots[0]},{tots[1]},{tots[2]},{tots[3]},{tots[4]}\n")
        print(f"total  {tots[0]:9d}  {tots[1]:9d}  {tots[2]:9d}  "
              f"{tots[3]:10d}  {tots[4]:12d}")
    print("table ->", out)


if __name__ == "__main__":
    main()
