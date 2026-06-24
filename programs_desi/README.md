# programs_desi

A **DESI (Mayall 4-m, Kitt Peak) case** for the Phase III host-galaxy
spectroscopic follow-up of the Roman HLTDS supernova program in ELAIS-N1 —
the DESI counterpart to the Subaru/PFS study in `../programs_pfs/`.

**Start here:** read [`CLAUDE.md`](CLAUDE.md) — it has the goal, the data to
reuse (`../programs_pfs/outputs/csv/07_program_sne_ELAIS-N1.csv`, the 2270
remaining hosts), DESI vs PFS trade-offs (4.2× aperture penalty, ~3–4 fields vs
12–28, transit airmass 1.08 from KPNO), the observing model, and the suggested
scripts (`01_desi_visibility` … `05_desi_vs_pfs`).

Goal: a quantitative comparison — coverage, # pointings, integration time, and
total host spec-z yield — of DESI vs the PFS "Cases 1–6".

Outputs in `outputs/{png,csv}/`.
