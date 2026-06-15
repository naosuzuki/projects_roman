# projects_roman — Roman–Subaru SNe Ia follow-up planning

Working repo for planning **Subaru spectroscopic + imaging follow-up of the
Nancy Grace Roman HLTDS Type Ia supernova program**, focused on the northern
**ELAIS-N1** deep field. Most work lives in `programs_pfs/`. A companion paper
is drafted in `~/github/papers/26_roman_subaru_study/ms.tex` (Overleaf-linked).

Git remote: `origin = https://github.com/naosuzuki/projects_roman.git`, branch `main`.

---

## The big picture / goal

A **three-phase** observing program for ELAIS-N1, with a PFS **fiber-budget**
analysis as the central deliverable:

- **Phase I — Subaru/HSC imaging** (g, r2, i2): detect + photometrically
  monitor transients, measure host photometry, define target lists and the
  Roman Z-band magnitude cuts.
- **Phase II — PFS live-SN spectroscopy**: observe SNe Ia + core-collapse while
  their Roman Z087 light curve is brighter than Z = 24.
- **Phase III — PFS host-galaxy spectroscopy**: once a SN fades past Z = 24,
  observe its host if Z < 25, integrating to S/N = 5 at 10-pixel binning.

Adopted survey: 7-pointing hexagonal PFS "flower"; **2 visits/month** (config A
= 0°, config B = 30° rotation) in Subaru-observable months (Feb–Aug); **1 h per
pointing = 7 h per configuration**; **2028-06-01 → 2030-05-31** (29 visits; sim
timeline shifted −214 d so the main survey starts 2028-06-01, reflecting the
advanced 2026-08-30 Roman launch).

PFS fiber budget per pointing: **2394 total = 400 sky + 200 FluxStd + demand +
spare ⇒ 1794 available** for science targets.

---

## Data (external disk — must be mounted)

HOURGLASS2 SNANA simulation of the Roman HLTDS (PIPPIN / R. Kessler group),
on the external disk **`/Volumes/exdisk1/data/Roman/ltcv_sim/`**:

- `ELAIS-N1/` (north, Subaru field, gunzipped `.FITS`), `EDF-S/` (south,
  `.FITS.gz`), `original/` (un-split). astropy reads both; **skip macOS `._`
  AppleDouble files**.
- SNANA HEAD/PHOT/SPEC triplets, 9 transient models. Key models:
  `SNIaMODEL00` (Ia), `NONIaMODEL06` (CC = II/Ib/Ic), `NONIaMODEL02` (TDE).
  Roman bands `RZYJHFK`; Z-band = **Z087** (~877 nm).
- If `/Volumes/exdisk1` is missing: it's a USB WD "My Book"; replug the cable
  (bypass any USB hub), then `diskutil list` → it auto-mounts (was `/dev/disk4`).
  Note: the **per-field catalog CSVs are committed** in `programs_pfs/outputs/csv/`,
  so most analysis can proceed even with the disk offline (only re-reading the
  raw FITS needs it).

Full sim details: see memory file `project_roman_hourglass2_sim` and the data's
own `.README`.

---

## PFS ETC tool (external dependency)

Installed at **`~/github/spt_ExposureTimeCalculator`** (Subaru-PFS/spt_ExposureTimeCalculator),
built with **Homebrew gcc-15** for the OpenMP C core (Apple clang has no OpenMP):

```bash
git clone https://github.com/Subaru-PFS/spt_ExposureTimeCalculator.git
cd spt_ExposureTimeCalculator
CC=gcc-15 CXX=g++-15 python3 -m pip install .     # provides `pfs-run-etc`
```

`01_pfs_etc.py` wraps it (`run_etc`, `band_sn`); ~3 s/run (16 threads). Override
its location with the `PFS_ETC_DIR` env var. Python: conda base (`python3`),
astropy 7.2. Fonts: Times New Roman (System/Library/Fonts/Supplemental).

---

## Scripts (`programs_pfs/`, numbered; `02`+ load `01` by file path)

| File | What it does |
|---|---|
| `01_pfs_etc.py` | ETC wrapper: `run_etc()`, `band_sn(res,band)`, `band_resolution()`, `BANDS` (g/r/r2/i/i2/z/y), `ARM_DISP_NM`. |
| `02_plot_time_vs_sn.py` | PFS continuum S/N vs exposure time. `--band {g,r2,i,i2,z,…} --loglog --from-csv --bin M --vlines …`. √M binning with effective R. |
| `03_snia_radec.py` | RA/Dec of transients or hosts. `--model SNIaMODEL00/NONIaMODEL06/NONIaMODEL02 --tag --target sn/host --zmag-cut --color-by redshift/zmag --msize --alpha`. |
| `04_obs_calendar.py` | Roman cadence as calendar dates. `--shift-to-date 2028-06-01 --main-start 2029-01-01 --shade-visible --vis-hours 2 --label`. |
| `05_visibility.py` | ELAIS-N1 visibility from Subaru (elev>30, astro night). Monthly bars split by airmass (2.0–1.8 red / 1.8–1.5 green / ≤1.5 blue). |
| `06_pfs_pointings.py` | 7-pointing hex flower A (`--rotate 0`) / B (`--rotate 30`); per-pointing counts + A/B/A∪B coverage for SNe/CC/TDE/Ia-host/CC-host. |
| `07_fiber_budget.py` | **End-to-end fiber-budget sim** (the main result). Reads Z087 light curves, simulates visits, SN→host lifecycle, per-pointing fibers needed vs 1794, spare, remaining hosts. |
| `08_hsc_etc.py` | **Subaru/HSC imaging ETC** (g/r2/i2): point-source S/N vs exposure time for mag 24–28. CCD-equation model tuned to HSC-SSP depths. |

Outputs in `programs_pfs/outputs/{png,csv}/`. README in `programs_pfs/README.md`.

**Plot style convention** (match across figures): Times New Roman; axis labels
~17–18 pt with **first letter of each word capitalized**; `rainbow_r` colormap
(bright = red → faint = violet); log-log with plain-number ticks; 60/90/120-min
dotted guides; √t reference dashed.

---

## Key results

- **ETC**: PFS sky-bkg-limited (S/N ∝ √t). z-band (Z087) 1 h S/N/pix ≈ 2.0/1.3/0.8
  at z = 22/22.5/23. Host integration to S/N=5 at 10-pix: ~3/15/38/94 one-hour
  visits at Z = 23/24/24.5/25 → faint hosts are expensive.
- **HSC imaging** 5σ point-source depth (20 min): g≈27.2, r2≈26.7, i2≈26.1.
- **Targets** (ELAIS-N1): 13,094 Ia + 18,774 CC; Z<24 → 3,864 Ia, 2,862 CC; TDE tiny (21).
- **Pointing coverage** A∪B ≈ 76% (all SNe), ~70% for the Z<24 observable subset.
- **Fiber budget** (07): 5,860 program SNe; demand peaks ~130–162 fibers/pointing
  (host-dominated, season-opening backlog), mean ~75–86; **min spare 1632, mean
  ~1715 ⇒ ≳90% spare** for collaborator filler targets. 4,754 Z<25 hosts: 2,713
  completed, **2,041 remain incomplete** (mostly faint Z>24). Bottleneck =
  faint-host integration time, not fibers.
- ELAIS-N1 transits at elev 55.3° = **airmass 1.22** from Subaru (never <1.2);
  observable Feb–Aug, prime Apr–Jul (~4–8 h/night).

---

## Paper (`~/github/papers/26_roman_subaru_study/ms.tex`)

AASTeX, **Overleaf-linked** (remote `git.overleaf.com/...`, branch `master`).
Workflow: edit `ms.tex`, copy figures into `figures/`, `git add`,
`git pull --no-edit origin master` (merge Overleaf edits), `git push`.
Section order: Intro (3-phase framing) → **Phase I: HSC Imaging** → PFS ETC
(z-band; Phases II–III) → Spatial distribution → Schedule (+ shifted-to-2028
fig) → Pointing config A/B + coverage table → Visibility → Fiber budget →
Summary → Appendix A (i-band ETC) [→ Appendix B (HSC binning) — TO BE REMOVED].

---

## ⚠️ PENDING — pick up here

The **HSC ETC was redone** as a true *imaging* ETC (`08_hsc_etc.py` →
`08_hsc_etc_{g,r2,i2}band.png`, mag 24–28, PFS-aligned style, capitalized
labels/title). The paper still references the **old (mistaken)
PFS-in-HSC-band** figures (`time_vs_sn_{g,r2,i2}band_loglog*.png`). To finish:

1. Copy `08_hsc_etc_{g,r2,i2}band.png` into the paper `figures/`.
2. In `ms.tex` Phase-I section (`fig:hsc`): swap the 3 panels to the new
   `08_hsc_etc_*` figures; update the Phase-I text to imaging depths (5σ:
   g≈27.2, r2≈26.7, i2≈26.1 at 20 min; point source, native 1×1 pixels).
3. **Remove Appendix B** ("HSC-band binning", labels `app:hscbin`,
   `fig:hscbin_{g,r2,i2}`) — no binning for imaging.
4. Optionally `git rm` the unused `figures/time_vs_sn_{g,r2,i2}band*.png` from
   the paper.
5. Commit/push both repos.

Then the likely next direction is feeding **real SN Ia / host SEDs** into the
ETC (replacing the flat-fν continuum) and refining the cadence-aware yields.

---

## Resume checklist (on the new machine)

```bash
cd ~/github/projects_roman && git pull            # latest scripts + outputs
# mount /Volumes/exdisk1 (USB) if you need the raw FITS
# ensure ~/github/spt_ExposureTimeCalculator is built (see above) for ETC runs
cd programs_pfs && python3 08_hsc_etc.py           # e.g. regenerate HSC ETC
```

Paper: `cd ~/github/papers/26_roman_subaru_study && git pull origin master`.
