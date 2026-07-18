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
- **Phase II — PFS main survey (2028-06-01 → 2030-05-31)**: live SNe Ia + CC +
  TDE while Z087 < 24; when a SN fades, its host (Z < 25.5) gets fibers in the
  same visits, integrating to S/N = 5 at 10-pixel binning. (Paper convention:
  Phase II = the main-survey campaign incl. in-survey host obs.)
- **Phase III — dedicated host follow-up after the main survey**: the 2,270
  remaining hosts; baseline tiling **Configuration α** = E ring (12 pts,
  r=1.5°) + G central square (4), alternated with **β** (α rotated 15°).

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
| `08_hsc_etc.py` | **Subaru/HSC imaging ETC** (g/r2/i2): empirical sky-limited, anchored m5(1h)≈26. |
| `09–12` | Redshift histograms: Ia/CC SNe and hosts by Z-mag bins (blue/limegreen/red). |
| `13–17` | Successful spec-z N(z): Ia live/host (13/14), z-reach table (15), CC pair (16/17), 100/80/70% weather. |
| `18` | Remaining-hosts RA/Dec map colored by Z. |
| `19–25` | Phase III configs: C mosaic (19), E ring (20), E/F (21), radius opt →1.5° (22), G central fill (23), E+G=α (24), case summary (25). |
| `26–27` | Combined SN+host spec-z histogram (26), + [OII] line-z estimate (27). |
| `28` | HSC 7-pointing A/B (circular 1.77 deg², gap-free ring d=√3R=1.30°; 88%/88%/96%). |
| `29–31` | Phase III ledgers: remaining table (29), α visits (30), 80%-weather + α/β (31). |
| `32` | Configuration β (α+15°): coverage 92→97%, alternating recovery. |

Outputs in `programs_pfs/outputs/{png,csv}/`. README in `programs_pfs/README.md`.

**Plot style convention** (match across figures): Times New Roman; capitalized
axis labels, recently restyled figures use **24–28 pt labels/titles, 19–24 pt
ticks (inward), legends 13–20 pt placed to avoid data**; plain-number ticks on
log axes ("10", "100"); 60/90/120-min dotted guides; √t dashed (label "∝√t").
Workflow: numbered script in programs_pfs → PNG in outputs/png → copy to paper
figures/ → commit BOTH repos (paper: add → pull --no-rebase --no-edit → push).

---

## Key results

- **ETC**: PFS sky-limited (S/N∝√t); z-band 1 h S/N/pix ≈ 2.0/1.3/0.8 at
  Z=22/22.5/23; host S/N=5@10-pix needs ~3/15/38/94 h at Z=23/24/24.5/25.
- **Model (07)**: SNe Ia+CC+TDE live at Z<24; hosts Z<25.5. Program 5,867
  (3,091 Ia / 2,769 CC / 7 TDE). Live spec-z: visible 1,416/1,476/5 →
  observed 971/1,018/3 (success ~69%; ×0.8 weather → 777/814/2).
  Hosts 4,986 → completed 2,716 (Ia 1,661) → **remaining 2,270**
  (1,698 unobserved + 572 below S/N=5). Fibers never binding (min spare 1,611).
- **Phase III cases** (Table 8): C 19pt/100%cov; E 12pt; E+F r=1.5° no-spill;
  α=E+G 16 fields 92% cov (v1/v2/v3 = 934/1072/1160); β=α+15° → α∪β 97%,
  α,β alternation v2 = 1,590. 80% weather: completed 2,173, remaining 2,813,
  α/β visits 1,271/1,590, cumulative II+III 3,444/3,763 (Ia 82.2%).
- **HSC**: imaging depth m5(1h)≈26; 7-pointing A/B circles (d=1.30°) cover
  88%/88%, A∪B 96% of SNe, gap-free inside r=1.5°.
- **[OII] line-z**: ~80% of hosts star-forming → ~1,200 Ia-host line-z in
  ~1 h each; faint-tail continuum cost mostly avoidable (not yet in cost model).
- ELAIS-N1 from Subaru: transit 55.3° (X=1.22), observable Feb–Aug.

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

- **Two-tier host cost model**: fold [OII] line-z (1 h for ~80% SF hosts) into
  07 and re-run Cases 1–6 — would slash Phase III costs (script 27 has f_SF).
- **Figs 3 & Appendix i-band** (02 outputs) still on old 17 pt style; Fig 4
  (z-band bins) already restyled (24 pt, split legends, "∝√t").
- **CC spec-z histograms** (16/17) still "Number of CC..." y-label; Ia set uses
  "Expected Number of Successful Spec-z".
- **programs_desi/**: DESI case handover in programs_desi/CLAUDE.md; a separate
  session drafted scripts 04–09 + outputs (committed as WIP).
- Paper section order/labels shifted a lot this cycle — trust 
ef, and always
  `git pull --no-rebase` the paper first (user edits on Overleaf; binary PNG
  conflicts have deleted figures before — regenerate from programs_pfs).

---

## Resume checklist (on the new machine)

```bash
cd ~/github/projects_roman && git pull            # latest scripts + outputs
# mount /Volumes/exdisk1 (USB) if you need the raw FITS
# ensure ~/github/spt_ExposureTimeCalculator is built (see above) for ETC runs
cd programs_pfs && python3 08_hsc_etc.py           # e.g. regenerate HSC ETC
```

Paper: `cd ~/github/papers/26_roman_subaru_study && git pull origin master`.
