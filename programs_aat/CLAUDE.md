# projects_roman/programs_aat — Roman HLTDS south (EDF-S) follow-up with AAT 2dF/AAOmega

Repeat the Subaru/PFS + ELAIS-N1 study (`../programs_pfs/`, scripts 01–32) for
the **AAT 2dF/AAOmega** spectroscopic follow-up of the **Euclid Deep Field South
(EDF-S)**, the southern HLTDS deep field. **No HSC/imaging phase** — only
Phase II (live SNe) and Phase III (hosts).

Paper: `~/github/papers/26_roman_aat_study/ms.tex` (Overleaf, branch **main**,
a copy of the Subaru paper being converted section by section; HSC section
removed; text below the new visibility section is still the PFS/ELAIS-N1
template with PFS→AAT, ELAIS-N1→EDF-S substituted — **numbers are ELAIS-N1's
until each section is redone**).

## Field / site / instrument

- EDF-S center (HOURGLASS2 footprint): RA ≈ 62.6°, Dec ≈ −48.2°; SNe span
  ~6.4° × 5° (two Roman pointings). 18,649 SNe Ia in the sim.
- AAT 3.9 m, Siding Spring: lat −31.2754°, lon 149.0672°, 1164 m. AEST = UTC+10.
- 2dF: 2° diameter FoV (~3.1 deg²), ~400 fibers (392 science + 8 guide),
  ~2.1″ fibers. AAOmega: dual beam, ~370–950 nm, R ~ 1300–1700 (low-res gratings).
  OzDES (DES SN program) is the direct precedent; `../data/AATmaglimit_v2.dat`
  is an AAT redshift-efficiency vs mag table used in the old white paper.
- ELAIS-N1 is invisible from the AAT (max elevation 4°).

## Data

- Committed: `../programs_pfs/outputs/csv/03_snia_radec.csv` (field, SNID, RA,
  Dec, z_cmb, peakmag_Z) — includes EDF-S. **peakmag_Z = −9 is a sentinel**
  (no Z data) — exclude before applying magnitude cuts.
- Light curves, hosts, CC/TDE need the HOURGLASS2 FITS on
  `/Volumes/exdisk1/data/Roman/ltcv_sim/EDF-S/` (`.FITS.gz`).

## Script map (PFS → AAT)

| AAT | PFS source | Status |
|---|---|---|
| `01_aat_visibility.py` | `05_visibility.py` | done: 319 nights, 1327 h/yr, best Sep–Jan, transit airmass 1.05 |
| `02_aat_etc.py` | `01/02` (PFS ETC) | todo — AAOmega S/N model (no official ETC wrapper) |
| `03_snia_radec.py` | `03` | todo |
| `04_obs_calendar.py` | `04` | todo — EDF-S cadence, AAT season Jul–Mar |
| `05_2df_pointings.py` | `06` | todo — 2° circular FoV tiling |
| `06_fiber_budget.py` | `07` | todo — 392 fibers; needs EDF-S light curves |
| later | `09–32` | z-histograms, spec-z yields, remaining hosts, configs, cases |

Style: match `../programs_pfs` (Times New Roman, capitalized labels, rainbow_r).
