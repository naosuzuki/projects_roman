# projects_roman/programs_aat — Roman HLTDS south (EDF-S) follow-up with AAT 2dF/AAOmega

Repeat the Subaru/PFS + ELAIS-N1 study (`../programs_pfs/`, scripts 01–32) for
the **AAT 2dF/AAOmega** spectroscopic follow-up of the **Euclid Deep Field South
(EDF-S)**, the southern HLTDS deep field. **No HSC/imaging phase** — only
Phase II (live SNe) and Phase III (hosts).

Paper: `~/github/papers/26_roman_aat_study/ms.tex` (Overleaf, branch **main**,
a copy of the Subaru paper being converted section by section; HSC section
removed; text below the new visibility section is still the PFS/ELAIS-N1
template). Sections now: Intro → AAT visibility → AAT success model → EDF-S
targets → 2dF pointings → host campaign → Summary. The old PFS/ELAIS-N1
sections sit in a `\begin{comment}` block at the end for reference. Local TeX
lacks revtex4-1, so compile on Overleaf.

## Field / site / instrument

- EDF-S center (HOURGLASS2 footprint): RA ≈ 62.6°, Dec ≈ −48.2°; SNe span
  ~6.4° × 5° (two Roman pointings). 18,649 SNe Ia in the sim.
- AAT 3.9 m, Siding Spring: lat −31.2754°, lon 149.0672°, 1164 m. AEST = UTC+10.
- 2dF: 2° diameter FoV (~3.1 deg²), ~400 fibers (392 science + 8 guide),
  ~2.1″ fibers. AAOmega: dual beam, ~370–950 nm, R ~ 1300–1700 (low-res gratings).
  OzDES (DES SN program) is the direct precedent; `../data/AATmaglimit_v2.dat`
  is an AAT redshift-efficiency vs mag table used in the old white paper.
- ELAIS-N1 is invisible from the AAT (max elevation 4°).

## Data (all local — no external disk used)

- `../programs_pfs/outputs/csv/03_snia_radec.csv` (field, SNID, RA, Dec, z_cmb,
  peakmag_Z) — includes the 18,649 EDF-S SNe Ia. **peakmag_Z = −9 is a sentinel**
  (no Z peak; 2011 SNe) — treated as faint.
- EDF-S has **no host mags / light curves / CC** locally. Host Z087 is assigned
  by drawing an ELAIS-N1 SN Ia host at the same z (|dz|<0.05) from
  `03_snia_host_radec_ELAIS-N1.csv` (same HOURGLASS2 Ia model; hostless carried
  over). Replace with real EDF-S hosts if the FITS
  (`/Volumes/exdisk1/data/Roman/ltcv_sim/EDF-S/`) become available.
- AAT success calibration: `~/github/projects_hsc/specz/csvfiles_aat/{g,r,i,z,y}_magfraction.csv`
  (HSC-AAT host campaign, paper `~/github/papers/24_hscsnhost_aat` Fig. 2).

## Scripts

| AAT | PFS analog | Result |
|---|---|---|
| `01_aat_visibility.py` | `05` | EDF-S transit airmass 1.05; 319 nights, 1327 dark h/yr; best Sep–Jan, none in May |
| `02_aat_success.py` | `01/02` ETC | logistic fit to HSC-AAT success at 2 h: Z-band m50=21.30, w=0.97; depth ∝ 1.25 log10(t/2h). t50 = 1.2/7.3/46 h at Z=21/22/23 |
| `03_edfs_targets.py` | `03/09` | DEEP disc (63.596,−47.599; r 1.42°; 11,346 SNe; z_med 1.42) + WIDE (58.915,−49.317; r 1.81°; 7,303; z_med 1.17). Peak Z<22/23/24: 337/1454/5027 |
| `04_2df_pointings.py` | `06` | 2dF R=1.05°: DEEP 4 + WIDE 5 = 9 fields, 100% coverage, 1.9 fields/SN |
| `05_aat_host_yield.py` | `07/24/25` | 9,230 hosts Z<24; 1-h visits, 360 fibers, retire 20 h: 3,230/4,050/4,650 secure at 90/180/360 h (51% of targets, 25% of all SNe; 5% at z>1.5) |

Not yet done (need EDF-S light curves): obs calendar (`04`), Phase II live-SN
fiber budget with time windows (`07`), CC SNe/hosts.

Style: match `../programs_pfs` (Times New Roman, capitalized labels, rainbow_r).
