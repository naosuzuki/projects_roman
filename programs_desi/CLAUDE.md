# projects_roman/programs_desi — A DESI case for the Roman HLTDS host follow-up

Working sub-project to **make the case for using DESI (Mayall 4-m, Kitt Peak) as
an alternative to Subaru/PFS** for the **Phase III host-galaxy spectroscopic
follow-up** of the Nancy Grace Roman HLTDS Type Ia (and CC) supernova program in
the northern **ELAIS-N1** deep field.

The companion PFS study lives in `../programs_pfs/` (scripts 01–25) and the paper
`~/github/papers/26_roman_subaru_study/ms.tex`. This sub-project mirrors that
analysis for DESI and compares the two instruments.

Git: same repo, `origin = https://github.com/naosuzuki/projects_roman.git`,
branch `main`. Outputs go in `programs_desi/outputs/{png,csv}/`.

---

## The goal / the case to make

After the main HSC imaging + PFS live-SN survey, **2270 SN-host galaxies remain
incomplete** (not integrated to S/N = 5) and define the Phase III target list.
The PFS study showed these are mostly **faint (Z > 24) hosts in the outer annulus
(r = 1.6–2.2°)**, expensive in integration time on a 1.3° FoV.

**Question:** can DESI — with its huge 8 deg² field and 5000 fibers — do this
host follow-up more efficiently than PFS, despite its smaller 4-m aperture?
Deliver a quantitative comparison: coverage, # pointings, integration time, and
total host spec-z yield, DESI vs the PFS "cases" (Cases 1–6 in the paper).

---

## Reuse this data (no external disk needed)

The single key input is **committed** in the sibling project:

```
../programs_pfs/outputs/csv/07_program_sne_ELAIS-N1.csv
```

One row per program SN. Relevant columns:
- `host_ra, host_dec` — host position (deg)
- `host_Z` — host Roman Z087 AB mag (proxy for brightness; see ETC caveat)
- `host_target` — 1 if host qualifies (Z < 25.5) for Phase III
- `host_visits_needed` — **PFS** 1-h visits to reach S/N = 5 at 10-pix binning
- `host_visits_done` — PFS hours already accumulated in the main survey
- `host_completed`, `host_started`
- `z` — host/SN redshift (SIM_REDSHIFT_CMB)

**Remaining hosts** = `host_target==1 & host_completed==0` → **2270**.
Their **PFS remaining requirement** `rneed = ceil(host_visits_needed -
host_visits_done)` (hours, ≥ 1). Convert to DESI hours via the aperture scaling
below (and refine with a real DESI ETC).

Field center: RA0 = 242.498°, Dec0 = +54.497°; use `cosd = cos(Dec0)` for the
tangent plane (`xi = (RA−RA0)·cosd`, `eta = Dec−Dec0`), as in `../programs_pfs`.

---

## DESI / Mayall — instrument facts (the proposal)

- **Telescope:** Mayall **4.0-m**, Kitt Peak National Observatory.
  Site: lat **+31.9633°**, lon **−111.5997°**, elev **2120 m**.
- **Focal plane:** **5000** robotic fiber positioners.
- **Field of view:** **3.2° diameter (≈8 deg²)**. Fiber aperture **1.5″** diameter.
- **Spectrograph:** 360–980 nm in three arms (B/R/Z); resolution **R ≈ 2000
  (blue) → 5100 (red)**. Good for host emission/absorption-line redshifts.
- **Design philosophy:** massively-parallel survey instrument (rapid re-pointing,
  thousands of targets per pointing).

DESI software for a real ETC (optional): `desimodel`, `specsim`, `desispec`
(desihub on GitHub). If not installed, scale from the PFS ETC (next section).

---

## DESI vs PFS — the trade-offs (grounding numbers already computed)

| | Subaru/PFS | DESI/Mayall |
|---|---|---|
| Aperture | 8.2 m | 4.0 m → **4.2× longer** integration for equal S/N |
| FoV | 1.3° (~1.3 deg²) | **3.2° (~8 deg²)** |
| Fibers | 2394 (1794 science) | **5000** |
| Fields to tile ELAIS-N1 | 12–28 | **≈3–4** (see below) |
| ELAIS-N1 transit | 55.3°, airmass 1.22 | **67.5°, airmass 1.08** (better) |

Computed facts to build on:
- ELAIS-N1 (Dec +54.5°) transits at **67.5° elevation (airmass 1.08)** from KPNO
  — *higher/better than from Subaru* (lat closer to the field's Dec).
- The remaining-host footprint is **~4.3° across** (r₉₉ = 2.15°). DESI's 3.2° FoV
  (r < 1.6°) centered on the field covers only **40%** of the remaining hosts —
  because the leftovers sit in the **outer annulus** (the inner ones were already
  completed by the main survey). So expect **≈3–4 DESI pointings** for the field,
  vs 12–28 PFS fields. **Still far fewer pointings, and all 2270 hosts fit in
  5000 fibers.**
- Aperture scaling: **DESI hours ≈ 4.2 × PFS hours** for the same host continuum
  S/N (first-order; the 1.5″ vs ~1.05″ fiber and throughput differences partly
  offset for extended hosts — refine with a real ETC).

**The crux of the case:** DESI trades *more integration time per host* (4× the
aperture penalty) for *far fewer pointings* (one DESI field observes thousands of
hosts at once). Whether that nets out favorably is exactly the analysis to do.

---

## The DESI observing model to build (mirror ../programs_pfs 07/24/25)

1. **Tile** the ELAIS-N1 host field with N DESI pointings (3.2° FoV). Find the
   minimum N for ~full coverage; record which hosts each covers.
2. One DESI pointing's 5000 fibers observe **all** hosts in its FoV at once.
   Tiles barely overlap, so per host n_cov ≈ 1 (unlike the overlapping PFS rings).
3. **Integrate** each tile to depth T (hours). A host completes when
   `T ≥ DESI_rneed` (≈ 4.2 × PFS rneed). Total time ≈ T × N (no re-pointing within
   a tile — every host integrates in parallel).
4. **Compare** to the PFS cases: host spec-z completeness vs total exposure hours
   (the analog of paper Table `tab:cases` / Fig `fig:casesummary`).

Note the structural difference: PFS gains from **overlap** (a host in n fields
gets n hours/visit); DESI gains from **multiplexing the whole field in one shot**.
The comparison should make this explicit.

---

## Suggested scripts (new numbering, this directory)

- `01_desi_visibility.py` — ELAIS-N1 from KPNO: elevation/airmass vs date,
  observable months, dark hours. Mirror `../programs_pfs/05_visibility.py`
  (swap the site to Kitt Peak: `EarthLocation.of_site("kpno")` or lat/lon above).
- `02_desi_pointings.py` — tile the field with the 3.2° DESI FoV; minimum
  pointings for coverage; per-host coverage. Mirror `06_pfs_pointings.py`
  (circular FoV instead of hexagon).
- `03_desi_etc.py` — host integration time to S/N = 5 at 10-pix: either run a
  real DESI ETC (specsim/desimodel) or scale the PFS ETC by the 4.2× aperture
  factor. Tabulate hours vs host Z (cf. PFS ~3/15/38/94 h at Z=23/24/24.5/25 →
  DESI ~13/63/160/395 h — these are the numbers to check/replace).
- `04_desi_recovery.py` — completeness (hosts to S/N=5) vs total exposure hours
  for the DESI tiling; overlay the PFS cases. Mirror `24/25`.
- `05_desi_vs_pfs.py` — the summary comparison figure/table for the paper.

---

## Key open questions to resolve

- **DESI host ETC:** what continuum S/N = 5 (10-pix) integration time vs host
  magnitude? Validate the 4.2× scaling against a real DESI ETC, and decide the
  reference band (DESI is optical 360–980 nm; we have Roman Z087 = host_Z, so
  adopt an optical-r proxy, e.g. r ≈ Z + color, or use Z directly as PFS did with
  i-band — state the assumption).
- **Faint-host reach:** on a 4-m, hosts at Z ≳ 24.5 need ~hundreds of hours — is
  the faint tail simply out of reach for DESI? (Likely yes; the case may be that
  DESI cheaply mops up the *bright* hosts the whole field at once.)
- **DESI operations / feasibility:** DESI is a survey instrument. Frame the case
  as a **dedicated deep-field campaign or ancillary/secondary program**, and note
  scheduling realities (dark time, the main DESI survey's footprint/cadence).
- **Sky coverage:** confirm Dec +54.5° is within DESI's accessible declination
  range and observable for a useful window from KPNO.

---

## Environment

- Python: conda base `python3`, `astropy` (coordinates/Time for KPNO visibility),
  `numpy`, `matplotlib`. Fonts: Times New Roman
  (`/System/Library/Fonts/Supplemental/`) — match the PFS plot style.
- Optional DESI tooling: `pip install specsim desimodel` (desihub) for a real ETC.
- No external disk needed — the host catalog CSV is committed in `../programs_pfs`.

Plot-style convention (match the PFS figures): Times New Roman; capitalized
axis labels; clear colors; per-hour / per-visit comparisons; cap exposure axes
at the affordable budget where relevant.

---

## Resume checklist (start the new session here)

```bash
cd ~/github/projects_roman/programs_desi
# 1. read this CLAUDE.md
# 2. sanity-check the input catalog:
python3 -c "import numpy as np; d=np.genfromtxt('../programs_pfs/outputs/csv/07_program_sne_ELAIS-N1.csv',delimiter=',',names=True,dtype=None,encoding='utf-8'); m=(d['host_target']==1)&(d['host_completed']==0); print('remaining hosts:',m.sum())"
# 3. build 01_desi_visibility.py (KPNO) and 02_desi_pointings.py (3.2 deg FoV tiling)
```

First deliverable: **how many DESI pointings cover ELAIS-N1, and the host-spec-z
completeness vs exposure hours, compared to the PFS Cases 1–6.**
