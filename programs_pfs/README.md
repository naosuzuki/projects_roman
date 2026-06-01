# programs_pfs

Code for Subaru **PFS** (Prime Focus Spectrograph) exposure-time / S/N studies,
in support of the Roman–Subaru SNIa program (Roman HLTDS transients followed up
with PFS spectroscopy, primarily in the northern **ELAIS-N1** field).

## Dependencies: the PFS ETC

These scripts drive the official **PFS Exposure Time Calculator**
([Subaru-PFS/spt_ExposureTimeCalculator](https://github.com/Subaru-PFS/spt_ExposureTimeCalculator)),
installed separately. On this machine it lives at
`~/github/spt_ExposureTimeCalculator` and was built with Homebrew `gcc-15`
(needed for the OpenMP C core):

```bash
git clone https://github.com/Subaru-PFS/spt_ExposureTimeCalculator.git
cd spt_ExposureTimeCalculator
CC=gcc-15 CXX=g++-15 python3 -m pip install .   # provides the `pfs-run-etc` CLI
```

Point elsewhere via the `PFS_ETC_DIR` environment variable if you relocate it.

## Files

| File | Purpose |
|---|---|
| `01_pfs_etc.py` | Reusable wrapper: `run_etc()` runs the ETC; `band_sn()` summarizes continuum S/N within a g/r/i/z/y filter; `band_resolution()` gives R for a pixel binning. Loaded by path (leading digit isn't import-safe). |
| `02_plot_time_vs_sn.py` | Driver: continuum S/N vs total exposure time, per band, over a magnitude grid. |
| `03_snia_radec.py` | Sky positions (RA/Dec) of all HOURGLASS2 Type Ia SNe (two fields), colored by redshift; also writes a positions catalog CSV. |
| `outputs/png/`, `outputs/csv/` | Generated figures and tables. |

## Usage

```bash
cd ~/github/projects_roman/programs_pfs
python 02_plot_time_vs_sn.py                       # i-band, linear, default grid
python 02_plot_time_vs_sn.py --loglog              # log-log version (-> *_loglog.png)
python 02_plot_time_vs_sn.py --loglog --from-csv   # restyle from cached grid (no ETC re-run)
python 02_plot_time_vs_sn.py --loglog --from-csv --bin 5   # 5-pixel binning (S/N x sqrt5, R~1840)
python 02_plot_time_vs_sn.py --band z --mags 23 24 25
python 03_snia_radec.py                            # SNIa sky positions, both fields
```

Curves are color-coded as a **red→blue rainbow from bright (red) to faint (blue)**
via `--cmap` (default `rainbow_r`). `--from-csv` replots from the cached
`outputs/csv/time_vs_sn_<band>.csv` so you can iterate on styling instantly.
The linear plot keeps a small log-log orientation inset; `--loglog` makes the
whole panel log-log and writes a separate `*_loglog.png` (so both are kept).

## Notes / conventions

- "S/N in i-band" = **median continuum S/N per PFS pixel** in 700–830 nm
  (red arm; effective wavelength ~770 nm). Bin to a resolution element
  (~2–3 pix) or use `band_sn(..., stat="quad")` for the broadband value.
- PFS is sky-background-limited in the optical, so S/N ∝ √t.
- Default conditions = ETC defaults: seeing 0.8″, zenith 45°, new moon,
  field edge (0.675°), point-like source `REFF=0.3″`. A single frame is 450 s.
- Current model uses a **flat fν continuum** at the given AB mag; feeding real
  SN Ia / host-galaxy SEDs (e.g. from the HOURGLASS2 sim) is a planned next step.
