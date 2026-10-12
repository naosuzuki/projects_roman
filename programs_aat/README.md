# programs_aat

AAT 2dF/AAOmega spectroscopic follow-up of the Roman HLTDS southern deep field
(EDF-S): the AAT counterpart of the Subaru/PFS ELAIS-N1 study in `../programs_pfs/`.
Paper: `~/github/papers/26_roman_aat_study/ms.tex`. Details and results: [`CLAUDE.md`](CLAUDE.md).

```bash
python3 01_aat_visibility.py    # EDF-S from Siding Spring
python3 02_aat_success.py       # redshift success vs Z and exposure (HSC-AAT calibrated)
python3 03_edfs_targets.py      # EDF-S SN Ia map / z-hist
python3 04_2df_pointings.py     # 9-field 2dF layout
python3 05_aat_host_yield.py    # host-redshift campaign yield vs hours
```
Outputs in `outputs/{png,csv}/`.
