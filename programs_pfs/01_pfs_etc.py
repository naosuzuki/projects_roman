"""
pfs_etc.py -- thin reusable wrapper around the Subaru PFS Exposure Time
Calculator (Subaru-PFS/spt_ExposureTimeCalculator).

The ETC itself is a C+Python tool installed separately (see README.md). This
module runs it for a given configuration, parses the continuum-S/N output, and
provides helpers to summarize S/N within a broadband filter (g/r/i/z/y).

Typical use
-----------
    from pfs_etc import run_etc, band_sn, BANDS
    res = run_etc(mag=22.5, exp_time=450, exp_num=8)   # -> dict of arrays
    sn_i = band_sn(res, "i")                            # median S/N/pix in i-band
"""
import os
import subprocess
import numpy as np

# ----------------------------------------------------------------------------
# Configuration: where the ETC lives and its default parameter file.
# Override with the PFS_ETC_DIR environment variable if you move the ETC.
# ----------------------------------------------------------------------------
ETC_DIR = os.environ.get(
    "PFS_ETC_DIR", os.path.expanduser("~/github/spt_ExposureTimeCalculator")
)
DEFAULTS = os.path.join(ETC_DIR, "scripts", "run_etc.defaults")
OUT_SNC = os.path.join(ETC_DIR, "out", "ref.snc.dat")

# PFS spectrograph arms (vacuum wavelength coverage, nm)
ARMS = {0: ("blue", 380, 650), 1: ("red", 630, 970), 2: ("nir", 940, 1260)}

# Per-arm pixel dispersion [nm/pixel] = coverage / 4096 px (measured from ETC).
ARM_DISP_NM = {0: 0.0659, 1: 0.0830, 2: 0.0781}

# Broadband filters: (lambda_min_nm, lambda_max_nm, preferred_arm).
# Ranges approximate HSC/SDSS; "arm" disambiguates the small arm overlaps.
BANDS = {
    "g": (400.0, 550.0, 0),   # HSC-g, eff ~478 nm (blue arm)
    "r": (550.0, 690.0, 0),
    "r2": (560.0, 650.0, 0),  # HSC-r2, eff ~620 nm (blue-arm portion)
    "i": (700.0, 830.0, 1),   # effective wavelength ~770 nm (red arm)
    "i2": (700.0, 830.0, 1),  # HSC-i2, eff ~771 nm (red arm; ~ i)
    "z": (830.0, 925.0, 1),
    "y": (950.0, 1060.0, 2),
}


def run_etc(mag, exp_time=450, exp_num=8, defaults=DEFAULTS, **overrides):
    """Run the PFS ETC once and return parsed continuum-S/N arrays.

    Parameters
    ----------
    mag : float or str
        Target AB magnitude (flat continuum) or path to an input spectrum.
    exp_time : float
        Single-exposure time [sec].
    exp_num : int
        Number of exposures (total time = exp_time * exp_num).
    overrides : dict
        Any other ETC keyword, e.g. SEEING=0.6, ZENITH_ANG=30, MOON_PHASE=0.25.

    Returns
    -------
    dict with numpy arrays: arm, pixel, wl (nm), snc (S/N per pixel),
    plus 'meta' (the resolved parameters).
    """
    cmd = [
        "pfs-run-etc", f"@{defaults}",
        f"--MAG_FILE={mag}", f"--EXP_TIME={exp_time:.0f}", f"--EXP_NUM={int(exp_num)}",
        "--OUTFILE_SNL=-", "--OUTFILE_OII=-",   # skip line-S/N outputs (faster)
    ]
    for k, v in overrides.items():
        cmd.append(f"--{k}={v}")
    subprocess.run(cmd, cwd=ETC_DIR, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    d = np.genfromtxt(OUT_SNC)
    return {
        "arm": d[:, 0].astype(int),
        "pixel": d[:, 1].astype(int),
        "wl": d[:, 2],
        "snc": d[:, 3],
        "meta": dict(mag=mag, exp_time=exp_time, exp_num=exp_num,
                     total_time_s=exp_time * exp_num, **overrides),
    }


def band_mask(res, band):
    """Boolean mask selecting pixels of `res` inside `band`'s window+arm."""
    lo, hi, arm = BANDS[band]
    return (res["arm"] == arm) & (res["wl"] >= lo) & (res["wl"] <= hi)


def band_resolution(band, npix=1):
    """Effective spectral resolution R = lambda / (npix * dispersion) for a
    binned spectrum, evaluated at the band's central wavelength.

    Returns (lambda_eff_nm, R). Binning npix pixels also boosts continuum S/N
    by sqrt(npix) for a flat continuum (handled by the caller).
    """
    lo, hi, arm = BANDS[band]
    lam = 0.5 * (lo + hi)
    R = lam / (npix * ARM_DISP_NM[arm])
    return lam, R


def band_sn(res, band, stat="median"):
    """Summarize continuum S/N per pixel within a broadband filter.

    stat : 'median' (default), 'mean', or 'quad' (sqrt(sum S/N^2),
           i.e. the effective broadband-continuum S/N).
    """
    s = res["snc"][band_mask(res, band)]
    if s.size == 0:
        return float("nan")
    if stat == "median":
        return float(np.median(s))
    if stat == "mean":
        return float(np.mean(s))
    if stat == "quad":
        return float(np.sqrt(np.sum(s ** 2)))
    raise ValueError(f"unknown stat: {stat}")
