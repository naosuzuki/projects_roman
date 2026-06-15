#!/usr/bin/env python3
"""
08_hsc_etc.py -- Subaru/HSC imaging exposure-time calculator.

Point-source S/N (aperture photometry, native 1x1 pixels) as a function of
total exposure time, for HSC g / r2 / i2, at source AB magnitudes 24-27.
One figure per band, with 60/90/120-min guides.

Sky-background-limited CCD equation:
    S/N = S / sqrt(S + B_sky + N_read),
    S      = src_rate(mag) * t * f_enc        [source e- in aperture]
    B_sky  = sky_rate_pix  * t * n_pix         [sky e- in aperture]
    N_read = (t/t_exp) * n_pix * RN^2

Throughputs/sky are tuned to reproduce HSC-SSP Wide-like depths
(~5 sigma at i~26 in ~20 min). Parameters at top are easy to adjust.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV_DIR = os.path.join(HERE, "outputs", "csv")

# --- telescope / detector ---
AREA_CM2 = 52.81e4           # Subaru 8.2 m effective collecting area [cm^2]
PIXSCALE = 0.168             # HSC pixel scale [arcsec/pix]
RN = 4.5                     # read noise [e-/pix]
SEEING = 0.70               # PSF FWHM [arcsec]
T_EXP = 300.0                # single-frame time [s] (for read-noise accounting)
NPH0 = 5.48e6                # photons/s/cm^2 per (dlam/lam) for AB=0

# HSC bands: (lambda_nm, dlam/lam, total throughput, sky AB mag/arcsec^2)
HSC = {
    "g":  dict(lam=477.0, dll=0.29, thru=0.42, sky=22.0, color="#d62728"),
    "r2": dict(lam=623.0, dll=0.23, thru=0.48, sky=21.1, color="#1f77b4"),
    "i2": dict(lam=771.0, dll=0.19, thru=0.40, sky=20.2, color="#2ca02c"),
}
MAGS = [24.0, 25.0, 26.0, 27.0, 28.0]
VLINES = [60, 90, 120]       # minutes


def src_rate(mag, b):        # source electrons / s
    return NPH0 * 10 ** (-0.4 * mag) * b["dll"] * AREA_CM2 * b["thru"]


def sky_rate_pix(b):         # sky electrons / s / pixel
    return (NPH0 * 10 ** (-0.4 * b["sky"]) * b["dll"] * AREA_CM2 * b["thru"]
            * PIXSCALE ** 2)


def aperture():
    """n_pix in a radius=FWHM aperture and the enclosed Gaussian flux fraction."""
    r = SEEING                                   # aperture radius = FWHM [arcsec]
    n_pix = np.pi * r ** 2 / PIXSCALE ** 2
    sigma = SEEING / 2.3548
    f_enc = 1.0 - np.exp(-r ** 2 / (2 * sigma ** 2))
    return n_pix, f_enc


def sn(mag, t_min, b, n_pix, f_enc):
    t = np.asarray(t_min) * 60.0
    S = src_rate(mag, b) * t * f_enc
    Bsky = sky_rate_pix(b) * t * n_pix
    Nread = (t / T_EXP) * n_pix * RN ** 2
    return S / np.sqrt(S + Bsky + Nread)


def main():
    os.makedirs(PNG_DIR, exist_ok=True); os.makedirs(CSV_DIR, exist_ok=True)
    n_pix, f_enc = aperture()
    t_min = np.logspace(np.log10(2), np.log10(360), 200)   # 2 min .. 6 h

    print(f"HSC imaging ETC  (seeing {SEEING}\", aperture r={SEEING}\", "
          f"n_pix={n_pix:.0f}, f_enc={f_enc:.2f})")
    # depth check + CSV
    csv = os.path.join(CSV_DIR, "08_hsc_etc.csv")
    with open(csv, "w") as fo:
        fo.write("band,mag,t_min,SN\n")
        for bn, b in HSC.items():
            s60 = sn(26.0, 60.0, b, n_pix, f_enc)
            print(f"  {bn}: mag 26 in 60 min -> S/N {s60:.1f};  "
                  f"5sigma depth at 20min = mag "
                  f"{26 + 2.5*np.log10(sn(26.0,20.0,b,n_pix,f_enc)/5.0):.2f}")
            for mag in MAGS:
                for tm in t_min[::20]:
                    fo.write(f"{bn},{mag},{tm:.1f},{sn(mag,tm,b,n_pix,f_enc):.3f}\n")
    print("table ->", csv)

    # --- plots: one per band ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    for fnt in ("Times New Roman.ttf", "Times New Roman Bold.ttf", "Times New Roman Italic.ttf"):
        fp = f"/System/Library/Fonts/Supplemental/{fnt}"
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    from matplotlib.ticker import FixedLocator, FuncFormatter
    cmap = plt.get_cmap("rainbow_r")                          # same as PFS plots
    cols = [cmap(x) for x in np.linspace(0, 1, len(MAGS))]    # 24 red -> 28 violet
    tmk = np.array([5, 7.5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 240, 300, 360])
    tref = np.logspace(np.log10(tmk.min()), np.log10(tmk.max()), 100)
    yticks = [0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500]

    for bn, b in HSC.items():
        fig, ax = plt.subplots(figsize=(8, 6))
        for mag, c in zip(MAGS, cols):
            ax.plot(tmk, sn(mag, tmk, b, n_pix, f_enc), "o-", color=c, lw=1.8,
                    ms=5, label=f"{bn} = {mag:.0f} AB")
        # sqrt(t) reference anchored to the middle magnitude
        mref = MAGS[len(MAGS) // 2]
        s0 = sn(mref, tmk[0], b, n_pix, f_enc)
        ax.plot(tref, s0 * np.sqrt(tref / tmk[0]), "k--", lw=1.2, alpha=0.7,
                label=r"$\propto\sqrt{t}$ (bkg-limited)")
        for tv in VLINES:
            ax.axvline(tv, color="0.35", ls=":", lw=1.3)
            ax.text(tv, 0.985, f"{tv} min", transform=ax.get_xaxis_transform(),
                    rotation=90, va="top", ha="right", fontsize=13, color="0.3")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.yaxis.set_major_locator(FixedLocator(yticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_xlabel("Total Exposure Time  [min]", fontsize=17)
        ax.set_ylabel("Imaging S/N (Point Source)", fontsize=17)
        ax.set_title(f"Subaru/HSC Imaging ETC — {bn}-band S/N vs Exposure Time (log–log)\n"
                     f"Seeing {SEEING}″, Dark Sky, "
                     f"$\\lambda\\approx{b['lam']:.0f}$ nm", fontsize=15)
        ax.tick_params(labelsize=13)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(frameon=True, fontsize=12, loc="upper left")
        png = os.path.join(PNG_DIR, f"08_hsc_etc_{bn}band.png")
        fig.tight_layout(); fig.savefig(png, dpi=140)
        print("plot ->", png)


if __name__ == "__main__":
    main()
