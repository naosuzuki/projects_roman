#!/usr/bin/env python3
"""
29_survey_depths.py -- Survey landscape: depth vs. area for deep extragalactic fields.

Reproduces the "Survey Area vs. Total Stacked Depth" diagram used in the
Roman-Subaru proposal / roman2026 (Figures/survey_depths.pdf).  Each field is a
marker keyed by facility (HST = circle, JWST = triangle, Roman = square) and
coloured by wavelength coverage (Optical = blue, NIR = red, Both = magenta).
Roman HLTDS Deep/Wide and the proposed XDF ultra-deep field are highlighted.
Dashed grey lines are iso-"survey-power" contours (constant depth + area grasp);
power increases toward the upper right.

Usage:
    python3 29_survey_depths.py            # writes outputs/survey_depths.{pdf,png}

No external data or the PFS ETC are needed -- the field list is tabulated below.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter

# ---------------------------------------------------------------- Times-Roman fonts
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"   # Times-compatible math ($\sigma$, arcmin$^2$)

# ---------------------------------------------------------------- style / colors
COL = {"Optical": "b", "NIR": "#e12b25", "Both": "#e024c0"}
MARK = {"HST": "o", "JWST": "^", "Roman": "s", "Subaru": "h", "Euclid": "D"}

# ---------------------------------------------------------------- field table
# name, area[arcmin^2], depth[AB, 5sigma point source], facility, wavelength,
# label offset (dex in x, mag in y), leader line?, highlight?, marker size
FIELDS = [
    ("NGDEEP (HUDF)",          6,   29.68, "JWST",  "NIR",     (-0.046, 0.12, "left"), True,  False, 120),
    ("HUDF",                   11,  29.00, "HST",   "Optical", ( 0.12,  0.02), False, False, 130),
    ("UNCOVER (Abell 2744)",   45,  29.72, "JWST",  "NIR",     ( 0.06,  0.00), True,  False, 120),
    ("CEERS (EGS)",            100, 29.25, "JWST",  "NIR",     ( 0.06,  0.00), True,  False, 120),
    ("GLASS",                  18,  28.75, "JWST",  "NIR",     ( 0.06,  0.00), False, False, 120),
    ("HDF-N+HDF-S",            18,  28.58, "HST",   "Optical", ( 0.06,  0.00), False, False, 120),
    ("Hubble Frontier Fields", 130, 28.85, "HST",   "Both",    ( 0.06,  0.00), False, False, 120),
    ("JADES/PRIMER UDS",       185, 28.75, "JWST",  "NIR",     ( 0.06, -0.07), True,  False, 120),
    ("GOODS+CANDELS",          400, 27.22, "HST",   "Both",    ( 0.15,  0.02), True,  False, 120),
    ("COSMOS-Web",             2100, 28.00, "JWST",  "NIR",    (-0.06,  0.00), False, False, 120),
    ("COSMOS",                 7000, 27.22, "HST",   "Optical", ( 0.06,  0.00), False, False, 130),
    ("COSMOS-DASH",            7000, 26.35, "HST",   "NIR",     (-0.08,  0.00), True,  False, 120),
    ("HLTDS Deep\nELAIS-N1, 2.4 deg$^2$",   8640,  29.15, "Roman", "NIR", ( 0.00,  0.42), False, True,  380),
    ("HLTDS Wide\nELAIS-N1, 11.6 deg$^2$",  41760, 28.25, "Roman", "NIR", ( 0.00,  0.34), False, True,  380),
    ("HLTDS HSC, g, r2\nThis Proposal",  50220, 28.00, "Subaru", "Optical", ( 0.11,  0.00), False, True,  380),
    ("HLTDS HSC, i2\nThis Proposal",     50220, 27.50, "Subaru", "Optical", ( 0.11,  0.00), False, True,  380),
    ("HLWAS Ultra-Deep",       18000, 28.20, "Roman", "NIR",   (-0.06,  0.00), False, False, 150),
    ("HLWAS Deep",             50000, 27.70, "Roman", "NIR",   (-0.06,  0.00), False, False, 150),
    ("H20 HSC",                72000, 27.00, "Subaru", "Optical", ( 0.08,  0.00), False, False, 150),
    ("H20 Euclid",             72000, 26.40, "Euclid", "NIR",     ( 0.12,  0.00), False, False, 150),
]

# ---------------------------------------------------------------- figure builder
def build_figure(hero_fs=17):
    """Draw the survey-depth diagram. hero_fs sets the font size of the
    highlighted labels (HLTDS Deep/Wide and HLTDS HSC, both Roman-red and
    Subaru-blue)."""
    fig, ax = plt.subplots(figsize=(11.5, 7.6))

    # iso-"survey power" contours: depth = b - 0.33 * log10(area)
    xline = np.array([3.0, 8e5])
    for b in (30.73, 29.83, 28.93, 28.03):
        ax.plot(xline, b - 0.33 * np.log10(xline), ls="--", color="0.55",
                lw=1.4, zorder=1)

    # points + labels
    for name, area, depth, fac, wl, off, leader, hi, size in FIELDS:
        color = COL[wl]
        marker = "*" if fac == "Star" else MARK[fac]
        ax.scatter(area, depth, marker=marker, s=size, color=color,
                   edgecolor="black", linewidth=1.0,
                   zorder=6 if fac == "Star" else 5)
        dx, dy = off[0], off[1]
        lx, ly = area * 10 ** dx, depth + dy
        # optional 3rd element of the offset overrides horizontal alignment
        ha = off[2] if len(off) > 2 else (
            "right" if dx < 0 else ("center" if dx == 0 else "left"))
        ax.text(lx, ly, name, ha=ha, va="center", color=color, zorder=7,
                fontsize=(hero_fs if hi else 15),
                weight=("bold" if hi else "normal"))

    # axes
    ax.set_xscale("log")
    ax.set_xlim(3, 8e5)
    ax.set_ylim(26, 30.2)
    ax.set_yticks([26, 27, 28, 29, 30])
    ax.set_xlabel(r"Survey Area  (Arcmin$^2$)", fontsize=28)
    ax.set_ylabel("Total Stacked Depth per Filter\n(AB mag, point source at 5$\\sigma$)",
                  fontsize=28)
    ax.tick_params(axis="both", which="both", direction="in", labelsize=28)
    ax.grid(True, which="major", color="0.85", lw=0.7)

    # secondary x-axis on TOP: square degrees (1 deg^2 = 3600 arcmin^2)
    axtop = ax.twiny()
    axtop.set_xscale("log")
    axtop.set_xlim(3 / 3600.0, 8e5 / 3600.0)
    axtop.set_xlabel(r"Survey Area  (Deg$^2$)", fontsize=28)
    axtop.tick_params(axis="x", which="both", direction="in", labelsize=28)
    axtop.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2])
    axtop.set_xticklabels(["0.001", "0.01", "0.1", "1.0", "10.0", "100.0"])
    axtop.xaxis.set_minor_formatter(NullFormatter())

    # two legends
    tel_handles = [
        Line2D([], [], marker="o", ls="", ms=11, mfc="0.6", mec="black", label="HST"),
        Line2D([], [], marker="^", ls="", ms=11, mfc="0.6", mec="black", label="JWST"),
        Line2D([], [], marker="s", ls="", ms=11, mfc="0.6", mec="black", label="Roman"),
        Line2D([], [], marker="h", ls="", ms=12, mfc="0.6", mec="black", label="Subaru"),
        Line2D([], [], marker="D", ls="", ms=10, mfc="0.6", mec="black", label="Euclid"),
    ]
    wl_handles = [
        Line2D([], [], marker="o", ls="", ms=11, mfc=COL["Optical"], mec="black", label="Optical"),
        Line2D([], [], marker="o", ls="", ms=11, mfc=COL["NIR"], mec="black", label="NIR"),
        Line2D([], [], marker="o", ls="", ms=11, mfc=COL["Both"], mec="black", label="Both"),
    ]
    leg1 = ax.legend(handles=tel_handles, title="Telescope", loc="lower left",
                     bbox_to_anchor=(0.02, 0.02), fontsize=13, title_fontsize=14,
                     framealpha=0.95)
    leg1._legend_box.align = "left"
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=wl_handles, title="Wavelength", loc="lower left",
                     bbox_to_anchor=(0.20, 0.02), fontsize=13, title_fontsize=14,
                     framealpha=0.95)
    leg2._legend_box.align = "left"

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- save default + v2
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
os.makedirs(os.path.join(OUTDIR, "png"), exist_ok=True)
# ("", 18) = default;  ("_v2", 20) = HLTDS / HLTDS-HSC labels at 20 pt
for suffix, hfs in (("", 18), ("_v2", 20)):
    fig = build_figure(hero_fs=hfs)
    pdf = os.path.join(OUTDIR, f"survey_depths{suffix}.pdf")
    png = os.path.join(OUTDIR, "png", f"survey_depths{suffix}.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=140)
    plt.close(fig)
    print("wrote", pdf)
    print("wrote", png)
