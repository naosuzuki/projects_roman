#!/usr/bin/env python3
"""
27_snia_specz_elg_zhist.py -- as 26 (combined live-SN + host spec-z redshift
histogram), plus a rough estimate of EMISSION-LINE redshifts: star-forming hosts
yield an [OII] 3727 line redshift in ~1 visit even without reaching the
continuum S/N=5 requirement, so every Ia host that receives >=1 h of integration
and is star-forming counts.

Assumptions (rough, stated in the figure):
  * star-forming (emission-line) fraction  f_SF(z) = min(0.85, 0.70 + 0.15 z)
    (magnitude-limited samples at z~1 are ~75-85% emission-line galaxies, and
    SN Ia hosts skew star-forming at high z)
  * a line redshift needs only >=1 h on target (host_started or completed)
  * the 80% weather factor applies as elsewhere

    python 27_snia_specz_elg_zhist.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_DIR = os.path.join(HERE, "outputs", "png")
CSV = os.path.join(HERE, "outputs", "csv", "07_program_sne_ELAIS-N1.csv")

CASES = [(1.0, "SNIa 100% (no weather loss)", "blue"),
         (0.8, "SNIa 80% weather",            "g"),
         (0.7, "SNIa 70% weather",            "red")]
W = 0.8                                        # weather factor for host/ELG


def f_sf(z):
    """Star-forming (emission-line) fraction vs redshift (rough)."""
    return np.minimum(0.85, 0.70 + 0.15 * np.asarray(z, float))


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    typ = np.asarray(d["type"]).astype(str)
    obs = np.asarray(d["observed"]).astype(int)
    hcomp = np.asarray(d["host_completed"]).astype(int)
    hstart = np.asarray(d["host_started"]).astype(int)
    z = np.asarray(d["z"]).astype(float)
    ia = typ == "Ia"
    zo = z[ia & (obs == 1)]                    # live SN spec-z
    zh = z[ia & (hcomp == 1)]                  # host continuum spec-z (S/N=5)
    m_elg = ia & ((hstart == 1) | (hcomp == 1))  # any integration -> line-z if SF
    zel, wel = z[m_elg], f_sf(z[m_elg])
    n_elg = W * wel.sum()
    print(f"live SN spec-z: {len(zo)};  host continuum spec-z: {len(zh)}")
    print(f"Ia hosts with >=1h: {m_elg.sum()};  expected ELG line-z "
          f"(x f_SF x {W} weather): {n_elg:.0f}")

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

    zmax = max(zo.max(), zh.max(), zel.max())
    bins = np.arange(0.0, np.ceil(zmax * 10) / 10 + 0.051, 0.05)
    counts, edges = np.histogram(zo, bins)
    hcounts, _ = np.histogram(zh, bins)
    ecounts, _ = np.histogram(zel, bins, weights=wel)   # expected ELG counts
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = (edges[1] - edges[0]) * 0.93
    fig, ax = plt.subplots(figsize=(9, 6))
    # ELG line-z estimate: outermost background, orange
    ax.bar(centers, ecounts * W, width=bw, color="orange", alpha=0.45,
           edgecolor="white", linewidth=0.4, zorder=0,
           label=f"Host [OII] line-$z$ 80% weather  ($N={n_elg:.0f}$)")
    # host continuum case in the background, light green alpha=0.6
    nh = W * len(zh)
    ax.bar(centers, hcounts * W, width=bw, color="lightgreen", alpha=0.6,
           edgecolor="white", linewidth=0.4, zorder=1,
           label=f"Host continuum 80% weather  ($N={nh:.0f}$)")
    # nested SN bars on top (pure colors)
    for i, (w, lab, col) in enumerate(CASES):
        n = w * len(zo)
        ax.bar(centers, counts * w, width=bw, color=col, edgecolor="white",
               linewidth=0.4, zorder=i + 2, label=f"{lab}  ($N={n:.0f}$)")

    ax.set_xlabel("Redshift", fontsize=21)
    ax.set_ylabel("Expected Number of Successful Spec-$z$", fontsize=21)
    ax.set_title("ELAIS-N1 SN Ia, Host Continuum, and Host Line-$z$ vs Redshift",
                 fontsize=17)
    ax.tick_params(labelsize=16)
    ax.set_xlim(0, bins[-1])
    ax.legend(fontsize=12.5, loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(0.985, 0.02, "$f_{\\rm SF}(z)=\\min(0.85,\\ 0.70+0.15z)$;  "
            "line-$z$ assumes $\\geq$1 h on target",
            transform=ax.transAxes, fontsize=10.5, color="0.25",
            ha="right", va="bottom", zorder=10,
            bbox=dict(fc="white", ec="0.7", alpha=0.85))

    png = os.path.join(PNG_DIR, "27_snia_specz_elg_zhist_ELAIS-N1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    print("plot ->", png)


if __name__ == "__main__":
    main()
