#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""Fig 5 (cross-algorithm count) in the Summer-26 house style (matches fig2b/3a/4a).
Colors: algorithms use a NON-metric triad (metrics are blue/red/green in the set) --
Classical=black(o), Deep-ULM=purple(s), LISTA=dark goldenrod(^).
jun23-led, combined tube, M3 (rung 8) excluded from points and fits.
"""
import os, sys, numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inference"))

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"], "font.size": 11.5,
    "axes.titlesize": 13, "axes.titleweight": "bold", "axes.linewidth": 0.9,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "grid.color": "0.85", "grid.linewidth": 0.6, "legend.fontsize": 10.5,
})

from paths import OUTPUT_ROOT as Z, SCRATCH_DIR as SC  # machine-local roots; see learned/README.md
M3 = 7
METH = [("Classical (LAT-ULM)", "#000000", "o", "classical"),
        ("Deep-ULM", "#7E2F8E", "s", "deepulm"),
        ("LISTA", "#B8860B", "^", "lista")]


def read_ro(path):
    with h5py.File(path, "r") as f:
        g = f["readout"]
        cell = lambda n: ["".join(chr(int(c)) for c in np.array(f[r]).ravel()) for r in np.array(g[n]).ravel()]
        return dict(domains=cell("domains"), roiNames=cell("roiNames"),
                    conc=np.array(g["rungConc"]).astype(float).ravel(),
                    loc_sub=np.array(g["loc_sub"]).astype(float).transpose(),
                    loc_sem=np.array(g["loc_sem"]).astype(float).transpose())


def block_sem(per_block, rung_tags):
    """Per-rung SEM from a {tag: rate} dict over each rung's block tags."""
    out = []
    for tags in rung_tags:
        v = np.array([per_block[t] for t in tags if t in per_block], float)
        out.append(np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan)
    return np.array(out)


def main():
    import datasets as DS, count_locrate as CL
    cfg = DS.get_dataset("jun23"); conc = np.asarray(cfg["conc"], float)
    ro = read_ro(os.path.join(Z, "jun23", "readout_jun23.mat"))
    d = ro["domains"].index("PI"); rC = ro["roiNames"].index("combinedTube")
    out = np.load(os.path.join(SC, "locrate_jun23.npy"), allow_pickle=True).item()
    stats = np.load(os.path.join(SC, "ceiling_stats.npy"), allow_pickle=True).item()["jun23"]["combinedTube"]
    rung_tags = DS.rung_blocks(cfg)

    keep = np.array([i for i in range(15) if i != M3])
    x = conc[keep]

    def series(key):
        if key == "classical":
            return ro["loc_sub"][d, :, rC], ro["loc_sem"][d, :, rC]
        r = out["results"]["PI"]["combinedTube"][key]
        return np.asarray(r["loc_sub"], float), block_sem(r["locRate"], rung_tags)

    bp = stats["classical"]["breakpoint"]           # physical ceiling onset (method-invariant)

    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    lo, hi = float(x.min()), float(x.max())
    # regime shading: below vs above the count ceiling (fold-change cancels the method-specific offset)
    ax.axvspan(lo, bp, color="#EAF3EA", zorder=0)                     # below ceiling (faint green)
    ax.axvspan(bp, hi, color="#FBECEC", zorder=0)                     # above ceiling (faint red)
    ax.axvline(bp, color="0.45", ls="--", lw=1.2, zorder=1)

    fmax = 0.0
    for label, col, mk, key in METH:
        y, sem = series(key)
        f = y / y[0]; ferr = sem / y[0]                              # fold-change from the lowest rung
        yk, ek = f[keep], ferr[keep]
        m = np.isfinite(yk) & (yk > 0)
        ax.errorbar(x[m], yk[m], yerr=np.where(np.isfinite(ek[m]), ek[m], 0.0),
                    color=col, marker=mk, ms=7, lw=1.5, capsize=2.5, mec="white", mew=0.6,
                    label=label, zorder=3)
        fmax = max(fmax, np.nanmax(yk[m]))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi)                                              # data edge to edge
    ax.set_xlabel("Concentration (MB/mL)")
    ax.set_ylabel("Localization count (fold-change from lowest rung)")
    ax.set_title("Learned localizers sustain the count above the classical ceiling,\n"
                 "extending the usable concentration range (15-Rung, PI, Combined Tube; M3 excluded)")
    ax.grid(True, which="both", alpha=0.6)
    ax.legend(loc="upper left", frameon=True, edgecolor="0.7")

    # region labels (top of each shaded band)
    yt = ax.get_ylim()[1]
    ax.text(np.sqrt(lo * bp), yt * 0.72, "below ceiling\ncount tracks concentration",
            ha="center", va="top", fontsize=10, color="#2E6B2E", style="italic")
    ax.text(np.sqrt(bp * hi), yt * 0.72, "above ceiling\nsaturation",
            ha="center", va="top", fontsize=10, color="#8B3A3A", style="italic")
    ax.text(bp, ax.get_ylim()[0] * 1.15, " ceiling ~%.0e MB/mL" % bp, ha="left", va="bottom",
            fontsize=8.5, color="0.35")

    fig.tight_layout()
    outp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "real", "fig5_house_v3.png")
    fig.savefig(outp, dpi=200, bbox_inches="tight"); print("saved ->", outp)


if __name__ == "__main__":
    main()
