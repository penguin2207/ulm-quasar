#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""Fig 5 mechanism panel: PI vs fundamental slopes (jun23 combined tube, M3 excluded).
Shows the supralinearity is a PI-contrast effect (linear-range slope relaxes in fundamental) while the
range-extension is domain-independent (above-ceiling divergence persists). House style; non-metric colors.
Slopes from ceiling_analysis / the PI-vs-fundamental count (both nets swept in both domains 2026-07-03)."""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 11.5,
                     "axes.titlesize": 12, "axes.titleweight": "bold", "axes.linewidth": 0.9,
                     "figure.facecolor": "white", "grid.color": "0.85", "legend.fontsize": 10})

# (PI, fundamental) count-vs-concentration slope, jun23 combined tube, M3 excluded
LINEAR = {"Classical (LAT-ULM)": (0.76, 0.60, "#000000", "o"),
          "Deep-ULM":            (1.14, 0.81, "#7E2F8E", "s"),
          "LISTA":               (1.22, 0.86, "#B8860B", "^")}
HIGH   = {"Classical (LAT-ULM)": (0.23, 0.24, "#000000", "o"),
          "Deep-ULM":            (0.31, 0.34, "#7E2F8E", "s"),
          "LISTA":               (0.41, 0.40, "#B8860B", "^")}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.2))
for ax, data, ttl in [(a1, LINEAR, "Linear range (L1-U1): steepness RELAXES in fundamental"),
                      (a2, HIGH, "Above ceiling (U1-U5): extension PERSISTS in both domains")]:
    for name, (pi, fu, col, mk) in data.items():
        ax.plot([0, 1], [pi, fu], color=col, marker=mk, ms=10, lw=2.2, mec="white", mew=0.9,
                label=name if ax is a1 else None)
    ax.axhline(1.0, color="0.55", ls=":", lw=1.3)
    ax.text(1.28, 1.0, "linear\n(slope 1)", ha="left", va="center", fontsize=8.5, color="0.4")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["PI\n(nonlinear contrast)", "fundamental\n(linear)"])
    ax.set_xlim(-0.32, 1.5); ax.set_ylim(0, 1.45)
    ax.set_ylabel("count-vs-concentration slope")
    ax.set_title(ttl, fontsize=11)
    ax.grid(True, axis="y", alpha=0.55)
a1.legend(loc="lower left", frameon=True, edgecolor="0.7")
fig.suptitle("Two effects separate by imaging domain: the below-ceiling supralinearity is PI contrast; "
             "the\nabove-ceiling range extension is super-resolution (jun23, Combined Tube, M3 excluded)",
             fontsize=12.5, fontweight="bold", y=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.93))
outp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "real", "fig5_pivsfund.png")
fig.savefig(outp, dpi=200, bbox_inches="tight"); print("saved ->", outp)
