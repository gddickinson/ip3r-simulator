"""Figures for the VUS stratification (Paper 5 §8).

``draw_fractions`` is the exhibit of ``P5.vus_stratification``: per gene and
layer, the share of VUS at or past each labelled median, on a fixed 0–1
scale. ``draw_strip`` is the Variants tab's picture of one gene on one layer:
every scored position of each class with the two medians drawn.
"""

from __future__ import annotations

import numpy as np

from .vus_strata import CLASS_COLORS, STRATUM_COLORS, Stratification

__all__ = ["draw_fractions", "draw_strip"]

_TEXT = "#d7dbe3"
_LAYERS = ("deep", "shallow", "vert", "family")


def draw_fractions(ax, d: dict) -> None:
    per = d["strata"]
    genes = sorted(per)
    x, labels = [], []
    for gi, gene in enumerate(genes):
        for li, layer in enumerate(_LAYERS):
            v = per[gene].get(layer)
            if v is None:
                continue
            pos = gi * (len(_LAYERS) + 1) + li
            hi, lo = v["frac_vus_above_pathogenic_median"], v["frac_vus_below_benign_median"]
            ax.bar(pos, hi, color=STRATUM_COLORS["pathogenic-like"], width=0.8)
            ax.bar(pos, 1 - hi - lo, bottom=hi, color=STRATUM_COLORS["between"],
                   width=0.8, alpha=0.5)
            ax.bar(pos, lo, bottom=1 - lo, color=STRATUM_COLORS["benign-like"], width=0.8)
            x.append(pos)
            labels.append(f"{gene[-1]}·{layer[:4]}")
    ax.set_xticks(x, labels, rotation=60, fontsize=6)
    ax.set_ylim(0, 1)
    ax.set_ylabel("share of scored VUS positions")
    ax.text(0.01, 0.98, "red: ≥ P/LP median   blue: ≤ B/LB median",
            transform=ax.transAxes, va="top", fontsize=7, color=_TEXT)


def draw_strip(ax, s: Stratification) -> None:
    """One gene, one layer: the three classes' scores and both medians."""
    rng = np.random.default_rng(0)                  # fixed jitter: same picture every time
    groups = ((s.pathogenic, CLASS_COLORS["P/LP"]),
              (s.vus, None), (s.benign, CLASS_COLORS["B/LB"]))
    for i, (vals, col) in enumerate(groups):
        y = np.fromiter(vals.values(), float)
        xj = i + rng.uniform(-0.25, 0.25, len(y))
        if col is None:
            col = [STRATUM_COLORS[s.strata[r]] for r in vals]
        ax.scatter(xj, y, s=5, color=col, alpha=0.8, linewidths=0)
    ax.axhline(s.median_pathogenic, color=CLASS_COLORS["P/LP"], lw=0.9, ls="--",
               label=f"P/LP median {s.median_pathogenic:.3f} (n={len(s.pathogenic)})")
    ax.axhline(s.median_benign, color=CLASS_COLORS["B/LB"], lw=0.9, ls="--",
               label=f"B/LB median {s.median_benign:.3f} (n={len(s.benign)})")
    ax.set_xticks(range(3), ["P/LP", f"VUS (n={s.n_vus})", "B/LB"])
    ax.set_ylim(0, 1)                               # fixed, never auto-ranged
    ax.set_ylabel(f"{s.layer} JSD")
    ax.set_title(f"{s.gene}, {s.layer}: {s.n_above_pathogenic} VUS ≥ P/LP median, "
                 f"{s.n_below_benign} ≤ B/LB median", fontsize=8, color=_TEXT)
    ax.legend(fontsize=7, frameon=False, labelcolor=_TEXT, loc="lower left")
