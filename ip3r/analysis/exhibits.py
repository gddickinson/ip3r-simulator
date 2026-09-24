"""Figures that illustrate a finding, drawn from a check's own outcome.

Each exhibit takes a matplotlib ``Axes`` and the :class:`Outcome` data the
check computed — never a separately loaded table — so the picture shown
beside a verdict is a picture of the numbers that verdict was reached on.
Used by the GUI's findings panel and by ``python -m ip3r findings --figures``.
"""

from __future__ import annotations

import numpy as np

__all__ = ["EXHIBITS", "has_exhibit", "draw"]

_C = ["#5b9ef2", "#f28c4d", "#72cc80", "#cc80e6", "#f2cc4d", "#66d9d9"]


def _pore(ax, d):
    ax.plot(d["z"], d["r_published"], color=_C[1], lw=1.4, label="published (S0)")
    ax.plot(d["z"], d["r_recomputed"], color=_C[0], lw=0.9, ls="--",
            label="recomputed here")
    ax.set_xlabel("z along the four-fold axis (Å)")
    ax.set_ylabel("min heavy-atom radius (Å)")
    ax.set_ylim(0, 20)
    ax.legend(fontsize=7, frameon=False, labelcolor="#d7dbe3")


def _aucs(ax, d):
    aucs = d["aucs"]
    order = sorted(aucs, key=lambda k: -aucs[k])
    ax.bar(order, [aucs[k] for k in order], color=[_C[i] for i in range(len(order))])
    ax.axhline(0.5, color="#8a8f99", lw=0.8, ls=":")
    ax.set_ylim(0.4, 1.0)
    ax.set_ylabel("AUC, P/LP vs B/LB (all-layers-scorable positions)")


def _shares(ax, d):
    s = d["share"]
    ax.bar(list(s), list(s.values()), color=_C[:len(s)])
    ax.set_ylim(0, 1)
    ax.set_ylabel("teleost genomes with ≥ 2 copies")


def _states(ax, d):
    s = dict(sorted(d["states"].items(), key=lambda x: -x[1]))
    ax.barh(list(s), list(s.values()), color=_C[0])
    ax.set_xscale("log")
    ax.set_xlabel("genome × paralog cells (log)")
    ax.invert_yaxis()


def _contacts(ax, d):
    per = d["contacts"]
    labels = sorted({r for v in per.values() for r in v},
                    key=lambda s: int("".join(c for c in s if c.isdigit())))
    m = np.array([[r in per[k] for r in labels] for k in sorted(per)], float)
    ax.imshow(m, cmap="Greens", vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(range(len(per)), [f"site {k}" for k in sorted(per)])
    ax.set_xticks(range(len(labels)), labels, rotation=60, fontsize=6)


def _modules(ax, d):
    per = d["paralogs"]
    for i, (gene, v) in enumerate(sorted(per.items())):
        ax.scatter(v["pore"], v["core"], s=6, alpha=0.6, color=_C[i],
                   label=f"{gene} (mean core − pore {v['mean_difference']:+.3f})")
    lo = min(min(min(v["core"]), min(v["pore"])) for v in per.values())
    ax.plot([lo, 1], [lo, 1], color="#8a8f99", lw=0.8, ls=":")
    ax.set_xlabel("pore-module identity to human (per orthologue)")
    ax.set_ylabel("ligand-core identity to human")
    ax.legend(fontsize=7, frameon=False, labelcolor="#d7dbe3")


def _shell_trend(ax, d):
    for i, (gene, v) in enumerate(sorted(d["trend"].items())):
        ax.scatter(v["distance"], v["jsd"], s=7, alpha=0.7, color=_C[i],
                   label=f"{gene} (ρ {v['rho']:+.2f})")
        ax.axhline(v["protein"], color=_C[i], lw=0.7, ls="--", alpha=0.6)
    for e in d["edges"]:
        ax.axvline(e, color="#8a8f99", lw=0.7, ls=":")
    ax.set_xlabel("all-atom distance to IP3, median of six depositions (Å)")
    ax.set_ylabel("deep-layer conservation (JSD)")
    ax.legend(fontsize=7, frameon=False, labelcolor="#d7dbe3")


def _shell_means(ax, d):
    names = ("contact", "second", "third", "fourth")
    for i, (gene, v) in enumerate(sorted(d["shells"].items())):
        ax.plot(names, [v[k] for k in names], marker="o", color=_C[i], label=gene)
        ax.axhline(v["protein"], color=_C[i], lw=0.7, ls="--", alpha=0.6)
    ax.set_ylabel("mean JSD per shell (dashed: whole protein)")
    ax.legend(fontsize=7, frameon=False, labelcolor="#d7dbe3")


def _tree(ax, d):
    from .newick import parse
    from .tree_figure import draw_tree
    draw_tree(ax, parse(d["newick"]))


def _grid(name):
    def draw(ax, d):
        from . import grid_figure
        getattr(grid_figure, name)(ax, d)
    return draw


def _range(name):
    def draw(ax, d):
        from . import range_figure
        getattr(range_figure, name)(ax, d)
    return draw


def _vus(ax, d):
    from .vus_figure import draw_fractions
    draw_fractions(ax, d)


EXHIBITS = {"P5.vus_stratification": _vus, "S0.pore_profile": _pore, "P5.deep_ranks_third": _aucs,
            "P2.teleost_itpr1": _shares, "P3.no_absent_cells": _states,
            "S0.ip3_contacts": _contacts, "P6.module_contrast": _modules,
            "P6.loop_reverses": _modules, "P6.shell_trend": _shell_trend,
            "P6.shell_constraint": _shell_means, "P2.sister_pair": _tree,
            "P2.paralog_clades": _tree, "P2.cyclostome_lineages": _tree,
            "P2.support_bar": _tree, "P3.miss_by_contiguity": _grid("draw_misses"),
            "P3.contiguity_tests": _grid("draw_logistic"),
            "P4.recovery_channels": _grid("draw_recovery"),
            "P1.presence_range": _range("draw_presence"),
            "P1.kingdom_absences": _range("draw_presence"),
            "P1.absence_targets": _range("draw_presence"),
            "P1.relaxed_controls": _range("draw_relaxed"),
            "P1.absences": _range("draw_absences"),
            "P1.copy_number": _range("draw_copies"),
            "P1.record_chase": _range("draw_chase")}


def has_exhibit(check_id: str, outcome) -> bool:
    return check_id in EXHIBITS and bool(outcome.data)


def draw(ax, check_id: str, outcome) -> bool:
    """Draw the exhibit for ``check_id``; False if it has none."""
    if not has_exhibit(check_id, outcome):
        return False
    EXHIBITS[check_id](ax, outcome.data)
    ax.set_title(check_id, fontsize=9)
    return True
