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


EXHIBITS = {"S0.pore_profile": _pore, "P5.deep_ranks_third": _aucs,
            "P2.teleost_itpr1": _shares, "P3.no_absent_cells": _states,
            "S0.ip3_contacts": _contacts}


def has_exhibit(check_id: str, outcome) -> bool:
    return check_id in EXHIBITS and bool(outcome.data)


def draw(ax, check_id: str, outcome) -> bool:
    """Draw the exhibit for ``check_id``; False if it has none."""
    if not has_exhibit(check_id, outcome):
        return False
    EXHIBITS[check_id](ax, outcome.data)
    ax.set_title(check_id, fontsize=9)
    return True
