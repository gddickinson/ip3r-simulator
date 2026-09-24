"""Conservation against distance to IP3, over S22's ligand shells.

The pocket is measured here on the six IP3-bound ITPR3 depositions
(:mod:`ip3r.structure.shells`), carried to ITPR1 and ITPR2 by this project's
own pairwise alignment (:mod:`ip3r.core.pairwise`), and joined to S17's
per-residue deep-layer conservation. Three questions, as Paper 6 asks them:

* per shell, is conservation above the whole protein? (one-sided
  Mann–Whitney, each shell against every scored residue)
* does conservation fall with distance? (Spearman over the pocket)
* is there a step at the contact radius? (contact shell against the second,
  and whether the 4.5 Å boundary carries the largest drop between shells)
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from ..core import genes_data as G
from ..core.pairwise import paralog_transfer
from ..io.loader import load
from ..parameters import PARAMETERS as _P
from ..structure.shells import SHELLS, ShellResidue, consensus_shells, deposit_distances
from .checks_constraint import per_residue
from .stats import mann_whitney_greater, spearman

__all__ = ["DEPOSITS", "STRUCTURE_PARALOG", "measured_shells", "Pocket",
           "pocket", "shell_rows", "trend", "contact_step", "clear_caches"]

#: The six IP3-bound human ITPR3 depositions S22 measured, S0's first.
DEPOSITS = ("6DQN", "8TKG", "8TKF", "8TKH", "7T3P", "8TLA")
STRUCTURE_PARALOG = "ITPR3"


@lru_cache(maxsize=1)
def measured_shells() -> tuple[ShellResidue, ...]:
    """The consensus pocket over :data:`DEPOSITS`, in ITPR3 numbering."""
    return tuple(consensus_shells({p: deposit_distances(load(p)) for p in DEPOSITS}))


def _column(gene: str, name: str) -> dict[int, float]:
    return {int(r["resi"]): G.as_float(r[name]) for r in per_residue(gene)}


@dataclass
class Pocket:
    """The measured pocket in one paralog's numbering, with its conservation."""
    paralog: str
    resi: np.ndarray            # residue numbers in this paralog
    structure_resi: np.ndarray  # the ITPR3 residue each was measured on
    distance: np.ndarray        # consensus median distance, Å
    shell: np.ndarray           # shell name per residue
    jsd: np.ndarray             # S17 deep JSD (NaN = not scored)
    frac_modal: np.ndarray
    whole_protein: np.ndarray   # every scored residue's deep JSD
    n_unaligned: int            # pocket residues with no aligned partner

    def in_shell(self, name: str) -> np.ndarray:
        return self.shell == name


@lru_cache(maxsize=None)
def pocket(gene: str) -> Pocket:
    shells = measured_shells()
    move = paralog_transfer(STRUCTURE_PARALOG, gene)
    kept = [(move[s.resi], s) for s in shells if s.resi in move]
    jsd, fm = _column(gene, "deep_jsd"), _column(gene, "deep_frac_modal")
    resi = np.array([r for r, _ in kept], int)
    whole = np.array([v for v in jsd.values() if np.isfinite(v)])
    return Pocket(gene, resi, np.array([s.resi for _, s in kept], int),
                  np.array([s.median for _, s in kept]),
                  np.array([s.shell for _, s in kept]),
                  np.array([jsd.get(r, np.nan) for r in resi]),
                  np.array([fm.get(r, np.nan) for r in resi]),
                  whole, len(shells) - len(kept))


def _mean(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else float("nan")


def shell_rows(gene: str) -> list[dict]:
    """One row per shell, the fields of S22's ``shell_constraint.tsv``."""
    pk = pocket(gene)
    rows = []
    for name in SHELLS:
        m = pk.in_shell(name)
        v = pk.jsd[m]
        rows.append({
            "paralog": gene, "shell": name, "n_residues": int(np.isfinite(v).sum()),
            "mean_jsd": _mean(v),
            "median_jsd": float(np.nanmedian(v)) if np.isfinite(v).any() else float("nan"),
            "mean_frac_modal": _mean(pk.frac_modal[m]),
            "whole_protein_mean_jsd": _mean(pk.whole_protein),
            "p_greater_than_protein": mann_whitney_greater(v, pk.whole_protein)["p"],
        })
    return rows


def trend(gene: str) -> dict:
    """Spearman of conservation (JSD and modal fraction) against distance."""
    pk = pocket(gene)
    j, f = spearman(pk.distance, pk.jsd), spearman(pk.distance, pk.frac_modal)
    return {"n_residues": j["n"], "rho_distance_vs_jsd": j["rho"],
            "p_distance_vs_jsd": j["p"], "rho_distance_vs_frac_modal": f["rho"],
            "p_distance_vs_frac_modal": f["p"]}


def contact_step(gene: str) -> dict:
    """The drop in mean conservation across each shell boundary, and a
    one-sided test of the contact shell against the second."""
    pk = pocket(gene)
    means = {s: _mean(pk.jsd[pk.in_shell(s)]) for s in SHELLS}
    drops = {f"{a}→{b}": means[a] - means[b] for a, b in zip(SHELLS, SHELLS[1:])}
    test = mann_whitney_greater(pk.jsd[pk.in_shell("contact")],
                                pk.jsd[pk.in_shell("second")])
    return {"means": means, "drops": drops, "p_contact_gt_second": test["p"],
            "largest_drop": max(drops, key=drops.get)}


def clear_caches() -> None:
    """Forget the joined tables (the measured structures are kept: they do
    not come from ip3r_genes)."""
    pocket.cache_clear()


def _parameters_changed() -> None:
    """The shells hang on registered edges: forget both levels."""
    measured_shells.cache_clear()
    pocket.cache_clear()


_P.subscribe(_parameters_changed)
