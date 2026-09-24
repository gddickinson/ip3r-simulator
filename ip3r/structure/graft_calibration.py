"""How good is a fill? Hide stretches a deposit resolves, fill them, compare.

The stretches tested are the ones that matter: residues some *other* ITPR3
deposit leaves unresolved, but the host deposit builds. Each is cut out of the
host, filled from the prediction by exactly the route :mod:`.graft` uses, and
scored against the host's own coordinates (C-alpha RMSD).

Two baselines use no prediction, and a fill must beat them to be worth
drawing:

* ``line``: C-alphas evenly spaced on the straight line between the two
  flanking C-alphas. It is the least a fill could claim.
* ``global``: the prediction superposed on the whole chain at once, not on
  the stretch's own anchors. It shows what the local fit buys.

The seam distances of fills whose truth is known also calibrate
``graft.join_tolerance``: a tolerance that true fills fail is too tight.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.structure import Structure
from .graft import Stretch, fill_stretches, placed_ca, unresolved
from .numbering import chain_residues
from .symmetry import kabsch

__all__ = ["Trial", "candidate_stretches", "hide", "trial", "calibrate"]


@dataclass(frozen=True)
class Trial:
    host: str
    stretch: Stretch
    n: int                        # C-alphas compared
    rmsd_fill: float
    rmsd_line: float
    rmsd_global: float
    plddt: float
    joins: tuple
    anchor_rmsd: float

    @property
    def beats_line(self) -> bool:
        return self.rmsd_fill < self.rmsd_line


def candidate_stretches(host: Structure, others: list[Structure],
                        chain: str = "A") -> list[tuple[int, int]]:
    """Gaps of the other deposits that ``host`` resolves entirely, unstubbed."""
    res = chain_residues(host, chain)
    out = set()
    for o in others:
        oc = max(o.chains, key=lambda c: int((o.mask_ca() & (o.chain == c)).sum()))
        for s in unresolved(o, oc, pred_last=0):
            span = range(s.first - 1, s.last + 2)       # the seams too
            if all(r in res and not res[r][1] for r in span):
                out.add((s.first, s.last))
    return sorted(out)


def hide(st: Structure, chain: str, first: int, last: int) -> Structure:
    """``st`` without one stretch of one chain."""
    drop = (st.chain == chain) & (st.res_seq >= first) & (st.res_seq <= last) & ~st.hetero
    out = st.subset(~drop, name=st.name)
    out.meta = dict(st.meta)
    return out


def _true_ca(st: Structure, chain: str) -> dict[int, np.ndarray]:
    m = st.mask_ca() & (st.chain == chain) & ~st.hetero
    return {int(r): x.astype(np.float64) for r, x in zip(st.res_seq[m], st.xyz[m])}


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def trial(host: Structure, pred: Structure, chain: str, first: int, last: int,
          global_fit: tuple | None = None) -> Trial | None:
    """One hidden stretch, filled and scored; None if it could not be filled."""
    truth = _true_ca(host, chain)
    model = fill_stretches(hide(host, chain, first, last), pred,
                           [Stretch(chain, first, last, "gap")])
    if not model.fills:
        return None
    f = model.fills[0]
    got = placed_ca(model, f)
    res = [r for r in range(first, last + 1) if r in got and r in truth]
    true = np.array([truth[r] for r in res])
    fill = np.array([got[r] for r in res])
    a, b = truth[first - 1], truth[last + 1]
    t = (np.arange(1, last - first + 2) / (last - first + 2))[:, None]
    line = (a + t * (b - a))[[r - first for r in res]]
    rmsd_global = float("nan")
    if global_fit is not None:
        r_, t_, pca = global_fit
        g = np.array([pca[r] for r in res]) @ r_.T + t_
        rmsd_global = _rmsd(g, true)
    return Trial(host.name, f.stretch, len(res), _rmsd(fill, true), _rmsd(line, true),
                 rmsd_global, f.plddt, f.joins, f.anchor_rmsd)


def _global_fit(host: Structure, pred: Structure, chain: str) -> tuple:
    truth = _true_ca(host, chain)
    pca = _true_ca(pred, pred.chains[0])
    dep = chain_residues(host, chain)
    shared = [r for r in truth if r in pca and not dep[r][1]]
    r_, t_ = kabsch(np.array([pca[r] for r in shared]), np.array([truth[r] for r in shared]))
    return r_, t_, pca


def calibrate(host: Structure, others: list[Structure], pred: Structure,
              chain: str = "A") -> list[Trial]:
    """Every candidate stretch of ``host``, hidden, filled and scored."""
    g = _global_fit(host, pred, chain)
    out = []
    for first, last in candidate_stretches(host, others, chain):
        t = trial(host, pred, chain, first, last, g)
        if t is not None:
            out.append(t)
    return out
