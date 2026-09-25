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

**Long stretches** (Round 7.9). The real gaps are 20-70 residues and mostly
below pLDDT 50, and the stretches above are 4-10 residues at pLDDT 53-77.
:func:`windows` cuts resolved windows of a chosen length out of the host.
Those are ordered residues (pLDDT ~80), so they test length and not the
real gaps' regime. No deposit resolves a run below pLDDT 50 with its own
anchors either side. What the maps do resolve are **islands**
(:func:`islands`): short resolved runs with an unresolved stretch on each
side. Hiding one and filling the whole span from the outer anchors is the
real case, a 30-140-residue fill scored where the experiment saw residues.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.structure import Structure
from .graft import Stretch, fill_stretches, placed_ca, unresolved
from .graft_numbering import BY_NUMBER, NumberMap, runs
from .numbering import chain_residues
from .symmetry import kabsch

__all__ = ["Trial", "candidate_stretches", "hide", "trial", "calibrate",
           "windows", "islands", "global_fit", "window_trials", "island_trials"]


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
    plddt_scored: float = float("nan")    # mean over the compared residues only

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
          global_fit: tuple | None = None, numbering: NumberMap = BY_NUMBER,
          span: tuple[int, int] | None = None) -> Trial | None:
    """One hidden stretch, filled and scored; None if it could not be filled.

    ``span`` is the stretch filled, when it is wider than the one hidden (an
    island, :func:`islands`). The score is always on the hidden residues.
    """
    truth = _true_ca(host, chain)
    lo, hi = span or (first, last)
    model = fill_stretches(hide(host, chain, first, last), pred,
                           [Stretch(chain, lo, hi, "gap")], numbering=numbering)
    if not model.fills:
        return None
    f = model.fills[0]
    got = placed_ca(model, f)
    res = [r for r in range(first, last + 1) if r in got and r in truth]
    true = np.array([truth[r] for r in res])
    fill = np.array([got[r] for r in res])
    a, b = truth[lo - 1], truth[hi + 1]
    t = (np.arange(1, hi - lo + 2) / (hi - lo + 2))[:, None]
    line = (a + t * (b - a))[[r - lo for r in res]]
    rmsd_global = float("nan")
    if global_fit is not None:
        r_, t_, pca = global_fit
        g = np.array([pca[r] for r in res]) @ r_.T + t_
        rmsd_global = _rmsd(g, true)
    ca = pred.mask_ca() & ~pred.hetero
    plddt = dict(zip(pred.res_seq[ca].tolist(), pred.b_factor[ca].tolist()))
    scored = float(np.mean([plddt[numbering(r)] for r in res]))
    return Trial(host.name, f.stretch, len(res), _rmsd(fill, true), _rmsd(line, true),
                 rmsd_global, f.plddt, f.joins, f.anchor_rmsd, scored)


def global_fit(host: Structure, pred: Structure, chain: str,
               numbering: NumberMap = BY_NUMBER) -> tuple:
    """The prediction superposed on the whole chain: ``(R, t, deposit number -> C-alpha)``."""
    truth = _true_ca(host, chain)
    own = _true_ca(pred, pred.chains[0])
    pca = {r: own[numbering(r)] for r in truth if numbering(r) in own}
    dep = chain_residues(host, chain)
    shared = [r for r in truth if r in pca and not dep[r][1]]
    r_, t_ = kabsch(np.array([pca[r] for r in shared]), np.array([truth[r] for r in shared]))
    return r_, t_, pca


def calibrate(host: Structure, others: list[Structure], pred: Structure,
              chain: str = "A") -> list[Trial]:
    """Every candidate stretch of ``host``, hidden, filled and scored."""
    g = global_fit(host, pred, chain)
    out = []
    for first, last in candidate_stretches(host, others, chain):
        t = trial(host, pred, chain, first, last, g)
        if t is not None:
            out.append(t)
    return out


def _clean(host: Structure, chain: str, numbering: NumberMap, pred: Structure) -> set:
    """Residues resolved, unstubbed and paired with the same amino acid."""
    dep = chain_residues(host, chain)
    pres = chain_residues(pred, pred.chains[0])
    return {r for r, (aa, stub) in dep.items()
            if not stub and numbering(r) in pres and pres[numbering(r)][0] == aa}


def windows(host: Structure, pred: Structure, length: int, chain: str = "A",
            numbering: NumberMap = BY_NUMBER) -> list[tuple[int, int]]:
    """Resolved windows of ``length`` residues to hide, spaced two lengths apart.

    A window needs the anchor window clean on each side
    (``graft.anchor_window``), so the fill is judged on its own stretch and
    never on anchors that are themselves missing.
    """
    from ..parameters import PARAMETERS as _P
    w = int(_P.value("graft.anchor_window"))
    ok = _clean(host, chain, numbering, pred)
    out, r, top = [], min(ok) + w, max(ok)
    while r + length + w <= top:
        if all(x in ok for x in range(r - w, r + length + w)):
            out.append((r, r + length - 1))
            r += 3 * length
        else:
            r += 1
    return out


def islands(host: Structure, pred: Structure, chain: str = "A",
            numbering: NumberMap = BY_NUMBER, min_length: int = 5,
            max_length: int = 60) -> list[tuple[int, int, int, int]]:
    """Resolved runs of a loop's length with an unresolved stretch on each side.

    Each is returned as ``(first, last, lo, hi)``: the island, and the stretch
    a fill must span once it is hidden (both neighbouring gaps and the island
    between them). Hiding the island and filling ``lo..hi`` tests a long fill
    in a real gap, scored where the experiment saw residues.
    """
    seen = sorted(chain_residues(host, chain))
    ok = _clean(host, chain, numbering, pred)
    out = []
    for i, (a, b) in enumerate(runs(seen)):
        if i == 0 or b == seen[-1] or not min_length <= b - a + 1 <= max_length:
            continue
        if not all(r in ok for r in range(a, b + 1)):
            continue
        lo = max(r for r in seen if r < a) + 1
        hi = min(r for r in seen if r > b) - 1
        out.append((a, b, lo, hi))
    return out


def window_trials(host: Structure, pred: Structure, lengths=(10, 20, 30, 45, 60),
                  chain: str = "A", numbering: NumberMap = BY_NUMBER) -> list[Trial]:
    """Every :func:`windows` stretch of each length, hidden, filled and scored."""
    g = global_fit(host, pred, chain, numbering)
    out = []
    for n in lengths:
        for first, last in windows(host, pred, n, chain, numbering):
            t = trial(host, pred, chain, first, last, g, numbering)
            if t is not None:
                out.append(t)
    return out


def island_trials(host: Structure, pred: Structure, chain: str = "A",
                  numbering: NumberMap = BY_NUMBER) -> list[tuple[tuple, Trial]]:
    """Every island hidden and its span filled: ``((first, last), trial)``."""
    g = global_fit(host, pred, chain, numbering)
    out = []
    for first, last, lo, hi in islands(host, pred, chain, numbering):
        t = trial(host, pred, chain, first, last, g, numbering, span=(lo, hi))
        if t is not None:
            out.append(((first, last), t))
    return out
