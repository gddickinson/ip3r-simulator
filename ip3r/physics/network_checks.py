"""Is the elastic network a model of the receptor, or of how it was sampled?

Round 2's overlap (``transition_modes``) left two questions open; this module
answers both by measuring the network against itself.

**Where do the local modes come from?** :func:`local_modes` names the residue
each low-collectivity mode sits on (its squared amplitude summed over the four
subunits) and any unresolved stretch between that residue's sampled
neighbours. On 8TKG at stride 3 all five (#11-15) are one piece, residue 86,
the first residue after the unresolved 77-85 loop: four copies of one flap
give an A, a B and an E pair. At strides 1 and 2 the flap keeps enough springs
and there are no local modes; RyR1 9R8O has its own (residue 896, even at
stride 1). Sequence-neighbour springs do not remove them at any strength,
because the flap swings as a body rather than stretching along its chain, so
the network is not bridged; the modes are named and skipped.

**Is the lowest-A-mode number a property of the receptor?** :func:`cutoff_scan`
re-solves the network over the registered cutoff grid. On 8TKG at stride 1
the lowest collective A mode's overlap falls from 0.48 (12 A) to 0.24
(21.6 A) while the three lowest A modes together hold 0.66-0.68: the cutoff
redistributes the move among near-degenerate A modes. The subspace is the
result; the single mode is a cutoff choice.

:func:`stride_agreement` scores a strided network against stride 1 by the
root-mean-square inner product of the first k modes (RMSIP, Amadei et al.
1999), restricted to the sites both share.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from ..structure.transition import Transition
from .anm import ModeSet

__all__ = ["LocalMode", "local_modes", "describe_local", "CutoffRow",
           "cutoff_grid", "cutoff_scan", "rmsip", "stride_agreement"]


@dataclass(frozen=True)
class LocalMode:
    index: int                      # 0-based mode index
    irrep: str
    kappa: float
    residue: int                    # the residue carrying most of the mode
    share: float                    # its share of |u|^2, all subunits
    gap: tuple[int, int] | None     # unresolved stretch beside it


def _gap_beside(residue: int, sampled: np.ndarray,
                resolved: np.ndarray | None) -> tuple[int, int] | None:
    """The unresolved run between ``residue``'s sampled neighbours, if any."""
    if resolved is None:
        return None
    j = int(np.searchsorted(sampled, residue))
    lo = int(sampled[j - 1]) if j > 0 else int(resolved[0]) - 1
    hi = int(sampled[j + 1]) if j + 1 < len(sampled) else int(resolved[-1]) + 1
    have = set(int(r) for r in resolved)
    missing = [r for r in range(lo + 1, hi) if r not in have]
    if not missing:
        return None
    # The run nearest the residue (a sampled window can hold two).
    runs, start = [], missing[0]
    for a, b in zip(missing, missing[1:] + [None]):
        if b != a + 1:
            runs.append((start, a))
            start = b
    return min(runs, key=lambda r: min(abs(r[0] - residue), abs(r[1] - residue)))


def local_modes(modes: ModeSet, residues: np.ndarray,
                resolved: np.ndarray | None = None) -> list[LocalMode]:
    """Every mode below the collectivity threshold, located on the sequence.

    ``residues`` are the sampled residues (one subunit block, the ANM's site
    order); ``resolved`` the unstrided residues all subunits resolve, used to
    find the missing stretch beside the flap.
    """
    residues = np.asarray(residues)
    per = len(residues)
    n_sub = modes.vectors.shape[1] // per
    kappa = modes.collectivity()
    out = []
    for i in np.flatnonzero(~modes.is_collective()):
        u2 = np.einsum("ij,ij->i", modes.vectors[i], modes.vectors[i])
        by_res = u2.reshape(n_sub, per).sum(0) / u2.sum()
        j = int(np.argmax(by_res))
        sym = modes.symmetry[i] if modes.symmetry is not None else "?"
        out.append(LocalMode(int(i), str(sym), float(kappa[i]), int(residues[j]),
                             float(by_res[j]),
                             _gap_beside(int(residues[j]), residues, resolved)))
    return out


def describe_local(local: list[LocalMode]) -> list[str]:
    """One phrase per residue: ``#11-15 (B E E A B) on residue 86, after ...``."""
    lines = []
    for res in dict.fromkeys(m.residue for m in local):
        group = [m for m in local if m.residue == res]
        idx = [m.index + 1 for m in group]
        span = f"#{idx[0]}" if len(idx) == 1 else (
            f"#{idx[0]}-{idx[-1]}" if idx == list(range(idx[0], idx[-1] + 1))
            else "#" + ", #".join(map(str, idx)))
        gap = group[0].gap
        where = "" if gap is None else (
            f", after unresolved {gap[0]}-{gap[1]}" if gap[1] < res
            else f", before unresolved {gap[0]}-{gap[1]}")
        lines.append(f"{span} ({' '.join(m.irrep for m in group)}) on residue {res}"
                     f"{where}")
    return lines


@dataclass(frozen=True)
class CutoffRow:
    cutoff: float
    n_local: int
    lowest_a: int | None            # 0-based index of the lowest collective A mode
    lowest_a_overlap: float         # NaN when there is none
    collective_a: float             # all collective A modes together
    cumulative: float               # every computed mode
    null: float                     # random direction of the same irreps


def cutoff_grid() -> np.ndarray:
    lo, hi = _P.value("anm.scan_cutoff_low"), _P.value("anm.scan_cutoff_high")
    step = _P.value("anm.scan_cutoff_step")
    return np.round(np.arange(lo, hi + step / 2, step), 6)


def cutoff_scan(tr: Transition, cutoffs=None, stride: int | None = None,
                reference: str = "start", n_modes: int | None = None
                ) -> list[CutoffRow]:
    """The transition overlap re-measured at each cutoff."""
    from .transition_modes import transition_overlap
    rows = []
    for c in (cutoff_grid() if cutoffs is None else cutoffs):
        ov = transition_overlap(tr, reference, stride=stride, n_modes=n_modes,
                                cutoff=float(c))
        coll = ov.modes.is_collective()
        a = ov.lowest_a
        mask = coll & (ov.modes.symmetry == "A")
        rows.append(CutoffRow(float(c), int((~coll).sum()), a,
                              float(ov.overlap[a]) if a is not None else float("nan"),
                              float(np.sqrt(np.sum(ov.overlap[mask] ** 2))),
                              float(ov.cumulative[-1]), float(ov.null_cumulative[-1])))
    return rows


def rmsip(a: np.ndarray, b: np.ndarray, k: int) -> float:
    """Root-mean-square inner product of two sets of ``k`` unit vectors."""
    a = a[:k] / np.linalg.norm(a[:k], axis=1, keepdims=True)
    b = b[:k] / np.linalg.norm(b[:k], axis=1, keepdims=True)
    return float(np.sqrt(np.sum((a @ b.T) ** 2) / k))


def stride_agreement(tr: Transition, strides=(2, 3, 4), reference_stride: int = 1,
                     reference: str = "start", n_modes: int | None = None
                     ) -> dict[int, float]:
    """RMSIP of each strided network's modes with the reference stride's.

    The reference modes are restricted to the coarse network's sites (which
    are a subset when the strides divide) and renormalised.
    """
    from .transition_modes import transition_overlap
    ref = transition_overlap(tr, reference, stride=reference_stride, n_modes=n_modes)
    out = {}
    for s in strides:
        ov = transition_overlap(tr, reference, stride=s, n_modes=n_modes)
        pos = np.searchsorted(ref.residues, ov.residues)
        if not np.array_equal(ref.residues[np.minimum(pos, len(ref.residues) - 1)],
                              ov.residues):
            raise ValueError(f"stride {s} sites are not a subset of stride "
                             f"{reference_stride}'s")
        per, n = len(ref.residues), tr.n_subunits
        idx = np.concatenate([k * per + pos for k in range(n)])
        a = ref.modes.vectors[:, idx].reshape(ref.modes.n_modes, -1)
        b = ov.modes.vectors.reshape(ov.modes.n_modes, -1)
        out[int(s)] = rmsip(a, b, min(len(a), len(b)))
    return out
