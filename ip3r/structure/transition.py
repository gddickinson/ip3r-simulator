"""Two states of one paralog, reduced to a common, matched, superposed basis.

Before any interpolation or displacement means anything, three things must
hold, and each fails silently if skipped:

1. **Same residue, same site.** Both deposits must be in the same human
   paralog numbering (the S24 rule, :mod:`ip3r.structure.numbering`), and the
   basis is restricted to residues whose C-alpha is resolved on all four
   subunits of *both* deposits, whose amino acid matches the reference, and
   which are not stubbed (backbone + CB with an unassigned sequence — 6DQN's
   1434-1546 segment is the example; its residue numbers are not evidence).
2. **Same subunit, same block.** Subunits are put in right-handed order about
   the cytosol-up axis of each deposit; the four cyclic correspondences are
   all tried and the best-fitting kept. For a C4-symmetric pair the four fit
   equally well (8TKG→8TKF: 14.18 Å each), which is a statement about the
   symmetry and not a failure; the tie is broken toward the deposited chain
   labels and all four RMSDs are reported.
3. **One frame.** The start is taken **as deposited** — the coordinates the
   viewer draws — and the end is superposed onto it, so a displacement field
   can be added to what is on screen (the PIEZO1 morph once built its path in
   one frame and drew it in another, and landed 36 Å from its own endpoint).

``fit="pore"`` (default) superposes on the PF00520 channel domain, holding
the membrane still as it is in a cell, so the displacement reads as what the
cytosolic assembly does relative to the pore. ``fit="global"`` minimises the
overall RMSD. Neither changes the *shape* change; they differ by a rigid
body, which the mode overlap removes explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.annotations import element_array, elements, reference_sequence
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .numbering import best_numbering, chain_residues, check_numbering
from .symmetry import Frame, kabsch, subunit_ca, tetramer_frame

__all__ = ["Transition", "TransitionUnavailable", "prepare_transition",
           "atom_site_index", "displaced_coords", "atom_displacement", "FITS"]

FITS = ("pore", "global")


class TransitionUnavailable(ValueError):
    """The two deposits cannot be put on a common residue basis."""


@dataclass
class Transition:
    start_id: str
    end_id: str
    paralog: str
    residues: np.ndarray            # (m,) canonical residue numbers
    chains_start: list              # subunit blocks, start deposit
    chains_end: list                # the matched end chains, same order
    start: np.ndarray               # (4m, 3) as deposited
    end: np.ndarray                 # (4m, 3) superposed onto start
    fit: str
    fit_rmsd: float                 # over the fitted sites
    rmsd: float                     # over all sites, after the fit
    frame: Frame                    # the start deposit's axis frame
    meta: dict = field(default_factory=dict)

    @property
    def n_subunits(self) -> int:
        return len(self.chains_start)

    @property
    def displacement(self) -> np.ndarray:
        return self.end - self.start

    def site_distance(self) -> np.ndarray:
        """|end - start| per site (Å)."""
        return np.linalg.norm(self.displacement, axis=1)

    def residue_distance(self) -> np.ndarray:
        """Per residue, the mean over the four subunits of |displacement|."""
        return self.site_distance().reshape(self.n_subunits, -1).mean(0)

    def element_means(self) -> dict[str, float]:
        """Mean residue displacement per functional element (Å)."""
        names = element_array(self.paralog, self.residues)
        d = self.residue_distance()
        return {n: float(d[names == n].mean()) for n in dict.fromkeys(names)}


def _usable(st: Structure, chain: str, seq: str) -> set[int]:
    """Residues of a chain that are resolved, unstubbed and match ``seq``."""
    return {r for r, (aa, stub) in chain_residues(st, chain).items()
            if 1 <= r <= len(seq) and not stub and aa == seq[r - 1]}


def _paralog_of(st: Structure) -> str:
    num = best_numbering(st)
    if num is None:
        raise TransitionUnavailable(
            f"{st.name} is in no human paralog numbering; residues cannot be "
            "matched by number")
    return num.paralog


def _fit_sites(fit: str, paralog: str, residues: np.ndarray, n: int) -> np.ndarray:
    if fit == "global":
        return np.ones(n * len(residues), bool)
    pore = [e for e in elements(paralog) if e.name == "channel"][0]
    return np.tile((residues >= pore.start) & (residues <= pore.end), n)


def _superpose(mobile, target, sel):
    r, t = kabsch(mobile[sel], target[sel])
    moved = mobile @ r.T + t
    rms = lambda m: float(np.sqrt(((moved[m] - target[m]) ** 2).sum(1).mean()))
    return moved, rms(sel), rms(np.ones(len(target), bool)), (r, t)


def prepare_transition(st_start: Structure, st_end: Structure,
                       fit: str = "pore") -> Transition:
    """The common basis of two deposits, end superposed onto start."""
    if fit not in FITS:
        raise ValueError(f"unknown fit {fit!r}; one of {FITS}")
    paralog = _paralog_of(st_start)
    if _paralog_of(st_end) != paralog:
        raise TransitionUnavailable(
            f"{st_start.name} and {st_end.name} are different paralogs")
    seq = reference_sequence(paralog)
    fa, fb = tetramer_frame(st_start), tetramer_frame(st_end)
    usable = [_usable(st_start, c, seq) for c in fa.chains] + \
             [_usable(st_end, c, seq) for c in fb.chains]
    residues = np.array(sorted(set.intersection(*usable)))
    if len(residues) < _P.value("transition.min_residues"):
        raise TransitionUnavailable(
            f"only {len(residues)} residues common to both deposits")
    ca_a, ca_b = subunit_ca(st_start), subunit_ca(st_end)
    start = np.vstack([[ca_a[c][r] for r in residues] for c in fa.chains])
    n = len(fa.chains)
    sel = _fit_sites(fit, paralog, residues, n)

    tol = _P.value("transition.tie_tolerance")
    trials = []
    for shift in range(n):
        order = fb.chains[shift:] + fb.chains[:shift]
        end = np.vstack([[ca_b[c][r] for r in residues] for c in order])
        moved, fit_rms, rms, rt = _superpose(end, start, sel)
        relabelled = sum(a != b for a, b in zip(fa.chains, order))
        trials.append((fit_rms, relabelled, order, moved, rms, rt))
    best_rms = min(t[0] for t in trials)
    fit_rms, _, order, moved, rms, rt = min(
        (t for t in trials if t[0] <= best_rms + tol), key=lambda t: t[1])

    excluded = {c: len(set(ca_a[c]) - set(residues)) for c in fa.chains}
    meta = {
        "fit_rmsd_by_shift": [round(t[0], 3) for t in trials],
        "correspondence_determined": sum(t[0] <= best_rms + tol for t in trials) == 1,
        "numbering_identity": (check_numbering(st_start, paralog).identity,
                               check_numbering(st_end, paralog).identity),
        "start_residues_outside_basis": excluded,
        "n_fit_sites": int(sel.sum()),
        "end_transform": rt,    # (R, t): end deposit xyz @ R.T + t = this frame
    }
    return Transition(st_start.name, st_end.name, paralog, residues, list(fa.chains),
                      list(order), start, moved, fit, fit_rms, rms, fa, meta)


def atom_site_index(st: Structure, tr: Transition) -> np.ndarray:
    """For every atom of the start deposit, the basis site it moves with.

    An atom of a basis residue takes its own residue's site, so whole
    residues move together. Every other atom — residues outside the basis,
    IP3, ions, lipids — takes the **nearest site in space**, never a site
    found by residue number (the PIEZO1 morph tied lipids to a C-alpha
    64.8 Å away that way).
    """
    from scipy.spatial import cKDTree
    if st.name != tr.start_id:
        raise ValueError(f"{st.name} is not the start of {tr.start_id}->{tr.end_id}")
    m = len(tr.residues)
    out = np.full(st.n_atoms, -1, np.int64)
    pos = {int(r): i for i, r in enumerate(tr.residues)}
    for k, ch in enumerate(tr.chains_start):
        sel = np.flatnonzero((st.chain == ch) & ~st.hetero)
        idx = np.array([pos.get(int(r), -1) for r in st.res_seq[sel]])
        ok = idx >= 0
        out[sel[ok]] = k * m + idx[ok]
    rest = np.flatnonzero(out < 0)
    if len(rest):
        _, near = cKDTree(tr.start).query(st.xyz[rest].astype(np.float64))
        out[rest] = near
    return out


def displaced_coords(base_xyz: np.ndarray, frames: np.ndarray, i: int,
                     site_index: np.ndarray) -> np.ndarray:
    """The start deposit's atoms carried to frame ``i`` of a site path.

    ``frames[0]`` is the start's own C-alphas, so frame 0 returns
    ``base_xyz`` exactly and every basis C-alpha lands exactly on
    ``frames[i]``. Side chains and ligands are carried rigidly with their
    site: only the C-alpha trace is the end deposit.
    """
    return np.asarray(base_xyz, np.float64) + (frames[i] - frames[0])[site_index]


def atom_displacement(st: Structure, tr: Transition, site_index: np.ndarray) -> np.ndarray:
    """Per atom, its residue's mean displacement (Å); NaN off the basis.

    An atom that only rides the nearest site was not measured, so it is
    painted grey rather than given its neighbour's value.
    """
    per_res = dict(zip(tr.residues.tolist(), tr.residue_distance()))
    out = np.full(st.n_atoms, np.nan)
    for ch in tr.chains_start:
        sel = np.flatnonzero((st.chain == ch) & ~st.hetero)
        out[sel] = [per_res.get(int(r), np.nan) for r in st.res_seq[sel]]
    return out
