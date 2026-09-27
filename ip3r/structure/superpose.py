"""A second deposit placed on the one shown, by the transition's own fit.

Round 7.20: any two states compared by eye without building a morph. The
fit is :func:`~ip3r.structure.transition.prepare_transition`'s, unchanged:
the residue-matched basis (same paralog numbering, resolved and unstubbed on
all eight chains), the best cyclic subunit correspondence, and the pore or
global superposition of the second deposit onto the first **as deposited**.
Its rigid transform (``meta["end_transform"]``) is then applied to *every*
atom of the second deposit — side chains, ligands, residues outside the
basis — so what is drawn is the deposit itself, moved, never a C-alpha trace
carried along.

The second deposit's subunits are relabelled to the chains they were matched
to, so "subunit A" means the same block in both and the Subunits toggles
hide the pair together. A chain outside the correspondence (none among the
registered deposits) keeps its label with a prime, so it can never be
mistaken for a matched subunit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.structure import Structure
from ..io.registry import load_registry
from .transition import Transition, prepare_transition

__all__ = ["Superposition", "superpose", "candidates"]


@dataclass
class Superposition:
    shown_id: str
    other_id: str
    structure: Structure            # the other deposit, moved, chains relabelled
    transition: Transition          # the fit it came from
    chain_map: dict                 # other deposit's chain -> the shown chain

    @property
    def fit(self) -> str:
        return self.transition.fit

    @property
    def n_residues(self) -> int:
        return len(self.transition.residues)

    def summary(self) -> str:
        tr = self.transition
        moved = {k: v for k, v in self.chain_map.items() if k != v}
        relabel = ("" if not moved else "; subunits matched "
                   + ", ".join(f"{k}→{v}" for k, v in moved.items()))
        return (f"{self.other_id} on {self.shown_id} ({tr.fit} fit, "
                f"{self.n_residues} residues × {tr.n_subunits}): RMSD "
                f"{tr.fit_rmsd:.2f} Å over the fitted sites, {tr.rmsd:.2f} Å "
                f"over all{relabel}")


def superpose(shown: Structure, other: Structure, fit: str = "pore") -> Superposition:
    """``other`` moved onto ``shown`` as deposited; raises
    :class:`~ip3r.structure.transition.TransitionUnavailable` when the two
    cannot share a residue basis (another paralog, no human numbering)."""
    tr = prepare_transition(shown, other, fit=fit)
    r, t = tr.meta["end_transform"]
    xyz = other.xyz.astype(np.float64) @ r.T + t
    chain_map = dict(zip(tr.chains_end, tr.chains_start))
    labels = np.array([chain_map.get(c, f"{c}'") for c in other.chain])
    moved = other.copy_with_coords(xyz)
    moved.chain = labels          # a bijection: residue boundaries unchanged
    return Superposition(shown.name, other.name, moved, tr, chain_map)


def candidates(pdb_id: str | None, paralog: str | None) -> list[tuple[str, str]]:
    """(pdb id, label) of every deposit that can be superposed on ``pdb_id``:
    the same paralog numbering, human or RyR1 — the Transition tab's rule."""
    if not pdb_id or not paralog:
        return []
    return [(e.pdb_id, f"{e.pdb_id} — {e.state}") for e in load_registry()
            if e.paralog == paralog and (e.human or e.family == "RyR")
            and e.pdb_id != pdb_id]
