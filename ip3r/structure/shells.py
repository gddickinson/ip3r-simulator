"""Ligand shells: how far every residue is from IP3, over all its atoms.

S22's measurement, rebuilt on this project's reader. Per deposit, a residue's
distance is its minimum all-atom distance (hydrogens included, as S22's
reader keeps them) to the IP3 bound on **its own subunit**, the best over the
four subunits; the consensus is the median over the IP3-bound depositions
that resolve it. Residues beyond ``ligand.shell_radius`` are absent rather
than pooled into an open bin.

The shells are ``contact`` (< ligand.contact_cutoff), ``second``, ``third``
and ``fourth`` (edges ``ligand.shell_second_edge``, ``shell_third_edge``,
``shell_radius``). A Cα trace would put an arginine several ångström
further from a phosphate than the guanidinium that binds it, which is why
all atoms are used.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .ligand import ligand_sites, residue_distances

__all__ = ["SHELLS", "ShellResidue", "shell_edges", "shell_of",
           "deposit_distances", "chain_distances", "consensus_shells",
           "atom_ligand_distance"]

SHELLS = ("contact", "second", "third", "fourth")


def shell_edges() -> tuple[float, float, float, float]:
    """Outer edge of each shell, in Å, read from the registry at call time."""
    return (_P.value("ligand.contact_cutoff"), _P.value("ligand.shell_second_edge"),
            _P.value("ligand.shell_third_edge"), _P.value("ligand.shell_radius"))


def shell_of(d: float) -> str | None:
    """The shell a distance falls in; ``None`` beyond the search radius."""
    for name, hi in zip(SHELLS, shell_edges()):
        if d < hi:
            return name
    return None


def chain_distances(st: Structure) -> dict[str, dict[int, float]]:
    """Per protein chain: residue → all-atom distance to that chain's own IP3.

    A chain whose site is empty is absent — the same residue beside an empty
    site is not evidence of anything. Used to paint the displayed structure.
    """
    radius = _P.value("ligand.shell_radius")
    out: dict[str, dict[int, float]] = {}
    for site in ligand_sites(st):
        d = residue_distances(st, site, radius=radius, heavy_only=False)
        cur = out.setdefault(site.subunit, {})
        for r, v in d.items():
            cur[r] = min(cur.get(r, np.inf), v)
    return out


def deposit_distances(st: Structure) -> dict[int, float]:
    """Residue → best distance over the subunits (S22's per-deposit value)."""
    best: dict[int, float] = {}
    for per in chain_distances(st).values():
        for r, v in per.items():
            best[r] = min(best.get(r, np.inf), v)
    return dict(sorted(best.items()))


@dataclass(frozen=True)
class ShellResidue:
    resi: int
    median: float                   # Å, over the depositions resolving it
    n_structures: int
    n_contact: int                  # depositions placing it inside the cutoff
    shell: str
    per_deposit: dict

    @property
    def spread(self) -> float:
        v = list(self.per_deposit.values())
        return max(v) - min(v)


def consensus_shells(per_deposit: dict[str, dict[int, float]]) -> list[ShellResidue]:
    """The consensus over several depositions (``pdb → deposit_distances``).

    Taking a dict of already-measured deposits keeps this pure; the caller
    loads and measures (on a worker, in the GUI).
    """
    cutoff = _P.value("ligand.contact_cutoff")
    out = []
    for r in sorted({r for m in per_deposit.values() for r in m}):
        ds = {p: m[r] for p, m in per_deposit.items() if r in m}
        med = float(np.median(list(ds.values())))
        shell = shell_of(med)
        if shell is None:
            continue
        out.append(ShellResidue(r, med, len(ds), sum(v <= cutoff for v in ds.values()),
                                shell, ds))
    return out


def atom_ligand_distance(st: Structure) -> np.ndarray:
    """Per atom: its residue's distance to its own subunit's IP3 (Å).

    NaN for hetero atoms, residues beyond the search radius and subunits
    with no bound IP3 — drawn grey, never as the far end of a scale.
    """
    out = np.full(st.n_atoms, np.nan)
    prot = ~st.hetero
    for chain, per in chain_distances(st).items():
        if not per:
            continue
        sel = np.flatnonzero(prot & (st.chain == chain))
        keys = np.fromiter(per.keys(), int)
        vals = np.fromiter(per.values(), float)
        order = np.argsort(keys)
        keys, vals = keys[order], vals[order]
        rs = st.res_seq[sel].astype(int)
        pos = np.clip(np.searchsorted(keys, rs), 0, len(keys) - 1)
        hit = keys[pos] == rs
        out[sel[hit]] = vals[pos[hit]]
    return out
