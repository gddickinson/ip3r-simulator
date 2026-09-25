"""The second pKa route: PROPKA 3 (Olsson et al. 2011) on the deposit.

PROPKA is an empirical predictor with its own desolvation, hydrogen-bond
and charge-charge terms, calibrated on measured protein pKas. It shares no
code and no constants with the network in :mod:`.pka`, so where the two
agree the answer does not rest on either model's choices. Each group gets
one pKa with its neighbours in their default states, so PROPKA cannot show
a ring of four titrating together. That is what the network adds.

The deposit is written out as PDB text: protein heavy atoms only, and only
whole residues with an atom within ``pka.propka_radius`` of a group asked
about. PROPKA's longest range is its desolvation count, so a context that
reaches well past it gives the same pKa as the whole tetramer (tested at
25 and 35 A). Chain IDs must be one character; a longer one is refused.

``propka`` is imported when a route is run, so the rest of the package does
not need it installed.
"""

from __future__ import annotations

import io
import logging

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["context_mask", "pdb_text", "propka_pka", "charge_at"]

_ACIDS = {"ASP", "GLU"}
_TYPES = {"ASP", "GLU", "HIS", "LYS", "ARG"}


def context_mask(st: Structure, around: np.ndarray,
                 radius: float | None = None) -> np.ndarray:
    """Every atom of each protein residue with a heavy atom within
    ``radius`` A of a point of ``around``."""
    radius = _P.value("pka.propka_radius") if radius is None else radius
    heavy = ~st.hetero & (st.element != "H")
    idx = np.flatnonzero(heavy)
    near = np.zeros(len(idx), dtype=bool)
    for p in np.atleast_2d(np.asarray(around, dtype=float)):
        near |= np.einsum("ij,ij->i", st.xyz[idx] - p, st.xyz[idx] - p) <= radius ** 2
    keys = set(zip(st.chain[idx[near]].tolist(), st.res_seq[idx[near]].tolist()))
    whole = np.array([k in keys for k in zip(st.chain.tolist(),
                                              st.res_seq.tolist())])
    return whole & heavy


def pdb_text(st: Structure, mask: np.ndarray) -> str:
    """PDB ATOM records for ``mask`` (serials renumbered)."""
    lines = []
    for n, i in enumerate(np.flatnonzero(mask), start=1):
        chain = str(st.chain[i])
        if len(chain) != 1:
            raise ValueError(f"chain ID {chain!r} does not fit a PDB record")
        name = str(st.atom_name[i])
        name = f" {name:<3s}" if len(name) < 4 else name[:4]
        x, y, z = (float(v) for v in st.xyz[i])
        lines.append(f"ATOM  {n % 100000:5d} {name:4s} {st.res_name[i]:>3s} "
                     f"{chain}{int(st.res_seq[i]):4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                     f"  1.00  0.00          {str(st.element[i]):>2s}")
    return "\n".join(lines) + "\nEND\n"


def propka_pka(st: Structure, around: np.ndarray,
               radius: float | None = None) -> dict[tuple[str, int], tuple[str, float]]:
    """``(chain, residue) -> (residue name, pKa)`` for every Asp, Glu, His,
    Lys and Arg side chain PROPKA reports in the context around ``around``."""
    import propka.run
    logging.getLogger("propka").setLevel(logging.ERROR)
    text = pdb_text(st, context_mask(st, around, radius))
    mol = propka.run.single("context.pdb", stream=io.StringIO(text),
                            write_pka=False)
    out = {}
    for g in mol.conformations["AVR"].groups:
        # side chains only: PROPKA types both carboxylates and the chain's
        # C-terminus "COO", so the termini are told apart by class
        if g.residue_type in _TYPES and g.atom is not None \
                and type(g).__name__ not in ("CtermGroup", "NtermGroup"):
            out[(str(g.atom.chain_id), int(g.atom.res_num))] = (
                g.residue_type, float(g.pka_value))
    return out


def charge_at(res_name: str, pka: float, ph: float) -> float:
    """Henderson-Hasselbalch mean charge of one group."""
    protonated = 1.0 / (1.0 + 10.0 ** (ph - pka))
    return protonated - 1.0 if res_name in _ACIDS else protonated
