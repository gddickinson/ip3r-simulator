"""Salt bridges: ion pairs whose charges cancel before any wall charge counts.

:mod:`ip3r.physics.pore_charge` admits a lining Asp, Glu, Lys or Arg as a
full formal charge. A carboxylate already paired with a guanidinium or an
ammonium is not a free charge acting on the lumen: the pair is net neutral
at the distance an ion in the pore sees it from. This module finds those
pairs so the lining rule can drop them.

**The rule is Barlow & Thornton's** (1983): an ion pair is an acidic and a
basic side chain with a charged oxygen (OD1/OD2, OE1/OE2) and a charged
nitrogen (Lys NZ; Arg NE, NH1, NH2) within ``pore_charge.salt_bridge_cutoff``
(4 Å). Histidine carries no charge here, as in the lining rule.

**One partner each.** An Arg can reach two carboxylates, but it has one
charge to cancel. Pairs are therefore matched one-to-one, closest first, so
every bridge removes exactly +1 and -1 and the charge the lumen loses is
always balanced. The partner is searched in the *whole* deposit, not only
among lining groups: D2478 of 8TKF is bridged by R2471 of the neighbouring
subunit, which does not line the pore itself.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["Bridge", "salt_bridges", "ACID_ATOMS", "BASE_ATOMS"]

#: Charged oxygens of the acidic side chains.
ACID_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}

#: Charged nitrogens of the basic side chains.
BASE_ATOMS = {"LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2")}

Key = tuple[str, int]                    # (chain, residue number)


@dataclass(frozen=True)
class Bridge:
    """One matched ion pair."""

    acid: Key
    acid_name: str
    base: Key
    base_name: str
    distance: float                      # A, closest charged N-O

    def members(self) -> tuple[Key, Key]:
        return self.acid, self.base

    def label(self) -> str:
        return (f"{self.acid_name}{self.acid[1]}/{self.acid[0]}-"
                f"{self.base_name}{self.base[1]}/{self.base[0]} "
                f"{self.distance:.2f} A")


def _charged_atoms(st: Structure, table: dict) -> tuple[np.ndarray, list, list]:
    """Coordinates of the charged atoms, with each atom's residue key and name."""
    m = np.zeros(len(st.xyz), bool)
    for name, atoms in table.items():
        m |= (st.res_name == name) & np.isin(st.atom_name, atoms)
    m &= ~st.hetero
    idx = np.flatnonzero(m)
    keys = [(str(st.chain[i]), int(st.res_seq[i])) for i in idx]
    return st.xyz[idx].astype(float), keys, [str(st.res_name[i]) for i in idx]


def salt_bridges(st: Structure, cutoff: float | None = None) -> list[Bridge]:
    """Every ion pair in the deposit, matched one-to-one, closest first."""
    cutoff = (_P.value("pore_charge.salt_bridge_cutoff") if cutoff is None
              else cutoff)
    ox, o_key, o_name = _charged_atoms(st, ACID_ATOMS)
    nx, n_key, n_name = _charged_atoms(st, BASE_ATOMS)
    if not len(ox) or not len(nx):
        return []
    best: dict[tuple[Key, Key], float] = {}
    names: dict[Key, str] = {}
    pairs = cKDTree(ox).sparse_distance_matrix(cKDTree(nx), cutoff,
                                               output_type="coo_matrix")
    for i, j, d in zip(pairs.row, pairs.col, pairs.data):
        k = (o_key[i], n_key[j])
        names[o_key[i]], names[n_key[j]] = o_name[i], n_name[j]
        best[k] = min(best.get(k, np.inf), float(d))
    taken: set[Key] = set()
    out = []
    for (acid, base), d in sorted(best.items(), key=lambda kv: (kv[1], kv[0])):
        if acid in taken or base in taken:
            continue
        taken.update((acid, base))
        out.append(Bridge(acid, names[acid], base, names[base], d))
    return out
