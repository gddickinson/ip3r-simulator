"""The IP3 binding site, measured on the atoms.

A contact is a protein residue with any heavy atom within the registered
cutoff (4.5 Å, S0's and S22's convention) of any heavy atom of a bound IP3
(PDB component ``I3P``). Contacts are split by whether the residue is on the
ligand's own subunit or a neighbour's — the IP3 site sits at the interface
of the β-trefoil and the armadillo domain of **one** subunit, and a
cross-subunit contact would be a finding worth seeing.

Which subunit a ligand "belongs to" is decided by proximity, not by its
chain letter: depositors do not all assign ligands to the chain they bind.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ..config import IP3_COMP_ID
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["LigandSite", "ligand_sites", "contacts", "residue_distances", "neighbourhood"]


@dataclass
class LigandSite:
    chain: str              # chain the ligand is filed under in the deposit
    residue: int
    subunit: str            # protein chain its atoms are nearest to
    atoms: np.ndarray       # atom indices into the structure

    def centre(self, st: Structure) -> np.ndarray:
        return st.xyz[self.atoms].mean(0)


def ligand_sites(st: Structure, comp: str = IP3_COMP_ID) -> list[LigandSite]:
    """Every bound copy of ``comp``, with the subunit it sits on."""
    lig = (st.res_name == comp) & (st.element != "H")
    if not lig.any():
        return []
    prot = (st.element != "H") & st.mask_protein() & ~st.hetero
    tree = cKDTree(st.xyz[prot])
    prot_chain = st.chain[prot]
    sites = []
    idx = np.flatnonzero(lig)
    keys = sorted({(str(st.chain[i]), int(st.res_seq[i])) for i in idx})
    for ch, rs in keys:
        atoms = idx[(st.chain[idx] == ch) & (st.res_seq[idx] == rs)]
        _, near = tree.query(st.xyz[atoms], k=8)
        chains, counts = np.unique(prot_chain[near.ravel()], return_counts=True)
        sites.append(LigandSite(ch, rs, str(chains[np.argmax(counts)]), atoms))
    return sites


def contacts(st: Structure, site: LigandSite,
             cutoff: float | None = None) -> dict[str, list[tuple[str, int]]]:
    """Residues within ``cutoff`` of the ligand: own subunit and others.

    Each residue is ``(resname, resnum)``; sorted by number.
    """
    cutoff = _P.value("ligand.contact_cutoff") if cutoff is None else cutoff
    prot = (st.element != "H") & st.mask_protein() & ~st.hetero
    pidx = np.flatnonzero(prot)
    tree = cKDTree(st.xyz[site.atoms])
    d, _ = tree.query(st.xyz[pidx], k=1, distance_upper_bound=cutoff + 1e-6)
    hit = pidx[np.isfinite(d) & (d <= cutoff)]
    own, other = set(), set()
    for i in hit:
        key = (str(st.res_name[i]), int(st.res_seq[i]))
        (own if st.chain[i] == site.subunit else other).add(key)
    return {"same_subunit": sorted(own, key=lambda k: k[1]),
            "other_subunit": sorted(other, key=lambda k: k[1])}


def neighbourhood(st: Structure, site: LigandSite,
                  radius: float | None = None) -> np.ndarray:
    """Atom indices of the ligand and every atom within ``radius`` of it
    (default ``ligand.shell_radius``, S22's pocket): what a camera centred
    on one site must keep in view."""
    radius = _P.value("ligand.shell_radius") if radius is None else radius
    d, _ = cKDTree(st.xyz[site.atoms]).query(st.xyz, k=1,
                                              distance_upper_bound=radius + 1e-6)
    return np.flatnonzero(np.isfinite(d))


def residue_distances(st: Structure, site: LigandSite, chain: str | None = None,
                      radius: float = 15.0, heavy_only: bool = True) -> dict[int, float]:
    """Minimum distance (Å) from each residue to the ligand.

    The S22 "shell" measurement: every residue with an atom within ``radius``
    of IP3, keyed by residue number on ``chain`` (default: the ligand's own
    subunit). ``heavy_only=False`` counts hydrogens on both sides, which is
    what S22 did (its reader keeps every atom); S0's contacts are heavy-atom.
    The difference is not cosmetic: in 8TKG and 8TKH, Arg503 is 4.78 and
    4.83 Å from IP3 by heavy atoms and inside 4.5 Å only through a hydrogen.
    """
    chain = site.subunit if chain is None else chain
    m = st.mask_protein() & ~st.hetero & (st.chain == chain)
    lig = site.atoms
    if heavy_only:
        m &= st.element != "H"
    else:
        lig = np.flatnonzero((st.res_name == st.res_name[site.atoms[0]])
                             & (st.chain == site.chain)
                             & (st.res_seq == site.residue))
    idx = np.flatnonzero(m)
    tree = cKDTree(st.xyz[lig])
    d, _ = tree.query(st.xyz[idx], k=1)
    out: dict[int, float] = {}
    for i, di in zip(idx, d):
        if di <= radius:
            r = int(st.res_seq[i])
            out[r] = min(out.get(r, np.inf), float(di))
    return dict(sorted(out.items()))
