"""Fixed charge on the pore wall, from the deposit's own side-chain atoms.

The method is PIEZO1's (``piezo1/physics/pore_charge.py``): ionisable groups
that line the conduction pathway become a signed space-charge density on the
pore profile's slices, each spread along the axis by a Gaussian and divided
by the lumen cross-section it sits in. Spreading wall charge across the lumen
is legitimate when the screening length exceeds the radius; the permeation
result reports both so a reader can see where that holds (the narrow pore)
and where it does not (the wide vestibules).

**What differs from PIEZO1, and why.** PIEZO1's one open structure had no
side chains, so every charge there was placed at its C-alpha and admitted by
a fully-extended reach test. The IP3R deposits have side chains, so the
charge is placed where it is: the carboxylate oxygens' midpoint (Asp, Glu),
NZ (Lys) or CZ (Arg). A group lines the pore when that centre lies within
``pore_charge.lining_margin`` (about one water layer) of the profile's
atom-centre radius at its own height. A residue whose side chain was not
modelled (stubbed) is **not** guessed at: it is counted as *unplaced* when
its C-alpha is near enough that it could have lined the pore, and the count
travels with the result.

There is no curated route. The ip3r_genes annotation names filter and gate
residues, none of them ionisable, so a curated charge set would be empty.
Histidine and chain termini carry no charge here, as in PIEZO1: His is
mostly neutral at pH 7.4, and deposited chain ends are construct boundaries.

With ``pair_bridges``, a lining group that is one side of a salt bridge
(:mod:`ip3r.physics.salt_bridges`) is dropped, since the pair is net
neutral; the dropped bridges travel with the result.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.symmetry import Frame
from .salt_bridges import Bridge, salt_bridges

__all__ = ["ChargedGroup", "PoreCharge", "charged_groups", "map_charge",
           "pore_charge", "CHARGE", "CENTRE_ATOMS", "AVOGADRO"]

#: Avogadro constant, 1/mol (definitional since 2019).
AVOGADRO = 6.02214076e23

#: Formal side-chain charge at physiological pH.
CHARGE = {"ASP": -1.0, "GLU": -1.0, "LYS": 1.0, "ARG": 1.0}

#: The atoms whose mean is the charge centre.
CENTRE_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
                "LYS": ("NZ",), "ARG": ("CZ",)}


@dataclass(frozen=True)
class ChargedGroup:
    """One ionisable side chain that lines the pore."""

    res_seq: int
    res_name: str
    chain: str
    z: float                 # A along the axis, the charge centre's height
    radial: float            # A from the axis, the charge centre
    pore_radius: float       # A, the profile's r_min at this height
    charge: float            # e

    @property
    def margin(self) -> float:
        return self.radial - self.pore_radius

    def label(self) -> str:
        return f"{self.res_name}{self.res_seq}"


@dataclass
class PoreCharge:
    """A fixed-charge density profile and the groups it was built from."""

    z: np.ndarray            # A
    density: np.ndarray      # mol/m^3, signed (rho_fixed / F)
    groups: list[ChargedGroup] = field(default_factory=list)
    unplaced: list[str] = field(default_factory=list)
    meta: dict = field(default_factory=dict)
    bridged: list[Bridge] = field(default_factory=list)   # dropped as pairs

    @property
    def net_charge(self) -> float:
        return float(sum(g.charge for g in self.groups))

    @property
    def peak_density(self) -> float:
        return float(np.max(np.abs(self.density))) if len(self.density) else 0.0

    def residues(self) -> list[tuple[str, int, float]]:
        """``(label, copies, mean z)`` per residue number, ordered along z."""
        by: dict[str, list[float]] = {}
        for g in self.groups:
            by.setdefault(g.label(), []).append(g.z)
        rows = [(k, len(v), float(np.mean(v))) for k, v in by.items()]
        return sorted(rows, key=lambda r: r[2])

    def summary(self) -> str:
        paired = (f"; {len(self.bridged)} salt-bridged lining groups dropped"
                  if self.meta.get("pair_bridges") else "")
        return (f"{len(self.groups)} charged groups, net {self.net_charge:+.0f} e, "
                f"peak {self.peak_density / 1000.0:.2f} M-equivalent; "
                f"{len(self.unplaced)} stubbed ionisable residues unplaced"
                + paired)


def _centres(st: Structure) -> dict[tuple[str, int], tuple[str, np.ndarray]]:
    """Charge centre per (chain, residue) for every modelled ionisable group."""
    protein = ~st.hetero
    out = {}
    for name, atoms in CENTRE_ATOMS.items():
        m = protein & (st.res_name == name) & np.isin(st.atom_name, atoms)
        idx = np.flatnonzero(m)
        keys = list(zip(st.chain[idx], st.res_seq[idx].tolist()))
        grouped: dict[tuple[str, int], list[int]] = {}
        for k, i in zip(keys, idx):
            grouped.setdefault((str(k[0]), int(k[1])), []).append(i)
        for k, ii in grouped.items():
            if len(ii) == len(atoms):
                out[k] = (name, st.xyz[ii].astype(float).mean(axis=0))
    return out


def charged_groups(st: Structure, frame: Frame, profile
                   ) -> tuple[list[ChargedGroup], list[str]]:
    """Ionisable groups whose charge centre lines the lumen, and the stubbed
    residues that could have (``RESNUM/chain`` labels)."""
    margin = _P.value("pore_charge.lining_margin")
    stub_reach = _P.value("pore_charge.stub_reach")
    zp = np.asarray(profile.z, dtype=float)
    rp = np.asarray(profile.r_min, dtype=float)
    lo, hi = float(zp.min()), float(zp.max())

    centres = _centres(st)
    groups = []
    for (chain, seq), (name, xyz) in centres.items():
        f = frame.to_frame(xyz[None, :])[0]
        if not lo <= f[2] <= hi:
            continue
        g = ChargedGroup(seq, name, chain, float(f[2]), float(np.hypot(f[0], f[1])),
                         float(np.interp(f[2], zp, rp)), CHARGE[name])
        if g.margin <= margin:
            groups.append(g)

    unplaced = []
    ca = st.mask_ca() & np.isin(st.res_name, list(CHARGE))
    for i in np.flatnonzero(ca):
        key = (str(st.chain[i]), int(st.res_seq[i]))
        if key in centres:
            continue
        f = frame.to_frame(st.xyz[i][None, :].astype(float))[0]
        if not lo <= f[2] <= hi:
            continue
        if np.hypot(f[0], f[1]) - np.interp(f[2], zp, rp) <= stub_reach + margin:
            unplaced.append(f"{st.res_name[i]}{key[1]}/{key[0]}")
    return sorted(groups, key=lambda g: g.z), sorted(unplaced)


def map_charge(groups: list[ChargedGroup], z_A: np.ndarray, radius_A: np.ndarray,
               smoothing: float | None = None) -> np.ndarray:
    """Signed molar-equivalent density (mol/m^3) on the slices ``z_A``.

    Each charge is a Gaussian along the axis, normalised over the grid used
    so the total is conserved; the area is floored at the K+ radius so one
    carboxylate over a vanishing lumen cannot manufacture an unbounded density.
    """
    smoothing = (_P.value("pore_charge.smoothing") if smoothing is None
                 else smoothing)
    z = np.asarray(z_A, dtype=float)
    per_length = np.zeros_like(z)                       # e per A
    for g in groups:
        kernel = np.exp(-0.5 * ((z - g.z) / smoothing) ** 2)
        total = float(np.trapezoid(kernel, z))
        if total > 0.0:
            per_length += g.charge * kernel / total
    floor = _P.value("permeation.radius_potassium")
    area = np.pi * np.maximum(np.asarray(radius_A, dtype=float), floor) ** 2
    return (per_length / area) * 1e30 / AVOGADRO        # e/A^3 -> mol/m^3


def pore_charge(st: Structure, frame: Frame, profile,
                pair_bridges: bool = False,
                neutralise: frozenset[int] = frozenset()) -> PoreCharge:
    """Find the lining charges of a deposit and map them onto its profile;
    ``pair_bridges`` drops every lining group that is half of a salt bridge,
    ``neutralise`` every group at those residue numbers (a charge mutant)."""
    groups, unplaced = charged_groups(st, frame, profile)
    groups = [g for g in groups if g.res_seq not in neutralise]
    bridged: list[Bridge] = []
    if pair_bridges:
        lining = {(g.chain, g.res_seq) for g in groups}
        bridged = [b for b in salt_bridges(st)
                   if lining.intersection(b.members())]
        paired = {k for b in bridged for k in b.members()}
        groups = [g for g in groups if (g.chain, g.res_seq) not in paired]
    density = map_charge(groups, profile.z, np.maximum(profile.r_free, 0.0))
    meta = {"lining_margin_A": _P.value("pore_charge.lining_margin"),
            "smoothing_A": _P.value("pore_charge.smoothing"),
            "pair_bridges": pair_bridges, "structure": st.name,
            "neutralised": sorted(neutralise)}
    if pair_bridges:
        meta["salt_bridge_cutoff_A"] = _P.value("pore_charge.salt_bridge_cutoff")
    return PoreCharge(np.asarray(profile.z, dtype=float), density, groups,
                      unplaced, meta=meta, bridged=bridged)
