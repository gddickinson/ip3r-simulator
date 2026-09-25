"""Protonation readings of the pore wall, through each family's own
selectivity experiment (Round 7.4).

Every charged reading so far gave each lining Asp/Glu -1 and each Lys/Arg
+1. This module replaces that assumption with measured-model estimates and
asks what the selectivity does:

* **formal**: the assumption, as before.
* **network**: the Tanford-Kirkwood titration (:mod:`.pka`) at the
  recording pH, sigmoidal permittivity.
* **network eps N**: the same with a uniform permittivity N, the bound on
  how strongly the pore might couple its charges (protein-like at 4).
* **PROPKA**: PROPKA 3's pKa per group (:mod:`.pka_propka`), then
  Henderson-Hasselbalch at the recording pH.

and, independently of any pKa, :func:`corners`: each lining residue number
(all four copies) either formal or neutral, every combination. The largest
P_Ca:P_K over the corners bounds what *any* protonation state could give,
because a mean charge between 0 and formal lies inside the box the corners
span. (The bound is a maximum over corners, not a proof of monotonicity;
the corner table is reported so the reader can see the shape.)

The experiment is the family's: Vais 2010 for IP3R (pH 7.3, 140 mM KCl, both
ratios), Xu 2006 for RyR1 (pH 7.4, 250 mM KCl, P_Ca:P_K by their Eq. 1).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from ..core.annotations import is_ryr
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from . import pka as _pka
from .pore_charge import CENTRE_ATOMS, pore_charge

__all__ = ["Wall", "WallReading", "family_conditions", "lining_wall",
           "network_charges", "propka_charges", "readings", "corners",
           "interior", "measure", "EPS_BOUNDS"]

#: The uniform permittivities of the network's bound.
EPS_BOUNDS = (10.0, 4.0)


@dataclass
class Wall:
    """A deposit's pore, ready to be charged in different ways."""

    structure: Structure
    summary: object
    profile: object
    radius: np.ndarray
    ryr: bool
    groups: list = field(default_factory=list)       # the formal lining groups

    def centres(self) -> np.ndarray:
        st = self.structure
        return np.array([st.xyz[(st.chain == g.chain) & (st.res_seq == g.res_seq)
                                & np.isin(st.atom_name, CENTRE_ATOMS[g.res_name])
                                & ~st.hetero].astype(float).mean(axis=0)
                         for g in self.groups])

    def density(self, charges: dict | None = None,
                neutralise: frozenset[int] = frozenset()):
        ch = pore_charge(self.structure, self.summary.frame, self.profile,
                         neutralise=neutralise, charges=charges)
        return ch.density, ch.net_charge


@dataclass
class WallReading:
    label: str
    net_charge: float
    by_residue: dict[str, float]         # mean charge per lining residue
    pca_pk: float
    pcl_pk: float                        # NaN for RyR1 (Xu's ruler has none)
    conductance: float                   # S, symmetric bath
    converged: bool

    def row(self) -> str:
        pcl = "" if np.isnan(self.pcl_pk) else f"   P_Cl:P_K {self.pcl_pk:5.2f}"
        flag = "" if self.converged else "  (n.c.)"
        return (f"{self.label:15s} wall {self.net_charge:+6.2f} e   "
                f"P_Ca:P_K {self.pca_pk:6.2f}{pcl}   "
                f"g {self.conductance * 1e12:6.1f} pS{flag}")


def family_conditions(ryr: bool) -> tuple[float, float]:
    """(pH, ionic strength M) of the family's selectivity recordings."""
    if ryr:
        return (_P.value("pka.ph_xu"),
                _P.value("permeation.ryr1_bath_concentration"))
    return _P.value("pka.ph_vais"), _P.value("permeation.bath_concentration")


def lining_wall(st: Structure, summary=None) -> Wall:
    from ..structure.channel import measure_channel
    from .unitary import permeation_profile
    summary = summary or measure_channel(st)
    prof = permeation_profile(st, summary)
    ryr = bool(summary.numbering and is_ryr(summary.numbering.paralog))
    groups = pore_charge(st, summary.frame, prof).groups
    return Wall(st, summary, prof, np.maximum(prof.r_free, 0.0), ryr, groups)


def network_charges(wall: Wall, eps: float | None = None
                    ) -> tuple[dict, _pka.Titration]:
    ph, ionic = family_conditions(wall.ryr)
    site_list, _ = _pka.sites(wall.structure, wall.centres())
    t = _pka.titrate(site_list, ph, ionic, eps=eps)
    return t.charges(), t


def propka_charges(wall: Wall) -> tuple[dict, dict]:
    from .pka_propka import charge_at, propka_pka
    ph, _ = family_conditions(wall.ryr)
    pk = propka_pka(wall.structure, wall.centres())
    lining = {(g.chain, g.res_seq) for g in wall.groups}
    missing = lining - set(pk)
    if missing:
        raise RuntimeError(f"PROPKA reported no pKa for {sorted(missing)}")
    return {k: charge_at(name, v, ph) for k, (name, v) in pk.items()}, pk


def measure(wall: Wall, label: str, charges: dict | None = None,
             neutralise: frozenset[int] = frozenset()) -> WallReading:
    from .permeation import potassium_species, solve_pnp
    from .selectivity import ryr1_calcium_ratio, selectivity
    fixed, net = wall.density(charges, neutralise)
    z = wall.profile.z
    if wall.ryr:
        _, pca, ok = ryr1_calcium_ratio(z, wall.radius, fixed)
        pcl = float("nan")
        bath = _P.value("permeation.ryr1_bath_concentration")
    else:
        sel = selectivity(z, wall.radius, fixed, label)
        pca, pcl, ok, bath = sel.pca_pk, sel.pcl_pk, sel.converged, None
    g = solve_pnp(z, wall.radius, species=potassium_species(bath=bath),
                  fixed_charge=fixed)
    by: dict[str, list[float]] = {}
    for grp in wall.groups:
        if grp.res_seq in neutralise:
            continue
        q = grp.charge if charges is None else charges.get(
            (grp.chain, grp.res_seq), grp.charge)
        by.setdefault(grp.label(), []).append(q)
    return WallReading(label, net, {k: float(np.mean(v)) for k, v in by.items()},
                       float(pca), float(pcl), g.conductance,
                       bool(ok and g.converged))


def readings(wall: Wall, with_propka: bool = True
             ) -> tuple[list[WallReading], dict]:
    """The wall under every protonation reading; and the titrations
    (``network``, ``network eps N``, ``propka`` pKas) they came from."""
    rows = [measure(wall, "formal")]
    detail: dict = {}
    for eps in (None, *EPS_BOUNDS):
        label = "network" if eps is None else f"network eps {eps:g}"
        charges, t = network_charges(wall, eps)
        detail[label] = t
        rows.append(measure(wall, label, charges))
    if with_propka:
        charges, pk = propka_charges(wall)
        detail["propka"] = pk
        rows.append(measure(wall, "PROPKA", charges))
    return rows, detail


def corners(wall: Wall, residues: list[int] | None = None,
            progress=None) -> list[tuple[frozenset[int], WallReading]]:
    """Every combination of lining residues (all copies) neutral, the rest
    formal: ``(neutral residues, reading)``, largest P_Ca:P_K first. On
    8TKF (6 rings, 64 corners) this takes about 2.5 min."""
    residues = sorted({g.res_seq for g in wall.groups}
                      if residues is None else residues)
    out = []
    combos = [frozenset(c) for n in range(len(residues) + 1)
              for c in itertools.combinations(residues, n)]
    for i, off in enumerate(combos):
        label = "-" + ",".join(map(str, sorted(off))) if off else "formal"
        out.append((off, measure(wall, label, neutralise=off)))
        if progress:
            progress(i + 1, len(combos))
    return sorted(out, key=lambda t: -t[1].pca_pk)


def interior(wall: Wall, n: int | None = None, seed: int | None = None
             ) -> list[tuple[dict[int, float], WallReading]]:
    """Random fractional states inside the corners' box (each residue's
    copies at one fraction of formal), the first at one half: the check
    that the corners' maximum is not beaten between them."""
    n = int(_P.value("protonation.interior_samples") if n is None else n)
    rng = np.random.default_rng(int(_P.value("pka.mc_seed") if seed is None
                                    else seed))
    residues = sorted({g.res_seq for g in wall.groups})
    out = []
    for k in range(n):
        f = {r: 0.5 if k == 0 else float(rng.random()) for r in residues}
        ch = {(g.chain, g.res_seq): g.charge * f[g.res_seq] for g in wall.groups}
        out.append((f, measure(wall, "interior", ch)))
    return out
