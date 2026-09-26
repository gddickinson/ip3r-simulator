"""The charged lumen: the wall's potential and where the K+ drop falls with
it (Round 7.12).

Round 7.10 coloured the lumen by the neutral Laplace potential, which is
only geometry. Round 7.11 solved the wall charge on the same voxels
(:mod:`.charge3d`); this module reads that solution on the drawn lumen,
for one placement (``slice`` / ``local`` / ``pb``, or Round 7.13's
``dielectric``: every charged group in the box at its own centre, the
protein at ``dielectric.eps_protein``), two ways:

- **the wall potential** ``u`` (kT/e) at equilibrium, zero applied voltage:
  where the lining charge draws cations (u < 0) or repels them. The 1-D
  model's counterpart is its local Donnan potential on the inscribed
  circle's charge per length (:func:`.permeation.solve_pnp`'s closure).
- **the K+ drop**: the cation's electrochemical potential in linear
  response (:mod:`.charged3d`), a Laplace solve with conductivity
  ``σ e^{−u}``, normalised 0 at the luminal bath and 1 at the cytosolic.
  Set beside the neutral drop, it says where the charge moves the K+
  resistance. The 1-D counterpart is ``∫ dz / (A e^{−u})``. Each species
  has its own electrochemical drop in linear response; the electrical
  potential itself would need Poisson at first order, which no reading
  here solves, so the panel names this the K+ drop and nothing else.

The cation is the smallest ion of the family's bath, so the lumen field's
volume *is* :func:`.charged3d.wall_3d`'s electrostatic volume, and ``g``
here equals that function's K+ reading for the same closure (tested).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..structure.channel import ChannelSummary, measure_channel
from ._pnp_kernels import _donnan_potential
from .charge3d import CLOSURES_3D, WallField, group_positions, wall_field
from .dielectric3d import DIELECTRIC, dielectric_field
from .lumen_field import LumenField, plane_means
from .ohmic3d import geometric_conductance
from .permeation import potassium_species
from .pore_charge import PoreCharge, pore_charge

__all__ = ["ChargedLumen", "charged_lumen", "donnan_1d", "CLOSURE_LABELS",
           "LUMEN_CLOSURES"]

#: The placements the viewer offers: Round 7.11's three, then Round 7.13's.
LUMEN_CLOSURES = (*CLOSURES_3D, DIELECTRIC)

#: What the panel calls each placement.
CLOSURE_LABELS = {"slice": "1-D charge per length over the real area",
                  "local": "each group at its own centre (Donnan)",
                  "pb": "each group at its own centre, Poisson–Boltzmann",
                  DIELECTRIC: "every group in the box, protein ε (PB)"}


@dataclass
class ChargedLumen:
    """One placement of the wall charge, read on a :class:`LumenField`."""

    neutral: LumenField = field(repr=False)
    closure: str
    pair_bridges: bool
    wall: WallField = field(repr=False)
    charge: PoreCharge = field(repr=False)
    u: np.ndarray = field(repr=False)     # grid, kT/e; NaN outside the lumen
    mu: np.ndarray = field(repr=False)    # grid, K+ drop 0..1; NaN outside
    g: float                              # m, K+ geometric g with the energy
    u_3d: np.ndarray                      # mean u over each window plane
    u_1d: np.ndarray                      # 1-D Donnan u on the same planes
    mu_3d: np.ndarray                     # mean K+ drop over each plane (raw)
    area_1d: np.ndarray                   # Å², the 1-D model's K+ area
    converged: bool = True

    @property
    def z(self) -> np.ndarray:
        return self.neutral.z

    @property
    def ratio(self) -> float:
        """K+ conductance charged / neutral, 3-D."""
        return self.g / self.neutral.g if self.neutral.g > 0 else float("nan")

    @property
    def drop_3d(self) -> np.ndarray:
        """The K+ drop in the window, 0 at its luminal end, 1 at its cytosolic."""
        return _normalised(self.mu_3d)

    @property
    def drop_1d(self) -> np.ndarray:
        """``∫ dz / (A e^{−u})`` across the window, normalised (NaN if shut)."""
        if self.neutral.shut_1d or len(self.z) < 2:
            return np.full(len(self.z), np.nan)
        step = np.gradient(self.z) / (self.area_1d * np.exp(-self.u_1d))
        r = np.cumsum(step) - step / 2 - step[0] / 2
        return r / r[-1]

    def half_z(self, route: str = "3d") -> float:
        f = self.drop_3d if route == "3d" else self.drop_1d
        if not np.all(np.isfinite(f)):
            return float("nan")
        return float(np.interp(0.5, f, self.z))

    def steepest_z(self) -> float:
        """z where the 3-D K+ drop falls fastest, Å: where the charged pore's
        K+ resistance sits. Robust where :meth:`half_z` is not: a cation
        well carries almost none of the drop, and a plateau near one half
        moves the half-point by the well's length."""
        f = self.drop_3d
        if not np.all(np.isfinite(f)) or len(f) < 2:
            return float("nan")
        return float(self.z[int(np.argmax(np.gradient(f, self.z)))])

    def drop_across(self, z0: float, half_width: float) -> float:
        """Share of the window's K+ drop within ``z0 ± half_width`` (3-D)."""
        f = self.drop_3d
        if not np.all(np.isfinite(f)):
            return float("nan")
        return float(np.interp(z0 + half_width, self.z, f)
                     - np.interp(z0 - half_width, self.z, f))

    def well(self) -> tuple[float, float]:
        """The deepest cation well on the drawn voxels: (u in kT/e, its z)."""
        vals = np.where(np.isfinite(self.u), self.u, np.inf)
        k = np.unravel_index(int(np.argmin(vals)), vals.shape)
        return float(vals[k]), float(self.neutral.volume.zs[k[2]])

    def summary(self) -> str:
        pair = ", salt bridges paired" if self.pair_bridges else ""
        box = (f" ({self.wall.placed:+.1f} e in the box)"
               if self.closure == DIELECTRIC else "")
        u, z = self.well()
        text = (f"{self.neutral.name}, wall {self.charge.net_charge:+.1f} e{box}"
                f"{pair}, {self.closure}: K+ g ×{self.ratio:.2f} of neutral; "
                f"deepest cation well {u:+.1f} kT/e at z = {z:+.1f} Å; "
                f"the K+ drop is steepest at z = {self.steepest_z():+.1f} Å")
        return text


def _normalised(f: np.ndarray) -> np.ndarray:
    lo, hi = f[0], f[-1]
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        return np.full(len(f), np.nan)
    return (f - lo) / (hi - lo)


def donnan_1d(charge: PoreCharge, species, z: np.ndarray) -> np.ndarray:
    """The 1-D model's local Donnan potential (kT/e) on ``z``, Å."""
    valences = [s.valence for s in species]
    bath = np.array([[s.concentration * 1000.0] * len(charge.z) for s in species])
    u = _donnan_potential(valences, bath, np.asarray(charge.density, float), 1.0)
    return np.interp(z, charge.z, u)


def charged_lumen(st: Structure, neutral: LumenField, closure: str,
                  summary: ChannelSummary | None = None,
                  pair_bridges: bool = False, species=None) -> ChargedLumen:
    """The wall charge of ``st`` placed by ``closure`` on ``neutral``'s
    volume, its equilibrium potential and the K+ drop through it. Under
    ``dielectric`` every modelled group in the box carries charge (scope
    ``all``): unpaired, a salt bridge is its two charges (Round 7.13's
    dipole); paired, both partners are left out (its "pair omitted")."""
    from .unitary import bath_for, permeation_profile
    if closure not in LUMEN_CLOSURES:
        raise ValueError(f"closure must be one of {LUMEN_CLOSURES}, not {closure!r}")
    summary = summary or measure_channel(st)
    if species is None:
        paralog = summary.numbering.paralog if summary.numbering else None
        species = potassium_species(bath=bath_for(paralog))
    cation = next(s for s in species if s.valence > 0)
    if abs(cation.radius - neutral.probe) > 1e-9 or \
            min(s.radius for s in species) < cation.radius:
        raise ValueError("the lumen field must be cut for the bath's smallest "
                         "ion, the cation")
    vol = neutral.volume
    charge = pore_charge(st, summary.frame, permeation_profile(st, summary),
                         pair_bridges=pair_bridges)
    if closure == DIELECTRIC:
        wf = dielectric_field(st, summary.frame, vol, species, charge, scope="all")
    else:
        positions = group_positions(st, summary.frame, charge.groups)
        wf = wall_field(vol, closure, species, charge, positions=positions)
    lap = geometric_conductance(vol, energy=wf.energy(cation.valence))
    u = np.full(vol.mask.shape, np.nan)
    u[vol.mask] = wf.potential[vol.mask]
    mu = np.full(vol.mask.shape, np.nan)
    mu[vol.mask] = lap.potential
    lo, hi = neutral.window
    keep = (vol.zs >= lo) & (vol.zs <= hi)
    return ChargedLumen(
        neutral=neutral, closure=closure, pair_bridges=pair_bridges, wall=wf,
        charge=charge, u=u, mu=mu, g=lap.g,
        u_3d=plane_means(vol, u, keep), u_1d=donnan_1d(charge, species, neutral.z),
        mu_3d=plane_means(vol, mu, keep), area_1d=neutral.area_1d,
        converged=wf.converged and lap.converged)
