"""The image cost on the open pores' conductance (Round 7.15).

Each open deposit is read with :mod:`.born3d`'s self-energy W on its
electrostatic volume, beside the same readings without it:

* the neutral pore, and the neutral pore with the image cost (a symmetric
  salt stays neutral: both ions pay z²W alike, so it is one Boltzmann-weighted
  Laplace solve per species);
* Round 7.13's dipole on its own map (boundary at the ion centres), the
  dipole on the swept map (boundary an ion radius out, where W is solved),
  and the dipole with W in both its Poisson–Boltzmann equilibrium and its
  conduction; then the pair's two limits with W.

The image is screened by the bath's ionic strength. Screening by the mean
field's own ions is not taken: against a bath reference it would fold the
Debye–Hückel activity of a dense counter-ion cloud into W, which is a
different term (charge–space's, not the wall's).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from ..structure.channel import measure_channel
from ..structure.pore_volume import pore_volume
from .born3d import BornField, born_field
from .bridge_charge import lining_bridges
from .charged3d import Wall3D, wall_3d
from .permeation import _species_conductivity, potassium_species

__all__ = ["READINGS", "BornReading", "born_reading", "born_mutants",
           "SCAN_KEYS", "born_scan"]

#: label -> (surface, image, drop acid, drop base).
READINGS = {"dipole (7.13)": ("centres", False, False, False),
            "dipole, swept": ("swept", False, False, False),
            "dipole + image": ("swept", True, False, False),
            "base omitted + image": ("swept", True, False, True),
            "pair omitted + image": ("swept", True, True, True)}


@dataclass
class BornReading:
    """One deposit's conductances with and without the image cost."""

    name: str
    neutral: float                      # S, no image
    neutral_image: float                # S
    walls: dict[str, Wall3D]
    born: BornField
    axis: dict[str, tuple[float, float]]   # name -> (z, W on axis), kT per z²
    bridge: tuple[int, int] | None
    measured: dict[str, float] = field(default_factory=dict)
    spacing: float = 1.0
    anion_share: dict[str, float] = field(default_factory=dict)

    def g(self, label: str) -> float:
        return self.walls[label].charged["dielectric"]

    def ratio(self, label: str) -> float:
        """Against the neutral pore *without* the image: one ruler."""
        return self.g(label) / self.neutral

    def rows(self) -> list[str]:
        out = [f"{'neutral':22s} {self.neutral * 1e12:6.1f} pS",
               f"{'neutral + image':22s} {self.neutral_image * 1e12:6.1f} pS"
               f"  x{self.neutral_image / self.neutral:4.2f}"]
        for label, w in self.walls.items():
            tag = "" if w.converged else "  [n.c.]"
            out.append(f"{label:22s} {self.g(label) * 1e12:6.1f} pS  "
                       f"x{self.ratio(label):4.2f}  Cl- carries "
                       f"{self.anion_share.get(label, float('nan')):4.0%}{tag}")
        return out


def _species(summary):
    from .unitary import bath_for
    paralog = summary.numbering.paralog if summary.numbering else None
    return potassium_species(bath=bath_for(paralog)), paralog


def _axis(born: BornField, vol, summary) -> dict[str, tuple[float, float]]:
    zs, w = born.on_axis(vol)
    out = {}
    for name, c in summary.constrictions.items():
        k = int(np.argmin(np.abs(zs - c.z)))
        out[name] = (float(zs[k]), float(w[k]))
    lo, hi = summary.span
    inside = (zs >= lo) & (zs <= hi) & np.isfinite(w)
    k = int(np.flatnonzero(inside)[np.argmax(w[inside])])
    out["highest"] = (float(zs[k]), float(w[k]))
    return out


def born_reading(st, spacing: float = 1.0, workers: int | None = None,
                 readings=READINGS, cache: bool = True) -> BornReading:
    """:class:`BornReading` of ``st`` at ``spacing`` (the dielectric
    readings' 1 Å by default)."""
    from .unitary import published
    summary = measure_channel(st)
    species, paralog = _species(summary)
    smallest = min(species, key=lambda s: s.radius)
    vol = pore_volume(st, summary.frame, summary.span, smallest.radius,
                      spacing=spacing)
    born = born_field(st, summary.frame, vol, species, workers=workers,
                      cache=cache)
    pairs = lining_bridges(st)
    bridge = pairs[0] if pairs else None
    walls = {}
    plain = None
    for label, (surface, image, drop_a, drop_b) in readings.items():
        off = frozenset()
        if bridge:
            off = frozenset(r for r, d in zip(bridge, (drop_a, drop_b)) if d)
        walls[label] = wall_3d(st, summary, closures=("dielectric",),
                               spacing=spacing, neutralise=off,
                               surface=surface,
                               self_energy=born.energy if image else None)
        if not image and plain is None:
            plain = walls[label]
        if image:
            with_image = walls[label]
    return BornReading(st.name, plain.neutral, with_image.neutral, walls, born,
                       _axis(born, vol, summary), bridge,
                       published(paralog or "ITPR3"), spacing,
                       {k: _anion_share(w, species) for k, w in walls.items()})


def _anion_share(wall: Wall3D, species) -> float:
    """Fraction of the dielectric reading's conductance the anion carries."""
    t = _P.value("permeation.temperature")
    g = {s.name: _species_conductivity(s, t) * wall.per_species["dielectric"][s.name]
         for s in species}
    return sum(v for s in species if s.valence < 0 for v in [g[s.name]]) / sum(g.values())


def born_mutants(spacing: float = 1.0, workers: int | None = None,
                 progress=None) -> dict[str, tuple]:
    """Xu 2006's RyR1 mutants on 9HEO under the dielectric closure: Round
    7.13's map, the swept map, and the swept map with the image cost.
    ``label -> (wild type Wall3D, rows)``."""
    from ..io import loader
    from .charged3d import mutant_panel_3d
    from .ryr_mutants import open_deposit
    st = loader.load(open_deposit())
    summary = measure_channel(st)
    species, _ = _species(summary)
    vol = pore_volume(st, summary.frame, summary.span,
                      min(s.radius for s in species), spacing=spacing)
    born = born_field(st, summary.frame, vol, species, workers=workers)
    out = {}
    for label, kw in (("dipole (7.13)", {}), ("dipole, swept", {"surface": "swept"}),
                      ("dipole + image", {"surface": "swept",
                                          "self_energy": born.energy})):
        out[label] = mutant_panel_3d(st, closures=("dielectric",),
                                     spacing=spacing, progress=progress, **kw)
    return out


#: What :func:`born_scan` can move.
SCAN_KEYS = ("born.box_half_width", "born.reach", "dielectric.eps_protein")


def born_scan(st, key: str, values, spacing: float = 1.0,
              workers: int | None = None) -> list[tuple[float, BornReading]]:
    """:func:`born_reading` (neutral and dipole, each with the image) with
    one registered constant moved (restored afterwards, whatever happens)."""
    if key not in SCAN_KEYS:
        raise ValueError(f"key must be one of {SCAN_KEYS}, not {key!r}")
    chosen = {k: READINGS[k] for k in ("dipole, swept", "dipole + image")}
    before = _P.overrides().get(key)
    out = []
    try:
        for v in values:
            _P.set_value(key, v)
            out.append((float(v), born_reading(st, spacing, workers,
                                               readings=chosen)))
    finally:
        if before is None:
            _P.reset(key)
        else:
            _P.set_value(key, before)
    return out
