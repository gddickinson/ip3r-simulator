"""Is a lining salt bridge charged, and what does its field do? (Round 7.13)

8TKF's D2478 carboxylates line the filter, each bridged (2.4–2.6 Å N–O) to
R2471 of the next subunit, whose guanidinium sits 3.6 Å further out at the
same height. Round 7.11's 3-D readings turn on how the pair is counted:
PB ×3.1 with D2478 alone, ×0.5 with the pair removed. Two questions, two
answers here:

* **Protonation** (:func:`titrate_pair`): the pair's pKas by the
  Tanford–Kirkwood network (sigmoidal ε, and uniform ε 10 and 4) with and
  without the base among the sites, and by PROPKA. If the acid titrates
  near the recording pH, the question is protonation; if both stay ionised,
  the pair is two charges, and the question is their field.
* **Field** (:func:`pair_readings`): the dielectric closure
  (:mod:`.dielectric3d`) with every group in the box at its own centre (the
  pair as a *dipole*), with the base omitted (the *full* wall's
  assumption), with both omitted (*paired*'s), and with the lining alone.
  Round 7.11's ``pb`` is shown beside them.

The homologous pairs: ITPR3 D2478–R2471′ (8TKF, 7T3T) and RyR1
D4899–R4892′ (9HEO), whose D4899Q is measured (Xu 2006: ×0.20).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from ..structure.channel import measure_channel
from .salt_bridges import salt_bridges

__all__ = ["PairTitration", "lining_bridges", "titrate_pair", "pair_readings",
           "READINGS", "SCAN_KEYS", "pair_scan"]

#: The field readings: label -> (closure, scope, drop acid, drop base).
READINGS = {"pb (7.11)": ("pb", "lining", False, False),
            "lining": ("dielectric", "lining", False, False),
            "dipole": ("dielectric", "all", False, False),
            "base omitted": ("dielectric", "all", False, True),
            "pair omitted": ("dielectric", "all", True, True)}

#: What :func:`pair_scan` can move.
SCAN_KEYS = ("dielectric.eps_protein", "dielectric.charge_width",
             "permeation.permittivity_pore", "pore3d.box_half_width",
             "pore3d.bath_margin")


def lining_bridges(st, wall=None) -> list[tuple[int, int]]:
    """``(acid, base)`` residue numbers of every salt bridge with a lining
    member (one entry per residue pair, all copies together)."""
    from .protonation import lining_wall
    wall = wall or lining_wall(st)
    lining = {(g.chain, g.res_seq) for g in wall.groups}
    return sorted({(b.acid[1], b.base[1]) for b in salt_bridges(st)
                   if lining.intersection(b.members())})


@dataclass
class PairTitration:
    """Apparent pKa ranges (over copies) of the pair under one reading."""

    reading: str
    acid: tuple[float, float]
    base: tuple[float, float] | None       # None: the base left out
    ph: float

    def acid_charge(self) -> tuple[float, float]:
        """Mean charge range of the acid at the recording pH."""
        return tuple(-1.0 / (1.0 + 10.0 ** (p - self.ph)) for p in self.acid)

    def row(self) -> str:
        base = ("base omitted" if self.base is None else
                f"R pKa {self.base[0]:5.1f}-{self.base[1]:5.1f}")
        q = self.acid_charge()
        return (f"{self.reading:22s} D pKa {self.acid[0]:5.2f}-{self.acid[1]:5.2f} "
                f"(charge {max(q):+.2f} to {min(q):+.2f})   {base}")


def _span(values) -> tuple[float, float]:
    return (float(np.min(values)), float(np.max(values)))


def titrate_pair(st, acid: int, base: int, wall=None, with_propka: bool = True
                 ) -> list[PairTitration]:
    """The pair under the network (sigmoidal, ε 10, ε 4), each with and
    without the base among the sites, and PROPKA (the base kept)."""
    from . import pka
    from .protonation import EPS_BOUNDS, family_conditions, lining_wall
    wall = wall or lining_wall(st)
    ph, ionic = family_conditions(wall.ryr)
    site_list, _ = pka.sites(st, wall.centres())
    out = []
    for eps in (None, *EPS_BOUNDS):
        name = "network" if eps is None else f"network eps {eps:g}"
        for drop in (False, True):
            use = [s for s in site_list if not (drop and s.res_seq == base)]
            ap = pka.titrate(use, ph, ionic, eps=eps).apparent_pka()
            a = [p for s, p in zip(use, ap) if s.res_seq == acid]
            b = [p for s, p in zip(use, ap) if s.res_seq == base]
            out.append(PairTitration(f"{name}{' -base' if drop else ''}",
                                     _span(a), None if drop else _span(b), ph))
    if with_propka:
        from .pka_propka import propka_pka
        pk = propka_pka(st, wall.centres())
        a = [v for (_, r), (_, v) in pk.items() if r == acid]
        b = [v for (_, r), (_, v) in pk.items() if r == base]
        out.append(PairTitration("PROPKA", _span(a), _span(b), ph))
    return out


def pair_readings(st, acid: int, base: int, summary=None,
                  spacing: float | None = None, readings=READINGS
                  ) -> dict[str, object]:
    """``label -> Wall3D`` for each field reading (K+ g, neutral beside)."""
    from .charged3d import wall_3d
    summary = summary or measure_channel(st)
    out = {}
    for label, (closure, scope, drop_acid, drop_base) in readings.items():
        off = frozenset(r for r, d in ((acid, drop_acid), (base, drop_base)) if d)
        out[label] = wall_3d(st, summary, closures=(closure,), spacing=spacing,
                             scope=scope, neutralise=off)
    return out


def pair_scan(st, acid: int, base: int, key: str, values,
              spacing: float | None = None,
              readings=("dipole", "base omitted", "pair omitted")
              ) -> list[tuple[float, dict[str, float]]]:
    """The dielectric readings' g ratios with one registered constant moved
    (restored afterwards, whatever happens)."""
    if key not in SCAN_KEYS:
        raise ValueError(f"key must be one of {SCAN_KEYS}, not {key!r}")
    summary = measure_channel(st)
    chosen = {k: READINGS[k] for k in readings}
    before = _P.overrides().get(key)
    out = []
    try:
        for v in values:
            _P.set_value(key, v)
            got = pair_readings(st, acid, base, summary, spacing, chosen)
            out.append((float(v), {k: w.ratio(READINGS[k][0])
                                   for k, w in got.items()}))
    finally:
        if before is None:
            _P.reset(key)
        else:
            _P.set_value(key, before)
    return out
