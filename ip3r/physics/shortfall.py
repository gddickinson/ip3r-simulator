"""Round 7.6: why the open pore conducts 2.4× less than measured.

The neutral 1-D model gives activated 8TKF 65 pS against 358 / 545 pS, and
RyR1's open 9HEO 136 pS against 801 pS. Three candidates, each tested here
by one measurement:

* **a substate** — 8TKF is one laboratory's open class. Schmitz et al.
  2022's active-state 7T3T (another laboratory, another preparation) is
  measured the same way (:func:`open_panel`).
* **the continuum at 3 Å** — the 1-D model's geometry: the inscribed circle
  of every slice, eroded over ±``pore.slab``. :func:`reading` measures the
  same electrolyte in the voxelised lumen, in 3-D
  (:mod:`.ohmic3d`), and as a slice-area integral of that volume.
* **the exit window** — the 1-D profile runs S0's 12 Å past the span and
  adds Hall's access at its ends. :func:`window_scan` moves the window; the
  3-D solve has no window at all (lateral exits count, the box sides are
  bath), and :func:`seal_scan` / :func:`grid_scan` show which of its own
  choices it depends on.

The in-pore diffusivity is the one constant no scan can settle, so every row
also says what fraction of bulk the measured value would need
(``needed_scale``): ohmic conductance is proportional to it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.annotations import is_ryr
from ..io import loader
from ..io.registry import load_registry
from ..parameters import PARAMETERS as _P
from ..structure.channel import measure_channel
from ..structure.pore import PROFILE_MARGIN, pore_profile
from .ohmic3d import Ohmic3D, conductance_3d, extrapolate
from .permeation import potassium_species, series_conductance
from .unitary import bath_for, published, unitary

__all__ = ["Reading", "reading", "open_entries", "open_panel", "window_scan", "seal_scan",
           "grid_scan", "OPEN_STATES"]

# A registered state whose first word is one of these (so not
# "preactivated", "inactivated"): the only
# ones the 1-D model finds conducting (every other state is shut, r_free ≤
# 1.03 Å; Rounds 4 and 6.1).
OPEN_STATES = ("activated", "active", "open")


@dataclass
class Reading:
    pdb_id: str
    state: str
    paralog: str
    measured: dict[str, float]      # pS
    gate_radius: float              # S0's r_min at the gate, Å
    one_d: float                    # pS, neutral 1-D (the Round 4 number)
    three_d: Ohmic3D                # at pore3d.spacing
    coarse: Ohmic3D | None          # at twice it (for the extrapolation)
    best: Ohmic3D | None = None     # the sweep's favourable corner
    best_coarse: Ohmic3D | None = None

    @property
    def three_d_pS(self) -> float:
        return self.three_d.conductance_pS

    @property
    def extrapolated_pS(self) -> float:
        if self.coarse is None:
            return float("nan")
        return extrapolate(self.coarse, self.three_d) * 1e12

    @property
    def best_pS(self) -> float:
        """The favourable corner (bulk diffusivity, smallest K+ radius),
        extrapolated to h 0."""
        if self.best is None or self.best_coarse is None:
            return float("nan")
        return extrapolate(self.best_coarse, self.best) * 1e12

    @property
    def gain(self) -> float:
        """3-D (extrapolated) over the 1-D reading."""
        return self.extrapolated_pS / self.one_d

    def short_by(self, value: float | None = None) -> float:
        """Measured over the extrapolated 3-D reading (the lowest measurement
        unless ``value`` is given)."""
        value = min(self.measured.values()) if value is None else value
        return value / self.extrapolated_pS

    def needed_scale(self, value: float | None = None) -> float:
        """In-pore diffusivity (fraction of bulk) at which the 3-D reading
        would meet the measurement; above 1 no admissible value will."""
        return _P.value("permeation.diffusion_scale") * self.short_by(value)

    def row(self) -> str:
        meas = " / ".join(f"{v:.0f}" for v in self.measured.values())
        return (f"{self.pdb_id:5s} {self.paralog:5s} {self.state[:22]:22s} gate "
                f"{self.gate_radius:4.2f} A   1-D {self.one_d:5.0f}   3-D "
                f"{self.three_d_pS:5.0f} (h {self.three_d.spacing:g}) -> "
                f"{self.extrapolated_pS:5.0f} (h 0)   slice-area "
                f"{self.three_d.slice_pS:5.0f}   measured {meas} pS")


def reading(pdb_id: str, state: str = "", corner: bool = True) -> Reading:
    """One deposit through all three geometries, in its family's bath."""
    st = loader.load(pdb_id)
    summary = measure_channel(st)
    paralog = summary.numbering.paralog if summary.numbering else ""
    h = _P.value("pore3d.spacing")
    u = unitary(st, summary)
    best = best_coarse = None
    if corner:
        sp = potassium_species(bath=bath_for(paralog),
                               diffusion_scale=_P.value("permeation.sweep_scale_high"),
                               ion_radius=_P.value("permeation.sweep_radius_low"))
        best = conductance_3d(st, summary, species=sp, spacing=h)
        best_coarse = conductance_3d(st, summary, species=sp, spacing=2 * h)
    return Reading(pdb_id=st.name, state=state, paralog=paralog,
                   measured=published(paralog),
                   gate_radius=u.gate_radius, one_d=u.neutral.conductance_pS,
                   three_d=conductance_3d(st, summary, spacing=h),
                   coarse=conductance_3d(st, summary, spacing=2 * h), best=best,
                   best_coarse=best_coarse)


def open_entries() -> list:
    """Every registered deposit whose state is open: ITPR3's 8TKF and the
    7T3T control, RyR1's 9HEO."""
    return [e for e in load_registry()
            if e.state.lower().split()[0] in OPEN_STATES
            and (e.human or is_ryr(e.paralog))]


def open_panel(corner: bool = True, progress=None) -> list[Reading]:
    """:func:`reading` of every :func:`open_entries` deposit held locally."""
    entries = open_entries()
    out = []
    for i, e in enumerate(entries):
        if progress:
            progress(i, len(entries), e.pdb_id)
        if loader.is_local(e.pdb_id) or loader.ALLOW_FETCH:
            out.append(reading(e.pdb_id, e.state, corner=corner))
    return out


def window_scan(pdb_id: str, factors=(0.5, 1.0, 2.0)) -> list[tuple[float, float, float]]:
    """``(margin Å, pS, access share)`` of the 1-D series reading with the
    profile window at ``factor`` × S0's margin past each end of the span."""
    st = loader.load(pdb_id)
    summary = measure_channel(st)
    paralog = summary.numbering.paralog if summary.numbering else None
    sp = potassium_species(bath=bath_for(paralog))
    out = []
    lo, hi = summary.span
    for f in factors:
        margin = PROFILE_MARGIN * f
        prof = pore_profile(st, summary.frame, lo - margin, hi + margin)
        s = series_conductance(prof.z, np.maximum(prof.r_free, 0.0), sp)
        share = s["access_ohm"] / (s["pore_ohm"] + s["access_ohm"])
        out.append((margin, s["conductance"] * 1e12, share))
    return out


def seal_scan(pdb_id: str, factors=(0.8, 1.0, 1.4)) -> list[tuple[float, float]]:
    """``(seal radius Å, pS)`` at the registered spacing doubled (a scan, so
    the coarse grid): flat while the seal lies inside the protein, a jump
    once it reaches the empty lipid space."""
    st = loader.load(pdb_id)
    summary = measure_channel(st)
    h = 2 * _P.value("pore3d.spacing")
    base = _P.value("pore3d.seal_radius")
    return [(base * f, conductance_3d(st, summary, spacing=h,
                                      seal=base * f).conductance_pS)
            for f in factors]


def grid_scan(pdb_id: str, factors=(2.0, 1.5, 1.0)) -> list[tuple[float, float]]:
    """``(spacing Å, pS)``: the convergence the extrapolation relies on."""
    st = loader.load(pdb_id)
    summary = measure_channel(st)
    h = _P.value("pore3d.spacing")
    return [(h * f, conductance_3d(st, summary, spacing=h * f).conductance_pS)
            for f in factors]
