"""Unitary conductance of a deposit, and across the ITPR3 state panel.

The pore profile a deposit measures as (:mod:`ip3r.structure.pore`) becomes a
K+ conductance in symmetric KCl — the condition the measured values were
recorded in — three ways:

* ``series``: the closed-form resistor sum over the free radius, no solver;
* ``neutral``: drift-diffusion with no wall charge (must equal ``series``);
* ``charged``: drift-diffusion with the deposit's own lining charges
  (:mod:`ip3r.physics.pore_charge`) partitioning the ions;
* ``paired``: the same, with every salt-bridged lining group dropped
  (:mod:`ip3r.physics.salt_bridges`), since an ion pair is net neutral.

The profile is **protein only**, not S0's HETATM-inclusive one: a bound lipid
or detergent near the axis is a property of the preparation, not a wall an
ion meets in a membrane. ``r_free`` (less van der Waals radii) is the radius
an ion's centre is measured against.

**Recording condition by family.** An IP3R deposit is compared in
symmetric 140 mM KCl with the ITPR3 measurements; a RyR1 deposit in
symmetric 250 mM KCl with Xu et al. 2006's recombinant RyR1 (801 pS), the
bath each measurement was made in.

The answer depends on two constants nobody has measured for IP3R — the
in-pore diffusivity and the effective ion radius — so :func:`sensitivity`
reports the range over the registered sweep bounds, not one tuned number.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.annotations import is_ryr
from ..core.structure import Structure
from ..io import loader
from ..parameters import PARAMETERS as _P
from ..structure.channel import ChannelSummary, measure_channel
from ..structure.pore import PoreProfile, pore_profile
from ..structure.states import state_panel
from .permeation import (PermeationResult, potassium_species,
                         series_conductance, solve_pnp)
from .pore_charge import PoreCharge, pore_charge

__all__ = ["Unitary", "permeation_profile", "unitary", "unitary_panel",
           "sensitivity", "published"]


@dataclass
class Unitary:
    name: str
    state: str
    profile: PoreProfile
    charge: PoreCharge
    series: dict
    neutral: PermeationResult
    charged: PermeationResult
    paired_charge: PoreCharge
    paired: PermeationResult
    gate_radius: float = float("nan")      # S0's r_min at the gate, A
    sweep: dict | None = None              # sensitivity(), when asked for
    bath: float | None = None              # M, symmetric KCl; None = IP3R default

    @property
    def min_free_radius(self) -> float:
        return float(np.min(self.profile.r_free))

    @property
    def series_pS(self) -> float:
        return self.series["conductance"] * 1e12

    def row(self) -> str:
        if not self.neutral.is_conducting:
            return (f"{self.name:5s} {self.state:24s} r_free {self.min_free_radius:5.2f} A"
                    f"   closed ({self.neutral.blocked_by.split(':')[0]})")
        return (f"{self.name:5s} {self.state:24s} r_free {self.min_free_radius:5.2f} A"
                f"   series {self.series_pS:6.1f}  neutral "
                f"{_pS(self.neutral)}  charged {_pS(self.charged)}  paired "
                f"{_pS(self.paired)} pS"
                f"   (wall {self.charge.net_charge:+.0f} e, "
                f"{self.paired_charge.net_charge:+.0f} e unpaired)")


def _pS(r: PermeationResult) -> str:
    """A conductance, or ``n.c.`` where the solver did not converge (e.g. a
    wall charge that K+ alone cannot neutralise once Cl- is excluded)."""
    return f"{r.conductance_pS:6.1f}" if r.converged else "  n.c."


def published(paralog: str = "ITPR3") -> dict[str, float]:
    """The measured conductances the model of ``paralog`` is compared with,
    pS: ITPR3 in symmetric 140 mM KCl, RyR1 in symmetric 250 mM KCl."""
    if is_ryr(paralog):
        return {"Xu 2006 (bilayer)": _P.value("permeation.published_ryr1")}
    return {"Vais 2010 (DT40)": _P.value("permeation.published_itpr3_dt40"),
            "Mak 2000 (oocyte)": _P.value("permeation.published_itpr3_oocyte")}


def bath_for(paralog: str | None) -> float | None:
    """Symmetric KCl (M) the measurements of ``paralog`` were made in;
    ``None`` means the registered IP3R bath."""
    return _P.value("permeation.ryr1_bath_concentration") if is_ryr(paralog) else None


def permeation_profile(st: Structure, summary: ChannelSummary) -> PoreProfile:
    """Protein-only profile over the same window S0's profile spans."""
    lo, hi = summary.span
    return pore_profile(st, summary.frame, lo - 12.0, hi + 12.0,
                        include_hetero=False)


def unitary(st: Structure, summary: ChannelSummary | None = None,
            state: str = "", species=None,
            neutralise: frozenset[int] = frozenset()) -> Unitary:
    """Measure one deposit's conductance, in its family's recording bath;
    ``neutralise`` drops the charges at those residue numbers."""
    summary = summary or measure_channel(st)
    prof = permeation_profile(st, summary)
    radius = np.maximum(prof.r_free, 0.0)
    bath = bath_for(summary.numbering.paralog if summary.numbering else None)
    species = species or potassium_species(bath=bath)
    charge = pore_charge(st, summary.frame, prof, neutralise=neutralise)
    paired = pore_charge(st, summary.frame, prof, pair_bridges=True,
                         neutralise=neutralise)
    gate = summary.constrictions.get("gate")
    return Unitary(
        name=st.name, state=state, profile=prof, charge=charge,
        series=series_conductance(prof.z, radius, species),
        neutral=solve_pnp(prof.z, radius, species=species),
        charged=solve_pnp(prof.z, radius, species=species,
                          fixed_charge=charge.density),
        paired_charge=paired,
        paired=solve_pnp(prof.z, radius, species=species,
                         fixed_charge=paired.density),
        gate_radius=float("nan") if gate is None else gate.radius, bath=bath)


def unitary_panel(paralog: str = "ITPR3", progress=None,
                  sweep: bool = False) -> list[Unitary]:
    """Every human deposit of ``paralog``, ordered by gate radius (S11's
    panel); ``sweep`` also fills each row's :func:`sensitivity`."""
    rows = state_panel(paralog, progress=progress)
    out = [unitary(loader.load(r.pdb_id), r.summary, r.state) for r in rows]
    if sweep:
        for u in out:
            u.sweep = sensitivity(u)
    return out


def sensitivity(u: Unitary) -> dict[str, tuple[float, float]]:
    """``(min, max)`` pS over the sweep corners of diffusivity x ion radius.

    The corners bound the range because conductance is monotone in both:
    rising with diffusivity, falling with radius. A corner where the solver
    did not converge is left out; if none converged the range is NaN.
    """
    radius = np.maximum(u.profile.r_free, 0.0)
    out: dict[str, list[float]] = {"neutral": [], "charged": [], "paired": []}
    for scale in (_P.value("permeation.sweep_scale_low"),
                  _P.value("permeation.sweep_scale_high")):
        for ion in (_P.value("permeation.sweep_radius_low"),
                    _P.value("permeation.sweep_radius_high")):
            sp = potassium_species(bath=u.bath, diffusion_scale=scale,
                                   ion_radius=ion)
            for key, fixed in (("neutral", None), ("charged", u.charge.density),
                               ("paired", u.paired_charge.density)):
                r = solve_pnp(u.profile.z, radius, species=sp, fixed_charge=fixed)
                if r.converged:
                    out[key].append(r.conductance_pS)
    nan = float("nan")
    return {k: (min(v), max(v)) if v else (nan, nan) for k, v in out.items()}
