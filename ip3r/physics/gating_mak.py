"""IP3 receptor gating as measured by Mak, McBride & Foskett (1998).

Single IP3R-1 channels in *Xenopus* oocyte nuclear patches, Ca2+ and IP3
clamped on the cytoplasmic side. The steady-state open probability fits a
biphasic Hill equation with **one** denominator (their Eq. 1):

    P_open = P_max / [1 + (K_act / c)^H_act + (c / K_inh)^H_inh]

and of its five constants only the inhibitory one depends on IP3 (Eq. 2):

    K_inh(p) = K_inf / [1 + (K_IP3 / p)^H_IP3]

"Kinh being the only IP3-concentration-sensitive parameter ... Pmax, Kact,
Hact, and Hinh all remain unchanged" — IP3 activates the channel by relieving
Ca2+ inhibition, nothing else. The De Young-Keizer scheme (``gating``) moves
both flanks; this one moves the inhibitory flank alone, except where the
bell narrows so far at very low IP3 (K_inh approaching K_act) that the two
terms of the denominator meet.

It is a description of steady-state data, not a kinetic scheme: there is no
inhibition time constant and it cannot drive the cell or puff models.
Concentrations are µM.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from .bell import Bell, flank_shifts, measure_bell

__all__ = ["MakParams", "k_inh", "open_probability", "bell_at",
           "compare_flanks"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class MakParams:
    p_max: float = _v("mak.p_max")
    k_act: float = _v("mak.k_act")
    h_act: float = _v("mak.h_act")
    h_inh: float = _v("mak.h_inh")
    k_inf: float = _v("mak.k_inf")
    k_ip3: float = _v("mak.k_ip3")
    h_ip3: float = _v("mak.h_ip3")


def k_inh(p, g: MakParams | None = None):
    """Inhibitory half-constant (µM) at IP3 ``p`` (Eq. 2); 0 at no IP3."""
    g = g or MakParams()
    p = np.asarray(p, float)
    with np.errstate(divide="ignore"):
        return g.k_inf / (1.0 + (g.k_ip3 / p) ** g.h_ip3)


def open_probability(c, p, g: MakParams | None = None):
    """Steady-state P_open at clamped Ca2+ ``c`` and IP3 ``p`` (Eq. 1)."""
    g = g or MakParams()
    c = np.asarray(c, float)
    ki = k_inh(p, g)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = 1.0 + (g.k_act / c) ** g.h_act + (c / ki) ** g.h_inh
    return np.where(np.isfinite(denom), g.p_max / denom, 0.0)


def bell_at(p: float, g: MakParams | None = None) -> Bell:
    g = g or MakParams()
    return measure_bell(lambda c: open_probability(c, p, g))


def compare_flanks(p_lo: float | None = None, p_hi: float | None = None
                   ) -> dict[str, tuple[float, float]]:
    """Fold shift of (activating, inhibitory) half-points, IP3 ``p_lo`` ->
    ``p_hi``, for the three IP3R models on one ruler:
    ``{"DYK": ..., "Mak": ..., "Park/drive": ...}``."""
    from . import gating, park_drive
    p_lo = _P.value("mak.compare_ip3_low") if p_lo is None else p_lo
    p_hi = _P.value("mak.compare_ip3_high") if p_hi is None else p_hi
    return {"DYK": flank_shifts(gating.bell_at, p_lo, p_hi),
            "Mak": flank_shifts(bell_at, p_lo, p_hi),
            "Park/drive": flank_shifts(park_drive.bell_at, p_lo, p_hi)}
