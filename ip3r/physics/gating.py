"""IP3 receptor gating: the De Young-Keizer model, as reduced by Li & Rinzel.

Each of the four subunits carries three independent sites (De Young & Keizer
1992): one for IP3 and two for Ca2+ — a fast **activating** site and a slow
**inhibitory** one. A subunit is *active* when IP3 and activating Ca2+ are
bound and the inhibitory site is empty. Li & Rinzel (1994) showed the fast
sites can be put at equilibrium, leaving a Hodgkin-Huxley-like model:

    m_inf = p / (p + d1)                IP3 site occupancy
    n_inf = c / (c + d5)                activating-Ca2+ occupancy
    dh/dt = a2 (Q2 (1 - h) - c h),      Q2 = d2 (p + d1) / (p + d3)

with ``h`` the fraction of inhibitory sites **free**. The channel conducts in
proportion to ``(m n h)^3`` — three active subunits of four.

Two consequences this module exists to make visible:

* **The bell.** At fixed IP3 the steady-state open probability rises with
  Ca2+ (activation) and falls again (inhibition): the bell-shaped curve
  Bezprozvanny et al. (1991) measured on cerebellar receptors in bilayers.
  ``bell_peak`` finds its maximum numerically.
* **IP3 relieves inhibition.** ``Q2`` grows with IP3, so the inhibitory
  site's apparent affinity for Ca2+ falls as IP3 rises and the right flank
  of the bell moves out. In this model the left flank moves too, less
  (over 0.1 -> 10 µM IP3: right flank 2.4x, half-activation 1.8x). Mak et
  al. (1998) measured IP3 tuning inhibition *alone*; the DYK scheme does not
  reproduce that, and a Hill-type model that does is on the roadmap.

Concentrations are µM, time is s.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize_scalar

from ..parameters import PARAMETERS as _P

__all__ = ["GatingParams", "m_inf", "n_inf", "q2", "h_inf", "tau_h",
           "subunit_activity", "open_probability", "bell_peak",
           "hill_fit_left_flank"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class GatingParams:
    a2: float = _v("gating.a2")
    d1: float = _v("gating.d1")
    d2: float = _v("gating.d2")
    d3: float = _v("gating.d3")
    d5: float = _v("gating.d5")
    subunits: float = _v("gating.subunits_required")


def m_inf(p, g: GatingParams | None = None):
    g = g or GatingParams()
    p = np.asarray(p, float)
    return p / (p + g.d1)


def n_inf(c, g: GatingParams | None = None):
    g = g or GatingParams()
    c = np.asarray(c, float)
    return c / (c + g.d5)


def q2(p, g: GatingParams | None = None):
    g = g or GatingParams()
    p = np.asarray(p, float)
    return g.d2 * (p + g.d1) / (p + g.d3)


def h_inf(c, p, g: GatingParams | None = None):
    """Steady-state fraction of inhibitory sites free."""
    g = g or GatingParams()
    qq = q2(p, g)
    return qq / (qq + np.asarray(c, float))


def tau_h(c, p, g: GatingParams | None = None):
    """Time constant of inhibition (s) — the slow variable of the model."""
    g = g or GatingParams()
    return 1.0 / (g.a2 * (q2(p, g) + np.asarray(c, float)))


def subunit_activity(c, p, h=None, g: GatingParams | None = None):
    """Probability that one subunit is active (IP3 on, act. on, inh. off)."""
    g = g or GatingParams()
    h = h_inf(c, p, g) if h is None else h
    return m_inf(p, g) * n_inf(c, g) * h


def open_probability(c, p, h=None, g: GatingParams | None = None):
    """Channel open probability, ``activity ** subunits`` (Li-Rinzel: cube).

    With ``h=None`` this is the **steady-state** curve at clamped Ca2+ and
    IP3 — the quantity a bilayer experiment reports.
    """
    g = g or GatingParams()
    return subunit_activity(c, p, h, g) ** g.subunits


def bell_peak(p: float, g: GatingParams | None = None) -> tuple[float, float]:
    """``(Ca2+ at the peak, peak open probability)`` at IP3 ``p``."""
    g = g or GatingParams()
    res = minimize_scalar(lambda lc: -float(open_probability(10 ** lc, p, g=g)),
                          bounds=(-4.0, 3.0), method="bounded",
                          options={"xatol": 1e-6})
    return float(10 ** res.x), float(-res.fun)


def hill_fit_left_flank(p: float, g: GatingParams | None = None) -> float:
    """Ca2+ at half-maximal activation on the rising flank of the bell (µM)."""
    g = g or GatingParams()
    c_peak, po_peak = bell_peak(p, g)
    grid = np.logspace(-4, np.log10(c_peak), 2000)
    po = open_probability(grid, p, g=g)
    return float(np.interp(po_peak / 2, po, grid))
