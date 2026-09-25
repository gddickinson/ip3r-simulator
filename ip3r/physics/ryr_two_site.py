"""Stern's C channel with a two-site inactivation gate, fitted to the bell's slope.

Stern 1997's inactivation gate closes on one Ca2+, so its Hill slope is 1.
Murayama 2015's bell falls with nI 1.5, which they **fixed** in their fit
("fixed values for n A (1.2) and n I (1.5)"). The slope is theirs, not a
free measurement. The one-site fit (:func:`ryr_gating.fit_to_bell`)
matches the two half-peak points and leaves the slope wrong.

**The gate.** Two Ca2+ bind in sequence and the channel is inactivated
only when both are bound:

    N0 <-> N1 <-> N2 (inactivated),  on k_on x for each step,
                                     off k1_off and k2_off,

with ``x = c + mg_i mg`` (Mg2+ at the same site, as in the one-site gate),
``K1 = k1_off/k_on`` and ``K2 = k2_off/k_on``. At steady state the gate is
uninactivated with probability ``(1 + x/K1) / (1 + x/K1 + x^2/(K1 K2))``.
Its Hill slope at half inhibition, ``2 - c_h/(K1 + c_h)``, runs from 1
(K2 >> K1: the first Ca2+ as good as inactivates) to 2 (K2 << K1: the
second binds far tighter, which is cooperativity), so K2/K1 sets the slope. The activation gate is Stern's, unchanged and
independent, so the scheme has six states: activation gate shut or open
times 0, 1 or 2 inactivating Ca2+ bound.

**Fitted** (:func:`fit_two_site`). Three constants (Ka, K1, K2) against
three numbers of Murayama's bell on one ruler: the two half-peak points
and the log-log slope of the bell at its half-inhibition point. Stern's on
rates are kept for the activation gate and for both binding steps, and
only the off rates move, as in the one-site fit.

State ``s = a + 2 n`` (``a`` = activation gate open, ``n`` = inactivating
Ca2+ bound): C0, O0, C1, O1, C2, O2. O0 and O1 conduct; C2 and O2 are
inactivated. Exit kinds: activation gate, bind, unbind.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import ClassVar

import numpy as np
from scipy.optimize import fsolve

from ..parameters import PARAMETERS as _P
from .bell import Bell
from .ryr_gating import SternParams, bell_at, murayama_bell, open_probability

__all__ = ["TwoSiteParams", "STATES", "two_site", "fit_two_site",
           "bell_slope", "gate_hill_slope"]

STATES = ("C0", "O0", "C1", "O1", "C2", "O2")


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class TwoSiteParams(SternParams):
    dest: ClassVar[np.ndarray] = np.array(
        [[1, 2, 0], [0, 3, 1], [3, 4, 0], [2, 5, 1], [5, 4, 2], [4, 5, 3]])
    open_mask: ClassVar[np.ndarray] = np.array([False, True, False, True,
                                                False, False])
    inact_mask: ClassVar[np.ndarray] = np.array([False] * 4 + [True] * 2)
    trigger_map: ClassVar[np.ndarray] = np.array([1, 1, 3, 3, 4, 5])

    #: Second step's off rate, s^-1. The default makes K2 = K1 (Stern's Ki);
    #: the scheme used is the fitted one.
    k_inact2_off: float = _v("ryr.k_inact_off")

    @property
    def k_i2(self) -> float:
        return self.k_inact2_off / self.k_inact_on

    def uninactivated(self, c) -> np.ndarray:
        x = np.asarray(c, float) + self.mg_inact
        a = x / self.k_i
        return (1.0 + a) / (1.0 + a + a * x / self.k_i2)

    def exit_rates(self, state: np.ndarray, c: np.ndarray) -> np.ndarray:
        act_open, n = state % 2, state // 2
        r_act = np.where(act_open == 0, self.k_act_on * self.mg_factor * c * c,
                         self.k_act_off)
        r_bind = np.where(n < 2, self.k_inact_on * (c + self.mg_inact), 0.0)
        r_unbind = np.select([n == 1, n == 2],
                             [self.k_inact_off, self.k_inact2_off], 0.0)
        return np.concatenate([r_act, r_bind, r_unbind * np.ones_like(c)])

    def generator(self, c: float) -> np.ndarray:
        q = np.zeros((6, 6))
        states = np.arange(6)
        rates = self.exit_rates(states, np.full(6, float(c))).reshape(3, 6).T
        for s in states:
            for kind in range(3):
                d = self.dest[s, kind]
                if d != s:
                    q[s, d] += rates[s, kind]
        q[np.diag_indices(6)] = -q.sum(axis=1)
        return q

    def open_probability(self, c) -> np.ndarray:
        c = np.asarray(c, float)
        f_a = c ** 2 / (c ** 2 + self.k_a_eff ** 2)
        return f_a * self.uninactivated(c)


def two_site(k_a: float, k1: float, k2: float, sp: SternParams | None = None
             ) -> TwoSiteParams:
    """The two-site scheme with half-activation ``k_a`` and binding
    constants ``k1``, ``k2`` µM, keeping ``sp``'s on rates and Mg2+."""
    sp = sp or SternParams()
    base = {f: getattr(sp, f) for f in ("k_act_on", "k_inact_on", "mg",
                                        "k_mg_a", "mg_i")}
    return TwoSiteParams(**base, k_act_off=float(sp.k_act_on * k_a * k_a),
                         k_inact_off=float(sp.k_inact_on * k1),
                         k_inact2_off=float(sp.k_inact_on * k2))


def bell_slope(f, c: float, h: float = 1e-4) -> float:
    """d ln f / d ln c at ``c`` (central difference in ln c)."""
    up, down = f(c * np.exp(h)), f(c * np.exp(-h))
    return float((np.log(up) - np.log(down)) / (2 * h))


def gate_hill_slope(sp: TwoSiteParams) -> float:
    """Hill slope of the inactivation gate alone at its half point."""
    from scipy.optimize import brentq
    f_i = lambda x: 1.0 - sp.uninactivated(x)
    half = brentq(lambda x: f_i(x) - 0.5, 1e-6, 1e9)
    return bell_slope(lambda x: f_i(x) / (1.0 - f_i(x)), half)


def fit_two_site(target: Bell | None = None, sp: SternParams | None = None,
                 target_f=None) -> TwoSiteParams:
    """Ka, K1 and K2 such that the scheme's half-peak points and its slope
    at half inhibition are those of ``target_f`` (default Murayama's 25 C
    bell, ``target`` its measurement). Made at ``mg = 0``; raises unless all
    three are met to 1e-6."""
    from .ryr_gating import murayama_activity
    target = target or murayama_bell()
    target_f = target_f or murayama_activity
    want = bell_slope(target_f, target.c_half_inh)
    sp = replace(sp or SternParams(), mg=0.0)

    def resid(x):
        s = two_site(*np.exp(x), sp)
        b = bell_at(s)
        slope = bell_slope(lambda c: open_probability(c, s), b.c_half_inh)
        return [np.log(b.c_half_act / target.c_half_act),
                np.log(b.c_half_inh / target.c_half_inh),
                slope / want - 1.0]

    k_i = target.c_half_inh
    x, _, ok, msg = fsolve(resid, np.log([target.c_half_act, k_i, 0.5 * k_i]),
                           full_output=True)
    if ok != 1 or max(abs(r) for r in resid(x)) > 1e-6:
        raise RuntimeError(f"two-site fit did not converge: {msg}")
    return two_site(*np.exp(x), sp)
