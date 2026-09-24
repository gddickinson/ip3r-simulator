"""RyR1 gating by cytosolic Ca2+: one measured bell, one kinetic scheme.

Two sources, kept apart because they answer different questions:

* **Murayama et al. 2015** (PLoS One 10:e0130606), the measured bell:
  Ca2+-dependent [3H]ryanodine binding of recombinant rabbit RyR1,
  ``A = Amax fA (1 - fI)`` with ``fA = c^nA/(c^nA + KA^nA)`` and
  ``fI = c^nI/(c^nI + KI^nI)`` (their Eqs. 1-3; wild type at 25 C, S1
  Table). Binding is an index of channel activity, not an open
  probability, so it is drawn normalised to its own peak. Its flanks are
  what can be compared.
* **Stern, Pizarro & Rios 1997** (J Gen Physiol 110:415), the kinetics:
  the skeletal "C channel" as two gates in series. Activation opens on the
  cooperative binding of two Ca2+ (on ``k_o c^2``, off ``k_o-``), and
  inactivation closes on one Ca2+ (on ``k_i c``, off ``k_i-``). Four states,
  C, O, CI and I, and only O conducts. The authors say the scheme "has not
  been objectively 'fitted' to data". It is used because it is the only
  RyR1-labelled Markov scheme whose constants could be read. The flank
  comparison with Murayama shows how far it is from a measured bell.

Because the two gates are independent, the stationary open probability is
the product ``fA (1 - fI)`` with ``Ka = sqrt(k_o-/k_o)`` and
``Ki = k_i-/k_i``. :func:`stationary` computes it from the generator's null
space instead, and the test holds the two to each other.

**Fitted to the measured bell** (:func:`fit_to_bell`). The steady state
fixes only the two ratios ``Ka`` and ``Ki``, so the fit moves the off
rates and keeps Stern's on rates. How fast inactivation is remains free:
``spark_termination`` scans it. A one-Ca2+ gate has a Hill slope of 1
against Murayama's fixed nI 1.5, so the fit matches the two half-peak
points, not the slopes.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import fsolve

from ..parameters import PARAMETERS as _P
from .bell import Bell, measure_bell

__all__ = ["STATES", "OPEN", "SternParams", "MurayamaParams", "generator",
           "stationary", "open_probability", "bell_at", "murayama_activity",
           "murayama_bell", "compare_bells", "fit_to_bell", "with_constants"]

#: C closed, O open, CI closed and inactivated, I open-gate but inactivated.
STATES = ("C", "O", "CI", "I")
OPEN = 1


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class SternParams:
    k_act_on: float = _v("ryr.k_act_on")          # µM^-2 s^-1 (x c^2)
    k_act_off: float = _v("ryr.k_act_off")        # s^-1
    k_inact_on: float = _v("ryr.k_inact_on")      # µM^-1 s^-1 (x c)
    k_inact_off: float = _v("ryr.k_inact_off")    # s^-1

    @property
    def k_a(self) -> float:
        """Half-activation, µM: sqrt(off/on) for the two-Ca2+ gate."""
        return float(np.sqrt(self.k_act_off / self.k_act_on))

    @property
    def k_i(self) -> float:
        return self.k_inact_off / self.k_inact_on


@dataclass
class MurayamaParams:
    a_max: float = _v("ryr.murayama_amax")
    k_a: float = _v("ryr.murayama_ka")            # µM
    n_a: float = _v("ryr.murayama_na")
    k_i: float = _v("ryr.murayama_ki")            # µM
    n_i: float = _v("ryr.murayama_ni")

    @classmethod
    def at_37(cls) -> "MurayamaParams":
        """The wild-type row at 37 C (the Hill coefficients are shared)."""
        return cls(a_max=_P.value("ryr.murayama_amax_37"),
                   k_a=_P.value("ryr.murayama_ka_37"),
                   k_i=_P.value("ryr.murayama_ki_37"))


def generator(c: float, sp: SternParams | None = None) -> np.ndarray:
    """Rate matrix Q (rows sum to zero) at cytosolic Ca2+ ``c`` µM."""
    sp = sp or SternParams()
    a, a_ = sp.k_act_on * c * c, sp.k_act_off
    i, i_ = sp.k_inact_on * c, sp.k_inact_off
    _C, _O, _CI, _I = range(4)
    q = np.zeros((4, 4))
    q[_C, _O], q[_O, _C] = a, a_       # activation gate
    q[_CI, _I], q[_I, _CI] = a, a_
    q[_C, _CI], q[_CI, _C] = i, i_     # inactivation gate
    q[_O, _I], q[_I, _O] = i, i_
    q[np.diag_indices(4)] = -q.sum(axis=1)
    return q


def stationary(c: float, sp: SternParams | None = None) -> np.ndarray:
    """Stationary occupancy of the four states (the null space of Q^T)."""
    v = null_space(generator(c, sp).T)[:, 0]
    return v / v.sum()


def open_probability(c, sp: SternParams | None = None) -> np.ndarray:
    """Closed form of :func:`stationary`'s O: fA (1 - fI)."""
    sp = sp or SternParams()
    c = np.asarray(c, float)
    f_a = c ** 2 / (c ** 2 + sp.k_a ** 2)
    return f_a * sp.k_i / (c + sp.k_i)


def bell_at(sp: SternParams | None = None) -> Bell:
    return measure_bell(lambda c: open_probability(c, sp))


def murayama_activity(c, mp: MurayamaParams | None = None) -> np.ndarray:
    """Murayama 2015 Eq. 1: [3H]ryanodine binding activity (B/Bmax)."""
    mp = mp or MurayamaParams()
    c = np.asarray(c, float)
    f_a = c ** mp.n_a / (c ** mp.n_a + mp.k_a ** mp.n_a)
    f_i = c ** mp.n_i / (c ** mp.n_i + mp.k_i ** mp.n_i)
    return mp.a_max * f_a * (1.0 - f_i)


def murayama_bell(mp: MurayamaParams | None = None) -> Bell:
    return measure_bell(lambda c: murayama_activity(c, mp))


def compare_bells() -> dict[str, Bell]:
    """Both bells on one ruler (peak and half-peak flanks)."""
    return {"Stern 1997 scheme": bell_at(), "Murayama 2015 fit": murayama_bell()}


def with_constants(k_a: float, k_i: float, sp: SternParams | None = None,
                   rate: float = 1.0) -> SternParams:
    """Stern's scheme with half-activation ``k_a`` and inactivation constant
    ``k_i`` µM. Only the off rates move, so the activation and inactivation
    on rates stay Stern's; ``rate`` scales both inactivation rates at once
    (Ki unchanged)."""
    sp = sp or SternParams()
    k_on = sp.k_inact_on * rate
    return SternParams(k_act_on=sp.k_act_on,
                       k_act_off=float(sp.k_act_on * k_a * k_a),
                       k_inact_on=float(k_on), k_inact_off=float(k_on * k_i))


def fit_to_bell(target: Bell | None = None, sp: SternParams | None = None
                ) -> SternParams:
    """Ka and Ki such that the scheme's half-peak flanks are ``target``'s
    (default: Murayama's 25 C bell). Solved in log space; raises if the
    solve does not reach both flanks to 1e-6."""
    target = target or murayama_bell()
    sp = sp or SternParams()

    def resid(x):
        b = bell_at(with_constants(*np.exp(x), sp))
        return [np.log(b.c_half_act / target.c_half_act),
                np.log(b.c_half_inh / target.c_half_inh)]

    x, _, ok, msg = fsolve(resid, np.log([sp.k_a, sp.k_i * target.c_half_inh
                                          / bell_at(sp).c_half_inh]),
                           full_output=True)
    if ok != 1 or max(abs(r) for r in resid(x)) > 1e-6:
        raise RuntimeError(f"bell fit did not converge: {msg}")
    return with_constants(*np.exp(x), sp)
