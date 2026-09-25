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

**Mg2+** (``mg`` µM free, default 0: Murayama's bell and Stern's scheme are
both Mg2+-free). Two sourced mechanisms, no fitted constant:

* the activation site: Mg2+ competes with Ca2+ in rapid equilibrium with
  affinity ``k_mg_a`` (Laver et al. 2004, Table I). The two-Ca2+ on rate is
  divided by ``(1 + mg/k_mg_a)^2``, so the half-activation becomes
  ``Ka (1 + mg/k_mg_a)``: the competitive form of their Eq. 6.
* the inactivation site: Ca2+ and Mg2+ inhibit "at the same low-affinity
  inhibitory site" with identical effect (Laver et al. 1997), so the gate's
  ligand is ``c + mg_i mg`` at the same Ki, with ``mg_i`` (the Mg2+/Ca2+
  affinity ratio, ``ryr.mg_i_relative``) 1. Setting it to 0 removes the
  inactivation-site action, which is how the two are dissected. Laver's Hill coefficient ~2 is not
  kept; the gate binds one ion, as it does for Ca2+.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import ClassVar

import numpy as np

from scipy.optimize import fsolve

from ..parameters import PARAMETERS as _P
from .bell import Bell, measure_bell

__all__ = ["STATES", "OPEN", "SternParams", "MurayamaParams", "generator",
           "stationary", "open_probability", "bell_at", "murayama_activity",
           "murayama_bell", "compare_bells", "fit_to_bell", "with_constants",
           "with_mg"]

#: C closed, O open, CI closed and inactivated, I open-gate but inactivated.
STATES = ("C", "O", "CI", "I")
OPEN = 1


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class SternParams:
    #: What a simulator needs to know of the scheme (subclasses override):
    #: each state's exits (one column per kind, the order of
    #: :meth:`exit_rates`), which states conduct, which are inactivated, and
    #: where the V-channel trigger sends each state.
    dest: ClassVar[np.ndarray] = np.array([[1, 2], [0, 3], [3, 0], [2, 1]])
    open_mask: ClassVar[np.ndarray] = np.array([False, True, False, False])
    inact_mask: ClassVar[np.ndarray] = np.array([False, False, True, True])
    trigger_map: ClassVar[np.ndarray] = np.array([1, 1, 2, 3])

    k_act_on: float = _v("ryr.k_act_on")          # µM^-2 s^-1 (x c^2)
    k_act_off: float = _v("ryr.k_act_off")        # s^-1
    k_inact_on: float = _v("ryr.k_inact_on")      # µM^-1 s^-1 (x c)
    k_inact_off: float = _v("ryr.k_inact_off")    # s^-1
    mg: float = 0.0                               # µM free Mg2+
    k_mg_a: float = _v("ryr.k_mg_a")              # µM, activation site
    mg_i: float = _v("ryr.mg_i_relative")         # Mg2+/Ca2+ at inactivation

    @property
    def mg_inact(self) -> float:
        """Mg2+ as seen by the inactivation gate, µM of Ca2+ equivalent."""
        return self.mg_i * self.mg

    @property
    def mg_factor(self) -> float:
        """Share of the two-Ca2+ on rate left by Mg2+ at the activation site."""
        return 1.0 / (1.0 + self.mg / self.k_mg_a) ** 2

    @property
    def k_a_eff(self) -> float:
        """Half-activation under ``mg``, µM: Ka (1 + mg/k_mg_a)."""
        return self.k_a * (1.0 + self.mg / self.k_mg_a)

    @property
    def k_a(self) -> float:
        """Half-activation, µM: sqrt(off/on) for the two-Ca2+ gate."""
        return float(np.sqrt(self.k_act_off / self.k_act_on))

    @property
    def k_i(self) -> float:
        """Inactivation constant, µM: off/on. Infinite when the gate is
        removed (``k_inact_on`` 0), which is how a scheme with no Ca2+
        inactivation is written without an off rate large enough to swamp
        the Gillespie step."""
        on = self.k_inact_on
        return self.k_inact_off / on if on > 0 else float("inf")

    def exit_rates(self, state: np.ndarray, c: np.ndarray) -> np.ndarray:
        """Exit rates of channels in ``state`` at their own Ca2+ ``c``:
        activation gates, then inactivation gates (``dest``'s columns)."""
        act_gate_shut = (state == 0) | (state == 2)
        r_act = np.where(act_gate_shut, self.k_act_on * self.mg_factor * c * c,
                         self.k_act_off)
        r_inact = np.where(state <= 1, self.k_inact_on * (c + self.mg_inact),
                           self.k_inact_off)
        return np.concatenate([r_act, r_inact])

    def generator(self, c: float) -> np.ndarray:
        """Rate matrix Q (rows sum to zero) at cytosolic Ca2+ ``c`` µM."""
        a, a_ = self.k_act_on * self.mg_factor * c * c, self.k_act_off
        i, i_ = self.k_inact_on * (c + self.mg_inact), self.k_inact_off
        _C, _O, _CI, _I = range(4)
        q = np.zeros((4, 4))
        q[_C, _O], q[_O, _C] = a, a_       # activation gate
        q[_CI, _I], q[_I, _CI] = a, a_
        q[_C, _CI], q[_CI, _C] = i, i_     # inactivation gate
        q[_O, _I], q[_I, _O] = i, i_
        q[np.diag_indices(4)] = -q.sum(axis=1)
        return q

    def open_probability(self, c) -> np.ndarray:
        """Closed form of :func:`stationary`'s O: fA (1 - fI). Written as
        ``f_a / (1 + x/Ki)`` rather than ``f_a Ki/(x + Ki)`` so that a
        removed gate (``Ki`` infinite) gives ``f_a`` instead of 0/0."""
        c = np.asarray(c, float)
        f_a = c ** 2 / (c ** 2 + self.k_a_eff ** 2)
        return f_a / (1.0 + (c + self.mg_inact) / self.k_i)


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
    return (sp or SternParams()).generator(c)


def stationary(c: float, sp: SternParams | None = None) -> np.ndarray:
    """Stationary occupancy of the scheme's states: ``Q^T pi = 0`` with
    ``sum(pi) = 1``.

    Solved as that constrained system rather than as the null space of
    ``Q^T``, because an SVD null space loses accuracy when the rates span
    many decades (the activation on rate goes as ``c^2``) and when some
    states are unreachable, as they are in a scheme with a gate removed.
    The normalisation is then exact rather than applied afterwards. Any
    residual negative occupancy is numerical and is clipped.
    """
    q = generator(c, sp).T
    n = q.shape[0]
    a = np.vstack([q, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    pi = np.linalg.lstsq(a, b, rcond=None)[0]
    total = np.clip(pi, 0.0, None).sum()
    if not np.isfinite(total) or total <= 0:      # a degenerate generator
        raise FloatingPointError(f"no stationary state at Ca2+ {c:g} µM")
    return np.clip(pi, 0.0, None) / total


def open_probability(c, sp: SternParams | None = None) -> np.ndarray:
    """Closed form of :func:`stationary`'s open occupancy."""
    return (sp or SternParams()).open_probability(c)


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
    return replace(sp, k_act_off=float(sp.k_act_on * k_a * k_a),
                   k_inact_on=float(k_on), k_inact_off=float(k_on * k_i))


def with_mg(sp: SternParams | None = None, mg: float | None = None,
            k_mg_a: float | None = None, mg_i: float | None = None
            ) -> SternParams:
    """``sp`` under free Mg2+ ``mg`` µM (default ``ryr.mg_free``), with the
    activation-site affinity ``k_mg_a`` and the inactivation-site ratio
    ``mg_i`` (defaults: ``sp``'s own; ``k_mg_a`` = inf or ``mg_i`` = 0
    removes that site's action)."""
    sp = sp or SternParams()
    return replace(sp, mg=_P.value("ryr.mg_free") if mg is None else float(mg),
                   k_mg_a=sp.k_mg_a if k_mg_a is None else float(k_mg_a),
                   mg_i=sp.mg_i if mg_i is None else float(mg_i))


def fit_to_bell(target: Bell | None = None, sp: SternParams | None = None
                ) -> SternParams:
    """Ka and Ki such that the scheme's half-peak flanks are ``target``'s
    (default: Murayama's 25 C bell). Solved in log space; raises if the
    solve does not reach both flanks to 1e-6. The bell was measured without
    Mg2+, so the fit is made at ``mg = 0`` and the result carries none."""
    target = target or murayama_bell()
    sp = replace(sp or SternParams(), mg=0.0)

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
