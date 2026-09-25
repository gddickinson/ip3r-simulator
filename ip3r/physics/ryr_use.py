"""Use-dependent (flux-driven) RyR1 inactivation: the gate a bell cannot see.

Rounds 6.3-6.8 asked one question four ways: what ends a cleft spark, and
what stops the couplon releasing after repolarisation? Every answer tried
was a *steady-state* Ca2+ inactivation — Stern's one-Ca2+ gate, that gate
refitted to Murayama's bell, a two-site gate fitted to the bell's slope,
each with and without Mg2+, and SR depletion beside them. The result was
always the same: only Stern's Ki of 10 µM terminates, and it is 30x below
the measured bell. Round 6.8 closed with the remaining possibility, that
the terminator is an inactivation **an equilibrium binding curve cannot
measure**.

Such a process is measured. Laver & Lamb 1998 stepped RyR1 (and RyR2) in
bilayers in voltage and in cytosolic Ca2+ and found an inactivation that

* takes hold with ``tau`` about 1-3 s,
* affects one-half to two-thirds of channels,
* proceeds at a rate set by how much the channel is *open*: "inactivation
  rates increased with intraburst open probability (Po) and in proportion
  to the probability of a long-lived, RyR open state", and
* "depended on P(OL) and not on the particular activator (Ca2+ (microM),
  ATP, caffeine, and ryanodine), inhibitor (mM Ca2+ and Mg2+), or gating
  mode".

The last point is the one that matters here. The rate depends on the
channel's own conduction, not on the ligand. Put in a Markov scheme, that
means inactivation is **entered only from the conducting state**, at a
Ca2+-independent rate.

**The scheme.** Stern's two gates, unchanged and independent, plus a third:

    activation gate   a in {0, 1}   two Ca2+, on ``k_o c^2``, off ``k_o-``
    Ca2+ inactivation i in {0, 1}   one Ca2+, on ``k_i x``,   off ``k_i-``
    use gate          u in {0, 1}   on ``k_use`` FROM THE OPEN STATE ONLY,
                                    off ``k_use-``

State ``s = a + 2 i + 4 u``, so eight states, and each gate's exit flips
its own bit (``dest[s] = [s^1, s^2, s^4]``). Only ``s = 1`` (a = 1, i = 0,
u = 0) conducts; every state with ``i`` or ``u`` set is unavailable.

**It is not an equilibrium.** Around the cycle C -> O -> OU -> CU -> C the
forward rates are ``a k_use a- k_use-`` and the reverse are zero, because
C -> CU cannot happen: the gate has no way in except through conduction.
There is therefore a net stationary flux around that cycle and no detailed
balance, and the stationary open probability is *not* a product of the
three gates' equilibria (:func:`ryr_gating.SternParams.open_probability`'s
form does not apply and is overridden by a null-space solve). This is the
formal content of "a bell cannot see it": a [3H]ryanodine binding curve is
an equilibrium measurement, and this gate is driven by flux.

The strict zero is an idealisation of Laver & Lamb's result, not a claim
that a shut channel never inactivates; a real channel presumably has a
small entry rate from the closed states, and the model takes the limit in
which that rate is negligible beside entry from the open state.

**What is sourced and what is not** (Round 6.10, from the full text).
The 1-3 s is at +40 mV. Laver & Lamb's Fig. 4 puts the rate at the SR's
resting ~0 mV at 0.056 s^-1 (tau ~18 s, an extrapolation below their
data): that is ``ryr.k_use_on``. They report **no recovery rate at a fixed
potential** -- inactivated channels recovered only when the voltage was
reversed. What the model needs is not a recovery rate anyway: both rates
are seconds, a spark is milliseconds, so the bell, the refitted Ca2+ gate
and spark termination depend only on the ratio ``rho = k_use-/k_use``
(``ryr.use_recovery_ratio``; tested). The paper bounds rho only at +40 mV,
through the residual activity of its Fig. 8 (:func:`residual_bound`). At
0 mV rho is unmeasured, so it is scanned (:func:`bell_panel`,
:func:`spark_termination.ratio_scan`), not trusted.

The per-channel heterogeneity — that only half to two-thirds of channels
inactivate at all, stably, so that a cluster keeps a subpopulation that
never inactivates — is :mod:`ryr_mixed` (Round 6.11). This module's
schemes are the fraction-1 limit of it.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import ClassVar

import numpy as np

from ..parameters import PARAMETERS as _P
from .ryr_gating import SternParams, stationary

__all__ = ["STATES", "OPEN", "UseParams", "with_use", "no_ca_gate",
           "open_probability_no_ca_gate", "cycle_flux", "fit_with_use",
           "BellRow", "bell_panel", "tau_values", "ratio_values",
           "use_rate_at", "residual_bound"]

#: ``s = a + 2 i + 4 u``: activation gate open, Ca2+-inactivated, used up.
STATES = ("C", "O", "CI", "I", "CU", "OU", "CIU", "IU")
OPEN = 1


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class UseParams(SternParams):
    #: Each gate's exit flips its own bit of ``s``: activation, Ca2+
    #: inactivation, use gate.
    dest: ClassVar[np.ndarray] = np.array(
        [[s ^ 1, s ^ 2, s ^ 4] for s in range(8)])
    open_mask: ClassVar[np.ndarray] = np.arange(8) == OPEN
    inact_mask: ClassVar[np.ndarray] = np.arange(8) > OPEN
    #: The trigger opens the activation gate of an available channel only.
    trigger_map: ClassVar[np.ndarray] = np.array([1, 1, 2, 3, 4, 5, 6, 7])

    k_use_on: float = _v("ryr.k_use_on")      # s^-1, from the open state
    k_use_off: float = field(default_factory=lambda: _P.value(
        "ryr.use_recovery_ratio") * _P.value("ryr.k_use_on"))   # s^-1

    @property
    def tau_use(self) -> float:
        """Mean time to inactivate while held open, s."""
        return 1.0 / self.k_use_on if self.k_use_on > 0 else float("inf")

    @property
    def tau_use_recover(self) -> float:
        return 1.0 / self.k_use_off if self.k_use_off > 0 else float("inf")

    @property
    def available(self) -> float:
        """Share of channels the use gate leaves available to a channel held
        open indefinitely: ``k_use-/(k_use + k_use-)``. The ceiling on the
        open probability that no Ca2+ concentration can beat."""
        total = self.k_use_on + self.k_use_off
        return self.k_use_off / total if total > 0 else 1.0

    def exit_rates(self, state: np.ndarray, c: np.ndarray) -> np.ndarray:
        a, i, u = state % 2, (state // 2) % 2, state // 4
        r_act = np.where(a == 0, self.k_act_on * self.mg_factor * c * c,
                         self.k_act_off)
        r_ca = np.where(i == 0, self.k_inact_on * (c + self.mg_inact),
                        self.k_inact_off)
        # In only from the conducting state; out from anywhere it is set.
        r_use = np.where(u == 1, self.k_use_off,
                         np.where(state == OPEN, self.k_use_on, 0.0))
        return np.concatenate([r_act, r_ca, r_use * np.ones_like(c)])

    def generator(self, c: float) -> np.ndarray:
        q = np.zeros((8, 8))
        states = np.arange(8)
        rates = self.exit_rates(states, np.full(8, float(c))).reshape(3, 8).T
        for s in states:
            for kind in range(3):
                q[s, self.dest[s, kind]] += rates[s, kind]
        q[np.diag_indices(8)] = -q.sum(axis=1)
        return q

    def open_probability(self, c) -> np.ndarray:
        """Stationary occupancy of ``O``, from the generator's null space.

        There is no product form to fall back on: the use gate's entry rate
        depends on the other two gates, and the cycle it closes carries net
        flux (see the module docstring), so ``SternParams``'s closed form
        would be wrong rather than merely slower.
        """
        c = np.asarray(c, float)
        po = np.array([stationary(float(x), self)[OPEN]
                       for x in np.atleast_1d(c).ravel()])
        return po.reshape(c.shape) if c.ndim else float(po[0])


def with_use(sp: SternParams | None = None, k_use_on: float | None = None,
             k_use_off: float | None = None,
             ratio: float | None = None) -> UseParams:
    """``sp``'s activation and Ca2+ inactivation gates (Stern's, or a gate
    fitted to a bell) with the use gate added. Defaults are the registered
    rate and ratio; without ``k_use_off`` the recovery is ``ratio`` (default
    ``ryr.use_recovery_ratio``) times the entry rate, so changing the speed
    keeps the steady state. ``sp`` must be a one-site scheme: a
    ``TwoSiteParams``'s second site is not carried, so passing one is
    refused rather than silently dropped."""
    sp = sp or SternParams()
    if type(sp) is not SternParams:
        raise TypeError("with_use takes the one-site Ca2+ gate "
                        f"(SternParams), not {type(sp).__name__}")
    base = {f: getattr(sp, f) for f in ("k_act_on", "k_act_off", "k_inact_on",
                                        "k_inact_off", "mg", "k_mg_a", "mg_i")}
    on = _P.value("ryr.k_use_on") if k_use_on is None else float(k_use_on)
    if k_use_off is None:
        rho = _P.value("ryr.use_recovery_ratio") if ratio is None else ratio
        k_use_off = rho * on
    return UseParams(**base, k_use_on=on, k_use_off=float(k_use_off))


def use_rate_at(v: float) -> float:
    """Laver & Lamb's Fig. 4 line at bilayer potential ``v`` (volts,
    cytosol relative to lumen, positive limb): the macroscopic inactivation
    rate, s^-1. Exists so the tests can hold the reading of the figure to
    the paper's abstract (tau 1-3 s at +40 mV) and charge (z delta 1.14)."""
    return _P.value("ryr.k_use_on") * 10.0 ** (_P.value("ryr.use_rate_slope")
                                               * v)


def residual_bound(residual: float) -> float:
    """Upper bound on the recovery ratio from a residual activity ``R``
    (activity long after the step over its peak): a channel held at open
    probability ``Po`` keeps ``rho / (rho + Po)`` available, so
    ``rho = R Po / (1 - R) <= R / (1 - R)``; and a residual read before
    steady state is itself an upper bound."""
    return residual / (1.0 - residual)


def no_ca_gate(sp: SternParams) -> SternParams:
    """``sp`` with the Ca2+ inactivation gate removed (``k_inact_on`` 0, so
    ``Ki`` is infinite). The activation gate and any use gate are kept: for
    a plain ``SternParams`` this is the "no inactivation at all" control, in
    which only the cleft's geometry can end a spark (Laver et al. 2013's
    induction decay)."""
    return replace(sp, k_inact_on=0.0)


def open_probability_no_ca_gate(c: float, sp: UseParams) -> float:
    """Closed form of the stationary open probability when the Ca2+ gate is
    removed, so that only C, O, OU and CU are reachable.

    Solved by hand from the four balance equations, independently of
    :meth:`UseParams.open_probability`'s null space, which the tests hold it
    to. With ``a`` the activation on rate at ``c``, ``a-`` its off rate,
    ``u`` and ``u-`` the use gate's rates and ``S = a + a- + u-``:

        C/O = a- (S + u) / (a S),  OU/O = u (a + u-) / (u- S),
        CU/O = u a- / (u- S).
    """
    a = sp.k_act_on * sp.mg_factor * c * c
    a_, u, u_ = sp.k_act_off, sp.k_use_on, sp.k_use_off
    s = a + a_ + u_
    ratios = (a_ * (s + u) / (a * s), u * (a + u_) / (u_ * s), u * a_ / (u_ * s))
    return 1.0 / (1.0 + sum(ratios))


def cycle_flux(c: float, sp: UseParams) -> float:
    """Net stationary flux around C -> O -> OU -> CU -> C, per second, for a
    scheme whose Ca2+ gate is removed (:func:`no_ca_gate`).

    Those four states are then the whole reachable chain and the cycle is
    isolated, so the net flux is the same across each of its edges; it is
    read off the one irreversible edge, CU -> C. Zero for any scheme in
    detailed balance, and positive here because there is no way into the
    use gate except through conduction. This is the measure of how far from
    equilibrium the gate drives the channel, and hence of why an
    equilibrium binding curve cannot report it.
    """
    pi = stationary(c, sp)
    _C, _CU = 0, 4
    return float(pi[_CU] * sp.generator(c)[_CU, _C])


# --------------------------------------------------------------------------
# The bell as a constraint on the use gate.
#
# Fitting the Ca2+ gate to the measured bell and *then* adding a use gate
# counts the measurement twice: if the channel has both gates, Murayama's
# bell is already the composite of them. The right question is whether any
# Ca2+ gate, put beside a use gate of a given speed, reproduces the measured
# bell at all. :func:`fit_with_use` asks it and :func:`bell_panel` scans it.

def fit_with_use(k_use_on: float | None = None, target=None,
                 sp: SternParams | None = None,
                 k_use_off: float | None = None,
                 ratio: float | None = None,
                 fraction: float | None = None) -> UseParams:
    """Ka and Ki of the Ca2+ gate such that the *composite* scheme's
    half-peak points are ``target``'s (default Murayama's 25 C bell), with a
    use gate of rate ``k_use_on`` already present. With ``fraction``, only
    that share of channels carries the use gate and the fit is to the
    population's mean bell (:mod:`ryr_mixed`).

    The same solve as :func:`ryr_gating.fit_to_bell`, in log space, with the
    use gate inside the loop. Raises ``RuntimeError`` if no Ca2+ gate
    reproduces the measured bell beside this use gate — which is the
    informative outcome, not a failure of the solver.
    """
    from scipy.optimize import fsolve

    from .bell import measure_bell
    from .ryr_gating import murayama_bell, with_constants
    target = target or murayama_bell()
    sp = replace(sp or SternParams(), mg=0.0)

    def scheme(ka, ki):
        s = with_use(with_constants(ka, ki, sp), k_use_on, k_use_off, ratio)
        if fraction is None:
            return s
        from .ryr_mixed import mixed
        return mixed(s, fraction)

    def resid(x):
        with np.errstate(over="ignore"):
            ka, ki = np.exp(x)
        if not (np.isfinite(ka) and np.isfinite(ki)) or min(ka, ki) <= 0:
            return [1e3, 1e3]                    # the solve has wandered off
        try:
            b = measure_bell(scheme(ka, ki).open_probability)
        except (ValueError, FloatingPointError):
            # No half-peak crossing in the searched range, or a generator
            # too degenerate to have a stationary state.
            return [1e3, 1e3]
        if not (b.po_peak > 0):
            return [1e3, 1e3]
        if not np.isfinite(b.c_half_inh):
            # No descending limb within the measurable range: report the
            # distance to the edge of it so the solve has a gradient.
            return [np.log(b.c_half_act / target.c_half_act), 1e3]
        return [np.log(b.c_half_act / target.c_half_act),
                np.log(b.c_half_inh / target.c_half_inh)]

    # Several starts: one stalled solve is not evidence that no gate fits
    # (a single start stalled at rho 0.18 between two that fitted).
    x0 = np.log([target.c_half_act, target.c_half_inh])
    msg = ""
    for shift in _STARTS:
        x, _, ok, msg = fsolve(resid, x0 + np.asarray(shift), full_output=True)
        if ok == 1 and max(abs(r) for r in resid(x)) <= 1e-6:
            return scheme(*np.exp(x))
    raise RuntimeError(
        f"no Ca2+ gate reproduces the bell beside this use gate "
        f"(k_use {k_use_on}, k_use- {k_use_off}, ratio {ratio}): {msg}")


#: fsolve starting offsets in (log Ka, log Ki) around the measured half-peak
#: points, nearest first. A small recovery ratio drives the fitted Ka up and
#: Ki down by decades (ratio 0.002: Ka x8, Ki /900), so the grid reaches
#: there; the first start that lands is kept, so a fit is usually one solve.
_STARTS = tuple((a, i) for i in (0.0, -1.0, -2.0, -3.5, -5.0, -7.0, 1.0)
                for a in (0.0, 1.0, 2.0, -0.5))


@dataclass(frozen=True)
class BellRow:
    """One recovery ratio against the measured bell."""
    ratio: float                # k_use- / k_use
    available: float            # ceiling the gate leaves on P_open
    po_peak: float
    c_half_act: float           # µM
    c_half_inh: float           # µM, nan if the curve never falls to half
    width_decades: float
    fitted: bool                # a Ca2+ gate reproduces the measured bell
    k_a: float                  # µM, of that fit (nan if none)
    k_i: float                  # µM

    def row(self) -> str:
        fit = (f"Ka {self.k_a:6.2f} Ki {self.k_i:7.1f}" if self.fitted
               else "no Ca2+ gate fits  ")
        return (f"ratio {self.ratio:7.4g}  ceiling {self.available:7.4f}"
                f"  peak {self.po_peak:7.4f}  half {self.c_half_act:6.2f} -"
                f" {self.c_half_inh:9.1f} µM  width {self.width_decades:5.2f}"
                f" dec  {fit}")


def tau_values() -> np.ndarray:
    """The registered grid of use-gate time constants, s."""
    return np.geomspace(_P.value("spark.use_scan_min"),
                        _P.value("spark.use_scan_max"),
                        int(round(_P.value("spark.use_scan_points"))))


def ratio_values() -> np.ndarray:
    """The registered grid of recovery ratios around ``ryr.use_recovery_ratio``."""
    hi = _P.value("spark.use_ratio_scan_max")
    n = int(round(_P.value("spark.use_ratio_scan_points")))
    return _P.value("ryr.use_recovery_ratio") * np.geomspace(1.0 / hi, hi, n)


def bell_panel(ratios=None, sp: SternParams | None = None) -> list[BellRow]:
    """Each recovery ratio measured two ways: the bell it gives beside the
    Ca2+ gate fitted *without* it (the double-counting Round 6.9 avoids), and
    whether any Ca2+ gate reproduces the measured bell beside it. The speed
    is the registered one; it does not enter (tested)."""
    from .bell import Bell, measure_bell
    from .ryr_gating import fit_to_bell
    base = sp or fit_to_bell()
    rows = []
    for rho in (ratio_values() if ratios is None else ratios):
        scheme = with_use(base, ratio=rho)
        try:
            b = measure_bell(scheme.open_probability)
        except (ValueError, FloatingPointError):
            b = Bell(float("nan"), float("nan"), float("nan"), float("nan"))
        try:
            fit = fit_with_use(ratio=rho)
            k_a, k_i, ok = fit.k_a, fit.k_i, True
        except (RuntimeError, FloatingPointError):
            k_a = k_i = float("nan")
            ok = False
        rows.append(BellRow(float(rho), scheme.available, b.po_peak,
                            b.c_half_act, b.c_half_inh,
                            b.width_decades, ok, k_a, k_i))
    return rows
