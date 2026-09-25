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

**What is sourced and what is not.** ``ryr.k_use_on`` comes from the
measured 1-3 s, and is a *lower* bound on the rate from the open state
because the measured tau is macroscopic (rate ~ Po). ``ryr.k_use_off``
could not be read: the paper's full text is a paywalled page scan and
recovery is reported only qualitatively. Neither number is trusted here.
The question is answered by scanning the time constant over four decades
(:func:`spark_termination.use_scan`) and asking where termination lives,
not by running the model at one value.

The per-channel heterogeneity — that only half to two-thirds of channels
inactivate at all, stably, so that a cluster keeps a subpopulation that
never inactivates — is *not* modelled: every simulator here treats the
channels as identical. It is recorded as emergent, and it can only make
termination harder than this module reports.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import ClassVar

import numpy as np

from ..parameters import PARAMETERS as _P
from .ryr_gating import SternParams, stationary

__all__ = ["STATES", "OPEN", "UseParams", "with_use", "no_ca_gate",
           "open_probability_no_ca_gate", "cycle_flux", "fit_with_use",
           "BellRow", "bell_panel", "tau_values"]

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
    k_use_off: float = _v("ryr.k_use_off")    # s^-1

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
             k_use_off: float | None = None) -> UseParams:
    """``sp``'s activation and Ca2+ inactivation gates (Stern's, or a gate
    fitted to a bell) with the use gate added. Defaults are the registered
    rates. ``sp`` must be a one-site scheme: a ``TwoSiteParams``'s second
    site is not carried, so passing one is refused rather than silently
    dropped."""
    sp = sp or SternParams()
    if type(sp) is not SternParams:
        raise TypeError("with_use takes the one-site Ca2+ gate "
                        f"(SternParams), not {type(sp).__name__}")
    base = {f: getattr(sp, f) for f in ("k_act_on", "k_act_off", "k_inact_on",
                                        "k_inact_off", "mg", "k_mg_a", "mg_i")}
    kw = {}
    if k_use_on is not None:
        kw["k_use_on"] = float(k_use_on)
    if k_use_off is not None:
        kw["k_use_off"] = float(k_use_off)
    return UseParams(**base, **kw)


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
                 k_use_off: float | None = None) -> UseParams:
    """Ka and Ki of the Ca2+ gate such that the *composite* scheme's
    half-peak points are ``target``'s (default Murayama's 25 C bell), with a
    use gate of rate ``k_use_on`` already present.

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
    scheme = lambda ka, ki: with_use(with_constants(ka, ki, sp), k_use_on,
                                     k_use_off)

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

    x, _, ok, msg = fsolve(resid, np.log([target.c_half_act,
                                          target.c_half_inh]),
                           full_output=True)
    if ok != 1 or max(abs(r) for r in resid(x)) > 1e-6:
        raise RuntimeError(
            f"no Ca2+ gate reproduces the bell beside a use gate of "
            f"{k_use_on:g} s^-1: {msg}")
    return scheme(*np.exp(x))


@dataclass(frozen=True)
class BellRow:
    """One use-gate speed against the measured bell."""
    tau_use: float              # s, from the open state
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
        return (f"use tau {self.tau_use:8.4g} s  ceiling {self.available:7.4f}"
                f"  peak {self.po_peak:7.4f}  half {self.c_half_act:6.2f} -"
                f" {self.c_half_inh:9.1f} µM  width {self.width_decades:5.2f}"
                f" dec  {fit}")


def tau_values() -> np.ndarray:
    """The registered grid of use-gate time constants, s."""
    return np.geomspace(_P.value("spark.use_scan_min"),
                        _P.value("spark.use_scan_max"),
                        int(round(_P.value("spark.use_scan_points"))))


def bell_panel(taus=None, sp: SternParams | None = None) -> list[BellRow]:
    """Each use-gate speed measured two ways: the bell it gives beside the
    Ca2+ gate fitted *without* it (the double-counting Round 6.9 avoids), and
    whether any Ca2+ gate reproduces the measured bell beside it."""
    from .bell import Bell, measure_bell
    from .ryr_gating import fit_to_bell
    base = sp or fit_to_bell()
    rows = []
    for tau in (tau_values() if taus is None else taus):
        scheme = with_use(base, 1.0 / tau)
        try:
            b = measure_bell(scheme.open_probability)
        except (ValueError, FloatingPointError):
            b = Bell(float("nan"), float("nan"), float("nan"), float("nan"))
        try:
            fit = fit_with_use(1.0 / tau)
            k_a, k_i, ok = fit.k_a, fit.k_i, True
        except (RuntimeError, FloatingPointError):
            k_a = k_i = float("nan")
            ok = False
        rows.append(BellRow(float(tau), scheme.available, b.po_peak,
                            b.c_half_act, b.c_half_inh,
                            b.width_decades, ok, k_a, k_i))
    return rows
