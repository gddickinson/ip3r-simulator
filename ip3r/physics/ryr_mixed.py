"""A mixed RyR1 cluster: only some channels carry the use gate.

Every simulator up to Round 6.10 treated a cluster's channels as identical.
Laver & Lamb 1998 found they are not: of RyRs at Po > 0.2, 80 % of skeletal
(12 of 15) and 56 % of cardiac channels inactivated after voltage steps,
"one-half to two-thirds" overall, and "only those that showed inactivation
after voltage steps inactivated after [Ca2+] steps, and vice versa" -- a
stable property of the channel, not a chance of the trial. Laver & Curtis
1996 saw the same split after Ca2+ steps (70 % of 25 channels declined), and
Ma 1995 50-70 % of skeletal channels. So a cluster keeps a subpopulation
that never use-inactivates, and Round 6.9 argued that could only make
termination harder. This module measures it.

**The scheme.** :class:`MixedUseParams` is a :class:`ryr_use.UseParams`
with a fraction ``f`` (``ryr.use_inactivating_fraction``) and, once bound to
a cluster, a per-channel mask. A masked channel has the use gate; the rest
have its entry rate at zero, so they never leave ``u = 0``. Both keep the
same eight states, so the simulators' tables are unchanged. The Ca2+ gate
is shared: nothing in the papers says the two populations differ in Ca2+
sensitivity, and giving them different gates would be a second, unsourced
heterogeneity.

**The bell is a population measurement.** Murayama's bell is [3H]ryanodine
binding by a membrane preparation, the mean over every channel in it. So
the mixed scheme's :meth:`MixedUseParams.open_probability` is the mixture
``f Po_use + (1 - f) Po_plain``, and :func:`ryr_use.fit_with_use` with a
``fraction`` fits the shared Ca2+ gate to *that*. Fitting it to the use
population alone would put the use gate's share of the descending limb on
every channel, when a third of them do not have it.

**Which channels.** The mask is drawn for each run from its own random
generator (``round(f n)`` channels, without replacement), so a seed fixes
both the mask and the trajectory, and the spatial arrangement is averaged
over seeds rather than chosen.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from ..parameters import PARAMETERS as _P
from .ryr_gating import SternParams, stationary
from .ryr_use import UseParams, with_use

__all__ = ["MixedUseParams", "mixed", "bind", "fraction_values",
           "fit_mixed", "low_activity", "LA_READINGS"]


@dataclass
class MixedUseParams(UseParams):
    #: Share of channels that carry the use gate.
    fraction: float = field(default_factory=lambda: _P.value(
        "ryr.use_inactivating_fraction"))
    #: Per channel, True where the use gate is present; None until bound.
    mask: np.ndarray | None = None

    def inactivating(self) -> UseParams:
        """The population with the use gate, as a plain scheme."""
        return with_use(self._stern(), self.k_use_on, self.k_use_off)

    def non_inactivating(self) -> UseParams:
        """The population without it: entry rate zero, so ``u`` stays 0."""
        return with_use(self._stern(), 0.0, self.k_use_off)

    def _stern(self) -> SternParams:
        return SternParams(**{f: getattr(self, f) for f in (
            "k_act_on", "k_act_off", "k_inact_on", "k_inact_off", "mg",
            "k_mg_a", "mg_i")})

    def for_cluster(self, n: int, rng: np.random.Generator
                    ) -> "MixedUseParams":
        """This scheme bound to an ``n``-channel cluster: ``round(f n)``
        channels, drawn without replacement, carry the use gate."""
        mask = np.zeros(n, bool)
        mask[rng.choice(n, int(round(self.fraction * n)), replace=False)] = True
        return replace(self, mask=mask)

    def exit_rates(self, state: np.ndarray, c: np.ndarray) -> np.ndarray:
        if self.mask is None:
            raise ValueError("a mixed scheme must be bound to a cluster "
                             "(for_cluster) before it is simulated")
        if len(self.mask) != len(state):
            raise ValueError(f"mask for {len(self.mask)} channels, "
                             f"state for {len(state)}")
        r = super().exit_rates(state, c)
        n = len(state)
        # The unmasked channels cannot enter; one already in (never, from
        # a stationary start) would still leave, so only entry is removed.
        r[2 * n:] = np.where(self.mask | (state >= 4), r[2 * n:], 0.0)
        return r

    def draw(self, c: float, rng: np.random.Generator) -> np.ndarray:
        """Each bound channel's state drawn from its own population's
        stationary distribution at Ca2+ ``c``."""
        state = np.empty(len(self.mask), int)
        for keep, scheme in ((True, self.inactivating()),
                             (False, self.non_inactivating())):
            idx = np.flatnonzero(self.mask == keep)
            pi = np.cumsum(stationary(c, scheme))
            state[idx] = np.searchsorted(pi, rng.random(idx.size) * pi[-1])
        return state

    def open_probability(self, c) -> np.ndarray:
        """The population's mean open probability, as a binding bell sees
        it: ``f Po_use + (1 - f) Po_plain``."""
        f = self.fraction
        return (f * np.asarray(self.inactivating().open_probability(c))
                + (1.0 - f) * np.asarray(
                    self.non_inactivating().open_probability(c)))

    def population_available(self) -> float:
        """Mean ceiling the use gate leaves over the population."""
        return self.fraction * self.available + (1.0 - self.fraction)


def mixed(sp: UseParams, fraction: float | None = None) -> MixedUseParams:
    """``sp``'s gates on a cluster where only ``fraction`` (default the
    registered one) of channels carry its use gate. Unbound."""
    f = _P.value("ryr.use_inactivating_fraction") if fraction is None \
        else float(fraction)
    if not 0.0 <= f <= 1.0:
        raise ValueError(f"fraction {f} outside [0, 1]")
    base = {k: getattr(sp, k) for k in (
        "k_act_on", "k_act_off", "k_inact_on", "k_inact_off", "mg", "k_mg_a",
        "mg_i", "k_use_on", "k_use_off")}
    return MixedUseParams(**base, fraction=f)


def bind(gating, n: int, rng: np.random.Generator):
    """``gating`` bound to an ``n``-channel cluster if it needs to be (a
    mixed scheme), else unchanged. Every simulator calls it once per run,
    before drawing initial states, so any scheme runs anywhere."""
    return gating.for_cluster(n, rng) if hasattr(gating, "for_cluster") \
        else gating


def fraction_values() -> np.ndarray:
    """The registered grid of inactivating fractions, 0 to 1 inclusive."""
    return np.linspace(0.0, 1.0, int(round(_P.value(
        "spark.use_fraction_scan_points"))))


def fit_mixed(fraction: float | None = None, ratio: float | None = None,
              k_use_on: float | None = None,
              low_activity_reading: str | None = None) -> MixedUseParams:
    """The shared Ca2+ gate fitted so that the *population* bell has
    Murayama's half-peak points, with ``fraction`` of channels carrying a
    use gate of recovery ratio ``ratio``. With ``low_activity_reading`` (one
    of :data:`LA_READINGS`), ``ryr.la_fraction`` of the population is
    Copello's low-activity channels and the fraction applies to the rest.
    Raises ``RuntimeError`` if no gate fits (see
    :func:`ryr_use.fit_with_use`)."""
    from .ryr_use import fit_with_use
    f = _P.value("ryr.use_inactivating_fraction") if fraction is None \
        else float(fraction)
    bg = None
    if low_activity_reading is not None:
        bg = (_P.value("ryr.la_fraction"),
              lambda c: low_activity(c, low_activity_reading))
    return fit_with_use(k_use_on, ratio=ratio, fraction=f, background=bg)


# --------------------------------------------------------------------------
# Copello et al. 1997's low-activity channels (Round 6.12).
#
# About a third of skeletal RyRs gate in a low-activity mode: Po below ~0.1
# at every Ca2+, half-activated at 70-150 µM and half-inhibited at
# 100-300 µM, a narrow bump that sits under the population bell's
# descending half-point. They are not simulated in the cleft; they enter
# only the bell the high-activity channels' Ca2+ gate is fitted to.

#: How the low-activity bump is placed within Copello's ranges: both half
#: points at their low ends, at the geometric middles, or at their high
#: ends (nearest Murayama's half-inhibition, so the largest effect).
LA_READINGS = ("low", "mid", "high")


def low_activity(c, reading: str = "mid") -> np.ndarray:
    """Open probability of a low-activity channel at Ca2+ ``c`` (µM): a
    Hill activation times a Hill inhibition, peak scale ``ryr.la_po_max``,
    half points placed by ``reading``."""
    if reading not in LA_READINGS:
        raise ValueError(f"reading {reading!r} not in {LA_READINGS}")

    def at(key):
        lo, hi = (_P.value(f"ryr.la_{key}_{e}") for e in ("min", "max"))
        return {"low": lo, "mid": float(np.sqrt(lo * hi)), "high": hi}[reading]

    c = np.asarray(c, float)
    act = 1.0 / (1.0 + (at("ec50") / c) ** _P.value("ryr.la_hill_act"))
    inh = 1.0 / (1.0 + (c / at("ic50")) ** _P.value("ryr.la_hill_inh"))
    return _P.value("ryr.la_po_max") * act * inh
