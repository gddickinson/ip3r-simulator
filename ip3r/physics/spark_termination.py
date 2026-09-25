"""What ends a spark in the cleft: the inactivation gate, scanned.

Round 6.3 left cleft sparks 3x longer than the measured release and blamed
the gating scheme, whose steady-state inactivation is 6.7x more sensitive
than Murayama's measured bell. This module tests that suspicion by
measurement rather than by one refit:

* :func:`refit` fits Stern's Ka and Ki to Murayama's bell (25 and 37 C) and
  runs the cleft array with each.
* :func:`ki_scan` holds Stern's rates' scale and moves Ki from below his
  10 µM to beyond the measured KI.
* :func:`rate_scan` holds Ka and Ki and scales both inactivation rates together:
  the one thing a steady-state bell cannot fix.
* :func:`use_scan` (Round 6.9) puts a use-dependent gate
  (:mod:`ip3r.physics.ryr_use`) beside the Ca2+ gate and fits *both* to the
  measured bell, so the bell is not counted twice. It scans the gate's
  speed at the registered recovery ratio; :func:`ratio_scan` (Round 6.10)
  scans the ratio, the quantity the result actually rests on.
  :func:`no_inactivation` is the control beneath them all: neither gate, so
  only the cleft's geometry can end a spark.

Each row reads the sparks with :func:`puff_compare.spark_ends`. A spark
that is still running when the trace ends is counted as ``unterminated``,
and its duration is only a lower bound.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from .puff_compare import recruitment, spark_ends
from .ryr_gating import (MurayamaParams, SternParams, fit_to_bell,
                         murayama_bell, with_constants)
from .ryr_use import fit_with_use, no_ca_gate, ratio_values, tau_values
from .sparks_cleft import CleftSparkParams, simulate_sparks_cleft

__all__ = ["Termination", "measure", "refit", "ki_scan", "rate_scan",
           "ki_values", "rate_values", "use_scan", "no_inactivation",
           "ratio_scan"]


@dataclass(frozen=True)
class Termination:
    label: str
    k_a: float                 # µM
    k_i: float                 # µM
    k_inact_on: float          # µM^-1 s^-1
    n_sparks: int
    per_s: float               # sparks reaching half the cluster, per second
    duration_ms: float         # median; nan without sparks
    unterminated: int          # sparks still running at the trace's end
    inactivated_end: float     # median channels inactivated at the end
    open_fraction: float       # of channel-time, over the whole run

    def row(self) -> str:
        return (f"{self.label:24s} Ka {self.k_a:6.2f} Ki {self.k_i:7.1f} "
                f"k_i {self.k_inact_on:6.2f}  {self.n_sparks:3d} sparks "
                f"({self.per_s:4.2f}/s)  median {self.duration_ms:7.1f} ms  "
                f"unended {self.unterminated:2d}  inactivated at end "
                f"{self.inactivated_end:4.1f}  open {self.open_fraction:.3f}")


def measure(sp: SternParams, label: str, duration: float = 10.0,
            seeds: int = 4) -> Termination:
    """The cleft array under gating ``sp``, pooled over ``seeds`` runs."""
    ends, per_s, open_frac = [], [], []
    for seed in range(seeds):
        tr = simulate_sparks_cleft(0.0, duration, seed, CleftSparkParams(gating=sp))
        r = recruitment(tr)
        ends += spark_ends(tr)
        per_s.append(r["large_per_s"])
        open_frac.append(r["open_fraction"])
    med = lambda k: float(np.median([e[k] for e in ends])) if ends else float("nan")
    return Termination(label, sp.k_a, sp.k_i, sp.k_inact_on, len(ends),
                       float(np.mean(per_s)), 1e3 * med("duration"),
                       sum(e["unterminated"] for e in ends),
                       med("inactivated_end"), float(np.mean(open_frac)))


def refit(duration: float = 10.0, seeds: int = 4) -> list[Termination]:
    """Stern as published, then fitted to Murayama at 25 and 37 C."""
    return [measure(SternParams(), "Stern 1997", duration, seeds)] + [
        measure(fit_to_bell(murayama_bell(mp)), f"fitted to Murayama {t}",
                duration, seeds)
        for t, mp in (("25 C", MurayamaParams()), ("37 C", MurayamaParams.at_37()))]


def ki_values() -> np.ndarray:
    return np.geomspace(_P.value("spark.ki_scan_min"), _P.value("spark.ki_scan_max"),
                        int(round(_P.value("spark.ki_scan_points"))))


def rate_values() -> np.ndarray:
    hi = _P.value("spark.rate_scan_max")
    return np.geomspace(1.0 / hi, hi, int(round(_P.value("spark.rate_scan_points"))))


def ki_scan(duration: float = 10.0, seeds: int = 4, kis=None) -> list[Termination]:
    """Stern's scheme with Ki moved (his Ka and inactivation on rate kept)."""
    sp = SternParams()
    return [measure(with_constants(sp.k_a, ki, sp), f"Ki {ki:.3g} µM", duration, seeds)
            for ki in (ki_values() if kis is None else kis)]


def rate_scan(base: SternParams | None = None, duration: float = 10.0,
              seeds: int = 4, rates=None) -> list[Termination]:
    """Both inactivation rates of ``base`` (default Stern's scheme) scaled
    together, its Ka and Ki fixed: how fast the gate is, which the bell
    leaves free."""
    sp = base or SternParams()
    return [measure(with_constants(sp.k_a, sp.k_i, sp, r), f"rate x{r:.3g}",
                    duration, seeds)
            for r in (rate_values() if rates is None else rates)]


def use_scan(duration: float = 10.0, seeds: int = 4, taus=None
             ) -> list[Termination]:
    """The cleft array under a use-dependent gate at each speed, recovery
    held at the registered ratio to it.

    For each time constant the Ca2+ gate is refitted *with the use gate
    present* (:func:`ryr_use.fit_with_use`), so the scheme reproduces
    Murayama's half-peak points as the fitted one-site scheme does and the
    two rows are comparable. At a fixed ratio every speed gives the same
    fit; the rows differ only where the gate is fast enough to act within a
    spark. A speed with no fit is skipped.
    """
    rows = []
    for tau in (tau_values() if taus is None else taus):
        try:
            sp = fit_with_use(1.0 / tau)
        except (RuntimeError, FloatingPointError):
            continue
        rows.append(measure(sp, f"use tau {tau:.3g} s", duration, seeds))
    return rows


def no_inactivation(duration: float = 10.0, seeds: int = 4) -> Termination:
    """Neither inactivation gate: does the cleft's geometry alone end a
    spark? This is Laver et al. 2013's "induction decay" as this model can
    put it, and the control that says whether the rows above are measuring
    inactivation or measuring the cleft."""
    return measure(no_ca_gate(fit_to_bell()), "no inactivation", duration,
                   seeds)


def ratio_scan(k_use_on: float | None = None, duration: float = 10.0,
               seeds: int = 4, ratios=None) -> list[Termination]:
    """The cleft array against the use gate's recovery ratio
    ``k_use- / k_use`` -- the one number Round 6.9's result rests on, and
    one Laver & Lamb 1998 bound only at +40 mV.

    The ratio is load-bearing and not a free knob: it sets how much of the
    measured bell's descending limb the use gate accounts for, so the Ca2+
    gate is refitted at every point and ``Ki`` moves with it. ``k_use_on``
    defaults to the registered rate at 0 mV; it does not change the fit.
    """
    rows = []
    for rho in (ratio_values() if ratios is None else ratios):
        try:
            sp = fit_with_use(k_use_on, ratio=rho)
        except (RuntimeError, FloatingPointError):
            continue
        rows.append(measure(sp, f"ratio {rho:.3g}", duration, seeds))
    return rows
