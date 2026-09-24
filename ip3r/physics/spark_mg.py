"""Mg2+ and the cleft spark: does the missing terminator end it?

Round 6.4 fitted Stern's scheme to Murayama's Mg2+-free bell and found that
the array, once lit, never shuts. This module adds cytosolic Mg2+ through
the two sourced actions of :class:`ryr_gating.SternParams` (competition at
the activation site, ``ryr.k_mg_a``; equal binding at the inactivation
site, ``ryr.mg_i_relative``) and measures sparks two ways:

* **spontaneous** (:func:`spark_termination.measure`): does a spark start by
  itself, and does it end?
* **triggered** (:func:`triggered`): every channel that is closed but not
  inactivated is opened at t = 0, a stand-in for the V channels that start
  a spark in the fibre (they are not simulated). The spark is timed to the
  first moment no channel is open. Under Mg2+ nothing starts by itself, so
  this is the only way to ask whether a spark would end.

:func:`dissect` separates the two Mg2+ sites at the fibre's free Mg2+, and
:func:`mg_scan` runs from 0 to it. The activation-site affinity has two
readings: Laver 2004's absolute 54 µM (measured with ATP), their
Mg2+/Ca2+ selectivity carried onto the fitted Ka (:func:`k_mg_a_by_ratio`),
and Meissner 1997's constant, measured in the [3H]ryanodine assay itself and
carried into Murayama's 0.17 M NaCl with Na+ competing at the same site
(:mod:`mg_competition`). All three are reported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from .mg_competition import equivalent_k_mg_a
from .ryr_gating import SternParams, fit_to_bell, with_mg
from .spark_termination import Termination, measure
from .sparks_cleft import CleftSparkParams, simulate_sparks_cleft

__all__ = ["Triggered", "triggered", "mg_values", "k_mg_a_by_ratio",
           "READINGS", "k_mg_a_reading",
           "mg_scan", "dissect", "spontaneous"]


@dataclass(frozen=True)
class Triggered:
    label: str
    mg: float                  # µM free
    k_a_eff: float             # µM, half-activation under this Mg2+
    trials: int
    opened: float              # median channels open at t = 0
    ended: int                 # trials in which every channel closed
    duration_ms: float         # median time to all closed (ended trials)
    inactivated_start: float   # median channels inactivated at t = 0
    inactivated_end: float     # median at the moment the last one closed

    def row(self) -> str:
        return (f"{self.label:30s} Mg {self.mg:6.0f}  Ka' {self.k_a_eff:6.1f}  "
                f"opened {self.opened:4.1f}  ended {self.ended:2d}/{self.trials}"
                f"  median {self.duration_ms:7.1f} ms  inactivated "
                f"{self.inactivated_start:4.1f} -> {self.inactivated_end:4.1f}")


def triggered(sp: SternParams, label: str, trials: int | None = None,
              window: float | None = None) -> Triggered:
    """Open the available channels at t = 0 and time the array to shut."""
    trials = int(round(_P.value("spark.trigger_trials"))) if trials is None else trials
    window = _P.value("spark.trigger_window") if window is None else window
    opened, ends, i0, i1 = [], [], [], []
    for seed in range(trials):
        tr = simulate_sparks_cleft(0.0, window, seed, CleftSparkParams(gating=sp),
                                   trigger=True)
        opened.append(tr.n_open[0])
        i0.append(tr.n_inactivated[0])
        shut = np.flatnonzero(tr.peaks == 0)
        if len(shut):
            ends.append(tr.t[shut[0]])
            i1.append(tr.n_inactivated[shut[0]])
    med = lambda v: float(np.median(v)) if len(v) else float("nan")
    return Triggered(label, sp.mg, sp.k_a_eff, trials, med(opened), len(ends),
                     1e3 * med(ends), med(i0), med(i1))


def mg_values() -> np.ndarray:
    """0, then ``spark.mg_scan_points`` geometric steps up to ``ryr.mg_free``."""
    return np.concatenate([[0.0], np.geomspace(
        _P.value("spark.mg_scan_min"), _P.value("ryr.mg_free"),
        int(round(_P.value("spark.mg_scan_points"))))])


def k_mg_a_by_ratio(base: SternParams) -> float:
    """The activation-site Mg2+ affinity if Laver 2004's Mg2+/Ca2+
    selectivity (54 / 0.51) holds and only Ca2+'s affinity is Murayama's."""
    return base.k_a * _P.value("ryr.k_mg_a") / _P.value("ryr.k_ca_a_laver")


#: The two readings of the activation-site Mg2+ affinity (the largest
#: uncertainty left: Laver's absolute value was measured with ATP, which
#: Murayama's bell was not).
READINGS = {"measured": "measured (Laver 2004, with ATP)",
            "selectivity": "Laver 2004's Mg²⁺/Ca²⁺ selectivity on this Ka",
            "meissner": "Meissner 1997 in Murayama's 0.17 M NaCl (Na⁺ competing)"}


def k_mg_a_reading(base: SternParams, reading: str) -> float:
    """K_Mg,A under one of :data:`READINGS` for the scheme ``base``."""
    if reading == "measured":
        return _P.value("ryr.k_mg_a")
    if reading == "selectivity":
        return k_mg_a_by_ratio(base)
    if reading == "meissner":
        return equivalent_k_mg_a()
    raise ValueError(f"unknown K_Mg,A reading {reading!r}; one of {list(READINGS)}")


def _base(base: SternParams | None) -> SternParams:
    return base if base is not None else fit_to_bell()


def mg_scan(base: SternParams | None = None, k_mg_a: float | None = None,
            trials: int | None = None, mgs=None) -> list[Triggered]:
    """Triggered sparks of ``base`` (default: fitted to Murayama 25 C) over
    free Mg2+, at activation-site affinity ``k_mg_a`` (default registered)."""
    sp = _base(base)
    return [triggered(with_mg(sp, mg, k_mg_a), f"Mg {mg:.3g} µM", trials)
            for mg in (mg_values() if mgs is None else mgs)]


def dissect(base: SternParams | None = None, mg: float | None = None,
            trials: int | None = None) -> list[Triggered]:
    """At ``mg`` (default the fibre's): no Mg2+, each site alone, both."""
    sp = _base(base)
    mg = _P.value("ryr.mg_free") if mg is None else mg
    rows = [("no Mg2+", with_mg(sp, 0.0)),
            ("activation site only", with_mg(sp, mg, mg_i=0.0)),
            ("inactivation site only", with_mg(sp, mg, k_mg_a=np.inf)),
            ("both sites", with_mg(sp, mg))]
    return [triggered(s, label, trials) for label, s in rows]


def spontaneous(base: SternParams | None = None, k_mg_a: float | None = None,
                duration: float = 10.0, seeds: int = 4, mgs=None
                ) -> list[Termination]:
    """Untriggered runs over free Mg2+ (Round 6.4's ruler)."""
    sp = _base(base)
    return [measure(with_mg(sp, mg, k_mg_a), f"Mg {mg:.3g} µM", duration, seeds)
            for mg in (mg_values() if mgs is None else mgs)]
