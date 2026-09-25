"""The use gate in the millisecond regime (Round 6.13).

Rounds 6.9-6.12 ran Laver & Lamb 1998's use-dependent gate at bilayer
speed: entry ``ryr.k_use_on`` 0.056 s^-1 at 0 mV (tau ~18 s), recovery a
ratio ``rho`` of it. Both rates are then far slower than a spark, and Round
6.10 showed that the bell, the refitted Ca2+ gate and termination depend
only on ``rho``.

Rios & Pizarro 2026 (J Gen Physiol, PMC13387315) build the same gate --
entered only from the open state, "every channel that opens inactivates",
no explicit Ca2+ role -- into the same couplon, fitted to the cell-level
flux of mammalian fibres. Their rates are **milliseconds**: O -> I1 half
time 3.5 ms, I1 -> C 20 ms (``ryr.rios_i1_half_time``,
``ryr.rios_r1_half_time``), a ratio 0.175 that falls inside the band Round
6.10 needed, at a speed ~3500x ours. This module asks what the speed adds:

* :func:`speed_panel` -- the gate's speed scanned at Rios's ratio, from the
  registered bilayer rate to Rios's, the Ca2+ gate refitted to Murayama's
  bell beside it at each (does "only rho matters" survive?), each read by
  spontaneous sparks (:func:`spark_termination.measure`) *and* by a
  triggered array (:func:`spark_mg.triggered`), because a gate this fast
  can stop a spark from igniting as well as end one.
* :func:`controls` -- the Ca2+ gate removed from each fit (what the use
  gate ends by itself), and Rios's own claim: his gate beside the fitted
  activation gate and *no* Ca2+ gate at all.
* :func:`fraction_panel` -- Round 6.11's mixed cluster (the measured 0.8 of
  channels carrying the gate) at both speeds.

Two differences from Rios's scheme are deliberate and stated: the deep state
I2 is not carried (it slows recovery between pulses, it cannot shorten a
spark), and our use gate recovers without first closing the activation
gate, where his I1 returns to C. Both are second order at spark time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from .ryr_gating import SternParams, fit_to_bell
from .ryr_use import UseParams, fit_with_use, no_ca_gate, with_use
from .spark_mg import Triggered, triggered
from .spark_termination import Termination, measure

__all__ = ["rios_rates", "rios_ratio", "speed_values", "SpeedRow",
           "read", "speed_panel", "controls", "fraction_panel"]


def rios_rates() -> tuple[float, float]:
    """Rios & Pizarro's use-gate entry and recovery rates, s^-1 (ln2 over
    each half time)."""
    return (np.log(2.0) / (1e-3 * _P.value("ryr.rios_i1_half_time")),
            np.log(2.0) / (1e-3 * _P.value("ryr.rios_r1_half_time")))


def rios_ratio() -> float:
    on, off = rios_rates()
    return off / on


def speed_values() -> np.ndarray:
    """Entry rates from the registered bilayer rate to Rios's, s^-1, on the
    use-scan's point count (``spark.use_scan_points``), both ends exact."""
    n = int(round(_P.value("spark.use_scan_points")))
    return np.geomspace(_P.value("ryr.k_use_on"), rios_rates()[0], n)


@dataclass(frozen=True)
class SpeedRow:
    """One scheme read both ways: sparks arising by themselves, and the
    whole array opened at t = 0."""
    label: str
    k_use_on: float             # s^-1
    k_use_off: float            # s^-1
    spontaneous: Termination
    trig: Triggered

    def row(self) -> str:
        s, t = self.spontaneous, self.trig
        tau = (f"{1e3 / self.k_use_on:9.2f} ms" if self.k_use_on > 0
               else "   no gate")
        return (f"{self.label:30s} tau {tau}  "
                f"Ka {s.k_a:6.2f} Ki {s.k_i:7.1f}  spont {s.n_sparks:3d} "
                f"({s.per_s:4.2f}/s) {s.duration_ms:7.1f} ms unended "
                f"{s.unterminated:2d}  | triggered ended {t.ended:2d}/"
                f"{t.trials} {t.duration_ms:6.1f} ms  inactivated "
                f"{t.inactivated_end:4.1f}")


def read(sp: SternParams, label: str, duration: float = 10.0,
         seeds: int = 4, trials: int | None = None) -> SpeedRow:
    """``sp`` in the cleft, spontaneous and triggered."""
    return SpeedRow(label, getattr(sp, "k_use_on", 0.0),
                    getattr(sp, "k_use_off", 0.0),
                    measure(sp, label, duration, seeds),
                    triggered(sp, label, trials))


def _fit(k_on: float, ratio: float, fraction: float | None = None
         ) -> UseParams:
    return fit_with_use(k_on, k_use_off=ratio * k_on, fraction=fraction)


def speed_panel(ratio: float | None = None, speeds=None,
                duration: float = 10.0, seeds: int = 4,
                trials: int | None = None) -> list[SpeedRow]:
    """The fitted Ca2+ gate alone, then the use gate at each entry rate with
    recovery ``ratio`` x entry (default Rios's 0.175) and the Ca2+ gate
    refitted beside it. A speed with no fit is skipped."""
    rho = rios_ratio() if ratio is None else float(ratio)
    rows = [read(fit_to_bell(), "fitted, Ca2+ gate only", duration, seeds,
                 trials)]
    for k in (speed_values() if speeds is None else speeds):
        try:
            sp = _fit(float(k), rho)
        except (RuntimeError, FloatingPointError):
            continue
        rows.append(read(sp, f"use {k:.3g}/s, rho {rho:.3g}", duration,
                         seeds, trials))
    return rows


def controls(duration: float = 10.0, seeds: int = 4,
             trials: int | None = None) -> list[SpeedRow]:
    """What the use gate ends by itself: each speed's fit with its Ca2+ gate
    removed, and Rios's gate beside the fitted activation gate with no Ca2+
    gate ever (his model; it leaves the bell no descending limb)."""
    on, off = rios_rates()
    rho = off / on
    slow = _fit(_P.value("ryr.k_use_on"), rho)
    fast = _fit(on, rho)
    return [read(no_ca_gate(slow), "bilayer fit, Ca2+ gate removed",
                 duration, seeds, trials),
            read(no_ca_gate(fast), "Rios-speed fit, Ca2+ gate removed",
                 duration, seeds, trials),
            read(no_ca_gate(with_use(fit_to_bell(), on, off)),
                 "Rios alone (no Ca2+ gate)", duration, seeds, trials),
            read(no_ca_gate(fit_to_bell()), "no inactivation", duration,
                 seeds, trials)]


def fraction_panel(fractions=None, duration: float = 10.0, seeds: int = 4,
                   trials: int | None = None) -> list[SpeedRow]:
    """Round 6.11's mixed cluster at Rios's ratio, bilayer speed against
    Rios's speed; the shared Ca2+ gate refitted to the population bell at
    each (default fractions: the registered one and 0.6)."""
    on, off = rios_rates()
    rho = off / on
    fr = ((_P.value("ryr.use_inactivating_fraction"), 0.6)
          if fractions is None else fractions)
    rows = []
    for f in fr:
        for name, k in (("bilayer", _P.value("ryr.k_use_on")), ("Rios", on)):
            try:
                sp = _fit(k, rho, fraction=float(f))
            except (RuntimeError, FloatingPointError):
                continue
            rows.append(read(sp, f"{name} speed, {f:.2f} carry", duration,
                             seeds, trials))
    return rows
