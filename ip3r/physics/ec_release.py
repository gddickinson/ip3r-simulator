"""Release under voltage clamp: the couplon's ensemble, and what it says.

Stern et al. 1997 judged their couplon against whole-fibre release: a step
from ``ec.holding`` for ``ec.pulse`` gives a peak then a plateau, and
release stops when the membrane repolarises (their "control" property).
:func:`ensemble` averages ``ec.trials`` independent couplons, as they did,
and :func:`summarise` reads the numbers their Figs. 1, 2 and 12 show.

:func:`configurations` then asks the question Round 6.5 left: with the V
channels doing what they do in the fibre (opened by voltage, and free of
Mg2+, which is Stern's assumption and Laver 2018's lifted block), which C
channel scheme gives release that starts, stops at repolarisation, and
comes in events as brief as the measured 6.3 ms? The candidates are
Stern's as published, the scheme fitted to Murayama's bell, and that fit
under the fibre's 1 mM Mg2+ at the activation site alone or at both sites.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from ..parameters import PARAMETERS as _P
from .couplon import CouplonParams, CouplonTrace, simulate_couplon, step_protocol
from .puff_compare import spark_ends
from .puffs import PuffTrace, detect_events
from .ryr_gating import SternParams, fit_to_bell, with_mg
from .ryr_two_site import fit_two_site
from .ryr_use import fit_with_use, tau_values
from .spark_mg import k_mg_a_reading

__all__ = ["Ensemble", "Summary", "ensemble", "summarise", "voltages",
           "configurations", "with_gating", "as_puff_trace", "events", "EventStats",
           "event_stats"]


@dataclass(frozen=True)
class Ensemble:
    v: float                    # mV, the step
    t: np.ndarray               # s
    v_po: np.ndarray            # V-channel open probability
    c_po: np.ndarray            # C-channel open probability
    c_inactivated: np.ndarray   # fraction of C channels in CI or I
    flux: np.ndarray            # mean release current per couplon, pA
    trials: int


@dataclass(frozen=True)
class Summary:
    label: str
    v: float
    v_plateau: float
    c_peak: float
    c_peak_ms: float
    c_plateau: float
    c_after: float              # C open probability after repolarisation
    inactivated_end: float      # fraction at the end of the pulse
    peak_to_plateau: float      # of the total flux

    def row(self) -> str:
        return (f"{self.label:34s} {self.v:5.0f} mV  V {self.v_plateau:.4f}  C "
                f"peak {self.c_peak:.3f} at {self.c_peak_ms:4.0f} ms  plateau "
                f"{self.c_plateau:.3f}  after {self.c_after:.4f}  inact "
                f"{self.inactivated_end:.2f}  flux peak/plateau "
                f"{self.peak_to_plateau:5.2f}")


def voltages() -> np.ndarray:
    return np.arange(_P.value("ec.scan_v_min"),
                     _P.value("ec.scan_v_max") + 1e-9, _P.value("ec.scan_v_step"))


def ensemble(v: float, pp: CouplonParams | None = None,
             trials: int | None = None, seed0: int = 0,
             scale: tuple | None = None) -> Ensemble:
    """Mean of ``trials`` couplons stepped to ``v`` mV (then back); ``scale``
    is the SR-content path of :func:`couplon.simulate_couplon`."""
    pp = pp or CouplonParams()
    trials = int(round(_P.value("ec.trials"))) if trials is None else trials
    duration = _P.value("ec.pulse") + _P.value("ec.after")
    runs = [simulate_couplon(step_protocol(v), duration, seed0 + s, pp, scale)
            for s in range(trials)]
    n = min(len(r.t) for r in runs)
    mean = lambda f: np.mean([f(r)[:n] for r in runs], axis=0)
    n_c = int(round(pp.cleft.n_channels))
    return Ensemble(v, runs[0].t[:n], mean(lambda r: r.v_open) / n_c,
                    mean(lambda r: r.c_open) / n_c,
                    mean(lambda r: r.c_inactivated) / n_c,
                    mean(lambda r: r.flux()), trials)


def summarise(e: Ensemble, label: str = "") -> Summary:
    """Peak, plateau (the last ``ec.plateau_window`` of the pulse) and the
    C open probability over the last half of the time after it."""
    pulse, win = _P.value("ec.pulse"), _P.value("ec.plateau_window")
    during = e.t < pulse
    plateau = (e.t >= pulse - win) & during
    after = e.t >= pulse + 0.5 * _P.value("ec.after")
    i = int(np.argmax(np.where(during, e.c_po, -1.0)))
    flux_plateau = float(e.flux[plateau].mean())
    return Summary(label, e.v, float(e.v_po[plateau].mean()), float(e.c_po[i]),
                   1e3 * float(e.t[i]), float(e.c_po[plateau].mean()),
                   float(e.c_po[after].mean()) if after.any() else float("nan"),
                   float(e.c_inactivated[plateau][-1]),
                   float(e.flux[during].max() / flux_plateau)
                   if flux_plateau > 0 else float("inf"))


def configurations(reading: str = "meissner", mg: float | None = None,
                   two_site: bool = False, use: bool = False
                   ) -> dict[str, SternParams]:
    """The C-channel schemes compared (V channels are the same in all).
    ``two_site``: the fit is :func:`ryr_two_site.fit_two_site` (slope
    matched too) instead of the one-site :func:`fit_to_bell`.

    ``use``: the use-dependent gate of :mod:`ip3r.physics.ryr_use`, at every
    speed the measured bell permits, with the Ca2+ gate refitted beside it
    at each one so that every row still reproduces Murayama's half-peak
    points. Mg2+ plays no part in these rows (the gate is Ca2+-independent
    by construction), so they are compared against Stern's scheme and the
    Ca2+-gate-only fit rather than against the Mg2+ rows.
    """
    if use:
        if two_site:
            raise ValueError("--use and --two-site are different schemes for "
                             "the same gate; run them separately")
        rows = {"Stern 1997": SternParams(),
                "fitted, Ca2+ gate only": fit_to_bell()}
        for tau in tau_values():
            try:
                rows[f"fitted + use tau {tau:.3g} s"] = fit_with_use(1.0 / tau)
            except (RuntimeError, FloatingPointError):
                continue          # no Ca2+ gate fits beside this one
        return rows
    fit = fit_two_site() if two_site else fit_to_bell()
    name = "two-site" if two_site else "fitted"
    mg = _P.value("ryr.mg_free") if mg is None else mg
    k = k_mg_a_reading(fit, reading)
    return {"Stern 1997": SternParams(),
            f"{name}, no Mg2+": fit,
            f"{name}, Mg2+ activation site": with_mg(fit, mg, k_mg_a=k, mg_i=0.0),
            f"{name}, Mg2+ both sites": with_mg(fit, mg, k_mg_a=k)}


def with_gating(sp: SternParams, pp: CouplonParams | None = None
                ) -> CouplonParams:
    pp = pp or CouplonParams()
    return replace(pp, cleft=replace(pp.cleft, gating=sp))


def as_puff_trace(tr: CouplonTrace) -> PuffTrace:
    """The C channels as the puff ruler reads them."""
    return PuffTrace(tr.t, tr.c_open, np.zeros(len(tr.t)), tr.params.cleft,
                     0.0, tr.c_peak, tr.c_inactivated)


def events(tr: CouplonTrace) -> list[dict]:
    """C-channel events (any C channel open) that start during the pulse."""
    pulse = _P.value("ec.pulse")
    return [e for e in detect_events(as_puff_trace(tr)) if e["start"] < pulse]


@dataclass(frozen=True)
class EventStats:
    label: str
    v: float
    trials: int
    n_events: int
    median_ms: float            # all C events
    median_peak: float          # C channels open at once
    large: int                  # events reaching half the C channels
    large_median_ms: float
    unended: int                # still running when the record ends

    def row(self) -> str:
        return (f"{self.label:34s} {self.v:5.0f} mV  {self.n_events:4d} events "
                f"({self.n_events / self.trials:.2f}/couplon)  median "
                f"{self.median_ms:5.1f} ms, peak {self.median_peak:4.1f}  "
                f"large {self.large:3d} median {self.large_median_ms:6.1f} ms"
                f"  unended {self.unended}")


def event_stats(sp: SternParams, label: str, v: float | None = None,
                trials: int | None = None) -> EventStats:
    """Single C events in couplons stepped to ``v`` (default
    ``ec.spark_voltage``, where V openings are sparse)."""
    v = _P.value("ec.spark_voltage") if v is None else v
    trials = int(round(_P.value("ec.trials"))) if trials is None else trials
    pp = with_gating(sp)
    duration = _P.value("ec.pulse") + _P.value("ec.after")
    ev, large = [], []
    for s in range(trials):
        tr = simulate_couplon(step_protocol(v), duration, s, pp)
        ev += events(tr)
        pulse = _P.value("ec.pulse")
        large += [e for e in spark_ends(as_puff_trace(tr)) if e["start"] < pulse]
    med = lambda xs: float(np.median(xs)) if xs else float("nan")
    return EventStats(label, v, trials, len(ev),
                      1e3 * med([e["duration"] for e in ev]),
                      med([e["peak_open"] for e in ev]), len(large),
                      1e3 * med([e["duration"] for e in large]),
                      sum(e["unterminated"] for e in large))
