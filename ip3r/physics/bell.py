"""Model-agnostic measurements of a bell-shaped P_open(Ca2+) curve.

Both gating schemes (``gating`` — De Young-Keizer; ``gating_mak`` — Mak et
al. 1998) are measured here the same way, so a statement such as "IP3 moves
the inhibitory flank and not the activating one" is tested with one ruler.

A curve is any ``f(c) -> P_open`` over Ca2+ in µM. The flanks are the Ca2+
at which the curve crosses half its own peak, on either side of it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import brentq, minimize_scalar

__all__ = ["Bell", "measure_bell", "flank_shifts"]

_LOG_C = (-4.0, 4.0)          # search range, log10 µM — wider than any bell


@dataclass(frozen=True)
class Bell:
    c_peak: float             # µM
    po_peak: float
    c_half_act: float         # rising flank, µM
    c_half_inh: float         # falling flank, µM (nan if it never falls)

    @property
    def width_decades(self) -> float:
        return float(np.log10(self.c_half_inh / self.c_half_act))


def measure_bell(f: Callable) -> Bell:
    """Peak and half-peak crossings of ``f`` (log-Ca2+ search)."""
    g = lambda lc: float(f(10.0 ** lc))
    grid = np.linspace(*_LOG_C, 801)
    vals = np.array([g(x) for x in grid])
    i = int(np.argmax(vals))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    res = minimize_scalar(lambda x: -g(x), bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-8})
    x_pk, pk = float(res.x), float(-res.fun)
    half = lambda x: g(x) - pk / 2
    x_act = brentq(half, _LOG_C[0], x_pk, xtol=1e-10)
    x_inh = (brentq(half, x_pk, _LOG_C[1], xtol=1e-10)
             if half(_LOG_C[1]) < 0 else np.nan)
    return Bell(10 ** x_pk, pk, 10 ** x_act, 10 ** x_inh)


def flank_shifts(bell_at: Callable[[float], Bell], p_lo: float, p_hi: float
                 ) -> tuple[float, float]:
    """Fold change of (activating, inhibitory) half-points from IP3 ``p_lo``
    to ``p_hi``. ``bell_at(p)`` returns the model's ``Bell`` at IP3 ``p``."""
    a, b = bell_at(p_lo), bell_at(p_hi)
    return b.c_half_act / a.c_half_act, b.c_half_inh / a.c_half_inh
