"""Whole-cell Ca2+ dynamics: the closed-cell Li-Rinzel model.

Two variables — cytosolic Ca2+ ``c`` and the inhibition gate ``h`` — and
three fluxes between cytosol and ER (De Young & Keizer 1992; Li & Rinzel
1994). The cell is closed, so total Ca2+ is conserved and the ER
concentration follows from ``c``:

    c_ER   = (c0 - c) / c1
    J_chan = c1 v1 (m_inf n_inf h)^3 (c_ER - c)     IP3R release
    J_leak = c1 v2 (c_ER - c)                        passive leak
    J_pump = v3 c^2 / (k3^2 + c^2)                   SERCA uptake
    dc/dt  = J_chan + J_leak - J_pump
    dh/dt  = a2 (Q2 (1 - h) - c h)

What it demonstrates: at low IP3 the cell rests; in an intermediate window
Ca2+ **oscillates** without any oscillating input, because fast Ca2+-induced
Ca2+ release (activation) is followed by slow Ca2+ inhibition (``h``); at high
IP3 it settles again, high. :func:`oscillation_window` finds the window's
edges by simulation, which is the measured version of the model's two Hopf
bifurcations.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from ..parameters import PARAMETERS as _P
from .gating import GatingParams, h_inf, m_inf, n_inf, q2

__all__ = ["CellParams", "fluxes", "simulate", "Trace", "oscillation_metrics",
           "oscillation_window", "steady_state"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class CellParams:
    v1: float = _v("cell.v1")
    v2: float = _v("cell.v2")
    v3: float = _v("cell.v3")
    k3: float = _v("cell.k3")
    c0: float = _v("cell.c0")
    c1: float = _v("cell.c1")
    gating: GatingParams = field(default_factory=GatingParams)


def fluxes(c, h, p, cp: CellParams | None = None) -> dict[str, np.ndarray]:
    cp = cp or CellParams()
    g = cp.gating
    c = np.asarray(c, float)
    c_er = (cp.c0 - c) / cp.c1
    po = (m_inf(p, g) * n_inf(c, g) * h) ** g.subunits
    return {"channel": cp.c1 * cp.v1 * po * (c_er - c),
            "leak": cp.c1 * cp.v2 * (c_er - c),
            "pump": cp.v3 * c ** 2 / (cp.k3 ** 2 + c ** 2),
            "c_er": c_er, "po": po}


def _rhs(t, y, p, cp):
    c, h = y
    f = fluxes(c, h, p, cp)
    dc = f["channel"] + f["leak"] - f["pump"]
    dh = cp.gating.a2 * (q2(p, cp.gating) * (1.0 - h) - c * h)
    return [dc, dh]


@dataclass
class Trace:
    t: np.ndarray
    c: np.ndarray
    h: np.ndarray
    p: float
    c_er: np.ndarray
    po: np.ndarray


def simulate(p: float, t_end: float = 100.0, c_init: float = 0.1,
             h_init: float | None = None, dt_out: float = 0.05,
             cp: CellParams | None = None, p_of_t=None) -> Trace:
    """Integrate the model at IP3 ``p`` (µM), or ``p_of_t(t)`` if given."""
    cp = cp or CellParams()
    h0 = float(h_inf(c_init, p, cp.gating)) if h_init is None else h_init
    if p_of_t is None:
        fun = lambda t, y: _rhs(t, y, p, cp)  # noqa: E731
    else:
        fun = lambda t, y: _rhs(t, y, p_of_t(t), cp)  # noqa: E731
    t_eval = np.arange(0.0, t_end + 1e-12, dt_out)
    sol = solve_ivp(fun, (0.0, t_end), [c_init, h0], method="LSODA",
                    t_eval=t_eval, rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    ps = np.full_like(sol.t, p) if p_of_t is None else np.array([p_of_t(t) for t in sol.t])
    f = fluxes(sol.y[0], sol.y[1], ps, cp)
    return Trace(sol.t, sol.y[0], sol.y[1], p, f["c_er"], f["po"])


def oscillation_metrics(tr: Trace, discard: float = 0.5) -> dict:
    """Amplitude and period over the last ``1 - discard`` of a trace.

    Period is the mean interval between upward crossings of the midline;
    amplitude is peak-to-trough. A trace is *oscillating* only if its
    amplitude exceeds 1 nM **and is not decaying**: the second half of the
    window must keep at least 90 % of the first half's amplitude. Without the
    second clause a slowly damped spiral just past a Hopf point counts as an
    oscillation, and the first version put the upper edge of the window at
    0.73 µM for exactly that reason.
    """
    keep = tr.t >= tr.t[-1] * discard
    c, t = tr.c[keep], tr.t[keep]
    amp = float(c.max() - c.min())
    half = len(c) // 2
    amp1 = float(c[:half].max() - c[:half].min()) if half else 0.0
    amp2 = float(c[half:].max() - c[half:].min()) if half else 0.0
    if amp < 1e-3 or amp2 < 0.9 * amp1:
        return {"oscillates": False, "amplitude": amp, "period": float("nan"),
                "mean": float(c.mean())}
    mid = 0.5 * (c.max() + c.min())
    up = np.flatnonzero((c[:-1] < mid) & (c[1:] >= mid))
    period = float(np.diff(t[up]).mean()) if len(up) >= 2 else float("nan")
    return {"oscillates": len(up) >= 2, "amplitude": amp, "period": period,
            "mean": float(c.mean())}


def steady_state(p: float, cp: CellParams | None = None) -> float:
    """Cytosolic Ca2+ at the fixed point (µM), found by bisection on dc/dt."""
    cp = cp or CellParams()
    lo, hi = 1e-6, cp.c0 * 0.999
    def f(c):
        return (_rhs(0, [c, float(h_inf(c, p, cp.gating))], p, cp))[0]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def oscillation_window(p_grid: np.ndarray | None = None, t_end: float = 400.0,
                       cp: CellParams | None = None) -> tuple[float, float, list]:
    """IP3 range (µM) over which the model oscillates, found by simulation.

    Returns ``(p_low, p_high, rows)``; the edges are the first and last grid
    points that oscillate, so their precision is the grid spacing.
    """
    p_grid = np.round(np.arange(0.2, 1.0001, 0.01), 4) if p_grid is None else p_grid
    rows = []
    for p in p_grid:
        m = oscillation_metrics(simulate(float(p), t_end=t_end, cp=cp))
        rows.append({"p": float(p), **m})
    osc = [r["p"] for r in rows if r["oscillates"]]
    if not osc:
        return float("nan"), float("nan"), rows
    return min(osc), max(osc), rows
