"""The park/drive IP3 receptor: Siekmann et al. (2012), with the gating
variables of Cao et al. (2013).

Six states in two modes, fitted by MCMC to stationary single-channel records
of IP3R-1:

    drive:  C1 -- C2 -- C3        park:  C4 -- O5
                  |                       |
                  O6                      (C4 -- C2 joins the modes)

Inside each mode the rates are constant. A drive-mode channel is open about
70 % of the time (the fast C2 <-> O6 pair); a park-mode channel almost never
is (q54 / q45 ~ 300). Ca2+ and IP3 act *only* on the switch between modes:

    q24 = a24 + V24 (1 - m24 h24)     drive -> park
    q42 = a42 + V42 m42 h42           park  -> drive

Each gating variable ``G`` relaxes towards a Ca2+-dependent equilibrium,
``dG/dt = lambda_G (G_inf(c) - G)``. The ``a``, ``V`` and half-constants
depend on IP3. These are heuristic fits, "with no biophysical basis" (Cao et
al.). ``h42`` recovers slowly while the channel is closed and falls fast
while it is open. That asymmetry is what ends a puff and spaces the next.

Why it is here: at resting Ca2+ nearly every receptor is parked, so a
cluster is quiet until one receptor enters drive mode. The De Young-Keizer
subunit scheme (``puffs``) has no such mode.

The stationary functions below take the gating variables at equilibrium at a
clamped Ca2+ (the single-channel condition). Concentrations are µM, rates
per second.
"""

from __future__ import annotations

from dataclasses import field, make_dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from .bell import Bell, measure_bell

__all__ = ["ParkDriveParams", "STATES", "OPEN", "PARK", "ip3_functions",
           "gate_inf", "mode_rates", "constant_generator", "stationary",
           "open_probability", "park_fraction", "drive_open_probability",
           "bell_at"]

# Index order C1 C2 C3 C4 O5 O6, as in the source code.
STATES = ("C1", "C2", "C3", "C4", "O5", "O6")
OPEN = np.array([False, False, False, False, True, True])
PARK = np.array([False, False, False, True, True, False])
C2, C4 = 1, 3


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


_KEYS = ["q12", "q21", "q23", "q32", "q26", "q62", "q45", "q54",
         "v42_amp", "v42_k", "v42_h",
         "k42_base", "k42_amp", "k42_k", "k42_h", "n42",
         "kn42_base", "kn42_amp", "kn42_k", "kn42_h", "nn42",
         "a42_amp", "a42_k", "a42_h",
         "v24_base", "v24_amp", "v24_k", "v24_h",
         "k24", "n24", "kn24", "nn24",
         "a24_base", "a24_amp", "a24_k", "a24_h",
         "lam_m24", "lam_h24", "lam_m42", "lam_h42_closed", "lam_h42_open"]

#: Every constant of the receptor, read from ``pd.<name>`` at construction.
ParkDriveParams = make_dataclass(
    "ParkDriveParams", [(k, float, _v(f"pd.{k}")) for k in _KEYS])
ParkDriveParams.__module__ = __name__     # so worker processes can pickle it


def _rise(p, k, h):
    p = np.asarray(p, float)
    return p ** h / (p ** h + k ** h)


def ip3_functions(p: float, g: ParkDriveParams | None = None) -> dict:
    """The IP3-dependent constants of the mode switch at IP3 ``p``."""
    g = g or ParkDriveParams()
    return {
        "V42": g.v42_amp * _rise(p, g.v42_k, g.v42_h),
        "k42": g.k42_base + g.k42_amp * _rise(p, g.k42_k, g.k42_h),
        "kn42": g.kn42_base + g.kn42_amp * _rise(p, g.kn42_k, g.kn42_h),
        "a42": g.a42_amp * _rise(p, g.a42_k, g.a42_h),
        "V24": g.v24_base + g.v24_amp * (1 - _rise(p, g.v24_k, g.v24_h)),
        "a24": g.a24_base + g.a24_amp * (1 - _rise(p, g.a24_k, g.a24_h)),
    }


def gate_inf(c, f: dict, g: ParkDriveParams | None = None):
    """Equilibria ``(m24, h24, m42, h42)`` at Ca2+ ``c`` (Cao 2013 Eqs. 4-7)."""
    g = g or ParkDriveParams()
    c = np.asarray(c, float)
    m24 = _rise(c, g.k24, g.n24)
    h24 = 1 - _rise(c, g.kn24, g.nn24)
    m42 = _rise(c, f["k42"], g.n42)
    h42 = 1 - _rise(c, f["kn42"], g.nn42)
    return m24, h24, m42, h42


def mode_rates(gates, f: dict):
    """``(q24, q42)`` from the gating variables (Cao 2013 Eqs. 1-2)."""
    m24, h24, m42, h42 = gates
    return f["a24"] + f["V24"] * (1 - m24 * h24), f["a42"] + f["V42"] * m42 * h42


def constant_generator(g: ParkDriveParams | None = None) -> np.ndarray:
    """The 6x6 generator of the constant-rate transitions, without the
    mode switch (rows sum to zero)."""
    g = g or ParkDriveParams()
    q = np.zeros((6, 6))
    for i, j, r in ((0, 1, g.q12), (1, 0, g.q21), (1, 2, g.q23), (2, 1, g.q32),
                    (1, 5, g.q26), (5, 1, g.q62), (3, 4, g.q45), (4, 3, g.q54)):
        q[i, j] = r
    return q - np.diag(q.sum(axis=1))


def stationary(c: float, p: float, g: ParkDriveParams | None = None) -> np.ndarray:
    """Stationary occupancy of the six states at clamped ``c`` and ``p``.

    The scheme is a tree, so detailed balance gives it directly: every state
    relative to C2.
    """
    g = g or ParkDriveParams()
    f = ip3_functions(p, g)
    q24, q42 = mode_rates(gate_inf(c, f, g), f)
    with np.errstate(divide="ignore", invalid="ignore"):
        c4 = np.where(q42 > 0, q24 / q42, np.inf)
    w = np.array([g.q21 / g.q12, 1.0, g.q23 / g.q32, c4,
                  c4 * g.q45 / g.q54, g.q26 / g.q62], dtype=float)
    if not np.isfinite(c4):          # never leaves park: all weight there
        park = np.array([0, 0, 0, 1.0, g.q45 / g.q54, 0])
        return park / park.sum()
    return w / w.sum()


def open_probability(c, p: float, g: ParkDriveParams | None = None):
    """Stationary P_open at clamped Ca2+ ``c`` (scalar or array)."""
    c = np.asarray(c, float)
    out = [stationary(ci, p, g)[OPEN].sum() for ci in c.ravel()]
    return np.array(out).reshape(c.shape)


def park_fraction(c: float, p: float, g: ParkDriveParams | None = None) -> float:
    return float(stationary(c, p, g)[PARK].sum())


def drive_open_probability(g: ParkDriveParams | None = None) -> float:
    """P_open inside drive mode alone (the source says 'around 70 %')."""
    g = g or ParkDriveParams()
    w = np.array([g.q21 / g.q12, 1.0, g.q23 / g.q32, g.q26 / g.q62])
    return float(w[3] / w.sum())


def bell_at(p: float, g: ParkDriveParams | None = None) -> Bell:
    return measure_bell(lambda c: open_probability(c, p, g))

