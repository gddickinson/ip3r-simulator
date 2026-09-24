"""The V channel: Rios et al. 1993's allosteric voltage-sensor model.

Rios, Karhanek, Ma & Gonzalez 1993 (J Gen Physiol 102:449) put the release
channel opposite a tetrad of voltage sensors (DHPRs) into the form of
Monod, Wyman & Changeux. The channel is closed (C) or open (O) as a whole,
and each of its four sensors is resting or activating, independently and
driven by voltage. That gives ten states, C0-C4 and O0-O4, j = sensors
activated (their Fig. 10):

* a sensor activates at ``k_c = alpha/2 exp((V - V_bar)/8K)`` and returns
  at ``k_-c = alpha/2 exp(-(V - V_bar)/8K)`` (Eqs. 2-3), times the
  degeneracy (4 - j, or j);
* the channel opens at ``k_L / f^j`` and closes at ``k_-L f^j``, so each
  activated sensor favours opening by ``f^-2``. On an open channel the
  sensor rates are ``k_c / f`` and ``f k_-c`` (microscopic reversibility).

The steady state is MWC's (their Eqs. 6-7, with x = exp((V - V_bar)/4K)
and L = k_-L/k_L):

    Po    = (1 + x/f^2)^4 / [(1 + x/f^2)^4 + L (1 + x)^4]
    Q/Qmax = [x/f^2 (1 + x/f^2)^3 + L x (1 + x)^3] / [same denominator]

:func:`stationary` computes both from the generator's null space instead,
and the tests hold the two to each other. The constants are fiber 827 of
their Table I. Stern et al. 1997 multiplied every time-dimensioned rate by
2 (``ec.stern_rate_scale``) for the V channels of their couplon. Their
couplon figure (Fig. 12) follows the factor; their stand-alone check of the
model (Fig. 11) follows the printed rates, so it was made before it. The model has no Ca2+ or Mg2+ in it: Stern's V channels are
"neither activated nor inactivated by Ca2+".
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import expm, null_space

from ..parameters import PARAMETERS as _P

__all__ = ["RiosParams", "N_SENSORS", "STATES", "OPEN_STATES", "generator",
           "stationary", "open_probability", "charge", "step_response",
           "exits"]

N_SENSORS = 4
#: C0..C4 then O0..O4; state ``j`` or ``5 + j`` has j sensors activated.
STATES = tuple(f"C{j}" for j in range(5)) + tuple(f"O{j}" for j in range(5))
OPEN_STATES = np.arange(5, 10)


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class RiosParams:
    k: float = _v("ec.rios_k")                    # mV (steepness 4K)
    v_bar: float = _v("ec.rios_v_bar")            # mV
    k_l: float = _v("ec.rios_k_l")                # 1/s
    k_minus_l: float = _v("ec.rios_k_minus_l")    # 1/s
    f: float = _v("ec.rios_f")
    alpha: float = _v("ec.rios_alpha")            # 1/s
    rate_scale: float = _v("ec.stern_rate_scale")

    @property
    def big_l(self) -> float:
        """L = k_-L / k_L, the closing equilibrium with no sensor active."""
        return self.k_minus_l / self.k_l

    def sensor_rates(self, v: float) -> tuple[float, float]:
        """(k_c, k_-c) of one sensor on a closed channel at ``v`` mV, 1/s."""
        a = 0.5 * self.alpha * self.rate_scale
        e = (v - self.v_bar) / (8.0 * self.k)
        return a * np.exp(e), a * np.exp(-e)


def generator(v: float, rp: RiosParams | None = None) -> np.ndarray:
    """Rate matrix (rows sum to zero, 1/s) of the ten states at ``v`` mV."""
    rp = rp or RiosParams()
    kc, kmc = rp.sensor_rates(v)
    s, f = rp.rate_scale, rp.f
    q = np.zeros((10, 10))
    for j in range(5):
        c, o = j, 5 + j
        if j < N_SENSORS:
            q[c, c + 1] = (N_SENSORS - j) * kc
            q[o, o + 1] = (N_SENSORS - j) * kc / f
        if j > 0:
            q[c, c - 1] = j * kmc
            q[o, o - 1] = j * kmc * f
        q[c, o] = s * rp.k_l / f ** j
        q[o, c] = s * rp.k_minus_l * f ** j
    q[np.diag_indices(10)] = -q.sum(axis=1)
    return q


def exits(v: float, rp: RiosParams | None = None
          ) -> tuple[np.ndarray, np.ndarray]:
    """Per state, the destinations and rates of its (up to three) exits,
    padded with rate 0: the table the couplon simulation draws from."""
    q = generator(v, rp)
    dest = np.zeros((10, 3), dtype=int)
    rate = np.zeros((10, 3))
    for i in range(10):
        js = [j for j in range(10) if j != i and q[i, j] > 0]
        dest[i, :len(js)] = js
        rate[i, :len(js)] = q[i, js]
    return dest, rate


def stationary(v: float, rp: RiosParams | None = None) -> np.ndarray:
    """Stationary occupancy of the ten states (null space of Q^T)."""
    p = null_space(generator(v, rp).T)[:, 0]
    return p / p.sum()


def _mwc(v, rp: RiosParams):
    x = np.exp((np.asarray(v, float) - rp.v_bar) / (4.0 * rp.k))
    y = x / rp.f ** 2
    o, c = (1.0 + y) ** 4, rp.big_l * (1.0 + x) ** 4
    return x, y, o, c


def open_probability(v, rp: RiosParams | None = None) -> np.ndarray:
    """Rios 1993 Eq. 6: steady open probability at ``v`` mV."""
    rp = rp or RiosParams()
    _, _, o, c = _mwc(v, rp)
    return o / (o + c)


def charge(v, rp: RiosParams | None = None) -> np.ndarray:
    """Rios 1993 Eq. 7: steady fraction of sensor charge moved."""
    rp = rp or RiosParams()
    x, y, o, c = _mwc(v, rp)
    num = y * (1.0 + y) ** 3 + rp.big_l * x * (1.0 + x) ** 3
    return num / (o + c)


def step_response(v: float, t, rp: RiosParams | None = None,
                  holding: float | None = None) -> np.ndarray:
    """Open probability at times ``t`` (s) after a step from the stationary
    state at ``holding`` (default ``ec.holding``) to ``v``: the master
    equation Stern et al. integrated to check their Monte Carlo (Fig. 11)."""
    rp = rp or RiosParams()
    holding = _P.value("ec.holding") if holding is None else holding
    p0, q = stationary(holding, rp), generator(v, rp)
    t = np.atleast_1d(np.asarray(t, float))
    return np.array([(p0 @ expm(q * ti))[OPEN_STATES].sum() for ti in t])
