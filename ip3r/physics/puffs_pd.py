"""A stochastic cluster of park/drive receptors (``park_drive``).

The cluster is the one in ``puffs``, with the De Young-Keizer receptor
replaced by the park/drive one and nothing else changed. Every closed
receptor sees the mean-field cluster Ca2+ ``c = ca_rest + ca_per_open x
(number open)``. An open receptor also sees its own mouth, ``c + ca_mouth``,
which is Cao et al.'s two-concentration scheme.

**Integration.** Drive-mode flicker (q26 = 10,500 s^-1) is far faster than
any Ca2+-dependent rate, so each step of ``dt`` is split in two:

1. the constant-rate transitions inside both modes are taken *exactly*,
   each receptor drawing its next state from a row of ``expm(Q dt)``;
2. the mode switch C2 -> C4 (``q24``) and C4 -> C2 (``q42``) fires with
   probability ``1 - exp(-q dt)`` at the step's gating variables, and each
   gating variable relaxes exactly (exponentially) toward its equilibrium at
   the Ca2+ the receptor saw during the step.

The splitting error is of order ``q24 x dt`` (a few per cent at the 0.1 ms
step Cao et al. used as their maximum). ``tests/test_puffs_pd.py`` checks
it: an uncoupled receptor at clamped Ca2+ must reproduce the stationary
occupancy of ``park_drive.stationary``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import expm

from ..parameters import PARAMETERS as _P
from . import park_drive as pdm
from .puffs import PuffTrace

__all__ = ["ParkDrivePuffParams", "simulate_cluster_pd"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class ParkDrivePuffParams:
    n_channels: float = _v("puff.n_channels")
    ca_rest: float = _v("puff.ca_rest")
    ca_per_open: float = _v("puff.pd_ca_per_open")
    ca_mouth: float = _v("pd.ca_mouth")
    dt: float = _v("puff.pd_dt")
    record_dt: float = _v("puff.record_dt")
    receptor: pdm.ParkDriveParams = field(default_factory=pdm.ParkDriveParams)


def simulate_cluster_pd(p: float, duration: float = 20.0, seed: int = 0,
                        pp: ParkDrivePuffParams | None = None,
                        clamp_ca: float | None = None) -> PuffTrace:
    """Simulate a park/drive cluster at IP3 ``p`` (µM) for ``duration`` s.

    Every receptor starts parked (C4), with its gates at equilibrium at
    resting Ca2+. With ``clamp_ca`` set, every receptor sees that Ca2+
    whatever the cluster does, open or closed. That is the single-channel
    condition the stationary functions describe, and the test uses it.
    """
    pp = pp or ParkDrivePuffParams()
    g = pp.receptor
    rng = np.random.default_rng(seed)
    n = int(round(pp.n_channels))
    dt = pp.dt
    f = pdm.ip3_functions(p, g)
    cum = np.cumsum(expm(pdm.constant_generator(g) * dt), axis=1)
    cum[:, -1] = 1.0                                  # guard round-off
    open_ = pdm.OPEN
    # The cluster Ca2+ takes one of n + 1 values, so the gate equilibria are
    # tabulated once: row k = (closed, open) receptor with k channels open.
    n_open_all = np.arange(n + 1)
    cluster = pp.ca_rest + pp.ca_per_open * n_open_all
    if clamp_ca is not None:
        cluster = np.full(n + 1, float(clamp_ca))
    closed_inf = np.stack(pdm.gate_inf(cluster, f, g), axis=1)          # (n+1, 4)
    mouth = 0.0 if clamp_ca is not None else pp.ca_mouth
    open_inf = np.stack(pdm.gate_inf(cluster + mouth, f, g), axis=1)
    # Exact exponential relaxation per step; only h42 depends on open/closed.
    lam_closed = np.array([g.lam_m24, g.lam_h24, g.lam_m42, g.lam_h42_closed])
    lam_open = np.array([g.lam_m24, g.lam_h24, g.lam_m42, g.lam_h42_open])
    decay_closed = np.exp(-lam_closed * dt)[:, None]
    decay_open = np.exp(-lam_open * dt)[:, None]
    state = np.full(n, pdm.C4)
    gates = np.repeat(closed_inf[0][:, None], n, axis=1)                # (4, n)

    record_every = max(1, int(round(pp.record_dt / dt)))
    steps = int(round(duration / dt))
    n_rec = steps // record_every + 1
    t_out = np.empty(n_rec)
    open_out = np.empty(n_rec, dtype=np.int32)
    peak_out = np.zeros(n_rec, dtype=np.int32)
    ca_out = np.empty(n_rec)
    k = 0
    for step in range(steps + 1):
        is_open = open_[state]
        n_open = int(is_open.sum())
        if step % record_every == 0:
            t_out[k], open_out[k], ca_out[k] = step * dt, n_open, cluster[n_open]
            k += 1
        peak_out[k - 1] = max(peak_out[k - 1], n_open)
        # Gates relax toward the Ca2+ each receptor sees over this step: the
        # cluster's, or the cluster's plus its own mouth while open.
        target = np.where(is_open, open_inf[n_open][:, None], closed_inf[n_open][:, None])
        gates = target + (gates - target) * np.where(is_open, decay_open, decay_closed)
        q24, q42 = pdm.mode_rates(gates, f)
        # 1. exact constant-rate step inside each mode
        state = (rng.random(n)[:, None] > cum[state]).sum(axis=1)
        # 2. mode switch
        u = rng.random(n)
        to_park = (state == pdm.C2) & (u < -np.expm1(-q24 * dt))
        to_drive = (state == pdm.C4) & (u < -np.expm1(-q42 * dt))
        state = np.where(to_park, pdm.C4, np.where(to_drive, pdm.C2, state))
    return PuffTrace(t_out[:k], open_out[:k], ca_out[:k], pp, p, peak_out[:k])
