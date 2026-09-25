"""A park/drive cluster in a microdomain (Round 7.2).

The receptors are ``puffs_pd.ReceptorCluster``, unchanged. What changes is
their Ca2+. In ``puffs_pd`` a closed receptor sees ``ca_rest + ca_per_open x
N_open`` at once. Here it sees the microdomain ``cb`` of
``microdomain.rhs``, which the open receptors fill and diffusion, pumps and
fluo-4 drain. That is Cao et al. 2013's point-domain model, with Cao et al.
2014's pools. An open receptor sees its own mouth,
``mouth_per_store x cs``.

Each step of ``puff.pd_dt`` integrates the pools by RK4 with the number open
at the step's start (Cao et al. do the same between channel events). The
receptors then step toward the gate equilibria at the new ``cb`` and mouth.

The trace records, per ``puff.record_dt`` bin, the number open (snapshot and
peak), the microdomain, cytosol and store Ca2+, and F/F0 of the bound
indicator. ``puff_stats`` reads puffs from F/F0, as the experiments do.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from . import microdomain as md
from . import park_drive as pdm
from .puffs import PuffTrace
from .puffs_pd import ReceptorCluster

__all__ = ["DomainPuffParams", "DomainTrace", "simulate_cluster_domain"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class DomainPuffParams:
    n_channels: float = _v("puff.n_channels")
    dt: float = _v("puff.pd_dt")
    record_dt: float = _v("puff.record_dt")
    clamp: str = "store"               # microdomain.CLAMPS
    domain: md.DomainParams = field(default_factory=md.DomainParams)
    receptor: pdm.ParkDriveParams = field(default_factory=pdm.ParkDriveParams)


@dataclass
class DomainTrace(PuffTrace):
    """``PuffTrace`` (``ca`` = microdomain) plus the other pools and dye."""
    cytosol: np.ndarray | None = None
    store: np.ndarray | None = None
    f_ratio: np.ndarray | None = None      # F/F0, snapshot per bin
    f_peak: np.ndarray | None = None       # F/F0, highest within the bin
    rest: md.Rest | None = None


def simulate_cluster_domain(p: float, duration: float = 20.0, seed: int = 0,
                            dp: DomainPuffParams | None = None) -> DomainTrace:
    """Simulate the cluster at IP3 ``p`` (µM) for ``duration`` s from rest.

    ``dp.clamp`` is ``"store"`` by default (Cao 2013's assumption: the store
    held at its starting value); ``"none"`` lets it deplete and ``"bath"``
    also holds the cytosol (``microdomain.CLAMPS``).
    """
    dp = dp or DomainPuffParams()
    d = dp.domain
    rest = md.rest_state(p, d)
    if dp.clamp != "none":
        d = md.clamped(d, rest, dp.clamp)
    rng = np.random.default_rng(seed)
    n = int(round(dp.n_channels))
    dt = dp.dt
    rc = ReceptorCluster(n, p, dt, rng, dp.receptor, rest.c)
    y = tuple(float(v) for v in rest.y)
    b0 = rest.b if rest.b > 0 else np.nan

    record_every = max(1, int(round(dp.record_dt / dt)))
    steps = int(round(duration / dt))
    n_rec = steps // record_every + 1
    t_out = np.empty(n_rec)
    open_out = np.empty(n_rec, dtype=np.int32)
    peak_out = np.zeros(n_rec, dtype=np.int32)
    cols = np.empty((4, n_rec))          # cb, c, cs, F/F0
    f_peak = np.zeros(n_rec)
    k = 0
    for step in range(steps + 1):
        n_open = int(rc.is_open.sum())
        cs = md.store(y, d)
        if step % record_every == 0:
            t_out[k], open_out[k] = step * dt, n_open
            cols[:, k] = y[1], y[0], cs, y[3] / b0
            k += 1
        peak_out[k - 1] = max(peak_out[k - 1], n_open)
        f_peak[k - 1] = max(f_peak[k - 1], y[3] / b0)
        y = md.rk4_step(y, n_open, p, d, dt)
        cs = md.store(y, d)
        rc.step(rc.equilibrium(y[1]), rc.equilibrium(d.mouth_per_store * cs))
    cb, c, cs, f = cols[:, :k]
    return DomainTrace(t_out[:k], open_out[:k], cb, dp, p, peak_out[:k],
                       cytosol=c, store=cs, f_ratio=f, f_peak=f_peak[:k],
                       rest=rest)
