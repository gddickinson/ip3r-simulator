"""Stochastic RyR1 clusters: Ca2+ sparks, read with the puff ruler.

A cluster of ``spark.n_channels`` RyR1s, each the four-state Stern 1997
scheme (:mod:`ip3r.physics.ryr_gating`). Every channel sees one cluster
Ca2+, ``c = ca_rest + ca_per_open x (number open)``, the same mean-field
coupling the IP3R clusters use (:mod:`ip3r.physics.puffs`). So a spark and
a puff differ only in the receptor, and :mod:`ip3r.physics.puff_compare`
measures both with one ruler.

**The coupling is derived, not chosen.** It is the steady free-diffusion
Ca2+ at one channel spacing from a point source carrying one channel's
current: ``c = i / (2 F 4 pi D r)``. All three inputs are Stern's Table I
(``spark.unitary_current`` 0.3 pA, ``spark.d_ca`` 5e-6 cm^2/s,
``spark.channel_spacing`` 30 nm), which gives about 8 µM per open
neighbour. No buffers are included, so it is an upper estimate. The spark
scan sweeps a band around it.

**The step is exact within itself.** The cluster Ca2+ depends only on how
many channels are open, so the transition matrix ``expm(Q(c) dt)`` is
tabulated once per possible number open, and each step draws every
channel's next state from its row. This is the split step of
``puffs_pd``. Rates up to ``k_o c^2`` of several 10^4 s^-1 are then handled
without a vanishing step. ``spark.dt`` only sets how often the cluster
Ca2+ is updated.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import expm

from ..parameters import PARAMETERS as _P
from .puffs import PuffTrace
from .ryr_gating import OPEN, SternParams, generator, stationary

__all__ = ["SparkParams", "diffusion_coupling", "simulate_sparks",
           "spark_couplings"]

F_FARADAY = 96485.33212          # C/mol
N_A = 6.02214076e23


def diffusion_coupling() -> float:
    """µM at one channel spacing from one open channel (free diffusion)."""
    i = _P.value("spark.unitary_current") * 1e-12             # A
    d = _P.value("spark.d_ca") * 1e-4                          # cm^2/s -> m^2/s
    r = _P.value("spark.channel_spacing") * 1e-9               # m
    flux = i / (2.0 * F_FARADAY)                               # mol/s
    return flux / (4.0 * np.pi * d * r) * 1e3                  # mol/m^3 -> µM


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class SparkParams:
    n_channels: float = _v("spark.n_channels")
    ca_rest: float = _v("puff.ca_rest")
    ca_per_open: float = field(default_factory=diffusion_coupling)
    dt: float = _v("spark.dt")
    record_dt: float = _v("puff.record_dt")
    gating: SternParams = field(default_factory=SternParams)


def _tables(pp: SparkParams, n: int) -> np.ndarray:
    """Cumulative transition rows per number open: shape (n+1, 4, 4)."""
    out = np.empty((n + 1, 4, 4))
    for k in range(n + 1):
        p = expm(generator(pp.ca_rest + pp.ca_per_open * k, pp.gating) * pp.dt)
        p = np.clip(p, 0.0, None)
        out[k] = np.cumsum(p / p.sum(axis=1, keepdims=True), axis=1)
    return out


def simulate_sparks(p: float = 0.0, duration: float = 5.0, seed: int = 0,
                    pp: SparkParams | None = None) -> PuffTrace:
    """Simulate a RyR1 cluster for ``duration`` s. ``p`` (IP3) is ignored;
    it is accepted so sparks and puffs share one call signature."""
    pp = pp or SparkParams()
    rng = np.random.default_rng(seed)
    n = int(round(pp.n_channels))
    cum = _tables(pp, n)
    # Start from each channel's stationary occupancy at resting Ca2+.
    pi = np.cumsum(stationary(pp.ca_rest, pp.gating))
    state = np.searchsorted(pi, rng.random(n) * pi[-1])
    record_every = max(1, int(round(pp.record_dt / pp.dt)))
    steps = int(round(duration / pp.dt))
    n_rec = steps // record_every + 1
    t_out, ca_out = np.empty(n_rec), np.empty(n_rec)
    open_out = np.empty(n_rec, dtype=np.int32)
    peak_out = np.zeros(n_rec, dtype=np.int32)
    inact_out = np.zeros(n_rec, dtype=np.int32)
    k = 0
    for step in range(steps + 1):
        n_open = int((state == OPEN).sum())
        if step % record_every == 0:
            t_out[k], open_out[k] = step * pp.dt, n_open
            ca_out[k] = pp.ca_rest + pp.ca_per_open * n_open
            inact_out[k] = int((state >= 2).sum())
            k += 1
        peak_out[k - 1] = max(peak_out[k - 1], n_open)
        rows = cum[n_open][state]                              # (n, 4)
        state = (rows < rng.random(n)[:, None]).sum(axis=1)
        state = np.minimum(state, 3)
    return PuffTrace(t_out[:k], open_out[:k], ca_out[:k], pp, p, peak_out[:k],
                     inact_out[:k])


def spark_couplings(base: float | None = None) -> np.ndarray:
    """0, then ``puff.scan_points`` geometric steps over the band
    ``spark.scan_low``-``spark.scan_high`` times ``base`` (default: the
    derived coupling)."""
    base = diffusion_coupling() if base is None else base
    return np.concatenate([[0.0], base * np.geomspace(
        _P.value("spark.scan_low"), _P.value("spark.scan_high"),
        int(round(_P.value("puff.scan_points"))))])
