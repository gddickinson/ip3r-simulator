"""Stochastic IP3R clusters: blips and puffs.

A cluster of ``N`` receptors, each of four De Young-Keizer subunits, each
subunit carrying three independent two-state sites (IP3, activating Ca2+,
inhibitory Ca2+). Every site flips as a Markov process with the DYK rates,
simulated on a fixed step (the scheme Shuai & Jung 2002 used for single
stochastic receptors):

    IP3 site          on  a1 p        off a1 d1      (d1: IP3 affinity)
    activating site   on  a5 c        off a5 d5
    inhibitory site   on  a2 c        off a2 Q2(p)   (Li-Rinzel reduction)

A subunit is active when IP3 and activating Ca2+ are bound and the inhibitory
site is empty; a channel is open when ``gating.subunits_required`` of its
subunits are active.

**Coupling.** Every channel in the cluster sees ``c = ca_rest + ca_per_open
x (number open)`` — a mean-field stand-in for the Ca2+ that open neighbours
deliver. That is the whole mechanism of a puff (Swillens et al. 1999): with
the coupling at zero, openings are independent *blips*; with it on, one
opening recruits others by Ca2+-induced Ca2+ release until the slow
inhibitory sites shut the cluster down. ``puff_compare.coupling_effect``
measures exactly that contrast.

**What it measures, and a limitation.** The signature of coupling is the
Fano factor (variance / mean) of the number of channels open at once: 1 for
independent channels (Poisson-like), above 1 when one opening recruits
others. With the DYK constants the activating site is already half occupied
at resting Ca2+ (``n_inf(0.1) = 0.55``), so the resting open probability is
high and puffs are modest (Fano ~1.4 at 0.2 µM IP3, against ~1.0 uncoupled).
The park/drive receptor (``puffs_pd``) is the contrast, and
``puff_compare`` measures both clusters with one ruler.

The simulation is seeded and reproducible: the same seed gives the same
trace.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from .gating import GatingParams, q2

__all__ = ["PuffParams", "PuffTrace", "simulate_cluster", "detect_events",
           "fano"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class PuffParams:
    n_channels: float = _v("puff.n_channels")
    ca_rest: float = _v("puff.ca_rest")
    ca_per_open: float = _v("puff.ca_per_open")
    dt: float = _v("puff.dt")
    record_dt: float = _v("puff.record_dt")
    a1: float = _v("gating.a1")
    a5: float = _v("gating.a5")
    gating: GatingParams = field(default_factory=GatingParams)


@dataclass
class PuffTrace:
    t: np.ndarray            # s
    n_open: np.ndarray       # open channels at the start of each bin
    ca: np.ndarray           # cluster Ca2+ seen by every channel, µM
    params: object           # the simulator's parameter set
    p: float
    n_peak: np.ndarray | None = None   # most open at once within each bin
    n_inactivated: np.ndarray | None = None  # RyR1 only: channels in CI or I

    @property
    def peaks(self) -> np.ndarray:
        return self.n_open if self.n_peak is None else self.n_peak


def _flip(state: np.ndarray, p_on: np.ndarray | float, p_off: np.ndarray | float,
          rng: np.random.Generator) -> np.ndarray:
    u = rng.random(state.shape)
    return np.where(state, u >= p_off, u < p_on)


def simulate_cluster(p: float, duration: float = 20.0, seed: int = 0,
                     pp: PuffParams | None = None) -> PuffTrace:
    """Simulate a cluster at IP3 ``p`` (µM) for ``duration`` seconds,
    recorded every ``record_dt``."""
    pp = pp or PuffParams()
    record_every = max(1, int(round(pp.record_dt / pp.dt)))
    g = pp.gating
    rng = np.random.default_rng(seed)
    n = int(round(pp.n_channels))
    need = int(round(g.subunits))
    dt = pp.dt
    shape = (n, 4)
    # Start from the resting equilibrium of each site at the resting Ca2+.
    c = pp.ca_rest
    ip3 = rng.random(shape) < p / (p + g.d1)
    act = rng.random(shape) < c / (c + g.d5)
    inh = rng.random(shape) < c / (c + q2(p, g))
    steps = int(round(duration / dt))
    n_rec = steps // record_every + 1
    t_out = np.empty(n_rec)
    open_out = np.empty(n_rec, dtype=np.int32)
    ca_out = np.empty(n_rec)
    peak_out = np.zeros(n_rec, dtype=np.int32)
    # Transition probabilities that do not depend on Ca2+ are fixed.
    ip3_on = pp.a1 * p * dt
    ip3_off = pp.a1 * g.d1 * dt
    act_off = pp.a5 * g.d5 * dt
    inh_off = g.a2 * float(q2(p, g)) * dt
    k = 0
    for step in range(steps + 1):
        active = ip3 & act & ~inh
        n_open = int(((active.sum(axis=1)) >= need).sum())
        c = pp.ca_rest + pp.ca_per_open * n_open
        if step % record_every == 0:
            t_out[k], open_out[k], ca_out[k] = step * dt, n_open, c
            k += 1
        peak_out[k - 1] = max(peak_out[k - 1], n_open)
        ip3 = _flip(ip3, ip3_on, ip3_off, rng)
        act = _flip(act, pp.a5 * c * dt, act_off, rng)
        inh = _flip(inh, g.a2 * c * dt, inh_off, rng)
    return PuffTrace(t_out[:k], open_out[:k], ca_out[:k], pp, p, peak_out[:k])


def detect_events(tr: PuffTrace, min_open: int = 1) -> list[dict]:
    """Contiguous runs of bins with at least ``min_open`` channels open at
    some moment in the bin.

    Each event reports start, duration and the peak number of channels open
    at once — a blip peaks at one, a puff at several.
    """
    peaks = tr.peaks
    on = peaks >= min_open
    edges = np.flatnonzero(np.diff(np.concatenate([[0], on.astype(int), [0]])))
    events = []
    for a, b in zip(edges[::2], edges[1::2]):
        events.append({"start": float(tr.t[a]),
                       "duration": float(tr.t[min(b, len(tr.t) - 1)] - tr.t[a]),
                       "peak_open": int(peaks[a:b].max())})
    return events


def fano(tr: PuffTrace) -> float:
    """Variance / mean of the number of channels open at once."""
    m = float(tr.n_open.mean())
    return float(tr.n_open.var() / m) if m > 0 else float("nan")
