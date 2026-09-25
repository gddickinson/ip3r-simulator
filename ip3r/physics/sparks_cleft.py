"""RyR1 sparks in the junctional cleft: each channel sees its own Ca2+.

The mean-field cluster (:mod:`ip3r.physics.sparks`) gives every channel one
cluster Ca2+ from the number open. Here each C channel ``i`` sees

    c_i = ca_rest + sum_j G[i, j] open_j

with ``G`` the steady cleft coupling of :mod:`ip3r.physics.cleft` (Stern
1997's geometry). Near neighbours matter and distant ones hardly at all,
and an open channel sees its own release (the diagonal).

**Exact, with no step.** Between gating events the field is constant (the
steady-state approximation), so the array is a true Markov process. It is
simulated the way Stern et al. did (their Appendix): Gillespie's method,
with every channel's two exit rates evaluated at its own Ca2+. There is no
``dt`` to converge.

**The coupling knob.** ``ca_per_open`` is the coupling between the two
nearest C channels, µM; the scan and the GUI's "Ca²⁺ coupling" set it. The
whole off-diagonal of ``G`` is scaled to match, so the shape of the field
is kept. The diagonal is a channel's own release and is never scaled, so
"uncoupled" (0) leaves independent channels that still see themselves.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from .cleft import CleftGeometry, coupling_matrix, nearest_coupling
from .puffs import PuffTrace
from .ryr_gating import SternParams, stationary

__all__ = ["CleftSparkParams", "native_coupling", "couplings_for",
           "simulate_sparks_cleft", "DEST", "c_rates", "initial_states"]

#: Stern's scheme's exits (C, O, CI, I): activation gate, inactivation gate.
#: A simulation reads ``sp.dest`` so that any scheme runs; this is Stern's.
DEST = SternParams.dest


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


def native_coupling(n_c: int | None = None) -> float:
    """Nearest-neighbour coupling of the cleft as solved, µM."""
    return nearest_coupling(coupling_matrix(CleftGeometry.from_parameters(n_c)))


@dataclass
class CleftSparkParams:
    n_channels: float = _v("spark.n_channels")
    ca_rest: float = _v("puff.ca_rest")
    ca_per_open: float = field(default_factory=native_coupling)
    record_dt: float = _v("puff.record_dt")
    gating: SternParams = field(default_factory=SternParams)


def c_rates(state: np.ndarray, c: np.ndarray, sp: SternParams) -> np.ndarray:
    """Exit rates of every C channel, one block of ``n`` per exit kind (the
    columns ``sp.dest`` indexes), at each channel's own Ca2+."""
    return sp.exit_rates(state, c)


def initial_states(pp: CleftSparkParams, n: int, rng: np.random.Generator,
                   trigger: bool = False) -> np.ndarray:
    """Each channel drawn from the stationary state at ``ca_rest``;
    ``trigger`` opens the closed, uninactivated ones."""
    pi = np.cumsum(stationary(pp.ca_rest, pp.gating))
    state = np.searchsorted(pi, rng.random(n) * pi[-1])
    if trigger:
        state = pp.gating.trigger_map[state]
    return state


def couplings_for(pp: CleftSparkParams) -> np.ndarray:
    """``G`` for this cluster, off-diagonal scaled to ``ca_per_open``."""
    g = coupling_matrix(CleftGeometry.from_parameters(int(round(pp.n_channels))))
    own = np.diag(np.diag(g))
    native = nearest_coupling(g)
    scale = pp.ca_per_open / native if native > 0 else 0.0
    return own + (g - own) * scale


def simulate_sparks_cleft(p: float = 0.0, duration: float = 5.0, seed: int = 0,
                          pp: CleftSparkParams | None = None,
                          trigger: bool = False) -> PuffTrace:
    """Simulate the cleft array for ``duration`` s (``p`` is ignored, as in
    :func:`sparks.simulate_sparks`). ``PuffTrace.ca`` is the mean Ca2+ over
    the C channels at each bin's start. ``trigger`` opens every channel that
    is closed but not inactivated at t = 0 (the resting draw is kept for the
    inactivated ones): a stand-in for the V channels' stimulus, so that a
    spark can be timed where none starts by itself."""
    pp = pp or CleftSparkParams()
    sp = pp.gating
    rng = np.random.default_rng(seed)
    g = couplings_for(pp)
    n = g.shape[0]
    state = initial_states(pp, n, rng, trigger)
    is_open = sp.open_mask[state].astype(float)
    c = pp.ca_rest + g @ is_open
    n_rec = int(np.floor(duration / pp.record_dt + 1e-9)) + 1
    t_out = np.arange(n_rec) * pp.record_dt
    open_out = np.zeros(n_rec, dtype=np.int32)
    peak_out = np.zeros(n_rec, dtype=np.int32)
    ca_out = np.zeros(n_rec)
    inact_out = np.zeros(n_rec, dtype=np.int32)
    t, k, n_open = 0.0, 0, int(is_open.sum())
    while True:
        rates = c_rates(state, c, sp)
        total = rates.sum()
        t_next = t + rng.exponential(1.0 / total)
        # Record every bin that starts before the next event.
        while k < n_rec and t_out[k] <= t_next:
            open_out[k], ca_out[k] = n_open, c.mean()
            inact_out[k] = int(sp.inact_mask[state].sum())
            peak_out[k] = max(peak_out[k], n_open)
            k += 1
        if t_next > duration:
            break
        t = t_next
        pick = int(np.searchsorted(np.cumsum(rates), rng.random() * total))
        ch, gate = pick % n, pick // n
        was_open = sp.open_mask[state[ch]]
        state[ch] = sp.dest[state[ch], gate]
        now_open = sp.open_mask[state[ch]]
        if was_open != now_open:
            sign = 1.0 if now_open else -1.0
            c += sign * g[:, ch]
            n_open += int(sign)
            if k > 0:
                peak_out[k - 1] = max(peak_out[k - 1], n_open)
    return PuffTrace(t_out[:k], open_out[:k], ca_out[:k], pp, p, peak_out[:k],
                     inact_out[:k])
