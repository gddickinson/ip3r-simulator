"""The whole couplon: V channels opened by voltage, C channels by their Ca2+.

Stern, Pizarro & Rios 1997's couplon is two rows of release channels in a
chessboard. Half face a tetrad of voltage sensors (V channels, Rios 1993's
allosteric model, :mod:`allosteric_v`) and half do not (C channels, the
Ca2+-gated scheme of :mod:`ryr_gating`). The cleft (:mod:`cleft`) turns
each open channel into Ca2+ at every C channel: ``G`` from C channels,
``H`` from V channels at the V current.

**Exact, in two passes.** A V channel is "neither activated nor inactivated
by Ca2+" (Stern), so the V channels are independent Markov chains driven by
the voltage alone. Each is simulated first, by Gillespie's method, through
the piecewise-constant voltage protocol (a step boundary is exact because
the process is memoryless). The C array is then simulated by Gillespie's
method with each V opening and closing as a scheduled event: between
events every rate is constant, so nothing is approximated and there is no
step to converge.

The V channels carry no Mg2+ or Ca2+ action of any kind. That is Stern's
assumption, and it is also the reading of Laver 2018 that the voltage
sensor lifts Mg2+ inhibition from the channels it touches. The C channels
carry whatever :class:`ryr_gating.SternParams` says, Mg2+ included.

**SR content.** ``scale`` makes every unitary current (V and C) a fraction
of its full value, piecewise constant on a given time grid: Stern's Fig. 20
made "the unitary currents ... proportional to SR calcium content". The
fraction changes the cleft Ca2+ at each change point, which are scheduled
events like the V openings, so the simulation stays exact given the path;
:mod:`lumen` finds the path.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from .allosteric_v import RiosParams, exits, stationary as v_stationary
from .cleft import CleftGeometry, v_coupling_matrix
from .ryr_mixed import bind
from .sparks_cleft import CleftSparkParams, c_rates, couplings_for, initial_states

__all__ = ["CouplonParams", "CouplonTrace", "step_protocol", "v_trajectory",
           "simulate_couplon"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class CouplonParams:
    cleft: CleftSparkParams = field(default_factory=CleftSparkParams)
    v: RiosParams = field(default_factory=RiosParams)
    v_current: float = _v("ec.v_unitary_current")     # pA


@dataclass
class CouplonTrace:
    t: np.ndarray               # s, bin starts
    c_open: np.ndarray          # C channels open at each bin's start
    v_open: np.ndarray          # V channels open at each bin's start
    c_peak: np.ndarray          # most C channels open at once in each bin
    c_inactivated: np.ndarray   # C channels in CI or I at each bin's start
    params: CouplonParams
    protocol: tuple
    scale: np.ndarray | None = None     # current fraction at each bin's start

    def flux(self) -> np.ndarray:
        """Release current, pA: C and V channels at their unitary currents
        (times the SR-content fraction, when there is one)."""
        i_c = _P.value("spark.unitary_current")
        full = self.c_open * i_c + self.v_open * self.params.v_current
        return full if self.scale is None else full * self.scale


def step_protocol(v: float, holding: float | None = None,
                  pulse: float | None = None) -> tuple:
    """((0, v), (pulse, holding)): a step from ``holding`` at t = 0."""
    holding = _P.value("ec.holding") if holding is None else holding
    pulse = _P.value("ec.pulse") if pulse is None else pulse
    return ((0.0, float(v)), (float(pulse), float(holding)))


def v_trajectory(protocol: tuple, duration: float, rp: RiosParams,
                 rng: np.random.Generator, holding: float | None = None
                 ) -> tuple[list[float], list[int]]:
    """One V channel through ``protocol``: the times it opens or closes and
    the sign (+1 opens). It starts in the stationary state at ``holding``."""
    holding = _P.value("ec.holding") if holding is None else holding
    tables = [exits(v, rp) for _, v in protocol]
    bounds = [t for t, _ in protocol[1:]] + [duration]
    pi = np.cumsum(v_stationary(holding, rp))
    s = int(np.searchsorted(pi, rng.random() * pi[-1]))
    t, seg, times, signs = 0.0, 0, [], []
    dest, rate = tables[0]
    while True:
        r = rate[s]
        total = r.sum()
        t_next = t + rng.exponential(1.0 / total) if total > 0 else np.inf
        if t_next > bounds[seg]:
            t = bounds[seg]
            seg += 1
            if seg >= len(tables):
                return times, signs
            dest, rate = tables[seg]
            continue
        t = t_next
        k = int(np.searchsorted(np.cumsum(r), rng.random() * total))
        new = int(dest[s, min(k, 2)])
        if (new >= 5) != (s >= 5):
            times.append(t)
            signs.append(1 if new >= 5 else -1)
        s = new


def simulate_couplon(protocol: tuple, duration: float, seed: int = 0,
                     pp: CouplonParams | None = None,
                     scale: tuple[np.ndarray, np.ndarray] | None = None
                     ) -> CouplonTrace:
    """One couplon through the voltage ``protocol`` ((t_start, mV), ...),
    recorded every ``pp.cleft.record_dt`` for ``duration`` s. ``scale`` =
    (times, fractions): the unitary currents are ``fractions[i]`` of their
    value from ``times[i]`` on (``times[0]`` must be 0)."""
    pp = pp or CouplonParams()
    rng = np.random.default_rng(seed)
    g = couplings_for(pp.cleft)
    n = g.shape[0]
    h = v_coupling_matrix(pp.v_current, CleftGeometry.from_parameters(n))
    # Pass 1: the V channels, merged into one schedule of (time, channel, sign).
    sched = []
    for k in range(h.shape[1]):
        ts, ss = v_trajectory(protocol, duration, pp.v, rng)
        sched += [(t, k, s) for t, s in zip(ts, ss)]
    f_t, f_v = (np.zeros(1), np.ones(1)) if scale is None else scale
    sched += [(float(t), -2, i) for i, t in enumerate(f_t)       # changes only
              if i > 0 and f_v[i] != f_v[i - 1]]
    sched.sort()
    sched.append((np.inf, -1, 0))
    # Pass 2: the C array, with the V events as scheduled changes of Ca2+.
    sp = bind(pp.cleft.gating, n, rng)
    state = initial_states(pp.cleft, n, rng, gating=sp)
    v_open = np.zeros(h.shape[1])
    f = float(f_v[0])
    c_rel = g @ sp.open_mask[state].astype(float)     # at full current
    c = pp.cleft.ca_rest + f * c_rel
    dt = pp.cleft.record_dt
    n_rec = int(np.floor(duration / dt + 1e-9)) + 1
    t_out = np.arange(n_rec) * dt
    c_out = np.zeros(n_rec, np.int32)
    v_out = np.zeros(n_rec, np.int32)
    peak = np.zeros(n_rec, np.int32)
    inact = np.zeros(n_rec, np.int32)
    f_out = np.ones(n_rec)
    t, k, i_s = 0.0, 0, 0
    n_open, n_v = int(sp.open_mask[state].sum()), 0
    rates = c_rates(state, c, sp)
    total = rates.sum()
    while True:
        t_c = t + rng.exponential(1.0 / total) if total > 0 else np.inf
        t_s = sched[i_s][0]
        t_next = min(t_c, t_s, duration + dt)
        while k < n_rec and t_out[k] <= t_next:
            c_out[k], v_out[k], f_out[k] = n_open, n_v, f
            inact[k] = int(sp.inact_mask[state].sum())
            peak[k] = max(peak[k], n_open)
            k += 1
        if t_next > duration:
            break
        t = t_next
        if t_s <= t_c and sched[i_s][1] == -2:   # the SR content steps
            f = float(f_v[sched[i_s][2]])
            i_s += 1
        elif t_s <= t_c:                       # a V channel opens or closes
            _, ch, sign = sched[i_s]
            i_s += 1
            v_open[ch] += sign
            n_v += sign
            c_rel += sign * h[:, ch]
        else:                                  # a C channel's gate moves
            pick = int(np.searchsorted(np.cumsum(rates), rng.random() * total))
            ch, gate = pick % n, pick // n
            was_open = sp.open_mask[state[ch]]
            state[ch] = sp.dest[state[ch], gate]
            if was_open != sp.open_mask[state[ch]]:
                sign = 1.0 if not was_open else -1.0
                c_rel += sign * g[:, ch]
                n_open += int(sign)
                if k > 0:
                    peak[k - 1] = max(peak[k - 1], n_open)
        c = pp.cleft.ca_rest + f * c_rel
        rates = c_rates(state, c, sp)
        total = rates.sum()
    return CouplonTrace(t_out[:k], c_out[:k], v_out[:k], peak[:k], inact[:k],
                        pp, tuple(protocol), None if scale is None else f_out[:k])
