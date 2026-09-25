"""Puffs read from fluorescence, and their statistics (Round 7.2).

The experiments see F/F0, not channels, so puffs are read the way Cao et
al. 2013 read their model:

1. *Events* are runs of bins whose dF/F0 exceeds ``domain.detect_fraction``
   of one open receptor's steady dF/F0 (``blip_df``). Each carries its
   peak dF/F0, its peak microdomain rise and the most channels open in it.
2. The mean *blip* is the mean peak dF/F0 of the events with one channel
   open. A *puff* is an event whose peak exceeds ``domain.puff_threshold``
   blips (Cao's dF/F0 > 3 against a mean blip of 1.6).
3. *Inter-puff intervals* run from one puff's onset to the next.

Thurley et al. 2011's IPI density (Cao 2013 Eq. 14) is a Poisson process
whose rate recovers after each puff:

    P(t) = lam (1 - e^{-xi t}) exp(-lam t + lam (1 - e^{-xi t}) / xi)

With ``xi >> lam`` it is the exponential ``lam e^{-lam t}`` (Eq. 15).
``fit_thurley`` is its maximum-likelihood fit. The test draws samples from
known ``(lam, xi)`` and recovers them, and an exponential sample must fit
with ``xi`` far above ``lam``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from ..parameters import PARAMETERS as _P
from . import microdomain as md

__all__ = ["blip_df", "fluorescence_events", "PuffStats", "puff_stats",
           "thurley_pdf", "thurley_sample", "fit_thurley", "refractory_lr",
           "cv"]


def blip_df(rest: md.Rest, d: md.DomainParams) -> float:
    """Steady dF/F0 with one receptor open and the store clamped."""
    cb = (d.k_diff * rest.c + d.k_ipr * rest.cs) / (d.k_diff + d.k_ipr)
    kd = d.fluo_koff / d.fluo_kon
    return (cb / (cb + kd)) / (rest.c / (rest.c + kd)) - 1.0


def fluorescence_events(tr, detect: float | None = None) -> list[dict]:
    """Runs of bins above the detection level in a ``DomainTrace``."""
    d = tr.params.domain
    frac = _P.value("domain.detect_fraction") if detect is None else detect
    level = frac * blip_df(tr.rest, d)
    on = (tr.f_peak - 1.0) > level
    edges = np.flatnonzero(np.diff(np.concatenate([[0], on.astype(int), [0]])))
    out = []
    for a, b in zip(edges[::2], edges[1::2]):
        out.append({"start": float(tr.t[a]), "end": float(tr.t[b - 1]),
                    "peak_df": float(tr.f_peak[a:b].max() - 1.0),
                    "peak_ca": float(tr.ca[a:b].max() - tr.rest.c),
                    "peak_open": int(tr.peaks[a:b].max())})
    return out


@dataclass
class PuffStats:
    duration: float
    n_events: int
    n_blips: int
    blip_mean: float          # mean peak dF/F0 of one-channel events
    n_puffs: int
    amp_mean: float           # mean puff peak dF/F0
    amp_max: float
    ca_mean: float            # mean puff peak microdomain rise, µM
    open_mean: float          # mean puff peak number open
    ipis: np.ndarray
    lam: float = math.nan     # Thurley fit
    xi: float = math.nan
    refractory_lr: float = math.nan   # 2 (logL Eq. 14 - logL exponential)

    @property
    def rate(self) -> float:
        return self.n_puffs / self.duration

    @property
    def cv(self) -> float:
        return cv(self.ipis)

    @property
    def amp_in_blips(self) -> float:
        return self.amp_mean / self.blip_mean if self.blip_mean > 0 else math.nan


def puff_stats(tr, fit: bool = True, detect: float | None = None,
               threshold: float | None = None) -> PuffStats:
    """Blips, puffs and inter-puff intervals of a ``DomainTrace``.

    ``detect`` and ``threshold`` default to ``domain.detect_fraction`` and
    ``domain.puff_threshold``.
    """
    ev = fluorescence_events(tr, detect)
    blips = [e["peak_df"] for e in ev if e["peak_open"] == 1]
    blip = float(np.mean(blips)) if blips else blip_df(tr.rest, tr.params.domain)
    if threshold is None:
        threshold = _P.value("domain.puff_threshold")
    cut = threshold * blip
    puffs = [e for e in ev if e["peak_df"] > cut]
    ipis = np.diff([e["start"] for e in puffs])
    duration = float(tr.t[-1] - tr.t[0])
    s = PuffStats(duration, len(ev), len(blips), blip, len(puffs),
                  _mean([e["peak_df"] for e in puffs]),
                  max((e["peak_df"] for e in puffs), default=math.nan),
                  _mean([e["peak_ca"] for e in puffs]),
                  _mean([e["peak_open"] for e in puffs]), ipis)
    if fit and len(ipis) >= 5:
        s.lam, s.xi = fit_thurley(ipis)
        s.refractory_lr = refractory_lr(ipis, s.lam, s.xi)
    return s


def _mean(v) -> float:
    return float(np.mean(v)) if len(v) else math.nan


def cv(ipis) -> float:
    ipis = np.asarray(ipis, float)
    return float(ipis.std(ddof=1) / ipis.mean()) if len(ipis) > 1 else math.nan


def thurley_pdf(t, lam: float, xi: float):
    """Cao 2013 Eq. 14 (Thurley et al. 2011)."""
    t = np.asarray(t, float)
    rec = -np.expm1(-xi * t)
    return lam * rec * np.exp(-lam * t + lam * rec / xi)


def _log_likelihood(t, lam, xi):
    rec = -np.expm1(-xi * t)
    with np.errstate(divide="ignore"):
        return float(np.sum(np.log(lam * rec) - lam * t + lam * rec / xi))


def fit_thurley(ipis) -> tuple[float, float]:
    """Maximum-likelihood ``(lam, xi)`` of Eq. 14, searched in log space.

    The exponential (Eq. 15) is the limit ``xi -> infinity``; the search is
    bounded at 1000x the rate, where Eq. 14 and Eq. 15 differ by < 0.1 %
    in likelihood per interval.
    """
    t = np.asarray(ipis, float)
    lam0 = 1.0 / t.mean()
    lo, hi = np.log(lam0 / 100), np.log(lam0 * 1000)

    def nll(x):
        return -_log_likelihood(t, math.exp(x[0]), math.exp(x[1]))

    best = None
    for xi0 in (0.3, 1.0, 3.0, 30.0):
        r = minimize(nll, [np.log(lam0), np.log(xi0 * lam0)],
                     method="L-BFGS-B", bounds=[(lo, hi), (lo, hi)])
        if best is None or r.fun < best.fun:
            best = r
    return float(math.exp(best.x[0])), float(math.exp(best.x[1]))


def refractory_lr(ipis, lam: float, xi: float) -> float:
    """Likelihood-ratio statistic of Eq. 14 against the exponential (one
    extra parameter): ~chi2(1) when the intervals are exponential, so above
    3.84 the refractory period is significant at 5 %."""
    t = np.asarray(ipis, float)
    lam_exp = 1.0 / t.mean()
    ll_exp = float(np.sum(np.log(lam_exp) - lam_exp * t))
    return 2.0 * (_log_likelihood(t, lam, xi) - ll_exp)


def thurley_sample(lam: float, xi: float, n: int, rng) -> np.ndarray:
    """``n`` intervals from Eq. 14, by inverting its survival function."""
    u = rng.random(n)
    out = np.empty(n)
    for i, ui in enumerate(u):
        target = -math.log(ui)            # integrated hazard to reach
        a, b = 0.0, target / lam + 1.0 / xi + 1.0
        for _ in range(200):              # bisection on H(t) = target
            m = 0.5 * (a + b)
            if lam * m - lam * (-math.expm1(-xi * m)) / xi < target:
                a = m
            else:
                b = m
        out[i] = 0.5 * (a + b)
    return out
