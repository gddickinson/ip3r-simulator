"""The Round 7.2 measurements on the microdomain cluster.

- ``ah42_scan``: the IPI distribution against the h42 recovery rate a_h42
  (``pd.lam_h42_closed``). Cao 2013 Fig. 4 (N = 10, 0.1 µM IP3): as a_h42
  rises from 0.1 to 5 /s the distribution turns from peaked to exponential,
  with lam 0.1-0.5 /s and xi 0.5-2.2 /s.
- ``n_scan``: puff amplitude against cluster size. Cao 2013 Fig. 8 and S9:
  the mean dF/F0 amplitude bends over near N = 12 while the mean Ca2+
  amplitude stays linear, because the dye is nonlinear.
- ``store_scan``: the Round 4 question. The mean-field cluster sits at a
  sustained ~9 % open above ~0.5 µM coupling. Here the release rate is
  raised to the same couplings under five conditions (``CONDITIONS``),
  each adding one thing to the last: the mean field; its mouth deepened to
  the microdomain's; the microdomain's kinetics and dye in a fixed bath;
  the cytosol free; the store free. Whichever step removes the sustained
  state is what the mean-field cluster lacks.

Every point runs in its own process with one seed (``domain.scan_seed``).
The parameter sets are built in the parent and passed whole, so an edited
registry reaches the workers.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace

import numpy as np

from ..parameters import PARAMETERS as _P
from . import microdomain as md
from .puff_stats import PuffStats, puff_stats
from .puffs_domain import DomainPuffParams, simulate_cluster_domain
from .puffs_pd import ParkDrivePuffParams, simulate_cluster_pd

__all__ = ["ah42_values", "n_values", "scale_values", "ah42_scan", "n_scan",
           "CONDITIONS", "CONDITION_NOTES", "Occupancy", "StoreRow",
           "store_scan", "run_all"]


def _geom(lo, hi, n):
    return np.geomspace(lo, hi, int(round(n)))


def ah42_values() -> np.ndarray:
    return _geom(_P.value("domain.ah42_min"), _P.value("domain.ah42_max"),
                 _P.value("domain.ah42_points"))


def n_values() -> list[int]:
    v = _geom(_P.value("domain.n_min"), _P.value("domain.n_max"),
              _P.value("domain.n_points"))
    return sorted({int(round(x)) for x in v})


def scale_values() -> np.ndarray:
    return _geom(1.0, _P.value("domain.store_scale_max"),
                 _P.value("domain.store_points"))


def run_all(fn, jobs, workers: int | None = None) -> list:
    """``fn(*job)`` for every job, in processes (in order)."""
    if workers == 1 or len(jobs) == 1:
        return [fn(*j) for j in jobs]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(fn, *zip(*jobs)))


def _method():
    return (_P.value("domain.detect_fraction"),
            _P.value("domain.puff_threshold"),
            int(_P.value("domain.scan_seed")))


def _stats_job(p, duration, dp, method) -> PuffStats:
    detect, threshold, seed = method
    tr = simulate_cluster_domain(p, duration, seed, dp)
    return puff_stats(tr, detect=detect, threshold=threshold)


def ah42_scan(p: float = 0.1, n: int = 10, values=None,
              duration: float | None = None, workers=None) -> list:
    """``[(a_h42, PuffStats)]`` at IP3 ``p`` for a cluster of ``n``."""
    values = ah42_values() if values is None else values
    duration = duration or _P.value("domain.ipi_duration")
    method = _method()
    jobs = []
    for a in values:
        dp = DomainPuffParams(n_channels=n)
        dp.receptor = replace(dp.receptor, lam_h42_closed=float(a))
        jobs.append((p, duration, dp, method))
    return list(zip(values, run_all(_stats_job, jobs, workers)))


def n_scan(p: float = 0.2, values=None, duration: float | None = None,
           scale: float = 1.0, workers=None) -> list:
    """``[(N, PuffStats)]`` at IP3 ``p``, per-receptor release fixed (and
    multiplied by ``scale``: the calibration raises it until puff Ca2+
    reaches the dye's K_d, where dF/F0 must bend while Ca2+ does not)."""
    values = n_values() if values is None else values
    duration = duration or _P.value("domain.ipi_duration")
    method = _method()
    jobs = []
    for n in values:
        dp = DomainPuffParams(n_channels=n)
        dp.domain = dp.domain.scaled(scale)
        jobs.append((p, duration, dp, method))
    return list(zip(values, run_all(_stats_job, jobs, workers)))


#: The conditions of the store scan, from most idealised to least.
CONDITIONS = ("mean-field", "mean-field, deep mouth", "bath", "store", "none")
CONDITION_NOTES = {
    "mean-field": "puffs_pd at the same coupling (mouth pd.ca_mouth)",
    "mean-field, deep mouth": "puffs_pd, mouth = the microdomain's",
    "bath": "microdomain, cytosol and store held",
    "store": "microdomain, store held (Cao 2013)",
    "none": "microdomain, store free to deplete",
}


@dataclass
class Occupancy:
    open_fraction: float     # mean share of the cluster open
    late_fraction: float     # ... over the run's second half
    active: float            # share of time with at least one open
    fano: float              # of the number open

    @classmethod
    def of(cls, n_open, n: float) -> "Occupancy":
        half = len(n_open) // 2
        m = float(n_open.mean())
        return cls(m / n, float(n_open[half:].mean()) / n,
                   float((n_open > 0).mean()),
                   float(n_open.var() / m) if m > 0 else float("nan"))


@dataclass
class StoreRow:
    scale: float           # k_ipr multiple
    coupling: float        # µM per open receptor, store clamped
    runs: dict             # condition -> Occupancy
    store_start: float     # µM
    store_min: float       # µM, free store
    store_end: float
    cytosol_max: dict      # condition -> peak cytosolic Ca2+, µM


def _store_job(p, duration, scale, dp_base, pp_base, seed):
    d = dp_base.domain.scaled(scale)
    rest = md.rest_state(p, d)
    coupling = md.coupling_per_open(rest, d)
    n = float(dp_base.n_channels)
    runs, cyto, traces = {}, {}, {}
    for clamp in ("bath", "store", "none"):
        tr = simulate_cluster_domain(p, duration, seed,
                                     replace(dp_base, domain=d, clamp=clamp))
        runs[clamp] = Occupancy.of(tr.n_open, n)
        cyto[clamp] = float(tr.cytosol.max())
        traces[clamp] = tr
    for name, mouth in (("mean-field", pp_base.ca_mouth),
                        ("mean-field, deep mouth", d.mouth_per_store * rest.cs)):
        pp = replace(pp_base, ca_per_open=coupling, ca_mouth=mouth)
        runs[name] = Occupancy.of(simulate_cluster_pd(p, duration, seed, pp).n_open, n)
    free = traces["none"]
    return StoreRow(scale, coupling, {c: runs[c] for c in CONDITIONS},
                    rest.cs, float(free.store.min()), float(free.store[-1]),
                    cyto)

def store_scan(p: float = 0.2, scales=None, duration: float | None = None,
               workers=None) -> list[StoreRow]:
    """Open fractions at raised release rates, three ways (module doc)."""
    scales = scale_values() if scales is None else scales
    duration = duration or _P.value("domain.store_duration")
    dp, pp = DomainPuffParams(), ParkDrivePuffParams()
    seed = int(_P.value("domain.scan_seed"))
    jobs = [(p, duration, float(s), dp, pp, seed) for s in scales]
    return run_all(_store_job, jobs, workers)
