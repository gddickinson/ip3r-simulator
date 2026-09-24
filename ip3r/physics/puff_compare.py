"""Blips and puffs, measured the same way on either receptor.

Two clusters share everything except the receptor: the De Young-Keizer
subunit scheme (``puffs``) and the park/drive scheme (``puffs_pd``). Both
record, per ``puff.record_dt`` bin, the number open at the bin's start and
the most open at once within it. This module reads both with one ruler:

* **Fano factor** of the number open at once (``puffs.fano``): 1 for
  independent channels, above 1 when openings cluster.
* **Event sizes**: an event is a run of bins with any channel open, and its
  size is its peak. A *blip* peaks at one channel.
* **Recruitment share**: of the events that involve more than one channel,
  the fraction whose peak reaches ``puff.large_fraction`` of the cluster,
  and the rate of such events per second.
  Independent coincidences almost never get that far. An all-or-none release
  does. Blips and puffs "separate cleanly" when multi-channel events are
  either small coincidences or near-whole-cluster puffs, with few between.

``coupling_scan`` runs both receptors over the same couplings, so the
comparison does not rest on either model's chosen coupling.

**Sparks.** A RyR1 cluster (:mod:`ip3r.physics.sparks`, key ``ryr1``) is
read with the same ruler. It is not in ``MODELS`` (the IP3R pair the
coupling scan compares at one IP3), because its natural coupling is ~8 µM,
outside the IP3R scan; ``spark_scan`` sweeps a band around it instead.
The cleft cluster (:mod:`ip3r.physics.sparks_cleft`, key ``ryr1-cleft``)
is the same receptor with a spatial Ca2+ field; its coupling is the
nearest-neighbour value, and ``spark_ends`` measures how each spark ends.
"""

from __future__ import annotations

import math

import numpy as np

from ..parameters import PARAMETERS as _P
from .puffs import PuffParams, PuffTrace, detect_events, fano, simulate_cluster
from .puffs_pd import ParkDrivePuffParams, simulate_cluster_pd
from .sparks import SparkParams, simulate_sparks, spark_couplings
from .sparks_cleft import CleftSparkParams, native_coupling, simulate_sparks_cleft

__all__ = ["MODELS", "ALL_MODELS", "SPARK", "SPARK_CLEFT", "SPARKS", "MODEL_LABELS",
           "spark_scan", "params_for", "simulate", "event_sizes", "recruitment",
           "spark_ends", "coupling_effect", "scan_couplings", "coupling_scan"]

MODELS = ("dyk", "park-drive")
SPARK = "ryr1"
SPARK_CLEFT = "ryr1-cleft"
SPARKS = (SPARK, SPARK_CLEFT)
ALL_MODELS = MODELS + SPARKS
MODEL_LABELS = {"dyk": "De Young–Keizer subunits",
                "park-drive": "Park/drive (Siekmann; Cao 2013)",
                SPARK: "RyR1 sparks (Stern 1997 scheme)",
                SPARK_CLEFT: "RyR1 sparks in the cleft (Stern 1997 geometry)"}
_SIM = {"dyk": (simulate_cluster, PuffParams),
        "park-drive": (simulate_cluster_pd, ParkDrivePuffParams),
        SPARK: (simulate_sparks, SparkParams),
        SPARK_CLEFT: (simulate_sparks_cleft, CleftSparkParams)}


def params_for(model: str, coupling: float | None = None,
               n_channels: float | None = None):
    """The model's default cluster parameters, optionally re-coupled/resized."""
    pp = _SIM[model][1]()
    if coupling is not None:
        pp.ca_per_open = coupling
    if n_channels is not None:
        pp.n_channels = n_channels
    return pp


def simulate(model: str, p: float, duration: float, seed: int = 0,
             pp=None) -> PuffTrace:
    return _SIM[model][0](p, duration, seed, pp or params_for(model))


def event_sizes(tr: PuffTrace) -> np.ndarray:
    """Peak channels open in each event (one entry per event)."""
    return np.array([e["peak_open"] for e in detect_events(tr)], dtype=int)


def recruitment(tr: PuffTrace, large_fraction: float | None = None) -> dict:
    """The ruler: Fano, resting activity and the event-size split."""
    frac = _P.value("puff.large_fraction") if large_fraction is None else large_fraction
    n = int(round(tr.params.n_channels))
    large_at = max(2, math.ceil(frac * n))
    sizes = event_sizes(tr)
    duration = float(tr.t[-1] - tr.t[0]) if len(tr.t) > 1 else 0.0
    multi = int((sizes >= 2).sum())
    large = int((sizes >= large_at).sum())
    return {"fano": fano(tr),
            "open_fraction": float(tr.n_open.mean() / n),
            "max_open": int(tr.peaks.max()) if len(tr.peaks) else 0,
            "n_events": int(len(sizes)), "blips": int((sizes == 1).sum()),
            "multi": multi, "large": large, "large_at": large_at,
            "large_share": large / multi if multi else float("nan"),
            "large_per_s": large / duration if duration > 0 else float("nan")}


def spark_ends(tr: PuffTrace, large_fraction: float | None = None) -> list[dict]:
    """Every event reaching the large size: duration, peak, and how many
    channels were inactivated (CI or I) at its first and after its last bin.
    Needs a RyR1 trace (``n_inactivated``)."""
    large_at = recruitment(tr, large_fraction)["large_at"]
    step = float(tr.t[1] - tr.t[0])
    out = []
    for e in detect_events(tr):
        if e["peak_open"] < large_at:
            continue
        a = int(round((e["start"] - tr.t[0]) / step))
        b = min(a + int(round(e["duration"] / step)), len(tr.t) - 1)
        out.append({**e, "inactivated_start": int(tr.n_inactivated[a]),
                    "inactivated_end": int(tr.n_inactivated[b])})
    return out


def coupling_effect(p: float = 0.2, duration: float = 20.0, seed: int = 0,
                    pp=None, model: str = "dyk") -> dict:
    """Coupling on versus off, same seed and cluster, for one receptor.

    The measured statement of the puff mechanism: coupling should raise the
    Fano factor above the ~1 of independent channels.
    """
    base = pp or params_for(model)
    out = {}
    for label, cpo in (("coupled", base.ca_per_open), ("uncoupled", 0.0)):
        trial = type(base)(**{**base.__dict__, "ca_per_open": cpo})
        out[label] = recruitment(simulate(model, p, duration, seed, trial))
    return out


def scan_couplings() -> np.ndarray:
    """0, then ``puff.scan_points`` geometric steps over the scan range."""
    return np.concatenate([[0.0], np.geomspace(
        _P.value("puff.scan_min"), _P.value("puff.scan_max"),
        int(round(_P.value("puff.scan_points"))))])


def coupling_scan(p: float = 0.2, duration: float = 10.0, seed: int = 0,
                  models=MODELS, couplings=None) -> dict[str, list[dict]]:
    """Both receptors over the same couplings: ``{model: [row, ...]}``, each
    row a :func:`recruitment` dict with its ``coupling``."""
    cs = scan_couplings() if couplings is None else np.asarray(couplings, float)
    return {m: [{"coupling": float(c),
                 **recruitment(simulate(m, p, duration, seed, params_for(m, c)))}
                for c in cs] for m in models}


def spark_scan(duration: float = 10.0, seed: int = 0, couplings=None,
               model: str = SPARK) -> list[dict]:
    """A RyR1 cluster over the registered spark band (for the cleft, the
    band's factors times its own nearest-neighbour coupling): one
    :func:`recruitment` row per coupling."""
    if couplings is None:
        couplings = spark_couplings(native_coupling() if model == SPARK_CLEFT else None)
    return [{"coupling": float(c),
             **recruitment(simulate(model, 0.0, duration, seed, params_for(model, c)))}
            for c in np.asarray(couplings, float)]
