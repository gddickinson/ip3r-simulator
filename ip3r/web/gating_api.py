"""The three IP3R gating models' constants and reference bells, for the web.

The browser's ``gating_engine.js`` evaluates the steady-state open
probability of De Young-Keizer / Li-Rinzel (``physics.gating``), Mak et al.
1998 (``physics.gating_mak``) and park/drive (``physics.park_drive``) from
these constants. Every constant is the registered parameter (value, unit,
citation, note), read from the registry at request time, so the web shows
what the desktop app would compute; ``reference_bells`` are the Python
measurements (``physics.bell.measure_bell``) the page and the node test
compare the JavaScript against.
"""

from __future__ import annotations

from dataclasses import asdict

from ..parameters import PARAMETERS as _P
from ..parameters import references
from ..physics import gating, gating_mak, park_drive

__all__ = ["MODELS", "model_constants", "reference_bells", "gating_record"]

#: model key -> (label, registry prefix, parameter dataclass, bell function)
MODELS = {
    "dyk": ("De Young–Keizer, as reduced by Li & Rinzel", "gating",
            gating.GatingParams, gating.bell_at),
    "mak": ("Mak, McBride & Foskett 1998", "mak", gating_mak.MakParams, gating_mak.bell_at),
    "pd": ("Park/drive (Siekmann 2012, Cao 2013)", "pd", park_drive.ParkDriveParams,
           park_drive.bell_at),
}

_REGISTRY_KEY = {"dyk": {"subunits": "gating.subunits_required"}}


def _meta(key: str) -> dict:
    p = _P.get(key)
    return {"key": key, "name": p.name, "unit": p.unit, "citation": p.citation,
            "source_note": p.source_note}


def model_constants(model: str) -> dict:
    label, prefix, cls, _ = MODELS[model]
    values = asdict(cls())
    rename = _REGISTRY_KEY.get(model, {})
    return {"label": label, "values": values,
            "meta": {k: _meta(rename.get(k, f"{prefix}.{k}")) for k in values}}


def _bell_dict(b) -> dict:
    nan = b.c_half_inh != b.c_half_inh
    return {"c_peak": b.c_peak, "po_peak": b.po_peak, "c_half_act": b.c_half_act,
            "c_half_inh": None if nan else b.c_half_inh}


def reference_bells(ip3_levels: tuple | None = None) -> dict:
    """Python's bell for each model at each IP3 level (µM)."""
    levels = ip3_levels or (_P.value("mak.compare_ip3_low"), 1.0, _P.value("mak.compare_ip3_high"))
    return {m: [{"ip3": p, **_bell_dict(MODELS[m][3](p))} for p in levels] for m in MODELS}


def gating_record() -> dict:
    cited = set()
    models = {m: model_constants(m) for m in MODELS}
    for rec in models.values():
        cited |= {v["citation"] for v in rec["meta"].values()}
    refs = references()
    return {
        "models": models,
        "compare_ip3": [_P.value("mak.compare_ip3_low"), _P.value("mak.compare_ip3_high")],
        "reference_bells": reference_bells(),
        "references": {k: refs[k] for k in sorted(cited) if k in refs},
    }
