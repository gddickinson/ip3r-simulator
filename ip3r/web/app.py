"""FastAPI web app for the IP3R Structural Simulator.

Serves the front end in ``static/`` and a read-only JSON API under ``/api``.
Run with ``python run_web.py`` (http://127.0.0.1:8194, local only).
``scripts/build_static.py`` exports every path in ``STATIC_ENDPOINTS`` to JSON
for the GitHub Pages site, so keep that list in step with the GET routes.

Heavy results (structures measured, the full-resolution elastic network,
the 8TKG -> 8TKF morph and the state panel) are precomputed by the repo's
own code into ``precomputed/`` (``scripts/make_web_data.py``) and served as
they are; the gating constants and reference bells are computed here at
request time from the parameter registry.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

HERE = Path(__file__).resolve().parent
STATIC = HERE / "static"
PRECOMPUTED = HERE / "precomputed"

#: Deposits with a precomputed viewer file (export_structures.VIEWER_DEPOSITS),
#: and the two whose full network is precomputed (export_dynamics.MODE_DEPOSITS).
#: Literal here so the API imports without numpy/scipy's structure stack.
DEPOSITS = ("6DQN", "6DQJ", "7T3P", "8TKG", "8TKH", "8TLA", "8TKF", "7T3T")
MODE_DEPOSITS = ("8TKG", "6DQN")

app = FastAPI(title="IP3R Structural Simulator")
app.mount("/static", StaticFiles(directory=STATIC), name="static")


@lru_cache(maxsize=None)
def _load(name: str) -> dict:
    path = PRECOMPUTED / name
    if not path.exists():
        raise HTTPException(404, f"{name} not precomputed (run scripts/make_web_data.py)")
    return json.loads(path.read_text())


def _check(pdb: str, allowed) -> str:
    pdb = pdb.upper()
    if pdb not in allowed:
        raise HTTPException(404, f"{pdb} is not one of {', '.join(allowed)}")
    return pdb


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC / "index.html")


def _summary(pdb: str) -> dict:
    s = _load(f"structure_{pdb}.json")
    keep = ("pdb_id", "title", "paralog", "organism", "state", "resolution", "is_control",
            "painted_in", "constrictions")
    out = {k: s[k] for k in keep}
    out["n_ca"] = sum(len(u["resi"]) for u in s["subunits"])
    out["n_ip3"] = len(s["ligands"])
    return out


@app.get("/api/info")
def info():
    """Project summary and counts for the home page."""
    from ..io.registry import load_registry
    from ..parameters import PARAMETERS
    reg = load_registry()
    refs = json.loads((HERE.parent / "resources" / "references.json").read_text())["references"]
    return {
        "title": "IP3R Structural Simulator",
        "description": "An interactive, physics-based 3-D model of the inositol 1,4,5-trisphosphate "
                       "receptor (IP3R), the channel that releases calcium from the cell's internal "
                       "store, and an independent checking instrument for the ip3r_genes "
                       "publication project.",
        "deposits_shown": len(DEPOSITS),
        "registry": {"ip3r": sum(e.family == "IP3R" for e in reg),
                     "ryr": sum(e.family == "RyR" for e in reg)},
        "parameters": len(PARAMETERS.parameters),
        "references": len(refs),
        "provenance": _load("provenance.json"),
    }


@app.get("/api/deposits")
def deposits():
    """The deposits the web viewer offers, with their measured constrictions."""
    return {"deposits": [_summary(p) for p in DEPOSITS], "modes": list(MODE_DEPOSITS)}


@app.get("/api/palette")
def palette():
    return _load("palette.json")


@app.get("/api/structure/{pdb}")
def structure(pdb: str):
    """C-alpha trace (axis frame, cytosol +z), annotation and channel summary."""
    return _load(f"structure_{_check(pdb, DEPOSITS)}.json")


@app.get("/api/modes/{pdb}")
def modes(pdb: str):
    """The full-resolution ANM (registered stride, cutoff, springs) with C4 irreps."""
    return _load(f"modes_{_check(pdb, MODE_DEPOSITS)}.json")


@app.get("/api/transition")
def transition():
    """8TKG -> 8TKF: matched basis, restrained morph, gate path, mode overlap."""
    return _load("transition_8TKG_8TKF.json")


@app.get("/api/states")
def states():
    """The ITPR3 state panel's pore profiles and the 7T3T control."""
    return _load("states_ITPR3.json")


@app.get("/api/gating")
def gating():
    """Registered constants of the three gating models and Python's reference bells."""
    from .gating_api import gating_record
    return gating_record()


@app.get("/api/anm-parameters")
def anm_parameters():
    """The registered ANM parameters (the reduced model in the browser starts from these)."""
    from ..parameters import PARAMETERS as _P
    keys = ("anm.cutoff", "anm.gamma", "anm.d0", "anm.n_modes", "anm.stride",
            "anm.symmetry_tolerance", "anm.min_collectivity",
            "anm.scan_cutoff_low", "anm.scan_cutoff_high")
    return {k: {"value": _P.value(k), "name": _P.get(k).name, "unit": _P.get(k).unit,
                "citation": _P.get(k).citation, "source_note": _P.get(k).source_note}
            for k in keys}


# GET paths (without the /api prefix) written to the static site as api/<path>.json.
STATIC_ENDPOINTS = (["/info", "/deposits", "/palette", "/transition", "/states", "/gating",
                     "/anm-parameters"]
                    + [f"/structure/{p}" for p in DEPOSITS]
                    + [f"/modes/{p}" for p in MODE_DEPOSITS])
