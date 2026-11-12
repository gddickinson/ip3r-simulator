"""Precompute the web site's heavy results with the repo's own code.

* ``modes_record`` — ``python -m ip3r modes PDB`` (the ANM at the registered
  stride, cutoff and spring weights), with each mode's displacement field
  turned into the deposit's axis frame so the viewer can animate it.
* ``transition_record`` — ``python -m ip3r transition 8TKG 8TKF --gate``:
  the matched basis, the restrained morph, the gate along it, and the
  overlap of 8TKG's modes with the observed displacement.
* ``states_record`` — ``python -m ip3r states`` plus the 7T3T control: each
  deposit's pore profile aligned on its own filter, as in
  ``scripts/figure_states.py``.

Each needs the mmCIF files in ``ref/``; the browser only reads the output.
"""

from __future__ import annotations

import numpy as np

from ..io import loader
from ..io.registry import get_entry
from ..parameters import PARAMETERS as _P
from ..physics.anm import ANM, tetramer_sites
from ..physics.network_checks import describe_local, local_modes
from ..physics.transition_modes import transition_overlap
from ..structure.channel import measure_channel
from ..structure.morph import NOTE as MORPH_NOTE
from ..structure.morph import morph
from ..structure.morph_pore import atom_path, gate_path
from ..structure.states import state_panel
from ..structure.symmetry import tetramer_frame
from ..structure.transition import prepare_transition
from .codec import encode_coords, encode_i16

__all__ = ["MODE_DEPOSITS", "VECTOR_SCALE", "modes_record", "transition_record",
           "states_record", "anm_parameters"]

MODE_DEPOSITS = ("8TKG", "6DQN")
#: Stored mode vectors are unit-peak; one int16 step is 1/30000 of the peak.
VECTOR_SCALE = 1.0 / 30000.0
#: Morph frames kept for the browser (every other one, both ends kept).
MORPH_KEEP = 2


def _r(x, n=4):
    return [round(float(v), n) for v in x]


def anm_parameters() -> dict:
    keys = ("anm.cutoff", "anm.gamma", "anm.d0", "anm.n_modes", "anm.stride",
            "anm.symmetry_tolerance", "anm.min_collectivity")
    return {k: _P.value(k) for k in keys}


def _unit_peak(v: np.ndarray) -> np.ndarray:
    peak = np.linalg.norm(v, axis=1).max()
    return v / peak if peak > 0 else v


def modes_record(pdb_id: str) -> dict:
    st = loader.load(pdb_id)
    fr = tetramer_frame(st)
    coords, res = tetramer_sites(st, fr)
    anm = ANM(coords, axis=fr.axis)
    ms = anm.label_symmetry(anm.calc_modes())
    kappa = ms.collectivity()
    local = local_modes(ms, res, tetramer_sites(st, fr, stride=1)[1])
    vectors = [encode_i16(_unit_peak(ms.vectors[i] @ fr.basis).ravel(), VECTOR_SCALE)
               for i in range(ms.n_modes)]
    return {
        "pdb_id": pdb_id, "chains": list(fr.chains), "residues": res.tolist(),
        "n_sites": int(len(coords)), "sites": encode_coords(fr.to_frame(coords)),
        "eigenvalues": [float(f"{v:.6e}") for v in ms.eigenvalues],
        "symmetry": [str(s) for s in ms.symmetry], "character": _r(ms.character),
        "collectivity": _r(kappa, 3), "collective": [bool(x) for x in ms.is_collective()],
        "lowest_collective_A": ms.first("A"),
        "local_artefacts": describe_local(local),
        "zero_modes_dropped": int(ms.meta["zero_modes_dropped"]),
        "vector_scale": VECTOR_SCALE, "vectors": vectors,
        "parameters": anm_parameters(),
    }


def transition_record(start: str = "8TKG", end: str = "8TKF") -> dict:
    st0, st1 = loader.load(start), loader.load(end)
    tr = prepare_transition(st0, st1, "pore")
    mt = morph(tr.start, tr.end, "restrained")
    gp = gate_path(atom_path(st0, st1, tr), tr, mt)
    ov = transition_overlap(tr, "start")
    sub = ov.subspace("A")
    n = len(mt.frames)
    keep = sorted(set(range(0, n, MORPH_KEEP)) | {n - 1})
    fractions = np.linspace(0.0, 1.0, n)
    return {
        "start": start, "end": end, "paralog": tr.paralog, "fit": tr.fit,
        "chains": tr.chains_start, "residues": tr.residues.tolist(),
        "n_residues": int(len(tr.residues)), "rmsd": round(tr.rmsd, 3),
        "fit_rmsd": round(tr.fit_rmsd, 3),
        "residue_displacement": _r(tr.residue_distance(), 2),
        "element_means": {k: round(v, 3) for k, v in tr.element_means().items()},
        "displacement_max": _P.value("display.displacement_max"),
        "morph": {"method": mt.method, "n_frames": n, "kept": keep,
                  "fraction": _r(fractions[keep]),
                  "frames": [encode_coords(tr.frame.to_frame(mt.frames[i])) for i in keep],
                  "note": MORPH_NOTE},
        "gate_path": {"fraction": _r(gp.fraction, 3), "gate": _r(gp.gate, 3),
                      "rigid_gate": _r(gp.rigid_gate, 3), "filter": _r(gp.filter, 3),
                      "lining": [list(x) for x in gp.lining],
                      "half_open": round(gp.half_open(), 3),
                      "overshoot": round(gp.overshoot(), 3)},
        "overlap": {"reference": ov.reference, "stride": ov.stride,
                    "n_sites": int(ov.modes.meta["n_sites"]),
                    "symmetry": [str(s) for s in ov.modes.symmetry],
                    "collectivity": _r(ov.modes.collectivity(), 3),
                    "overlap": _r(ov.overlap), "cumulative": _r(ov.cumulative),
                    "null_cumulative": _r(ov.null_cumulative),
                    "irrep_fraction": {k: round(v, 4) for k, v in ov.irrep_fraction.items()},
                    "rigid_fraction": round(ov.rigid_fraction, 4),
                    "lowest_collective_A": ov.lowest_a, "best": ov.best,
                    "A_subspace": {"index": sub.index.tolist(), "cumulative": _r(sub.cumulative),
                                   "null": _r(sub.null), "ceiling": round(sub.ceiling, 4)},
                    "report": ov.report()},
    }


def _profile_row(pdb_id, state, resolution, ip3, summary, control=False) -> dict:
    f = summary.constrictions.get("filter")
    g = summary.constrictions.get("gate")
    z0 = f.z if f else 0.0
    p = summary.profile
    return {"pdb_id": pdb_id, "state": state, "resolution": resolution,
            "ip3_bound": bool(ip3), "control": control,
            "gate": None if g is None else {"radius": round(g.radius, 3), "z": round(g.z - z0, 2),
                                            "residues": g.residues},
            "filter": None if f is None else {"radius": round(f.radius, 3), "residues": f.residues},
            "profile": {"z": _r(p.z - z0, 2), "r_min": _r(p.r_min, 3)}}


def states_record() -> dict:
    rows = [_profile_row(r.pdb_id, r.state, r.resolution, r.ip3_bound, r.summary)
            for r in state_panel("ITPR3")]
    e = get_entry("7T3T")
    rows.append(_profile_row("7T3T", e.state, e.resolution, e.ip3_bound,
                             measure_channel(loader.load("7T3T")), control=True))
    return {"paralog": "ITPR3", "rows": rows,
            "definition": "r_min: the minimum distance from the four-fold axis to any heavy-atom "
                          "centre (HETATM included, as ip3r_genes S0); z relative to each "
                          "deposit's filter, cytosol positive"}
