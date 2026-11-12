"""The browser's reduced elastic network, computed by the repo's own ANM.

The web page runs ``anm_engine.js`` on a deposit's precomputed C-alpha trace
(``/api/structure/PDB``): the residues all four subunits resolve, every
``stride``-th one, springs within ``cutoff``. This module builds exactly
those sites from the same (quantised) JSON and solves them with
``physics.anm.ANM``, which is what the node test holds the JavaScript to
(``tests/fixtures/anm_reduced.json``, ``scripts/make_anm_fixture.py``).

A reduced model is an illustration, not the result: the repo's own check
(``network_checks.stride_agreement``) finds stride 2 agrees with stride 1
and stride 4 does not, so the published modes are the full ones
(``/api/modes``), and a coarse network needs a longer cutoff to stay in one
piece (the number of pieces is reported, each adds six zero modes).
"""

from __future__ import annotations

import numpy as np

from ..physics.anm import ANM
from .codec import decode_coords

__all__ = ["reduced_sites", "reduced_modes", "FIXTURE_CASES"]

#: (deposit, stride, cutoff Å, modes): the cases the fixture records.
FIXTURE_CASES = (("8TKG", 16, 25.0, 12), ("8TKF", 20, 30.0, 10))


def reduced_sites(record: dict, stride: int) -> tuple[np.ndarray, np.ndarray]:
    """``(coords (4m, 3), residues (m,))`` in subunit order, frame coordinates."""
    traces = []
    for sub in record["subunits"]:
        xyz = decode_coords(sub["xyz"])
        traces.append(dict(zip(sub["resi"], xyz)))
    shared = sorted(set.intersection(*(set(t) for t in traces)))[::max(int(stride), 1)]
    coords = np.vstack([[t[r] for r in shared] for t in traces])
    return coords, np.asarray(shared)


def reduced_modes(record: dict, stride: int, cutoff: float, n_modes: int) -> dict:
    coords, res = reduced_sites(record, stride)
    anm = ANM(coords, axis=np.array([0.0, 0.0, 1.0]), cutoff=float(cutoff))
    ms = anm.label_symmetry(anm.calc_modes(n_modes))
    return {"pdb_id": record["pdb_id"], "stride": int(stride), "cutoff": float(cutoff),
            "gamma": anm.gamma, "d0": anm.d0, "n_sites": int(len(coords)),
            "residues": res.tolist(), "n_components": anm.n_components(),
            "eigenvalues": ms.eigenvalues.tolist(), "symmetry": [str(s) for s in ms.symmetry],
            "character": ms.character.tolist(), "collectivity": ms.collectivity().tolist(),
            "vectors": [np.round(v.ravel(), 7).tolist() for v in ms.vectors]}
