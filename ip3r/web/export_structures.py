"""Precompute the web viewer's coarse-grained structures with the repo's own code.

Each deposit becomes a C-alpha trace of its four subunits in its own axis
frame (``structure.symmetry.tetramer_frame``: z on the four-fold axis, the
cytosol at +z), plus the per-residue annotation the desktop viewer paints
and the ``measure_channel`` summary. Nothing here is new science: the
numbers are the ones ``python -m ip3r info PDB`` prints.

Residue-keyed annotation (element, conservation) follows the desktop rule:
painted only when ``best_numbering`` places the deposit in a human paralog's
numbering, otherwise grey (``None`` here). Needs the mmCIF files in ``ref/``
(``python -m ip3r fetch``); the output is small JSON committed under
``ip3r/web/precomputed/``.
"""

from __future__ import annotations

import hashlib

import numpy as np

from ..core.annotations import (ELEMENT_COLORS, ELEMENT_LABELS, ELEMENT_ORDER,
                                constraint_at, element_array)
from ..io import loader
from ..io.registry import get_entry, local_path
from ..render.colormaps import CHAIN_PALETTE, CONSERVATION_RANGE, MISSING
from ..structure.channel import measure_channel
from ..structure.ligand import ligand_sites
from .codec import encode_coords

__all__ = ["VIEWER_DEPOSITS", "structure_record", "palette_record", "sha256"]

#: The deposits the web viewer offers: ip3r_genes S11's ITPR3 state panel
#: (the registry's ITPR3 entries) and the open-state control 7T3T.
VIEWER_DEPOSITS = ("6DQN", "6DQJ", "7T3P", "8TKG", "8TKH", "8TLA", "8TKF", "7T3T")


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _rgb(c) -> list[float]:
    return [round(float(x), 3) for x in c]


def palette_record() -> dict:
    """The desktop viewer's colour tables, so the web paints the same colours."""
    return {
        "elements": [{"name": n, "label": ELEMENT_LABELS[n], "rgb": _rgb(ELEMENT_COLORS[n])}
                     for n in ELEMENT_ORDER if n in ELEMENT_COLORS],
        "chains": [_rgb(c) for c in CHAIN_PALETTE[:4]],
        "missing": _rgb(MISSING),
        "conservation_range": list(CONSERVATION_RANGE),
        "conservation_layer": "deep",
    }


def _subunit_trace(st, chain: str, frame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m = st.mask_ca() & (st.chain == chain)
    res = st.res_seq[m].astype(int)
    order = np.argsort(res, kind="stable")
    res = res[order]
    keep = np.r_[True, np.diff(res) != 0]          # one CA per residue (alt locs)
    xyz = frame.to_frame(st.xyz[m][order][keep])
    return res[keep], xyz, st.b_factor[m][order][keep].astype(float)


def structure_record(pdb_id: str) -> dict:
    """The viewer's JSON for one deposit (see the module docstring)."""
    st = loader.load(pdb_id)
    entry = get_entry(pdb_id)
    s = measure_channel(st)
    frame = s.frame
    paralog = s.numbering.paralog if (s.numbering and s.numbering.passed) else None
    elem_index = {n: i for i, n in enumerate(n for n in ELEMENT_ORDER if n in ELEMENT_COLORS)}
    lo, hi = CONSERVATION_RANGE
    subunits = []
    for k, ch in enumerate(frame.chains):
        res, xyz, bfac = _subunit_trace(st, ch, frame)
        rec = {"chain": ch, "palette_index": k, "resi": res.tolist(),
               "xyz": encode_coords(xyz), "bfactor": [round(b, 1) for b in bfac]}
        if paralog:
            names = element_array(paralog, res)
            rec["element"] = [elem_index.get(n, -1) for n in names]
            jsd = constraint_at(paralog, res)
            rec["jsd"] = [None if not np.isfinite(v) else round(float(v), 3) for v in jsd]
        subunits.append(rec)
    contacts = {}
    for sub, c in s.ip3_contacts.items():
        contacts[sub] = {"same_subunit": [f"{n}{r}" for n, r in c["same_subunit"]],
                         "other_subunit": [f"{n}{r}" for n, r in c["other_subunit"]]}
    ligands = [{"subunit": site.subunit, "xyz": encode_coords(frame.to_frame(st.xyz[site.atoms]))}
               for site in ligand_sites(st)]
    return {
        "pdb_id": pdb_id,
        "title": entry.title, "paralog": entry.paralog, "organism": entry.organism,
        "state": entry.state, "resolution": entry.resolution,
        "is_control": entry.is_control,
        "numbering": None if s.numbering is None else {
            "paralog": s.numbering.paralog, "identity": round(s.numbering.identity, 4),
            "n_compared": s.numbering.n_compared, "passed": bool(s.numbering.passed)},
        "painted_in": paralog,
        "chains": list(frame.chains),
        "n_atoms": int(st.n_atoms),
        "axis_disagreement_deg": round(s.axis_disagreement_deg, 3),
        "c4_residual": round(float(s.c4_residual), 3),
        "span": [round(x, 2) for x in s.span], "span_source": s.span_source,
        "constrictions": {n: {"z": round(c.z, 2), "radius": round(c.radius, 3),
                              "residues": c.residues} for n, c in s.constrictions.items()},
        "ip3_contacts": contacts,
        "ligands": ligands,
        "subunits": subunits,
        "notes": s.notes,
        "source": {"file": local_path(pdb_id).name, "sha256": sha256(local_path(pdb_id)),
                   "url": f"https://www.rcsb.org/structure/{pdb_id}"},
    }
