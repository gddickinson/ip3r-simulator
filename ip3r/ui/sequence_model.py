"""What the sequence window shows, Qt-free (Round 7.18).

One chain of the loaded deposit as the *construct*: every residue the
mmCIF's poly-seq scheme lists, resolved or not, in the deposit's own author
numbers (the numbers a click on the model reports). Unresolved residues are
kept, dimmed, because a gap is information: an alignment by position that
dropped them would put every later residue on the wrong number.

Tracks paint annotation keyed by residue number only when the deposit is
in that paralog's numbering (``structure.numbering``), exactly as the
viewer's colourings do; otherwise they are grey and say why.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.annotations import (ELEMENT_COLORS, ELEMENT_LABELS, constraint_at,
                                element_array, functional_sites)
from ..core.structure import AA3TO1
from ..io.poly_seq import construct_sequence
from ..render.colormaps import constraint_colors

__all__ = ["ChainSequence", "chain_sequence", "TRACKS", "track", "Decoration",
           "SITE_COLOURS"]

#: Track key -> label.
TRACKS = {"chemistry": "chemistry (letters only)",
          "element": "functional element",
          "constraint": "deep JSD (fixed 0.50–0.95)",
          "resolved": "resolved in this deposit"}

#: Site classes underlined in every track, and their colour.
SITE_COLOURS = {"ip3_contact": "#ff5fd2", "filter_lining": "#ff7a4d",
                "gate_lining": "#ffd23f"}


@dataclass(frozen=True)
class ChainSequence:
    chain: str
    positions: tuple[int, ...]      # author residue numbers
    letters: str
    resolved: tuple[bool, ...]
    #: "construct" (poly-seq scheme) or "built" (the atoms alone).
    source: str

    def __len__(self) -> int:
        return len(self.positions)

    @property
    def n_resolved(self) -> int:
        return int(sum(self.resolved))


@dataclass(frozen=True)
class Decoration:
    background: str = ""
    underline: str = ""
    tooltip: str = ""


def chain_sequence(st, chain: str) -> ChainSequence:
    """The construct of ``chain`` (resolved or not); the built residues when
    the file carries no poly-seq scheme."""
    built_mask = (st.chain == chain) & ~st.hetero
    built = {int(r): str(n) for r, n in zip(st.res_seq[built_mask], st.res_name[built_mask])}
    construct = construct_sequence(st.source, chain) if st.source else {}
    source = "construct" if construct else "built"
    names = construct or built
    pos = tuple(sorted(names))
    return ChainSequence(chain, pos, "".join(AA3TO1.get(names[p], "X") for p in pos),
                         tuple(p in built for p in pos), source)


def _hex(rgb) -> str:
    return "#" + "".join(f"{int(round(255 * float(c))):02x}" for c in rgb[:3])


def track(seq: ChainSequence, key: str, paralog: str | None) -> dict[int, Decoration]:
    """Per-residue decoration for track ``key``. Residue-keyed tracks are
    empty (grey) unless ``paralog`` names the deposit's numbering."""
    pos = np.array(seq.positions, int)
    sites: dict[int, str] = {}
    if paralog is not None:
        for cls in SITE_COLOURS:
            for r in functional_sites(paralog).get(cls, ()):
                sites.setdefault(int(r), cls)
    out: dict[int, Decoration] = {}
    bg = [""] * len(pos)
    tips = [""] * len(pos)
    if key == "element" and paralog is not None:
        names = element_array(paralog, pos)
        for i, n in enumerate(names):
            if n:
                bg[i], tips[i] = _hex(ELEMENT_COLORS[n]), ELEMENT_LABELS.get(n, n)
    elif key == "constraint" and paralog is not None:
        jsd = constraint_at(paralog, pos, "deep")
        rgb = constraint_colors(jsd)
        for i, v in enumerate(jsd):
            if np.isfinite(v):
                bg[i], tips[i] = _hex(rgb[i]), f"deep JSD {v:.3f}"
            else:
                tips[i] = "deep JSD not scored"
    elif key == "resolved":
        for i, ok in enumerate(seq.resolved):
            bg[i] = "#4f8a5b" if ok else "#6b3a3a"
    for i, p in enumerate(seq.positions):
        cls = sites.get(p, "")
        state = "" if seq.resolved[i] else "unresolved"
        tip = "; ".join(t for t in (f"{seq.letters[i]}{p}", tips[i], cls, state) if t)
        out[p] = Decoration(bg[i], SITE_COLOURS.get(cls, ""), tip)
    return out
