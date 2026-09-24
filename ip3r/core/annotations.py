"""Per-paralog annotation: reference sequence, elements and functional sites.

All three are read from the resources ``scripts/sync_genes.py`` imports from
``ip3r_genes``, so a residue number here is the canonical human UniProt
number the S17 constraint tables use — no conversion happens in this module.

**Every residue has exactly one element.** Pfam domains and the structural
elements overlap (the filter and gate sit inside PF00520, the luminal loop
too), and a per-residue colour needs a single answer. The rule is the one S17
applied: the most specific element wins (filter, gate, luminal loop, then the
Pfam domain, then the named linker between two domains).

**Ryanodine receptors** (``config.RYR_ACC``) are a numbering, not a
paralog of the publication: their sequence and Pfam domains come from
``ryr1.json`` (``scripts/curate_ryr.py``). They have no conservation, sites
or variants, so those return "not scored" (NaN, empty), which paints grey.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from ..config import NUMBERINGS, RESOURCE_DIR, RYR_ACC

__all__ = ["Element", "is_ryr", "reference_sequence", "elements", "element_of",
           "residue_elements", "functional_sites", "ELEMENT_ORDER",
           "ELEMENT_COLORS", "ELEMENT_LABELS", "LAYERS", "LAYER_LABELS",
           "residue_constraint", "constraint_at", "variants", "element_array"]

#: Specificity order: earlier wins where elements overlap.
_SPECIFIC_FIRST = ("selectivity_filter", "gate", "luminal_loop")
#: Elements that are annotations over a span, not a place — never painted.
_NOT_PAINTED = {"ip3_contact_set"}

#: N- to C-terminal order, for legends and plots.
ELEMENT_ORDER = ("nterm_trefoil", "MIR", "RIH_N", "RIH_C", "RIH_assoc",
                 "channel", "luminal_loop", "selectivity_filter", "gate",
                 "linker", "RIH", "SPRY", "RyR_repeat", "junctional_solenoid",
                 "RyR_TM4_6")

ELEMENT_LABELS = {
    "nterm_trefoil": "IP3-binding core, β-trefoil (PF08709)",
    "MIR": "MIR domain (PF02815)",
    "RIH_N": "RIH domain, N (PF01365)",
    "RIH_C": "RIH domain, C (PF01365)",
    "RIH_assoc": "RIH-associated (PF08454)",
    "channel": "Pore domain (PF00520)",
    "luminal_loop": "Luminal loop (6DQN geometry)",
    "selectivity_filter": "Selectivity filter (GGGVGD)",
    "gate": "Gate (6DQN)",
    "linker": "Inter-domain linker",
    "RIH": "RIH domain (PF01365; RyR)",
    "SPRY": "SPRY domain (PF00622; RyR)",
    "RyR_repeat": "RyR repeat (PF02026; RyR)",
    "junctional_solenoid": "Junctional solenoid (PF21119; RyR)",
    "RyR_TM4_6": "RyR TM 4-6 (PF06459; RyR)",
}

#: A categorical palette, readable on the dark viewport. The gate and filter
#: take the warm accents because they are what the constraint result is about.
ELEMENT_COLORS = {
    "nterm_trefoil": (0.36, 0.62, 0.95),
    "MIR": (0.30, 0.78, 0.86),
    "RIH_N": (0.45, 0.80, 0.50),
    "RIH_C": (0.62, 0.84, 0.36),
    "RIH_assoc": (0.80, 0.75, 0.35),
    "channel": (0.72, 0.52, 0.92),
    "luminal_loop": (0.60, 0.60, 0.66),
    "selectivity_filter": (1.00, 0.45, 0.30),
    "gate": (1.00, 0.80, 0.20),
    "linker": (0.42, 0.45, 0.52),
    "RIH": (0.45, 0.80, 0.50),
    "SPRY": (0.95, 0.55, 0.75),
    "RyR_repeat": (0.55, 0.70, 0.95),
    "junctional_solenoid": (0.85, 0.65, 0.45),
    "RyR_TM4_6": (0.50, 0.55, 0.85),
}


@dataclass(frozen=True)
class Element:
    name: str
    start: int
    end: int
    kind: str
    source: str

    def __contains__(self, resi: int) -> bool:
        return self.start <= resi <= self.end


def _load(name: str) -> dict:
    return json.loads((RESOURCE_DIR / name).read_text())


def _check(paralog: str, allowed: tuple = NUMBERINGS) -> str:
    if paralog not in allowed:
        raise KeyError(f"unknown paralog {paralog!r}; expected one of "
                       f"{allowed}")
    return paralog


def is_ryr(paralog: str | None) -> bool:
    return paralog in RYR_ACC


@lru_cache(maxsize=None)
def reference_sequence(paralog: str) -> str:
    """Reference sequence (canonical human, or rabbit RyR1), residue 1 at index 0."""
    if is_ryr(_check(paralog)):
        return _load("ryr1.json")["paralogs"][paralog]["sequence"]
    return _load("sequences.json")["paralogs"][paralog]["sequence"]


@lru_cache(maxsize=None)
def elements(paralog: str) -> tuple[Element, ...]:
    rows = (_load("ryr1.json")["paralogs"][paralog]["domains"]
            if is_ryr(_check(paralog))
            else _load("domains.json")["paralogs"][paralog])
    return tuple(Element(e["element"], e["start"], e["end"], e["kind"],
                         e["source"]) for e in rows)


@lru_cache(maxsize=None)
def residue_elements(paralog: str) -> tuple[str, ...]:
    """One element name per residue (index 0 is residue 1)."""
    n = len(reference_sequence(paralog))
    out = ["linker"] * n
    els = [e for e in elements(paralog) if e.name not in _NOT_PAINTED]
    general = [e for e in els if e.name not in _SPECIFIC_FIRST]
    # Specific elements painted last, most specific last of all, so it wins.
    specific = sorted((e for e in els if e.name in _SPECIFIC_FIRST),
                      key=lambda e: -_SPECIFIC_FIRST.index(e.name))
    for e in general + specific:
        for r in range(max(e.start, 1), min(e.end, n) + 1):
            out[r - 1] = e.name
    return tuple(out)


def element_of(paralog: str, resi: int) -> str | None:
    """The element residue ``resi`` belongs to, or None if out of range."""
    table = residue_elements(paralog)
    return table[resi - 1] if 1 <= resi <= len(table) else None


@lru_cache(maxsize=None)
def functional_sites(paralog: str) -> dict[str, tuple[int, ...]]:
    """``site_class -> residue numbers`` (ip3_contact, filter_lining, ...)."""
    out: dict[str, list[int]] = {}
    if is_ryr(_check(paralog)):
        return {}
    for s in _load("sites.json")["paralogs"][paralog]:
        out.setdefault(s["site_class"], []).append(s["resi"])
    return {k: tuple(sorted(v)) for k, v in out.items()}


def element_array(paralog: str, resi: np.ndarray) -> np.ndarray:
    """Vectorised :func:`element_of`; ``""`` where out of range."""
    table = np.array(residue_elements(paralog) + ("",), dtype="U20")
    idx = np.where((resi >= 1) & (resi <= len(table) - 1), resi - 1,
                   len(table) - 1)
    return table[idx]


#: The four S17 conservation layers, most to least taxonomically local.
LAYERS = ("deep", "shallow", "vert", "family")
LAYER_LABELS = {
    "deep": "deep — 249-265 orthologues of this paralog",
    "shallow": "shallow — msa_v2 tips of this paralog (control)",
    "vert": "vertebrate — all three paralogs",
    "family": "family — every eukaryotic lineage",
}


@lru_cache(maxsize=None)
def residue_constraint(paralog: str, layer: str = "deep") -> np.ndarray:
    """Per-residue JSD (index 0 = residue 1); NaN where not scored."""
    if is_ryr(_check(paralog)):
        return np.full(len(reference_sequence(paralog)), np.nan)
    vals = _load("constraint.json")["paralogs"][paralog][layer]
    return np.array([np.nan if v is None else v for v in vals], dtype=float)


def constraint_at(paralog: str, resi: np.ndarray, layer: str = "deep") -> np.ndarray:
    """Vectorised lookup; NaN outside the sequence."""
    table = np.append(residue_constraint(paralog, layer), np.nan)
    resi = np.asarray(resi, int)
    idx = np.where((resi >= 1) & (resi <= len(table) - 1), resi - 1, len(table) - 1)
    return table[idx]


@lru_cache(maxsize=None)
def variants(paralog: str | None = None) -> tuple[dict, ...]:
    """The S17 variant harvest (optionally one paralog), as dicts."""
    rows = _load("variants.json")["variants"]
    return tuple(r for r in rows if paralog is None or r["gene"] == paralog)
