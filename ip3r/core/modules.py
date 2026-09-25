"""Paper 6's two functional modules, defined from this project's annotation.

The **ligand core** is the smallest contiguous span holding all ten measured
IP3 contacts; the **pore module** is the Pfam channel domain (PF00520) less
the luminal loop, the least conserved element in the receptor. A second pore
definition keeps the loop in (``channel_all``): the paired contrast reverses
between the two, so both are shown. A second core definition
(``ibc_literature``) is the published IP3-binding core, ITPR1 224-604
(Bosanac 2002; ``ligand.ibc_start``/``ligand.ibc_end``), carried to the other
paralogs by this project's own pairwise alignment.

The definitions are the ones S22 states in prose; the spans are rebuilt here
from ``functional_sites`` and ``elements`` (both imported from ip3r_genes
with source hashes) rather than read from S22's ``module_map.tsv``, so the
``P6.module_map`` check compares two constructions of the same rule.

Each module is validated against what it must and must not contain (all ten
contacts and no pore site; both filter and both gate residues and no
contact). A module that fails raises :class:`ModuleRefusal` rather than being
drawn wrong.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from ..parameters import PARAMETERS as _P
from .annotations import elements, functional_sites
from .pairwise import paralog_transfer

__all__ = ["Module", "ModuleRefusal", "module", "modules", "PORE_DEFINITIONS",
           "CORE_DEFINITIONS", "IBC_PARALOG", "MODULE_COLORS", "MODULE_LABELS",
           "MODULES_KEY"]

#: The key the GUI's site toggles use for "both primary modules".
MODULES_KEY = "p6_modules"

#: ``definition -> (module, excludes the luminal loop)``.
PORE_DEFINITIONS = {"channel_minus_luminal": True, "channel_all": False}
CORE_DEFINITIONS = ("contact_span", "ibc_literature")
#: The paralog whose numbering ``ligand.ibc_start``/``ligand.ibc_end`` are in.
IBC_PARALOG = "ITPR1"

MODULE_LABELS = {"contact_span": "Ligand core (span of the ten IP3 contacts)",
                 "channel_minus_luminal": "Pore module (PF00520 less the luminal loop)",
                 "channel_all": "PF00520 including the luminal loop",
                 "ibc_literature": "Published IP3-binding core (ITPR1 224-604)"}
# Outside the blue -> red conservation ramp, so the trace reads over it.
MODULE_COLORS = {"contact_span": (0.30, 0.95, 0.40),
                 "channel_minus_luminal": (0.95, 0.35, 0.95),
                 "channel_all": (0.80, 0.50, 0.85),
                 "ibc_literature": (0.45, 0.80, 0.45)}


class ModuleRefusal(RuntimeError):
    """A module definition failed its own content check."""


@dataclass(frozen=True)
class Module:
    paralog: str
    name: str                 # ligand_core | pore_module
    definition: str           # contact_span | channel_minus_luminal | channel_all
    start: int
    end: int
    excluded: tuple = ()      # ((lo, hi), ...) holes inside [start, end]

    @property
    def residues(self) -> tuple[int, ...]:
        return tuple(r for r in range(self.start, self.end + 1)
                     if not any(a <= r <= b for a, b in self.excluded))

    def __contains__(self, resi: int) -> bool:
        return (self.start <= resi <= self.end
                and not any(a <= resi <= b for a, b in self.excluded))


def _element(paralog: str, name: str):
    found = [e for e in elements(paralog) if e.name == name]
    if len(found) != 1:
        raise ModuleRefusal(f"{paralog}: expected one {name} element, "
                            f"found {len(found)}")
    return found[0]


def _validate(m: Module, sites: dict) -> Module:
    rs = set(m.residues)
    n = {k: len(rs & set(sites.get(k, ()))) for k in
         ("ip3_contact", "filter_lining", "gate_lining")}
    want = ({"ip3_contact": 10, "filter_lining": 0, "gate_lining": 0}
            if m.name == "ligand_core" else
            {"ip3_contact": 0, "filter_lining": 2, "gate_lining": 2})
    if n != want:
        raise ModuleRefusal(f"{m.paralog}/{m.definition}: holds {n}, "
                            f"expected {want}")
    return m


def _ibc_span(paralog: str) -> tuple[int, int]:
    lo, hi = int(_P.value("ligand.ibc_start")), int(_P.value("ligand.ibc_end"))
    if paralog == IBC_PARALOG:
        return lo, hi
    move = paralog_transfer(IBC_PARALOG, paralog)
    if lo not in move or hi not in move:
        raise ModuleRefusal(f"{paralog}: the IBC ends {lo}/{hi} do not align")
    return move[lo], move[hi]


@lru_cache(maxsize=None)
def module(paralog: str, definition: str) -> Module:
    """One module definition for one paralog, validated."""
    sites = functional_sites(paralog)
    if definition == "contact_span":
        c = sites.get("ip3_contact", ())
        if len(c) != 10:
            raise ModuleRefusal(f"{paralog}: {len(c)} IP3 contacts, expected 10")
        m = Module(paralog, "ligand_core", definition, min(c), max(c))
    elif definition == "ibc_literature":
        m = Module(paralog, "ligand_core", definition, *_ibc_span(paralog))
    elif definition in PORE_DEFINITIONS:
        chan = _element(paralog, "channel")
        holes = ()
        if PORE_DEFINITIONS[definition]:
            loop = _element(paralog, "luminal_loop")
            holes = ((loop.start, loop.end),)
        m = Module(paralog, "pore_module", definition, chan.start, chan.end,
                   holes)
    else:
        raise KeyError(f"unknown module definition {definition!r}")
    return _validate(m, sites)


def modules(paralog: str) -> tuple[Module, Module]:
    """The primary pair: (ligand core, pore module less the luminal loop)."""
    core = module(paralog, "contact_span")
    pore = module(paralog, "channel_minus_luminal")
    if set(core.residues) & set(pore.residues):
        raise ModuleRefusal(f"{paralog}: primary modules overlap")
    return core, pore


def _parameters_changed() -> None:
    """The IBC ends and the alignment's gap costs are registered."""
    module.cache_clear()


_P.subscribe(_parameters_changed)
