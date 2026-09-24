"""RyR1 charge mutants: the wall-charge model against measured conductances.

The ITPR3 comparison could only ask whether one number is the right size.
RyR1 has a stronger test. Xu et al. 2006 neutralised five acidic residues
one at a time (D4899Q, E4900N, D4938N, D4945N, E4955Q) and measured each
mutant's K+ conductance in the same bath as the wild type. The model makes
each mutant by dropping that residue's charge on all four subunits
(:func:`ip3r.physics.unitary.unitary` with ``neutralise``) and keeps
everything else. So the model's *ratio* mutant / wild type is compared
with the measured ratio. A ratio cancels the constants the absolute number
depends on (diffusivity, ion radius), which is what makes it a test of the
wall charge rather than of the transport constants.

The mutants are read from the registered ``permeation.published_ryr1_*``
parameters. The residue number is the name's digits, in P11716 numbering.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from ..core.structure import Structure
from ..io import loader
from ..io.registry import load_registry
from ..parameters import PARAMETERS as _P
from ..structure.channel import measure_channel
from .unitary import Unitary, unitary

__all__ = ["MutantRow", "mutants", "mutant_panel", "open_deposit",
           "open_mutant_panel"]

_PREFIX = "permeation.published_ryr1_"


@dataclass(frozen=True)
class MutantRow:
    name: str                     # e.g. "D4899Q"
    residue: int
    measured_pS: float
    measured_ratio: float         # mutant / wild type, measured
    charged_ratio: float          # model, every lining charge
    paired_ratio: float           # model, salt bridges cancelled
    lining: bool                  # the residue lines this deposit's pore
    bridged: bool                 # ... and is salt-bridged

    def row(self) -> str:
        where = ("not lining" if not self.lining else
                 "lining, salt-bridged" if self.bridged else "lining")
        return (f"{self.name:7s} measured {self.measured_pS:5.0f} pS "
                f"(x{self.measured_ratio:4.2f})   model charged "
                f"x{self.charged_ratio:4.2f}  paired x{self.paired_ratio:4.2f}"
                f"   [{where}]")


def mutants() -> dict[str, tuple[int, float]]:
    """``name -> (residue, measured pS)`` from the registered parameters."""
    out = {}
    for key in sorted(_P.as_dict(only_overrides=False)):
        if key.startswith(_PREFIX):
            name = key[len(_PREFIX):].upper()
            out[name] = (int(re.sub(r"\D", "", name)), _P.value(key))
    return out


def _ratio(mut: Unitary, wt: Unitary, attr: str) -> float:
    m, w = getattr(mut, attr), getattr(wt, attr)
    if not (m.converged and w.converged) or w.conductance_pS <= 0:
        return float("nan")
    return m.conductance_pS / w.conductance_pS


def mutant_panel(st: Structure) -> tuple[Unitary, list[MutantRow]]:
    """The wild type and every registered mutant, modelled on ``st``."""
    summary = measure_channel(st)
    wt = unitary(st, summary)
    measured_wt = _P.value("permeation.published_ryr1")
    lining = {g.res_seq for g in wt.charge.groups}
    bridged = {b.acid[1] for b in wt.paired_charge.bridged} | {
        b.base[1] for b in wt.paired_charge.bridged}
    rows = []
    for name, (res, pS) in mutants().items():
        mut = unitary(st, summary, neutralise=frozenset({res}))
        rows.append(MutantRow(name, res, pS, pS / measured_wt,
                              _ratio(mut, wt, "charged"), _ratio(mut, wt, "paired"),
                              res in lining, res in bridged))
    return wt, sorted(rows, key=lambda r: r.residue)


def open_deposit() -> str:
    """The curated RyR1 open deposit (the morph's end, by the curation rules)."""
    ids = [e.pdb_id for e in load_registry()
           if e.family == "RyR" and e.state == "open"]
    if not ids:
        raise LookupError("the RyR1 panel has no open deposit")
    return ids[0]


def open_mutant_panel() -> tuple[Unitary, list[MutantRow]]:
    return mutant_panel(loader.load(open_deposit()))
