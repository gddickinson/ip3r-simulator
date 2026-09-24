"""Which reference numbering a deposited structure is in.

A constraint value or a ClinVar variant is keyed on a canonical human UniProt
residue number. Painting it onto a structure is only valid if the structure's
residue ``n`` *is* that residue ``n``. The test is the one ``ip3r_genes`` S24
applied (D64): among the residues the structure resolves, the fraction whose
amino acid matches the reference at the same number must reach
``numbering.min_identity``. An offset still lines up most residues with
*something*, so the bar is deliberately strict — S24 found the rat ITPR1
reference (7LHF) and an AlphaFold isoform both fail it.

**Stubbed residues are not evidence.** 6DQN carries a 113-residue segment
(1434-1546 on every chain) built as backbone plus CB, with residue names that
do not match human ITPR3 at those numbers — a region the map could not assign
a sequence to. Scored naively it drops the deposit to 94.8 % and fails a
numbering it is actually in. A residue whose heavy atoms stop at CB (and that
is not named Gly or Ala) says nothing about numbering, so it is excluded and
counted, and every run of disagreement is reported as a segment, because a
local mis-assignment and a global offset are different findings.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import PARALOGS
from ..core.annotations import reference_sequence
from ..core.structure import AA3TO1, Structure
from ..parameters import PARAMETERS as _P

__all__ = ["NumberingCheck", "check_numbering", "best_numbering", "chain_residues"]


_BACKBONE_CB = {"N", "CA", "C", "O", "OXT", "CB"}


@dataclass
class NumberingCheck:
    paralog: str
    chain: str
    n_compared: int
    n_identical: int
    n_stubbed: int = 0
    #: Runs of >= 5 consecutive disagreeing residues, as (first, last).
    mismatch_segments: tuple = ()

    @property
    def identity(self) -> float:
        return self.n_identical / self.n_compared if self.n_compared else 0.0

    @property
    def passed(self) -> bool:
        return self.identity >= _P.value("numbering.min_identity")


def chain_residues(st: Structure, chain: str) -> dict[int, tuple[str, bool]]:
    """``resnum -> (one-letter, stubbed)`` for the protein residues of a chain."""
    m = (st.chain == chain) & st.mask_protein() & ~st.hetero & (st.element != "H")
    atoms: dict[int, set] = {}
    names: dict[int, str] = {}
    for r, n, a in zip(st.res_seq[m], st.res_name[m], st.atom_name[m]):
        atoms.setdefault(int(r), set()).add(str(a))
        names[int(r)] = str(n)
    return {r: (AA3TO1.get(names[r], "X"),
                atoms[r] <= _BACKBONE_CB and names[r] not in ("GLY", "ALA"))
            for r in sorted(atoms) if "CA" in atoms[r]}


def _segments(bad: list[int], min_len: int = 5) -> tuple:
    out, start, prev = [], None, None
    for r in bad + [None]:
        if r is not None and prev is not None and r == prev + 1:
            prev = r
            continue
        if start is not None and prev - start + 1 >= min_len:
            out.append((start, prev))
        start = prev = r
    return tuple(out)


def check_numbering(st: Structure, paralog: str, chain: str | None = None) -> NumberingCheck:
    chain = chain or max(st.chains, key=lambda c: int((st.mask_ca() & (st.chain == c)).sum()))
    seq = reference_sequence(paralog)
    res = chain_residues(st, chain)
    inside = {r: v for r, v in res.items() if 1 <= r <= len(seq)}
    bad = sorted(r for r, (aa, _) in inside.items() if aa != seq[r - 1])
    scored = {r: aa for r, (aa, stub) in inside.items() if not stub}
    same = int(np.sum([aa == seq[r - 1] for r, aa in scored.items()]))
    return NumberingCheck(paralog, chain, len(scored), same,
                          n_stubbed=len(inside) - len(scored),
                          mismatch_segments=_segments(bad))


def best_numbering(st: Structure) -> NumberingCheck | None:
    """The paralog whose numbering the structure is in, or None if none fits."""
    checks = [check_numbering(st, p) for p in PARALOGS]
    best = max(checks, key=lambda c: c.identity)
    return best if best.passed else None
