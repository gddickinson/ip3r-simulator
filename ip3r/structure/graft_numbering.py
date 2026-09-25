"""Which prediction residue fills which deposit residue.

Two routes, tried in this order:

* **By number.** The prediction and the deposit are in the same numbering:
  the residue with number *n* in one is residue *n* in the other. The rule
  is the one the variant painting uses (``numbering.min_identity`` over the
  unstubbed residues). All seven ITPR3 deposits and ``AF-Q14573-F1``.
* **By alignment.** AlphaFold DB models only isoform 8 of rat ITPR1, while
  7LHF is numbered in the canonical sequence. The two differ by splice
  segments, so by number they agree on 22 % of residues. Here the deposit's
  whole construct (:mod:`ip3r.io.poly_seq`, the unresolved residues
  included) is aligned to the prediction (:func:`ip3r.core.pairwise.align`),
  and the aligned pairs are the map. The prediction must agree with the
  construct over the aligned pairs to ``graft.align_min_identity``. That
  bar separates the same protein in another isoform from an ortholog or a
  paralog.

A deposit residue with no partner lies in a segment the prediction does
not have (a splice segment the isoform lacks). It is never filled. Such
segments are kept in :attr:`NumberMap.unmapped`, and each stretch that
touches one is skipped with its reason.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

from ..core.structure import AA3TO1, Structure
from ..parameters import PARAMETERS as _P

__all__ = ["NumberMap", "BY_NUMBER", "by_number", "by_alignment", "runs"]


@dataclass(frozen=True)
class NumberMap:
    route: str                          # "number" | "alignment"
    identity: float
    to_pred: dict | None = None         # None = the identity map
    unmapped: tuple = ()                # ((first, last), ...) in deposit numbers

    def __call__(self, r: int) -> int | None:
        """The prediction's residue number for deposit residue ``r``, or None."""
        return r if self.to_pred is None else self.to_pred.get(r)

    def span(self) -> tuple[int, int] | None:
        """The deposit numbers the prediction covers, first and last (alignment only)."""
        if not self.to_pred:
            return None
        return min(self.to_pred), max(self.to_pred)

    def lacking(self, first: int, last: int) -> list[tuple[int, int]]:
        """Unmapped segments that overlap deposit residues ``first..last``."""
        return [(a, b) for a, b in self.unmapped if a <= last and b >= first]

    def describe(self) -> str:
        how = "by number" if self.route == "number" else "by alignment"
        text = f"{self.identity:.1%} identical {how}"
        if self.unmapped:
            segs = ", ".join(f"{a}-{b}" for a, b in self.unmapped)
            text += f"; the model lacks {segs}"
        return text


BY_NUMBER = NumberMap("number", 1.0)


def runs(numbers) -> list[tuple[int, int]]:
    """Consecutive runs in sorted integers, as ``(first, last)``."""
    out = []
    for _, g in itertools.groupby(enumerate(sorted(numbers)), lambda x: x[1] - x[0]):
        g = [x[1] for x in g]
        out.append((g[0], g[-1]))
    return out


def by_number(dep: dict, pred: dict) -> NumberMap:
    """The identity map, scored over the deposit's unstubbed residues."""
    scored = [r for r, (aa, stub) in dep.items() if not stub and r in pred]
    same = sum(dep[r][0] == pred[r][0] for r in scored)
    return NumberMap("number", same / len(scored) if scored else 0.0)


def by_alignment(construct: dict[int, str], pred: dict) -> NumberMap:
    """Map a deposit's construct onto a prediction by global alignment.

    ``construct`` is ``{author number: three-letter name}``, and ``pred`` is
    ``chain_residues`` of the prediction. Identity is taken over the aligned
    pairs. Unmapped segments are the construct residues inside the aligned
    span that have no partner. Construct residues outside that span are
    beyond the prediction's ends, not missing from its middle.
    """
    from ..core.pairwise import align

    dn = sorted(construct)
    pn = sorted(pred)
    if not dn or not pn:
        return NumberMap("alignment", 0.0, {})
    dseq = "".join(AA3TO1.get(construct[n], "X") for n in dn)
    pseq = "".join(pred[n][0] for n in pn)
    _, pairs = align(dseq, pseq)
    if not pairs:
        return NumberMap("alignment", 0.0, {})
    to_pred = {dn[i - 1]: pn[j - 1] for i, j in pairs}
    same = sum(dseq[i - 1] == pseq[j - 1] for i, j in pairs)
    lo, hi = min(to_pred), max(to_pred)
    missing = [n for n in dn if lo <= n <= hi and n not in to_pred]
    return NumberMap("alignment", same / len(pairs), to_pred, tuple(runs(missing)))


def construct_of(st: Structure, chain: str) -> dict[int, str]:
    """The deposit chain's construct sequence, from its source file."""
    from ..io.poly_seq import construct_sequence
    return construct_sequence(st.source, chain) if st.source else {}


def align_bar() -> float:
    return _P.value("graft.align_min_identity")
