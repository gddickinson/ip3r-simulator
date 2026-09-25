"""Paper 1's family-call benchmark (S1), measured again from the sequences.

S1 benchmarked the discovery scorer on 25 positive controls (the three
vertebrate paralogues across species, plus fly, worm and *Dictyostelium*)
and 31 decoys, six of them ryanodine receptors. Two things are rebuilt here.

**The labelled-bait margin** (``margins``): each control's identity to the
nearest human IP3 receptor bait minus its identity to the nearest human
ryanodine receptor bait, a bait never scored against itself. S1 read
identity off one 56-sequence MAFFT alignment. Here every pair is aligned
separately with this project's Gotoh/BLOSUM62 aligner
(:mod:`ip3r.core.pairwise`), from the UniProt sequences S1 committed.
"Full" identity is identical columns over columns where either sequence
has a residue. That is S1's full-alignment metric, the one its scorer uses.
"Covered" identity counts only columns where both have a residue. A pairwise
optimum aligns unrelated long sequences at 0.25-0.28 covered identity, so
for twilight pairs the covered reading is not a re-measurement. It is
reported, and the verdict rests on the full metric.

**The scores** (``scores``): each control's score rebuilt from the
components S1 records as fired, with the registered ``bench.points_*``.
The two caps are decided here. The sister cap comes from the margin
measured above, and the evidence gate from which family-specific
components fired. S1's own ``SISTER``/``GATED`` flags are only compared.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache

from ..core import genes_data as G
from ..core.pairwise import align
from ..parameters import PARAMETERS as _P

__all__ = ["POSITIVES", "DECOYS", "MARGIN_TSV", "POSITIVE_TSV", "NEGATIVE_TSV",
           "SUMMARY", "ITPR_BAITS", "RYR_BAITS", "sequences", "identity",
           "Margin", "margins", "COMPONENTS", "FAMILY_SPECIFIC", "Score", "rescore",
           "clear_caches"]

_DIR = "benchmark_controls/"
POSITIVES, DECOYS = _DIR + "panel_positives.json", _DIR + "panel_decoys.json"
MARGIN_TSV, SUMMARY = _DIR + "bait_margin.tsv", _DIR + "summary.json"
POSITIVE_TSV, NEGATIVE_TSV = _DIR + "positive_controls.tsv", _DIR + "negative_controls.tsv"

#: The labelled baits: the human paralogues of each family (UniProt).
ITPR_BAITS = ("Q14643", "Q14571", "Q14573")
RYR_BAITS = ("P21817", "Q92736", "Q15413")

#: component tag in S1's tables → registered points key
COMPONENTS = {"outlier": "outlier", "fold": "fold", "pfam": "pfam", "sig": "signature",
              "size": "size", "breadth": "breadth", "cluster": "cluster",
              "homology": "homology", "split": "split"}
#: The components that can carry a candidate past the evidence gate.
FAMILY_SPECIFIC = frozenset(("outlier", "fold", "pfam", "sig", "split"))


@lru_cache(maxsize=None)
def sequences() -> dict[str, str]:
    """accession → sequence, over both panels."""
    out = {}
    for path in (POSITIVES, DECOYS):
        for r in json.loads(G.read_text(path)):
            out[r["accession"]] = r["sequence"]
    return out


@lru_cache(maxsize=None)
def _pair(a: str, b: str, go: float, ge: float) -> tuple[int, int]:
    """(identical columns, aligned pairs); keyed on content, not on path."""
    _, pairs = align(a, b, go, ge)
    return sum(a[i - 1] == b[j - 1] for i, j in pairs), len(pairs)


def identity(a: str, b: str) -> tuple[float, float]:
    """(full, covered) identity of two sequences under the registered gaps."""
    same, n = _pair(a, b, _P.value("align.gap_open"), _P.value("align.gap_extend"))
    return same / (len(a) + len(b) - n), (same / n if n else 0.0)


@dataclass(frozen=True)
class Margin:
    accession: str
    itpr: tuple        # (full, covered) identity to the nearest ITPR bait
    ryr: tuple         # the same to the nearest RyR bait

    @property
    def full(self) -> float:
        return self.itpr[0] - self.ryr[0]

    @property
    def covered(self) -> float:
        return self.itpr[1] - self.ryr[1]

    def call(self, metric: str = "full") -> str:
        m = self.full if metric == "full" else self.covered
        return "ITPR" if m > 0 else "RYR" if m < 0 else "tie"

    def band(self, metric: str = "full") -> str:
        """``ITPR`` / ``RYR`` beyond the D7 margin, else ``no call``."""
        m, d = (self.full if metric == "full" else self.covered), _P.value("bench.sister_margin")
        return "ITPR" if m > d else "RYR" if m < -d else "no call"


def _nearest(acc: str, baits) -> tuple[float, float]:
    seq = sequences()
    ids = [identity(seq[acc], seq[b]) for b in baits if b != acc]
    return max(i[0] for i in ids), max(i[1] for i in ids)


def margins() -> dict[str, Margin]:
    """Every control's labelled-bait margin, measured here."""
    return {acc: Margin(acc, _nearest(acc, ITPR_BAITS), _nearest(acc, RYR_BAITS))
            for acc in sequences()}


@dataclass(frozen=True)
class Score:
    accession: str
    raw: int
    score: int
    cap: str           # "" | "sister" | "gate"


def rescore(components: str, margin: Margin) -> Score:
    """A control's score from S1's fired components, with both caps decided
    here: the sister cap from ``margin`` (closer to the RyR baits by more
    than the D7 margin), the gate from the family-specific components."""
    fired = [c.split("+")[0] for c in components.split(",") if "+" in c]
    raw = min(sum(int(_P.value(f"bench.points_{COMPONENTS[c]}")) for c in fired),
              int(_P.value("bench.score_cap")))
    bar = int(_P.value("bench.promotion"))
    cap = ""
    if raw >= bar and margin.band("full") == "RYR":
        cap = "sister"
    elif raw >= bar and not FAMILY_SPECIFIC & set(fired):
        cap = "gate"
    return Score(margin.accession, raw, bar - 1 if cap else raw, cap)


def clear_caches() -> None:
    """Forget the panels (alignments are keyed on the sequences themselves)."""
    sequences.cache_clear()
