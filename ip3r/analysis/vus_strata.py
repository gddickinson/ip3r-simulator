"""Paper 5 §8: where the uncertain variants sit on a conservation layer.

S17 places each VUS against the labelled variants *of the same gene* on one
layer, using the labelled sets' own medians as thresholds so that no cut is
chosen to make a count: a VUS is *pathogenic-like* at or above the P/LP
median and *benign-like* at or below the B/LB median. It is a stratification,
not a call.

The rule is rebuilt here from S17's sentence and code comments, not imported:

* one score per **position** per class (two alleles at a residue are one
  observation of its constraint); a residue carrying both a VUS and a P/LP
  allele is counted in both classes;
* a layer scores a residue only at occupancy ≥ ``constraint.min_occupancy``;
* **every source** — the curated UniProt P/LP records are in the P/LP median
  here, although the classifier test (``P5.variant_auc``) is ClinVar only.
  ``sources=("clinvar",)`` gives the ClinVar-only variant for comparison.

A median from a few positions is itself uncertain (ITPR2's P/LP median is
one position's score). :func:`bands` draws an exact order-statistic
interval round each median (``vus.median_level``; unbounded when a class has
too few positions for that level), and :meth:`Bands.firm` says whether a
VUS's stratum survives any threshold inside them or sits *near a median*.

``score`` is any ``(gene, layer, resi) -> float`` (NaN = not scored), so the
check can use the publication's per-residue tables and the GUI the committed
resource without either importing the other.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field
from typing import Callable, Iterable

import numpy as np

from .stats import median_interval

__all__ = ["STRATA", "STRATUM_LABELS", "CLASS_COLORS", "STRATUM_COLORS",
           "Stratification", "position_scores",
           "stratify", "stratum_of", "Bands", "bands", "NEAR"]

#: What the Variants table appends to a stratum that a threshold inside its
#: median's interval could change (followed by which median).
NEAR = "near"

#: What a VUS can be on one layer; "not scored" is grey, never low.
STRATA = ("pathogenic-like", "between", "benign-like", "both", "not scored")
STRATUM_LABELS = {
    "pathogenic-like": "VUS ≥ P/LP median",
    "between": "VUS between the medians",
    "benign-like": "VUS ≤ B/LB median",
    "both": "VUS past both medians (medians inverted)",
    "not scored": "VUS not scored on this layer",
}

Score = Callable[[str, str, int], float]


def position_scores(rows: Iterable[dict], score: Score, gene: str, layer: str,
                    bucket: str, sources: tuple[str, ...] | None = None
                    ) -> dict[int, float]:
    """``resi -> score`` for one class of one gene, one entry per position."""
    out: dict[int, float] = {}
    for r in rows:
        if r["gene"] != gene or r["class_bucket"] != bucket:
            continue
        if sources is not None and r["source"] not in sources:
            continue
        resi = int(r["resi"])
        if resi in out:
            continue
        x = score(gene, layer, resi)
        if np.isfinite(x):
            out[resi] = float(x)
    return out


def stratum_of(x: float, p_median: float, b_median: float) -> str:
    if not np.isfinite(x):
        return "not scored"
    hi, lo = x >= p_median, x <= b_median
    if hi and lo:
        return "both"
    return "pathogenic-like" if hi else "benign-like" if lo else "between"


@dataclass(frozen=True)
class Stratification:
    gene: str
    layer: str
    pathogenic: dict[int, float]
    benign: dict[int, float]
    vus: dict[int, float]
    median_pathogenic: float
    median_benign: float
    median_vus: float
    strata: dict[int, str] = field(repr=False)

    @property
    def n_vus(self) -> int:
        return len(self.vus)

    def count(self, stratum: str) -> int:
        return sum(s == stratum or (s == "both" and stratum != "between")
                   for s in self.strata.values())

    @property
    def n_above_pathogenic(self) -> int:
        return self.count("pathogenic-like")

    @property
    def n_below_benign(self) -> int:
        return self.count("benign-like")

    def row(self) -> dict:
        """The fields of S17's ``vus_stratification.tsv``, unrounded."""
        n = self.n_vus
        return {"gene": self.gene, "layer": self.layer, "n_vus": n,
                "n_pathogenic": len(self.pathogenic), "n_benign": len(self.benign),
                "median_pathogenic": self.median_pathogenic,
                "median_benign": self.median_benign, "median_vus": self.median_vus,
                "n_vus_above_pathogenic_median": self.n_above_pathogenic,
                "frac_vus_above_pathogenic_median": self.n_above_pathogenic / n,
                "n_vus_below_benign_median": self.n_below_benign,
                "frac_vus_below_benign_median": self.n_below_benign / n}


def stratify(rows: Iterable[dict], score: Score, gene: str, layer: str,
             sources: tuple[str, ...] | None = None) -> Stratification | None:
    """Stratify one gene's VUS on one layer; None if a class is empty
    (S17 writes no row then)."""
    rows = list(rows)
    path, ben, vus = (position_scores(rows, score, gene, layer, b, sources)
                      for b in ("P/LP", "B/LB", "VUS"))
    if not (path and ben and vus):
        return None
    pm, bm = statistics.median(path.values()), statistics.median(ben.values())
    strata = {r: stratum_of(x, pm, bm) for r, x in vus.items()}
    return Stratification(gene, layer, path, ben, vus, pm, bm,
                          statistics.median(vus.values()), strata)


@dataclass(frozen=True)
class Bands:
    """Each labelled median's interval (± inf = unbounded) at ``level``."""
    pathogenic: tuple[float, float]
    benign: tuple[float, float]
    level: float

    @property
    def bounded(self) -> bool:
        return all(np.isfinite(self.pathogenic + self.benign))

    def unsure(self, x: float) -> tuple[str, ...]:
        """The thresholds ("P/LP", "B/LB") a value inside their interval
        could put ``x`` on either side of; empty = the stratum is firm."""
        if not np.isfinite(x):
            return ()                                # "not scored" is no threshold call
        (plo, phi), (blo, bhi) = self.pathogenic, self.benign
        out = []
        if (x >= phi) != (x >= plo):                 # pathogenic-like is x ≥ median
            out.append("P/LP")
        if (x <= blo) != (x <= bhi):                 # benign-like is x ≤ median
            out.append("B/LB")
        return tuple(out)

    def firm(self, x: float) -> bool:
        """Is the stratum of a score ``x`` the same for every pair of
        thresholds inside the two intervals?"""
        return not self.unsure(x)

    def near(self, s: "Stratification", threshold: str | None = None) -> list[int]:
        """VUS positions whose stratum is not firm (about ``threshold`` only,
        if given)."""
        return sorted(r for r, x in s.vus.items() if (
            threshold in self.unsure(x) if threshold else not self.firm(x)))


def bands(s: "Stratification", level: float | None = None) -> Bands:
    """The two thresholds' intervals for one stratification
    (``vus.median_level`` unless given)."""
    if level is None:
        from ..parameters import PARAMETERS
        level = PARAMETERS.value("vus.median_level")
    return Bands(median_interval(list(s.pathogenic.values()), level),
                 median_interval(list(s.benign.values()), level), level)


#: Class colours for the variant spheres and the figures. P/LP and B/LB take
#: the strong ends; a stratified VUS takes a paler shade of the end it sits
#: at, so "pathogenic-like" never reads as "pathogenic". Missing is grey.
CLASS_COLORS = {
    "P/LP": (0.90, 0.22, 0.22),
    "B/LB": (0.25, 0.50, 0.95),
    "VUS": (0.93, 0.80, 0.30),
    "conflicting": (0.75, 0.45, 0.90),
    "other": (0.60, 0.62, 0.66),
}
STRATUM_COLORS = {
    "pathogenic-like": (0.98, 0.60, 0.55),
    "between": (0.93, 0.80, 0.30),
    "benign-like": (0.60, 0.78, 0.98),
    "both": (0.85, 0.85, 0.85),
    "not scored": (0.42, 0.43, 0.46),
}
