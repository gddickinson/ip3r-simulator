"""Paper 3's lesion lead (S15b §8.2–8.3), rebuilt from the per-locus table.

S15a scored every located ITPR/RyR locus for frameshifts and internal stops
(``integrity_loci.tsv``) and compared each cell with *its own genome's*
other family loci: a positive difference is a paralogue carrying more
lesions per kilo-residue than its siblings in the same assembly. S15b
stratified that identity-matched sign test by vertebrate class and found
ITPR3's excess to be a bird result, mostly in assemblies below the
contiguity bar.

Everything here is rebuilt from the loci, not read from the pair or test
tables: :func:`best_loci` keeps each genome's best-covered scored locus per
cell, :func:`pairs` compares it with the median of its identity-matched
siblings, :func:`stratify` runs the exact sign test per cell × class with
Benjamini–Hochberg over the testable strata, :func:`controls` splits each
significant stratum at this project's registered contiguity bar, and
:func:`layer` paints the per-genome sign on the Genomes grid.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from .stats import benjamini_hochberg, sign_test

__all__ = ["LOCI", "PAIRS", "BY_CLASS", "CONTROLS", "Pair", "Stratum", "Split",
           "density", "is_scored", "best_loci", "pairs", "stratify", "controls",
           "layer", "LAYER_VALUES"]

LOCI = "loss_dynamics/integrity_loci.tsv"
PAIRS = "loss_dynamics/integrity_pairs.tsv"
BY_CLASS = "loss_counts/lesion_by_class.tsv"
CONTROLS = "loss_counts/lesion_class_controls.tsv"

EXCESS, DEFICIT, TIE = "excess", "deficit", "tie"
NO_SIBLING, NOT_SCORED = "no matched sibling", "not scored"
#: What a grid cell of the lesion layer can hold, in legend order.
LAYER_VALUES = (EXCESS, DEFICIT, TIE, NO_SIBLING, NOT_SCORED)


def density(row: dict) -> float:
    """Lesions (frameshifts + internal stops) per kilo-residue aligned, at
    the table's four decimals."""
    aa = float(row["aligned_aa"])
    lesions = int(row["frameshifts"]) + int(row["stop_codons"])
    return round(lesions / (aa / 1000.0), 4) if aa > 0 else 0.0


def is_scored(row: dict) -> bool:
    """A locus is scored at the sweep's full-coverage bar
    (``lesion.coverage_bar``)."""
    return float(row["coverage"] or 0.0) >= _P.value("lesion.coverage_bar")


def best_loci(loci) -> dict[str, dict[str, tuple[float, float, float]]]:
    """accession → cell → (coverage, density, identity) of the genome's
    best-covered scored locus (the first, on a tie)."""
    out: dict[str, dict[str, tuple]] = {}
    for r in loci:
        if not is_scored(r):
            continue
        cells = out.setdefault(r["accession"], {})
        cov = float(r["coverage"])
        if r["cell"] not in cells or cov > cells[r["cell"]][0]:
            cells[r["cell"]] = (cov, density(r), float(r["identity"]))
    return out


@dataclass(frozen=True)
class Pair:
    accession: str
    cell: str
    density: float
    identity: float
    others_median: float
    others_identity_median: float
    n_others: int

    @property
    def diff(self) -> float:
        return self.density - self.others_median

    @property
    def sign(self) -> str:
        return EXCESS if self.diff > 0 else DEFICIT if self.diff < 0 else TIE


def pairs(loci, matched: bool = True) -> list[Pair]:
    """Each cell against the median of its own genome's other scored loci;
    ``matched`` keeps only siblings within ``lesion.identity_window`` of its
    identity to the bait (the form the paper's statement is made on)."""
    window = _P.value("lesion.identity_window")
    out = []
    for acc, cells in best_loci(loci).items():
        for cell, (_, dens, ident) in cells.items():
            others = [(d, i) for c, (_, d, i) in cells.items() if c != cell]
            if matched:
                # rounded like the table's columns, so the window is exact
                others = [(d, i) for d, i in others
                          if round(abs(i - ident), 10) <= window]
            if not others:
                continue
            out.append(Pair(acc, cell, dens, ident,
                            float(np.median([d for d, _ in others])),
                            float(np.median([i for _, i in others])), len(others)))
    return out


@dataclass(frozen=True)
class Stratum:
    cell: str
    vclass: str
    n_genomes: int
    n_pos: int
    n_neg: int
    n_ties: int
    median_diff: float
    p: float                  # NaN when underpowered
    q: float = float("nan")

    @property
    def n(self) -> int:
        return self.n_pos + self.n_neg

    @property
    def direction(self) -> str:
        return (EXCESS if self.n_pos > self.n_neg else
                DEFICIT if self.n_neg > self.n_pos else "tied")

    @property
    def testable(self) -> bool:
        return math.isfinite(self.p)


def _test(diffs) -> tuple[int, int, int, float]:
    t = sign_test(diffs)
    return t["n_pos"], t["n_neg"], t["n_ties"], t["p"]


def _powered(n: int, p: float) -> float:
    return p if n >= _P.value("lesion.min_untied") else float("nan")


def stratify(prs, vclass_of: dict[str, str]) -> list[Stratum]:
    """The sign test per cell × class; a stratum with fewer than
    ``lesion.min_untied`` untied pairs gets no p and stays out of the BH
    family (one family, one correction)."""
    groups: dict[tuple[str, str], list[float]] = {}
    for p in prs:
        groups.setdefault((p.cell, vclass_of.get(p.accession, "unknown")),
                          []).append(p.diff)
    rows = []
    for (cell, vclass), diffs in sorted(groups.items()):
        pos, neg, ties, pv = _test(diffs)
        rows.append(Stratum(cell, vclass, len(diffs), pos, neg, ties,
                            float(np.median(diffs)), _powered(pos + neg, pv)))
    q = benjamini_hochberg([s.p for s in rows])
    return [Stratum(**{**s.__dict__, "q": float(qq)}) for s, qq in zip(rows, q)]


@dataclass(frozen=True)
class Split:
    """A significant stratum split at the contiguity bar."""
    stratum: Stratum
    above: tuple[int, int, float]     # (excess, deficit, p or NaN)
    below: tuple[int, int, float]
    median_identity: float
    median_others_identity: float
    siblings: tuple[str, ...]


def controls(strata, prs, vclass_of: dict[str, str], above_of: dict[str, bool],
             q_cut: float | None = None) -> list[Split]:
    """Each stratum with q ≤ ``q_cut`` (``check.alpha``) split above and
    below the bar, with the same class's other significant cells named as
    its siblings: one within-genome comparison read from both sides."""
    q_cut = _P.value("check.alpha") if q_cut is None else q_cut
    sig = [s for s in strata if s.testable and s.q <= q_cut]
    out = []
    for s in sig:
        g = [p for p in prs if p.cell == s.cell
             and vclass_of.get(p.accession, "unknown") == s.vclass]
        halves = []
        for side in (True, False):
            pos, neg, _, pv = _test([p.diff for p in g
                                     if above_of.get(p.accession, False) == side])
            halves.append((pos, neg, _powered(pos + neg, pv)))
        out.append(Split(s, halves[0], halves[1],
                         float(np.median([p.identity for p in g])),
                         float(np.median([p.others_identity_median for p in g])),
                         tuple(sorted(t.cell for t in sig
                                      if t.vclass == s.vclass and t.cell != s.cell))))
    return out


def layer(accessions, cells, loci=None) -> np.ndarray:
    """The Genomes grid's lesion layer: each cell's identity-matched sign
    (:data:`LAYER_VALUES`); ``""`` where the loci table has no row."""
    loci = G.read_tsv(LOCI) if loci is None else loci
    sign = {(p.accession, p.cell): p.sign for p in pairs(loci)}
    scored = {(a, c) for a, cs in best_loci(loci).items() for c in cs}
    seen = {(r["accession"], r["cell"]) for r in loci}
    m = np.full((len(accessions), len(cells)), "", dtype=object)
    for i, acc in enumerate(accessions):
        for j, cell in enumerate(cells):
            key = (acc, cell)
            if key in sign:
                m[i, j] = sign[key]
            elif key in scored:
                m[i, j] = NO_SIBLING
            elif key in seen:
                m[i, j] = NOT_SCORED
    return m
