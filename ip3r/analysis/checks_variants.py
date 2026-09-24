"""Paper 5 §8: the VUS stratification, re-derived from its inputs.

Every row of ``vus_stratification.tsv`` (3 genes × 4 layers × 10 numbers) is
rebuilt from ``variants.tsv`` and the per-residue tables with
:mod:`.vus_strata`, which implements the rule from S17's description, not
its code. The comparison is also run ClinVar-only, because the classifier test
excludes the curated UniProt records and the stratification does not.
"""

from __future__ import annotations

import numpy as np

from ..config import PARALOG_ACC
from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from .checks import agree, register
from .checks_constraint import LAYERS, VARIANTS, _score, _sites
from .vus_strata import stratify

VUS_STRATA = "constraint/vus_stratification.tsv"
_COUNTS = ("n_vus", "n_pathogenic", "n_benign", "n_vus_above_pathogenic_median",
           "n_vus_below_benign_median")
#: S17 writes medians and fractions to four decimals.
_ROUNDED = ("median_pathogenic", "median_benign", "median_vus",
            "frac_vus_above_pathogenic_median", "frac_vus_below_benign_median")


def table_score():
    """``(gene, layer, resi) -> JSD`` from the publication's per-residue tables."""
    sites = {g: _sites(g) for g in PARALOG_ACC}

    def score(gene: str, layer: str, resi: int) -> float:
        x = _score(sites[gene].get(resi), layer)
        return np.nan if x is None else x
    return score


def _compare(mine: dict, pub: dict, tol: float) -> list[str]:
    bad = [f"{k} {mine[k]} vs {pub[k]}" for k in _COUNTS if mine[k] != int(pub[k])]
    bad += [f"{k} {mine[k]:.4f} vs {pub[k]}" for k in _ROUNDED
            if abs(mine[k] - float(pub[k])) > max(tol, 5e-5)]
    return bad


@register("P5.vus_stratification", "constraint",
          "Each VUS placed against its own gene's labelled medians, per layer: "
          "every count, median and fraction of vus_stratification.tsv (on the "
          "family layer 10/6/11 % of VUS reach the P/LP median and 16/36/52 % "
          "sit at or below the B/LB median).",
          "Rebuilt from variants.tsv and the per-residue tables: one score per "
          "position per class, occupancy ≥ 0.5, thresholds = the labelled "
          "sets' medians; the VUS is pathogenic-like at ≥ the P/LP median and "
          "benign-like at ≤ the B/LB median. Also run on ClinVar records only.",
          "rederived", (VARIANTS, VUS_STRATA, "constraint/constraint_*.tsv"))
def vus_stratification():
    rows, score = G.read_tsv(VARIANTS), table_score()
    tol = _P.value("check.stat_tol")
    bad, n, clinvar_moves, out = [], 0, [], {}
    for pub in G.read_tsv(VUS_STRATA):
        gene, layer = pub["gene"], pub["layer"]
        s = stratify(rows, score, gene, layer)
        if s is None:
            bad.append(f"{gene} {layer}: a class is empty here")
            continue
        n += 1
        mine = s.row()
        bad += [f"{gene} {layer}: {b}" for b in _compare(mine, pub, tol)]
        cv = stratify(rows, score, gene, layer, sources=("clinvar",))
        if cv is None or _compare(cv.row(), pub, tol):
            clinvar_moves.append(f"{gene} {layer}")
        out.setdefault(gene, {})[layer] = {
            k: mine[k] for k in ("n_vus", "frac_vus_above_pathogenic_median",
                                 "frac_vus_below_benign_median", "median_pathogenic",
                                 "median_benign", "median_vus")}
    expected = len(PARALOG_ACC) * len(LAYERS)
    found = (f"{n}/{expected} rows rebuilt"
             + ("; every field agrees" if not bad else ""))
    detail = "; ".join(bad) if bad else (
        "ClinVar-only gives the same table: the curated UniProt records add no "
        "scored P/LP position." if not clinvar_moves else
        "ClinVar-only differs at " + ", ".join(clinvar_moves)
        + " — the curated UniProt records move the P/LP median there.")
    return agree(not bad and n == expected, f"{expected} rows", found, detail,
                 strata=out)
