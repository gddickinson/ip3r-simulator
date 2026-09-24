"""Where the S17 variants go on a deposit, and in what colour.

One sphere per variant residue per resolved subunit, on the residue's Cα
(the cartoon passes through it). A residue carrying alleles of several
classes takes the most decisive one: P/LP, B/LB, conflicting, VUS, other.
With a layer chosen, a VUS is coloured by where it sits on that layer
against its gene's labelled medians (``analysis.vus_strata``), from the
committed resources — so the viewer needs no ip3r_genes checkout.

Only called for a deposit in the variants' own numbering; the caller refuses
otherwise, as for every residue-keyed annotation.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from ..analysis.vus_strata import (CLASS_COLORS, STRATUM_COLORS, Stratification,
                                   stratify, stratum_of)
from ..core.annotations import constraint_at, variants

__all__ = ["CLASS_ORDER", "resource_stratification", "variant_classes",
           "variant_spheres"]

CLASS_ORDER = ("P/LP", "B/LB", "conflicting", "VUS", "other")
RADIUS = {"P/LP": 2.4, "B/LB": 2.0}
DEFAULT_RADIUS = 1.6


def _resource_score(gene: str, layer: str, resi: int) -> float:
    return float(constraint_at(gene, np.array([resi]), layer)[0])


@lru_cache(maxsize=None)
def resource_stratification(gene: str, layer: str) -> Stratification | None:
    return stratify(variants(gene), _resource_score, gene, layer)


def variant_classes(gene: str, buckets) -> dict[int, str]:
    """``resi -> class`` for the chosen classes, most decisive class first."""
    rank = {b: i for i, b in enumerate(CLASS_ORDER)}
    out: dict[int, str] = {}
    for v in variants(gene):
        b = v["class_bucket"]
        if b not in buckets:
            continue
        r = int(v["resi"])
        if r not in out or rank[b] < rank[out[r]]:
            out[r] = b
    return out


def variant_spheres(st, gene: str, buckets, layer: str | None = None,
                    chain_mask: np.ndarray | None = None):
    """``(atom_index, radius, rgb, labels)`` — one entry per sphere."""
    classes = variant_classes(gene, buckets)
    strat = resource_stratification(gene, layer) if layer else None
    ca = (st.atom_name == "CA") & ~st.hetero
    if chain_mask is not None:
        ca &= chain_mask
    idx = np.flatnonzero(ca & np.isin(st.res_seq, list(classes)))
    radius = np.empty(len(idx), np.float32)
    rgb = np.empty((len(idx), 3), np.float32)
    labels = []
    for k, i in enumerate(idx):
        r = int(st.res_seq[i])
        c = classes[r]
        radius[k] = RADIUS.get(c, DEFAULT_RADIUS)
        if c == "VUS" and strat is not None:
            s = stratum_of(strat.vus.get(r, np.nan), strat.median_pathogenic,
                           strat.median_benign)
            rgb[k] = STRATUM_COLORS[s]
            labels.append(s)
        else:
            rgb[k] = CLASS_COLORS[c]
            labels.append(c)
    return idx, radius, rgb, labels
