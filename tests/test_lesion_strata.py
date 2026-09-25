"""Paper 3's lesion strata (S15b §8.2-8.3): the rule on toy loci, then the
committed tables rebuilt from the per-locus one.

The toy cases pin what the code comments leave open: the best-covered
scored locus stands for a cell (the first on a tie), siblings enter only
within the identity window, the difference is against the siblings'
median, a stratum below the untied minimum gets no p and stays out of the
BH family, and a significant stratum's siblings are the same class's other
significant cells.
"""

import math

import numpy as np
import pytest

from conftest import needs_genes
from ip3r.analysis import lesion_strata as L


def _locus(acc, cell, cov=1.0, ident=0.95, fs=0, st=0, aa=1000):
    return {"accession": acc, "cell": cell, "coverage": str(cov), "identity": str(ident),
            "frameshifts": str(fs), "stop_codons": str(st), "aligned_aa": str(aa)}


def test_density_and_scoring():
    assert L.density(_locus("a", "X", fs=2, st=1, aa=2000)) == 1.5
    assert L.is_scored(_locus("a", "X", cov=0.70))
    assert not L.is_scored(_locus("a", "X", cov=0.6999))


def test_best_covered_locus_stands_for_the_cell():
    loci = [_locus("g", "ITPR1", cov=0.8, fs=5), _locus("g", "ITPR1", cov=0.95, fs=1),
            _locus("g", "ITPR1", cov=0.95, fs=9),                 # a tie keeps the first
            _locus("g", "ITPR2", cov=0.5, fs=7)]                  # not scored
    best = L.best_loci(loci)
    assert best == {"g": {"ITPR1": (0.95, 1.0, 0.95)}}


def test_pairs_by_hand():
    loci = [_locus("g", "A", ident=0.90, fs=4), _locus("g", "B", ident=0.91, fs=1),
            _locus("g", "C", ident=0.915, fs=3), _locus("g", "D", ident=0.80, fs=0)]
    prs = {p.cell: p for p in L.pairs(loci)}
    # A's matched siblings are B and C (D is 0.10 away): median 2, diff +2
    assert prs["A"].n_others == 2 and prs["A"].diff == pytest.approx(2.0)
    assert prs["A"].sign == L.EXCESS
    assert "D" not in prs                         # no sibling within the window
    unmatched = {p.cell: p for p in L.pairs(loci, matched=False)}
    assert unmatched["D"].n_others == 3 and unmatched["D"].sign == L.DEFICIT


def test_window_edge_is_inclusive():
    loci = [_locus("g", "A", ident=0.90), _locus("g", "B", ident=0.92)]
    assert len(L.pairs(loci)) == 2                # |0.92 - 0.90| = the window


def _pair(acc, cell, diff):
    return L.Pair(acc, cell, diff, 0.9, 0.0, 0.9, 1)


def test_stratify_underpowered_out_of_the_family():
    prs = ([_pair(f"b{i}", "X", 1.0) for i in range(10)]        # 10:0, testable
           + [_pair(f"f{i}", "X", 1.0) for i in range(3)]       # 3 untied: no p
           + [_pair(f"b{i}", "Y", -1.0) for i in range(10)])
    vclass = {**{f"b{i}": "Aves" for i in range(10)}, **{f"f{i}": "Fish" for i in range(3)}}
    s = {(x.cell, x.vclass): x for x in L.stratify(prs, vclass)}
    assert not s[("X", "Fish")].testable and math.isnan(s[("X", "Fish")].q)
    p = 2 * 0.5 ** 10
    assert s[("X", "Aves")].p == pytest.approx(p)
    assert s[("X", "Aves")].q == pytest.approx(p)          # m = 2, both equal
    assert s[("Y", "Aves")].direction == L.DEFICIT


def test_controls_split_and_siblings():
    prs = ([_pair(f"b{i}", "X", 1.0) for i in range(12)]
           + [_pair(f"b{i}", "Y", -1.0) for i in range(12)])
    vclass = {f"b{i}": "Aves" for i in range(12)}
    above = {f"b{i}": i < 3 for i in range(12)}
    strata = L.stratify(prs, vclass)
    splits = {c.stratum.cell: c for c in L.controls(strata, prs, vclass, above)}
    x = splits["X"]
    assert x.above[:2] == (3, 0) and math.isnan(x.above[2])     # 3 < 8 untied
    assert x.below[:2] == (9, 0) and x.below[2] == pytest.approx(2 * 0.5 ** 9)
    assert x.siblings == ("Y",) and splits["Y"].siblings == ("X",)


def test_layer_values():
    loci = [_locus("g", "A", ident=0.90, fs=4), _locus("g", "B", ident=0.905),
            _locus("g", "C", ident=0.50), _locus("g", "D", cov=0.3)]
    m = L.layer(["g", "h"], ["A", "B", "C", "D", "E"], loci)
    assert list(m[0]) == [L.EXCESS, L.DEFICIT, L.NO_SIBLING, L.NOT_SCORED, ""]
    assert list(m[1]) == [""] * 5


# ------------------------------------------------------------ real tables

@needs_genes
def test_pairs_rebuild_the_committed_table():
    from ip3r.core import genes_data as G
    loci = G.read_tsv(L.LOCI)
    assert all(L.density(r) == pytest.approx(float(r["lesion_density"]), abs=1e-9)
               for r in loci)
    for matched in (True, False):
        ours = {(p.accession, p.cell): p.diff for p in L.pairs(loci, matched)}
        pub = {(r["accession"], r["cell"]): float(r["diff"])
               for r in G.read_tsv(L.PAIRS) if r["matched"] == str(int(matched))}
        assert ours.keys() == pub.keys()
        assert max(abs(ours[k] - pub[k]) for k in ours) < 1e-6


@needs_genes
def test_grid_layer_is_the_pairs():
    from ip3r.analysis import genome_grid as GG
    g = GG.load_grid(("search", "lesion"))
    c = g.counts("lesion")
    assert c[L.EXCESS] + c[L.DEFICIT] + c[L.TIE] == 744
    birds = [i for i, x in enumerate(g.genomes) if x.vclass == "Aves"]
    col = g.layers["lesion"][birds, GG.CELLS.index("ITPR3")]
    assert (np.sum(col == L.EXCESS), np.sum(col == L.DEFICIT)) == (25, 2)
