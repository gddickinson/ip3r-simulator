"""Paper 5 §8's VUS stratification: the rule on toy data, then the real table.

The toy cases pin what the sentence leaves open: positions not alleles,
ties at a median count, a residue in two classes is in both, and a class
with no scored position gives no row. The real cases prove the viewer's
resource route gives S17's table, and that the spheres land where they
should.
"""

import numpy as np
import pytest

from conftest import needs_genes, needs_structure
from ip3r.analysis.vus_strata import stratify, stratum_of


def _v(resi, bucket, source="clinvar", gene="G"):
    return {"gene": gene, "resi": resi, "class_bucket": bucket, "source": source}


SCORES = {1: 0.9, 2: 0.8, 3: 0.2, 4: 0.3, 5: 0.9, 6: 0.5, 7: 0.2, 8: np.nan, 9: 0.95}


def _score(gene, layer, resi):
    return SCORES[resi]


def test_rule_by_hand():
    rows = [_v(1, "P/LP"), _v(1, "P/LP"), _v(2, "P/LP"),        # two alleles at 1
            _v(3, "B/LB"), _v(4, "B/LB"),
            _v(5, "VUS"), _v(6, "VUS"), _v(7, "VUS"), _v(8, "VUS"),
            _v(1, "VUS")]                                        # 1 is P/LP and VUS
    s = stratify(rows, _score, "G", "deep")
    assert (len(s.pathogenic), len(s.benign), s.n_vus) == (2, 2, 4)   # 8 unscored
    assert (s.median_pathogenic, s.median_benign) == pytest.approx((0.85, 0.25))
    assert s.strata == {5: "pathogenic-like", 6: "between", 7: "benign-like",
                        1: "pathogenic-like"}
    r = s.row()
    assert r["n_vus_above_pathogenic_median"] == 2
    assert r["frac_vus_below_benign_median"] == 0.25
    assert s.median_vus == pytest.approx(0.7)


def test_ties_count_and_inverted_medians():
    assert stratum_of(0.85, 0.85, 0.25) == "pathogenic-like"
    assert stratum_of(0.25, 0.85, 0.25) == "benign-like"
    assert stratum_of(0.5, 0.4, 0.6) == "both"
    assert stratum_of(np.nan, 0.85, 0.25) == "not scored"


def test_sources_and_empty_class():
    rows = [_v(1, "P/LP", "uniprot"), _v(2, "P/LP"), _v(3, "B/LB"), _v(5, "VUS")]
    assert len(stratify(rows, _score, "G", "deep").pathogenic) == 2
    assert len(stratify(rows, _score, "G", "deep", ("clinvar",)).pathogenic) == 1
    assert stratify([_v(2, "P/LP"), _v(5, "VUS")], _score, "G", "deep") is None
    assert stratify([_v(8, "B/LB"), _v(2, "P/LP"), _v(5, "VUS")],
                    _score, "G", "deep") is None       # only B/LB is unscored


@needs_genes
def test_resource_route_reproduces_s17():
    """What the viewer draws (committed resources) = S17's table."""
    from ip3r.analysis.checks_variants import VUS_STRATA, _compare
    from ip3r.core import genes_data as G
    from ip3r.render.variant_spheres import resource_stratification
    for pub in G.read_tsv(VUS_STRATA):
        s = resource_stratification(pub["gene"], pub["layer"])
        assert not _compare(s.row(), pub, 0.0), (pub["gene"], pub["layer"])


@needs_structure("6DQN")
def test_spheres_on_every_subunit():
    from ip3r.io.loader import load
    from ip3r.render.variant_spheres import variant_classes, variant_spheres
    st = load("6DQN")
    idx, radius, rgb, labels = variant_spheres(st, "ITPR3", ("P/LP",))
    assert (st.atom_name[idx] == "CA").all()
    resolved = {int(r) for r in st.res_seq[idx]}
    assert resolved <= set(variant_classes("ITPR3", ("P/LP",)))
    for r in resolved:                                 # one sphere per subunit
        assert (st.res_seq[idx] == r).sum() == len(st.chains)
    # A residue with a P/LP and a VUS allele is drawn P/LP when both are shown.
    both = variant_classes("ITPR3", ("P/LP", "VUS"))
    assert all(both[r] == "P/LP" for r in variant_classes("ITPR3", ("P/LP",)))
    # With a layer, VUS take a stratum; one chain hidden drops its spheres.
    idx2, _, _, lab2 = variant_spheres(st, "ITPR3", ("VUS",), "family",
                                       st.chain != st.chains[0])
    assert set(lab2) <= {"pathogenic-like", "between", "benign-like", "not scored"}
    assert st.chains[0] not in set(st.chain[idx2])
