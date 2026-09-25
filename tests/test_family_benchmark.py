"""The family-call benchmark's arithmetic, on inputs whose answer is known."""

import pytest

from ip3r.analysis.family_benchmark import Margin, identity, rescore


def test_identity_metrics_by_hand():
    # identical: both metrics 1
    assert identity("ACDEFGHIKL", "ACDEFGHIKL") == (1.0, 1.0)
    # a 10-mer inside a 14-mer: 10 of 14 columns shared (full), all covered identical
    full, cov = identity("ACDEFGHIKL", "WWACDEFGHIKLWW")
    assert cov == 1.0 and full == pytest.approx(10 / 14)


def _m(full_itpr, full_ryr):
    return Margin("X", (full_itpr, 0.0), (full_ryr, 0.0))


def test_band_and_call_follow_the_d7_margin():
    assert _m(0.60, 0.10).band() == "ITPR"
    assert _m(0.15, 0.10).band() == "no call" and _m(0.15, 0.10).call() == "ITPR"
    assert _m(0.10, 0.70).band() == "RYR"
    assert _m(0.2, 0.2).call() == "tie"


def test_rescore_points_and_caps():
    itpr, ryr = _m(0.7, 0.1), _m(0.1, 0.7)
    assert rescore("size+15,pfam+20,breadth+15", itpr).score == 50
    assert rescore("size+15,pfam+20", itpr).score == 35           # the missed fly
    # a RyR with the family's Pfams: 45 raw, capped by the margin, not a flag
    s = rescore("pfam+20,cluster+10,breadth+15", ryr)
    assert (s.raw, s.score, s.cap) == (45, 39, "sister")
    # generic components alone reach 40 but fail the evidence gate (TLN1)
    s = rescore("size+15,cluster+10,breadth+15", itpr)
    assert (s.raw, s.score, s.cap) == (40, 39, "gate")
    # S1's own flags are ignored: nothing to cap below the bar
    assert rescore("pfam+20,cluster+10,SISTER→39", ryr).score == 30
    # the ceiling
    assert rescore("outlier+20,fold+20,pfam+20,size+15,breadth+15,cluster+10,"
                   "homology+10,split+10", itpr).score == 100
