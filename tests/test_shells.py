"""Ligand shells: the edges, the consensus rule, and the real pocket."""

import pytest

from ip3r.structure.shells import consensus_shells, shell_edges, shell_of
from conftest import needs_genes, needs_structure

DEPOSITS = ("6DQN", "8TKG", "8TKF", "8TKH", "7T3P", "8TLA")


def test_edges_are_half_open_like_s22():
    assert shell_edges() == (4.5, 8.0, 11.5, 15.0)
    assert shell_of(0.0) == "contact"
    assert shell_of(4.49) == "contact" and shell_of(4.5) == "second"
    assert shell_of(11.5) == "fourth" and shell_of(14.99) == "fourth"
    assert shell_of(15.0) is None


def test_consensus_is_the_median_over_resolving_deposits():
    per = {"A": {1: 3.0, 2: 9.0, 3: 20.0},
           "B": {1: 5.0, 2: 10.0},
           "C": {1: 4.0}}
    got = {s.resi: s for s in consensus_shells(per)}
    assert got[1].median == 4.0 and got[1].n_structures == 3
    assert got[1].n_contact == 2 and got[1].shell == "contact"
    assert got[2].median == 9.5 and got[2].shell == "third"
    assert 3 not in got                        # beyond the radius: absent


@needs_genes
@needs_structure(*DEPOSITS)
def test_real_pocket_counts():
    from ip3r.analysis.shell_constraint import measured_shells
    shells = measured_shells()
    counts = {k: sum(s.shell == k for s in shells)
              for k in ("contact", "second", "third", "fourth")}
    assert counts == {"contact": 12, "second": 14, "third": 40, "fourth": 59}
    contact = {s.resi for s in shells if s.shell == "contact"}
    assert {276, 411} <= contact              # S22's two consensus extras


@needs_genes
@needs_structure(*DEPOSITS)
def test_pocket_carries_to_every_paralog():
    from ip3r.analysis.shell_constraint import pocket
    for gene in ("ITPR1", "ITPR2", "ITPR3"):
        pk = pocket(gene)
        assert pk.n_unaligned == 0 and len(pk.resi) == 125
        assert pk.distance.min() >= 0 and pk.distance.max() < 15.0
    with pytest.raises(KeyError):
        pocket("RYR1")
