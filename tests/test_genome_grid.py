"""The genome × paralog grid: the recovery-channel rule, the bar, ordering,
and the real grid's shape."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ip3r.analysis import genome_grid as GG
from ip3r.parameters import PARAMETERS
from conftest import needs_genes


def _row(present="1", proteome="1", family=3, any_cell=2, this_cell=1):
    return {"gene_present": present, "has_reference_proteome": proteome,
            "n_family_protein_records_species": str(family),
            "n_records_resolving_to_any_cell": str(any_cell),
            "n_records_resolving_to_cell": str(this_cell)}


@pytest.mark.parametrize("kw, want", [
    ({}, GG.REACHABLE),
    ({"present": "0"}, GG.NO_GENE),
    ({"proteome": "0", "family": 0}, GG.NO_PROTEOME),
    ({"family": 0, "any_cell": 0, "this_cell": 0}, GG.NO_RECORD),
    ({"any_cell": 0, "this_cell": 0}, GG.FRAGMENTS_ONLY),
    ({"this_cell": 0}, GG.OTHER_PARALOG),
])
def test_recovery_channel_rule(kw, want):
    assert GG.recovery_channel(_row(**kw)) == want


def test_bar_is_the_registered_parameter():
    bar = PARAMETERS.value("genomes.contiguity_bar_bp")
    assert GG.above_bar(bar) and not GG.above_bar(bar - 1)
    assert not GG.above_bar(float("nan"))
    try:
        PARAMETERS.set_value("genomes.contiguity_bar_bp", 10.0)
        assert GG.above_bar(bar - 1)
    finally:
        PARAMETERS.reset()


def test_order_puts_missing_n50_last_and_keeps_rows_together():
    gs = [GG.Genome("A", "a", "Aves", 5e4), GG.Genome("B", "b", "Aves", float("nan")),
          GG.Genome("C", "c", "Mammalia", 3e6)]
    lay = np.array([["x"] * 4, ["y"] * 4, ["z"] * 4], dtype=object)
    grid = GG.GenomeGrid(gs, {"search": lay})
    o = GG.order(grid, "contig N50")
    assert [g.accession for g in o.genomes] == ["C", "A", "B"]
    assert list(o.layers["search"][:, 0]) == ["z", "x", "y"]
    assert [g.accession for g in GG.order(grid, "class, then N50").genomes] == ["A", "B", "C"]


@needs_genes
def test_real_grid_shape_and_counts():
    grid = GG.load_grid()
    assert grid.n == 309 and grid.cells == GG.CELLS
    assert sum(grid.counts("miss").values()) == 1236
    assert grid.counts("miss", ["ITPR1", "ITPR2", "ITPR3"])["missed"] == 140
    assert grid.counts("miss", ["RYR"])["missed"] == 42
    assert sum(grid.counts("state").values()) == 927          # no RyR state
    assert grid.counts("recovery", ["ITPR1", "ITPR2", "ITPR3"])[GG.REACHABLE] == 179
    assert sum(g.above_bar for g in grid.genomes) == 189
    assert all(math.isfinite(g.contig_n50) for g in grid.genomes)
    row = grid.row(0)
    assert set(row) == set(GG.CELLS) and set(row["RYR"]) == set(GG.LAYERS)
