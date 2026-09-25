"""S23's genomes one by one: every genome placed without a guess, the copy
count equal to the ledger's, contiguity against each genome's own bar."""

import matplotlib

matplotlib.use("Agg")

import pytest

from conftest import needs_genes
from ip3r.analysis import range_genomes as RG


@pytest.fixture(scope="module")
def rows():
    return RG.genome_rows()


@needs_genes
def test_every_genome_placed(rows):
    from collections import Counter
    assert len(rows) == 194
    how = Counter(r.placed_by for r in rows)
    assert how["taxid"] == 109
    unplaced = [r for r in rows if r.clade == RG.NOT_SWEPT]
    assert len(unplaced) == 13 and all(r.placed_by == "" for r in unplaced)
    names = [c for c, _ in RG.clades_with_genomes(rows)]
    assert names[-1] == RG.NOT_SWEPT


@needs_genes
def test_copies_and_bar_agree_with_the_ledgers(rows):
    from ip3r.analysis import range_table as RT
    from ip3r.core import genes_data as G
    led = {r["accession"]: r for r in G.read_tsv(RT.COPY_LEDGER)}
    man = {r["accession"]: r for r in G.read_tsv(RT.MANIFEST)}
    assert all(r.copies == int(led[r.accession]["n_full"]) for r in rows)
    assert all(r.spans_gene == (man[r.accession]["spans_gene"] == "Y") for r in rows)
    assert max(r.copies for r in rows) == 18
    from ip3r.analysis.range_figure import COPIES_MAX
    assert max(r.copies for r in rows) <= COPIES_MAX      # the fixed scale holds all


@needs_genes
def test_in_clade_and_drawn(rows):
    import matplotlib.pyplot as plt
    from ip3r.analysis.range_figure import draw_genomes
    clade = RG.clades_with_genomes(rows)[0][0]
    shown = RG.in_clade(rows, clade)
    assert [r.copies for r in shown] == sorted((r.copies for r in shown), reverse=True)
    fig, ax = plt.subplots()
    info = draw_genomes(ax, shown)
    assert info["genomes"] == len(shown)
    assert info["copies"] == sum(r.copies for r in shown)
    plt.close(fig)
