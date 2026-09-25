"""Papers 3 and 4 as one picture: every genome × every cell of the sweep.

The retention sweep asked each of 309 vertebrate assemblies for ITPR1, ITPR2,
ITPR3 and a ryanodine-receptor control cell (RYR). Three committed tables
describe each of those 1,236 cells from different sides:

* ``methods/contiguity_cells.tsv`` — what the search found (``grade``), and
  whether a gene known to be present was missed (``false_negative``), beside
  the assembly's contig N50, level and annotation source;
* ``loss_dynamics/character_matrix.tsv`` — the S15a evidence state of each
  ITPR cell (the character Paper 3's parsimony is run on);
* ``methods/gene_recovery.tsv`` — whether any protein record reaches the gene
  (Paper 4). The *channel* is re-derived here from the table's count columns
  by :func:`recovery_channel`, not read from its label;
* ``loss_dynamics/integrity_loci.tsv`` — each located locus's lesions; the
  ``lesion`` layer is the cell's identity-matched within-genome sign
  (S15b §8.2), rebuilt by :func:`.lesion_strata.layer`.

:func:`load_grid` joins them into a :class:`GenomeGrid` of per-cell layers;
:func:`order` sorts its genomes. Nothing here draws (see ``grid_figure``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from . import lesion_strata as _L

__all__ = ["CONTIG", "MATRIX", "RECOVERY", "CELLS", "Genome", "GenomeGrid",
           "load_grid", "recovery_channel", "above_bar", "order", "ORDERS",
           "LAYERS"]

CONTIG = "methods/contiguity_cells.tsv"
MATRIX = "loss_dynamics/character_matrix.tsv"
RECOVERY = "methods/gene_recovery.tsv"
CELLS = ("ITPR1", "ITPR2", "ITPR3", "RYR")

#: Grid layers: key → (table, what a cell holds).
LAYERS = {
    "search": (CONTIG, "the sweep's grade of the cell"),
    "miss": (CONTIG, "a gene known to be present, missed by the search"),
    "state": (MATRIX, "S15a evidence state (ITPR cells only)"),
    "recovery": (RECOVERY, "how a protein-database search could reach the gene"),
    "lesion": (_L.LOCI, "lesions against the genome's identity-matched siblings"),
}

#: The recovery channels, in the order Paper 4 lists the reasons.
REACHABLE = "protein_database_and_genome"
NO_PROTEOME = "genome_only:no_reference_proteome"
FRAGMENTS_ONLY = "genome_only:records_exist_but_none_full_length"
OTHER_PARALOG = "genome_only:no_record_resolves_to_this_paralog"
NO_RECORD = "genome_only:no_family_record_for_species"
NO_GENE = "no_gene"


def recovery_channel(row: dict) -> str:
    """Paper 4's recovery channel of one gene_recovery.tsv row, from its
    counts alone: is there a gene, a reference proteome, any family record
    for the species, any full-length record resolving to a cell, and one
    resolving to *this* cell."""
    if row["gene_present"] != "1":
        return NO_GENE
    if row["has_reference_proteome"] != "1":
        return NO_PROTEOME
    if int(row["n_family_protein_records_species"]) == 0:
        return NO_RECORD
    if int(row["n_records_resolving_to_any_cell"]) == 0:
        return FRAGMENTS_ONLY
    if int(row["n_records_resolving_to_cell"]) == 0:
        return OTHER_PARALOG
    return REACHABLE


def above_bar(contig_n50: float) -> bool:
    """Can an assembly of this contig N50 hold the gene on one contig?
    (the registered ``genomes.contiguity_bar_bp``)."""
    return (math.isfinite(contig_n50)
            and contig_n50 >= _P.value("genomes.contiguity_bar_bp"))


@dataclass(frozen=True)
class Genome:
    accession: str
    organism: str
    vclass: str
    contig_n50: float = float("nan")
    level: str = ""
    source: str = ""

    @property
    def above_bar(self) -> bool:
        return above_bar(self.contig_n50)


@dataclass
class GenomeGrid:
    """Genomes (rows) × :data:`CELLS` (columns); each layer an object array
    of category strings, ``""`` where the table says nothing about a cell."""
    genomes: list[Genome]
    layers: dict[str, np.ndarray] = field(default_factory=dict)
    cells: tuple = CELLS

    @property
    def n(self) -> int:
        return len(self.genomes)

    def counts(self, layer: str, cells=None) -> dict[str, int]:
        cols = [self.cells.index(c) for c in (cells or self.cells)]
        vals, n = np.unique(self.layers[layer][:, cols], return_counts=True)
        return {str(v): int(k) for v, k in zip(vals, n) if v != ""}

    def row(self, i: int) -> dict[str, dict[str, str]]:
        """Everything the grid holds about genome ``i``, cell by cell."""
        return {c: {k: str(v[i, j]) for k, v in self.layers.items()}
                for j, c in enumerate(self.cells)}

    def subset(self, rows) -> "GenomeGrid":
        rows = list(rows)
        return GenomeGrid([self.genomes[i] for i in rows],
                          {k: v[rows] for k, v in self.layers.items()}, self.cells)


def _genome(r: dict) -> Genome:
    n50 = r.get("contig_n50", "")
    return Genome(r["accession"], r["organism"], r["vclass"],
                  float(n50) if n50 not in ("", None) else float("nan"),
                  r.get("level", ""), r.get("annotation_source", ""))


def load_grid(layers=tuple(LAYERS)) -> GenomeGrid:
    """Join the tables behind ``layers`` into one grid. Genomes come from
    the first table read (contiguity_cells.tsv when it is among them, for
    its assembly metadata); a genome absent from a later table is grey
    there."""
    tables: dict[str, list[dict]] = {}
    for key in layers:
        rel = LAYERS[key][0]
        if rel not in tables:
            tables[rel] = G.read_tsv(rel)
    first = tables.get(CONTIG) or next(iter(tables.values()))
    genomes, index = [], {}
    for r in first:
        if r["accession"] not in index:
            index[r["accession"]] = len(genomes)
            genomes.append(_genome(r))
    col = {c: j for j, c in enumerate(CELLS)}
    grid = GenomeGrid(genomes)
    for key in layers:
        if key == "lesion":
            grid.layers[key] = _L.layer([g.accession for g in genomes], CELLS,
                                        tables[_L.LOCI])
            continue
        m = np.full((len(genomes), len(CELLS)), "", dtype=object)
        for r in tables[LAYERS[key][0]]:
            i, j = index.get(r["accession"]), col.get(r["cell"])
            if i is None or j is None:
                continue
            m[i, j] = _cell_value(key, r)
        grid.layers[key] = m
    return grid


def _cell_value(key: str, r: dict) -> str:
    if key == "search":
        return r["grade"]
    if key == "miss":
        if r["control"] not in ("itpr_present", "ryr_sister"):
            return "not a control"
        return "missed" if r["false_negative"] == "1" else "found"
    if key == "state":
        return r["state"]
    return recovery_channel(r)


def _n50_key(g: Genome) -> float:
    return -g.contig_n50 if math.isfinite(g.contig_n50) else math.inf


#: Row orders: label → sort key over a Genome.
ORDERS = {
    "contig N50": lambda g: (_n50_key(g), g.organism),
    "class, then N50": lambda g: (g.vclass, _n50_key(g), g.organism),
    "organism": lambda g: (g.organism,),
}


def order(grid: GenomeGrid, by: str = "contig N50") -> GenomeGrid:
    """The grid with its genomes sorted by one of :data:`ORDERS`."""
    key = ORDERS[by]
    rows = sorted(range(grid.n), key=lambda i: key(grid.genomes[i]))
    return grid.subset(rows)
