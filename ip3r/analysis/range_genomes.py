"""Paper 1's genome sweep (S23) one genome at a time, placed in S20's clades.

The Range tab's clade bars say *where* the receptor is; S23 asked 194
non-vertebrate assemblies again, each with a positive control. This module
joins S23's per-genome ledgers into one :class:`GenomeRow` per assembly —
its control verdict, what the copy ledger found, its complete gene models
counted from ``copies.tsv``, and whether its contig N50 reaches the
*genome's own* contiguity bar (S23 sets one per group from measured gene
spans) — and files it under the S20 clade it belongs to, so a clicked clade
can be drawn genome by genome, the way the Genomes tab draws Paper 3.

A genome is placed by its taxon id in S20's taxonomy when S20 swept it, else
by its phylum, else its class, when either is an S20 clade; otherwise it is
in :data:`NOT_SWEPT` (a phylum S20 had no proteome of). Nothing is guessed.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from ..core import genes_data as G
from . import range_table as RT

__all__ = ["GenomeRow", "genome_rows", "in_clade", "clades_with_genomes",
           "NOT_SWEPT", "STATUS_ORDER", "VERDICT_ORDER"]

NOT_SWEPT = "not swept in S20"

#: The copy ledger's statuses, most to least of a gene.
STATUS_ORDER = ("found_annotated", "found_no_annotation", "found_unannotated",
                "fragment_only", "assembly_gap", "tblastn_trace", "no_locus")
#: The control ledger's verdicts; only the first two show the search works.
VERDICT_ORDER = RT.CONTROLLED + ("controlled_partial", "no_control_bait")


@dataclass(frozen=True)
class GenomeRow:
    accession: str
    organism: str
    group: str
    phylum: str
    klass: str
    clade: str               # the S20 clade it is filed under, or NOT_SWEPT
    placed_by: str           # taxid / phylum / class / ""
    verdict: str             # control ledger
    status: str              # copy-number ledger
    copies: int              # complete gene models (rows of copies.tsv)
    contig_n50: float
    bar_bp: float            # the genome's own contiguity bar

    @property
    def controlled(self) -> bool:
        return self.verdict in RT.CONTROLLED

    @property
    def spans_gene(self) -> bool:
        return self.contig_n50 >= self.bar_bp


def _place(m: dict, tax: dict, clades: set[str]) -> tuple[str, str]:
    t = tax.get(m["taxid"])
    if t is not None:
        return t["clade"], "taxid"
    for rank, key in (("phylum", "phylum"), ("class", "class")):
        if m[key] in clades:
            return m[key], rank
    return NOT_SWEPT, ""


def genome_rows() -> list[GenomeRow]:
    """Every S23 genome, from the manifest and the three per-genome ledgers."""
    tax = {r["taxid"]: r for r in G.read_tsv(RT.TAXONOMY)}
    clades = {r["clade"] for r in tax.values()}
    verdict = {r["accession"]: r["verdict"] for r in G.read_tsv(RT.CONTROLS)}
    status = {r["accession"]: r["status"] for r in G.read_tsv(RT.COPY_LEDGER)}
    copies = Counter(r["accession"] for r in G.read_tsv(RT.COPIES))
    out = []
    for m in G.read_tsv(RT.MANIFEST):
        clade, how = _place(m, tax, clades)
        out.append(GenomeRow(m["accession"], m["organism"], m["group"],
                             m["phylum"], m["class"], clade, how,
                             verdict.get(m["accession"], ""),
                             status.get(m["accession"], ""),
                             copies.get(m["accession"], 0),
                             float(m["contig_n50"] or "nan"),
                             float(m["contiguity_bar_bp"] or "nan")))
    return out


def in_clade(rows, clade: str) -> list[GenomeRow]:
    """The genomes filed under ``clade``: most copies first, then name."""
    return sorted((r for r in rows if r.clade == clade),
                  key=lambda r: (-r.copies, r.organism))


def clades_with_genomes(rows) -> list[tuple[str, int]]:
    """(clade, genomes) for every clade holding at least one S23 genome,
    most genomes first; :data:`NOT_SWEPT` last."""
    n = Counter(r.clade for r in rows)
    swept = sorted(((c, k) for c, k in n.items() if c != NOT_SWEPT),
                   key=lambda x: (-x[1], x[0]))
    return swept + ([(NOT_SWEPT, n[NOT_SWEPT])] if n.get(NOT_SWEPT) else [])
