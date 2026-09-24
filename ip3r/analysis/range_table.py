"""Paper 1's range: where in the tree of life an IP3 receptor gene exists.

Two sweeps answer it from different sides:

* **S20, proteomes** — 6,928 non-vertebrate reference proteomes searched with
  the two family profiles. ``s20_sweep/proteome_presence.tsv`` holds one row
  per proteome with its call count; ``taxonomy.tsv`` places each taxon; the
  per-record calls are in ``assignments_<group>.tsv`` (one per search group),
  and ``relaxed_hits.tsv`` holds the second-sensitivity search.
* **S23, genomes** — 194 non-vertebrate assemblies asked again, each with a
  positive control chosen for its clade. The per-genome ledgers
  (``genome_manifest_s23.tsv``, ``control_ledger.tsv``, ``copies.tsv``,
  ``copy_number_ledger.tsv``) are joined here into the per-clade absence rows
  that ``absence_at_genome.tsv`` publishes, so that table can be *checked*
  rather than read.

A clade is the taxonomy table's own ``clade`` (the phylum where UniProt gives
one, the next-deepest named group where it does not); an absence target is a
phylum or class, as S23 defines them. Nothing here draws (see
``range_figure``).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass

from ..core import genes_data as G
from ..parameters import PARAMETERS as _P

__all__ = ["PRESENCE", "TAXONOMY", "ASSIGNMENTS", "RELAXED", "CHASE", "ABSENCE",
           "MANIFEST", "CONTROLS", "COPIES", "COPY_LEDGER", "SUPERGROUPS",
           "CONTROLLED", "LINEAGES", "Proteome", "CladeRow", "GenomeAbsence",
           "load_proteomes", "clade_rows", "count_at", "taxon_calls",
           "absence_targets", "genome_absences", "copy_numbers",
           "substantial_table", "load_range", "Range"]

S20, S23 = "s20_sweep/", "s23_scope/"
PRESENCE = S20 + "proteome_presence.tsv"
TAXONOMY = S20 + "taxonomy.tsv"
ASSIGNMENTS = tuple(S20 + f"assignments_{g}.tsv" for g in
                    ("archaea", "bacteria_genus", "fungi", "metazoa_nonvert",
                     "protista_other", "viridiplantae"))
RELAXED = S20 + "relaxed_hits.tsv"
CHASE = S20 + "plant_fungal_verdicts.tsv"
ABSENCE = S23 + "absence_at_genome.tsv"
MANIFEST = S23 + "genome_manifest_s23.tsv"
CONTROLS = S23 + "control_ledger.tsv"
COPIES = S23 + "copies.tsv"
COPY_LEDGER = S23 + "copy_number_ledger.tsv"

#: Taxonomy supergroups, in the order the figure stacks them.
SUPERGROUPS = ("Metazoa (non-vertebrate)", "Fungi", "Viridiplantae", "Amoebozoa",
               "Discoba", "SAR", "Eukaryota (other)", "Archaea", "Bacteria")

#: Control-ledger verdicts under which a genome's search is shown to work.
#: ``controlled_partial`` and ``no_control_bait`` are not.
CONTROLLED = ("controlled_cross_kingdom", "controlled_by_target")

#: The paper's four lineages for the relaxed-sensitivity table: name →
#: (phyla, is it a claimed absence?).
LINEAGES = {"land plants": (("Streptophyta",), True),
            "Chlorophyta": (("Chlorophyta",), False),
            "Dikarya": (("Ascomycota", "Basidiomycota"), True),
            "Mucoromycota": (("Mucoromycota",), False)}


@dataclass(frozen=True)
class Proteome:
    upid: str
    organism: str
    taxid: str
    supergroup: str
    domain: str
    phylum: str
    klass: str
    clade: str
    n_itpr: int

    @property
    def present(self) -> bool:
        return self.n_itpr > 0

    def rank(self, rank: str) -> str:
        return {"phylum": self.phylum, "class": self.klass, "clade": self.clade}[rank]


@dataclass(frozen=True)
class CladeRow:
    clade: str
    supergroup: str
    domain: str
    n: int
    present: int

    @property
    def fraction(self) -> float:
        return self.present / self.n


@dataclass(frozen=True)
class GenomeAbsence:
    """One absence target asked of the genomes, rebuilt from the ledgers."""
    rank: str
    clade: str
    proteomes: int              # swept in S20
    genomes: int                # assemblies of this clade in the S23 sweep
    controlled: int             # of them, with a working positive control
    with_full: int              # of them, carrying a complete gene model
    trace_only: int             # of them, with translated-search traces only

    @property
    def holds(self) -> bool:
        return self.genomes > 0 and self.controlled == self.genomes \
            and self.with_full == 0


def load_proteomes() -> list[Proteome]:
    """Every swept proteome, placed by the taxonomy table."""
    tax = {r["taxid"]: r for r in G.read_tsv(TAXONOMY)}
    out = []
    for r in G.read_tsv(PRESENCE):
        t = tax[r["taxid"]]
        out.append(Proteome(r["upid"], r["organism"], r["taxid"], t["group"],
                            t["domain"], t["phylum"], t["class"], t["clade"],
                            int(r["n_itpr"])))
    return out


def clade_rows(proteomes: list[Proteome]) -> list[CladeRow]:
    """Proteomes and those with a call, per clade; supergroup order, then
    most-swept first."""
    acc: dict[str, list] = {}
    for p in proteomes:
        a = acc.setdefault(p.clade, [p.supergroup, p.domain, 0, 0])
        a[2] += 1
        a[3] += p.present
    rows = [CladeRow(c, s, d, n, k) for c, (s, d, n, k) in acc.items()]
    rank = {s: i for i, s in enumerate(SUPERGROUPS)}
    return sorted(rows, key=lambda r: (rank.get(r.supergroup, len(rank)), -r.n, r.clade))


def count_at(proteomes: list[Proteome], rank: str, name: str) -> tuple[int, int]:
    """(proteomes swept, proteomes with a call) in one named phylum or class."""
    ps = [p for p in proteomes if p.rank(rank) == name]
    return len(ps), sum(p.present for p in ps)


def taxon_calls() -> Counter:
    """ITPR calls per taxon id, counted record by record from the six
    assignment tables — a second route to presence, independent of the
    presence table's own count column."""
    calls: Counter = Counter()
    for rel in ASSIGNMENTS:
        for r in G.read_tsv(rel):
            if r["assignment"] == "ITPR":
                calls[r["taxon_id"]] += 1
    return calls


def absence_targets(proteomes: list[Proteome]) -> list[tuple[str, str]]:
    """S23's rule G3, from its one-sentence statement: every eukaryotic
    phylum or class with at least ``range.absence_min_proteomes`` swept
    proteomes and no call among them."""
    floor = _P.value("range.absence_min_proteomes")
    tally: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    for p in proteomes:
        if p.domain != "Eukaryota":
            continue
        for rank in ("phylum", "class"):
            if p.rank(rank):
                t = tally[(rank, p.rank(rank))]
                t[0] += 1
                t[1] += p.present
    return sorted(k for k, (n, k_) in tally.items() if n >= floor and k_ == 0)


def copy_numbers() -> dict[str, tuple[str, int]]:
    """accession → (organism, complete gene models) for every S23 genome,
    counted from the rows of ``copies.tsv`` (zero where it has none)."""
    n = Counter(r["accession"] for r in G.read_tsv(COPIES))
    return {m["accession"]: (m["organism"], n.get(m["accession"], 0))
            for m in G.read_tsv(MANIFEST)}


def genome_absences(proteomes: list[Proteome],
                    targets: list[tuple[str, str]]) -> list[GenomeAbsence]:
    """Each target's genomes, how many are controlled, how many carry a
    complete gene model and how many only a trace — from the per-genome
    ledgers, not from ``absence_at_genome.tsv``."""
    manifest = G.read_tsv(MANIFEST)
    control = {r["accession"]: r["verdict"] for r in G.read_tsv(CONTROLS)}
    status = {r["accession"]: r["status"] for r in G.read_tsv(COPY_LEDGER)}
    copies = copy_numbers()
    out = []
    for rank, name in targets:
        gs = [m["accession"] for m in manifest if m[rank] == name]
        out.append(GenomeAbsence(
            rank, name, count_at(proteomes, rank, name)[0], len(gs),
            sum(control.get(a) in CONTROLLED for a in gs),
            sum(copies[a][1] > 0 for a in gs),
            sum(status.get(a) == "tblastn_trace" for a in gs)))
    return out


def substantial_table(proteomes: list[Proteome] | None = None) -> dict:
    """The relaxed-sensitivity table: lineage → profile → (targets reported
    at the relaxed E-value, of those substantial). Substantial = full
    E-value ≤ ``range.substantial_evalue`` and model coverage ≥
    ``range.substantial_coverage``; targets are counted once per accession.
    Also returns the substantial non-MIR hits inside claimed absences."""
    tax = {r["taxid"]: r["phylum"] for r in G.read_tsv(TAXONOMY)}
    e_max = _P.value("range.substantial_evalue")
    cov = _P.value("range.substantial_coverage")
    phylum_of = {ph: name for name, (phyla, _) in LINEAGES.items() for ph in phyla}
    seen: dict = defaultdict(lambda: defaultdict(lambda: [set(), set()]))
    exceptions: dict[str, set[str]] = defaultdict(set)
    for r in G.read_tsv(RELAXED):
        lineage = phylum_of.get(tax.get(r["taxon_id"], ""))
        if lineage is None:
            continue
        cell = seen[lineage][r["profile"]]
        cell[0].add(r["accession"])
        if r["full_evalue"] and float(r["full_evalue"]) <= e_max \
                and float(r["hmm_coverage"]) >= cov:
            cell[1].add(r["accession"])
            if LINEAGES[lineage][1] and r["profile"] != "PF02815":
                exceptions[r["accession"]].add(r["profile"])
    table = {lin: {p: (len(v[0]), len(v[1])) for p, v in prof.items()}
             for lin, prof in seen.items()}
    return {"table": table, "exceptions": {a: sorted(p) for a, p in exceptions.items()}}


@dataclass
class Range:
    """Everything the Range tab draws, loaded in one pass on a worker."""
    proteomes: list[Proteome]
    clades: list[CladeRow]
    absences: list[GenomeAbsence]
    copies: dict[str, tuple[str, int]]

    def absences_in(self, clade: str) -> list[GenomeAbsence]:
        """Absence targets that are this clade or sit inside it."""
        inside = {p.klass for p in self.proteomes if p.clade == clade}
        return [a for a in self.absences if a.clade == clade
                or (a.rank == "class" and a.clade in inside)]


def load_range() -> Range:
    ps = load_proteomes()
    return Range(ps, clade_rows(ps), genome_absences(ps, absence_targets(ps)),
                 copy_numbers())
