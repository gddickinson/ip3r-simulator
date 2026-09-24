"""Evolutionary findings (papers 1-4) and the publication's claims ledgers.

These are ``rederived`` checks where the committed per-cell or per-genome
table allows the headline count to be recomputed, and ``read`` checks where
the only thing this project can do is confirm that a stated test outcome is
what the table says. The origin question is asked of the committed tree with
this project's own Newick reader.
"""

from __future__ import annotations

from ..core import genes_data as G
from .checks import agree, register
from .newick import mrca, leaves, parse

TREE = "phylogeny/rooted.nwk"
MEMBERS = "phylogeny/membership_audit.tsv"
AU = "phylogeny/au_test.tsv"
TELEOST = "duplication/teleost_copies.tsv"
MATRIX = "loss_dynamics/character_matrix.tsv"
CONTIG = "methods/contiguity_cells.tsv"
DOLLO = "loss_counts/dollo_counts.tsv"
RECOVERY = "methods/gene_recovery.tsv"
ABSENCE = "s23_scope/absence_at_genome.tsv"


def _core_tips() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for r in G.read_tsv(MEMBERS):
        if r["rule"] == "core":
            out.setdefault(r["assigned_to"], set()).add(r["label"])
    return out


def sister_pair() -> tuple[str, str] | None:
    """The two paralogs whose union is a clade excluding the third."""
    root = parse(G.read_text(TREE))
    core = _core_tips()
    names = sorted(core)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            third = next(n for n in names if n not in (a, b))
            under = {lf.label for lf in leaves(mrca(root, core[a] | core[b]))}
            if not (under & core[third]):
                return a, b
    return None


@register("P2.sister_pair", "origin",
          "ITPR2 and ITPR3 are sisters: ITPR1 split first, on the vertebrate "
          "stem, and ITPR2 from ITPR3 later, on the gnathostome stem.",
          "The committed RyR-rooted ML tree parsed with this project's Newick "
          "reader; for each pair of paralogs, the smallest clade holding "
          "both paralogs' core tips is tested for the third's.",
          "rederived", (TREE, MEMBERS))
def sister():
    pair = sister_pair()
    found = " + ".join(pair) if pair else "no pair forms a clade without the third"
    return agree(pair == ("ITPR2", "ITPR3"), "ITPR2 + ITPR3", found,
                 newick=G.read_text(TREE))


@register("P2.au_test", "origin",
          "The AU test retains only the ITPR2 + ITPR3 pairing; both "
          "alternatives are rejected.",
          "p_AU read for each constrained hypothesis in au_test.tsv.",
          "read", (AU,))
def au_test():
    rows = {r["tree"]: r for r in G.read_tsv(AU)}
    p = {k: float(v["p_AU"]) for k, v in rows.items() if k.startswith("H")}
    kept = sorted(k for k, v in p.items() if v >= 0.05)
    return agree(kept == ["H3_23"], "only H3_23 retained",
                 ", ".join(f"{k} p_AU = {v:.3g}" for k, v in sorted(p.items())))


@register("P2.teleost_itpr1", "origin",
          "In teleosts one ancestral duplication doubled ITPR1 alone: ITPR1 "
          "is two-copy across the 3R radiation, ITPR2 and ITPR3 are not.",
          "Per 3R genome, copy counts per paralog from teleost_copies.tsv; "
          "the share of genomes carrying two or more copies of each.",
          "rederived", (TELEOST,))
def teleost():
    rows = [r for r in G.read_tsv(TELEOST) if r["group"] == "teleost"]
    share = {g: sum(int(r[f"{g}_copies"]) >= 2 for r in rows) / len(rows)
             for g in ("ITPR1", "ITPR2", "ITPR3")}
    ok = share["ITPR1"] >= 0.5 and share["ITPR2"] < 0.5 and share["ITPR3"] < 0.5
    return agree(ok, "ITPR1 two-copy, ITPR2/3 single",
                 f"{len(rows)} teleost genomes; ≥ 2 copies: " + ", ".join(
                     f"{g} {v:.0%}" for g, v in share.items()), share=share)


@register("P3.no_absent_cells", "retention",
          "No genome × paralog cell of 927 reaches a state that licenses a "
          "loss.",
          "States counted in the S15a character matrix.",
          "rederived", (MATRIX,))
def no_absent():
    rows = G.read_tsv(MATRIX)
    states: dict[str, int] = {}
    for r in rows:
        states[r["state"]] = states.get(r["state"], 0) + 1
    ok = len(rows) == 927 and states.get("absent", 0) == 0
    return agree(ok, "927 cells, 0 absent", f"{len(rows)} cells; " + ", ".join(
        f"{k} {v}" for k, v in sorted(states.items(), key=lambda x: -x[1])),
        states=states)


@register("P3.dollo_zero", "retention",
          "Dollo parsimony places zero losses on the IP3 receptor family "
          "character.",
          "The base-setting family row of dollo_counts.tsv read.",
          "read", (DOLLO,))
def dollo_zero():
    base = [r for r in G.read_tsv(DOLLO) if r["coding"] == "family"
            and r["is_base"] == "1"]
    vals = {(r["dollo_losses_min"], r["dollo_losses_max"]) for r in base}
    return agree(vals == {("0", "0")}, "0", f"{len(base)} base row(s): {vals}")


@register("P3.false_negatives", "retention",
          "The search behind that zero misses 140 of 923 known-present cells "
          "overall but only 5 of 563 in assemblies contiguous enough to hold "
          "the gene.",
          "False negatives counted over the ITPR control cells, overall and "
          "among cells whose assembly spans the gene.",
          "rederived", (CONTIG,))
def false_negatives():
    rows = [r for r in G.read_tsv(CONTIG) if r["control"] == "itpr_present"]
    fn = sum(r["false_negative"] == "1" for r in rows)
    span = [r for r in rows if r["spans_gene"] == "1"]
    fn_span = sum(r["false_negative"] == "1" for r in span)
    ok = (fn, len(rows), fn_span, len(span)) == (140, 923, 5, 563)
    return agree(ok, "140/923 overall, 5/563 contiguous",
                 f"{fn}/{len(rows)} overall, {fn_span}/{len(span)} contiguous")


@register("P4.unreachable", "archive",
          "744 of the 923 IP3 receptor genes the sweep demonstrates are "
          "reachable by no protein-database search.",
          "Present ITPR genes counted in gene_recovery.tsv (the RyR sister "
          "cell excluded), and those whose only recovery channel is the "
          "genome.",
          "rederived", (RECOVERY,))
def unreachable():
    rows = [r for r in G.read_tsv(RECOVERY) if r["gene_present"] == "1"
            and r["cell"] in ("ITPR1", "ITPR2", "ITPR3")]
    genome_only = sum(r["recovery_channel"].startswith("genome_only") for r in rows)
    ok = (genome_only, len(rows)) == (744, 923)
    return agree(ok, "744 / 923", f"{genome_only} / {len(rows)} "
                 f"({genome_only / len(rows):.1%})")


@register("P1.absences", "range",
          "35 clade-level absences are confirmed in genome assemblies, each "
          "genome carrying a positive control chosen for its clade.",
          "Clades whose controlled genomes carry no full-length receptor "
          "counted in absence_at_genome.tsv, independently of its verdict "
          "column.",
          "rederived", (ABSENCE,))
def absences():
    rows = G.read_tsv(ABSENCE)
    held = [r for r in rows if int(r["genomes_controlled"]) > 0
            and r["genomes_with_full_itpr"] == "0"]
    agrees_verdict = all(r["verdict"].startswith("absence holds") for r in held)
    ok = len(held) == 35 and agrees_verdict
    return agree(ok, "35", f"{len(held)} clades (of {len(rows)}); verdict column "
                 f"{'agrees' if agrees_verdict else 'DISAGREES'}")


@register("LEDGER.claims", "ledger",
          "Every load-bearing number in the manuscript (287), the thesis (245) "
          "and the paper series (713) re-verified against its source table.",
          "The three claims_check.tsv ledgers read and their verdicts counted.",
          "read", ("../manuscript/claims_check.tsv", "../thesis/claims_check.tsv",
                   "../papers/claims_check.tsv"))
def ledgers():
    parts, ok = [], True
    for rel, want in (("../manuscript/claims_check.tsv", 287),
                      ("../thesis/claims_check.tsv", 245),
                      ("../papers/claims_check.tsv", 713)):
        rows = G.read_tsv(rel)
        bad = [r["id"] for r in rows if r["verdict"] != "ok"]
        ok &= len(rows) == want and not bad
        parts.append(f"{rel.split('/')[1]} {len(rows) - len(bad)}/{len(rows)} ok")
    return agree(ok, "287, 245 and 713, all ok", "; ".join(parts))
