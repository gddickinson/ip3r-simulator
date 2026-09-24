"""Papers 3 and 4 re-derived cell by cell from the genome × paralog grid.

Paper 3's method control (how often the search misses a gene known to be
present, and whether that tracks contiguity) is recomputed from
``contiguity_cells.tsv`` with this project's statistics and its own split at
the registered contiguity bar, and compared with the published test table
and summary. Paper 4's recovery channels are rebuilt from the count columns
of ``gene_recovery.tsv`` (:func:`.genome_grid.recovery_channel`), not read
from its label. Each outcome carries the per-cell arrays its exhibit draws.
"""

from __future__ import annotations

import math
import re

import numpy as np

from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from . import genome_grid as GG
from .checks import agree, register
from .stats import fisher_exact, logistic_fit, mann_whitney_greater, wilson

TESTS = "methods/contiguity_tests.tsv"
SUMMARY = "methods/contiguity_summary.json"
BY_CELL = "methods/gene_recovery_by_cell.tsv"
LEDGER4 = "../papers/archive/claims_check.tsv"
SERIES = ("itpr_present", "ryr_sister")


def _series(rows, name):
    return [r for r in rows if r["control"] == name]


def _missed(rows):
    return [r for r in rows if r["false_negative"] == "1"]


def _n50(rows):
    return np.array([float(r["contig_n50"]) for r in rows])


def _ints(text: str) -> list[int]:
    return [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]*", text)]


def _printed(text: str, name: str) -> str:
    return re.search(rf"{name}=([-\d.eE+]+)", text).group(1)


def _number(text: str, name: str) -> float:
    return float(_printed(text, name))


def _as_printed(value: float, printed: str) -> bool:
    """Does ``value`` round to ``printed`` at the precision it is printed?"""
    places = len(printed.split(".")[1]) if "." in printed else 0
    return round(value, places) == float(printed)


def _tests() -> dict[tuple[str, str], dict]:
    return {(r["test"], r["series"]): r for r in G.read_tsv(TESTS)}


def _grid_data(rows) -> dict:
    """Per-genome N50 and per-cell miss flags, for the exhibits."""
    by: dict[str, dict] = {}
    for r in rows:
        g = by.setdefault(r["accession"], {"n50": float(r["contig_n50"])})
        g[r["cell"]] = r["false_negative"] == "1" if r["control"] else None
    return {"n50": [g["n50"] for g in by.values()],
            "missed": {c: [g.get(c) for g in by.values()] for c in GG.CELLS},
            "bar": _P.value("genomes.contiguity_bar_bp")}


@register("P3.miss_by_contiguity", "retention",
          "A missed cell's assembly has a median contig N50 of 23,460 bp "
          "against 3,396,515 bp for a found one; chromosome-level assemblies "
          "miss 3 of 512 IP3R cells and 0 of 172 ryanodine cells; below the "
          "142,212 bp bar the miss rate is 0.3917/0.4333/0.3 for ITPR1/2/3, "
          "and above it 189 genomes carry 563 IP3R cells of which 5 are "
          "missed, and the ryanodine series misses none of 189.",
          "The two control series of contiguity_cells.tsv split by this "
          "project at the registered bar (contig N50 ≥ "
          "genomes.contiguity_bar_bp) rather than by the table's spans_gene "
          "column; medians, chromosome-level misses and per-paralogue rates "
          "counted, and compared with contiguity_tests.tsv and "
          "contiguity_summary.json.",
          "rederived", (GG.CONTIG, TESTS, SUMMARY))
def miss_by_contiguity():
    rows = G.read_tsv(GG.CONTIG)
    tests, summ = _tests(), G.read_json(SUMMARY)
    tol = _P.value("check.stat_tol")
    ok, pub, got = True, [], []
    for s in SERIES:
        cells = _series(rows, s)
        miss = _missed(cells)
        found = [r for r in cells if r["false_negative"] == "0"]
        med = (float(np.median(_n50(miss))), float(np.median(_n50(found))))
        want = _ints(tests[("mannwhitney_contig_n50_found_vs_missed", s)]["effect"])
        ok &= (round(med[1]), round(med[0])) == tuple(want)
        chrom = [r for r in cells if r["level"] == "Chromosome"]
        c_got = (len(_missed(chrom)), len(chrom))
        c_want = tuple(_ints(tests[("fisher_chromosome_vs_lower", s)]["n"])[:2])
        ok &= c_got == c_want
        above = [r for r in cells if GG.above_bar(float(r["contig_n50"]))]
        a_got = (len({r["accession"] for r in above}), len(above), len(_missed(above)))
        a_sum = summ["series"][s]["above_bar"]
        a_want = (summ["floor"][s]["at_d4_bar"]["genomes_retained"], a_sum["n"],
                  a_sum["n_false_negative"])
        ok &= a_got == a_want
        pub.append(f"{s}: median {want[1]:,} vs {want[0]:,}; chromosome "
                   f"{c_want[0]}/{c_want[1]}; above bar {a_want[0]} genomes, "
                   f"{a_want[2]}/{a_want[1]} missed")
        got.append(f"{s}: median {med[0]:,.0f} vs {med[1]:,.0f}; chromosome "
                   f"{c_got[0]}/{c_got[1]}; above bar {a_got[0]} genomes, "
                   f"{a_got[2]}/{a_got[1]} missed")
    rates_p, rates_g = [], []
    for cell in ("ITPR1", "ITPR2", "ITPR3"):
        below = [r for r in rows if r["cell"] == cell and r["control"] == "itpr_present"
                 and not GG.above_bar(float(r["contig_n50"]))]
        rate = len(_missed(below)) / len(below)
        want = float(tests[("false_negative_rate_below_bar", cell)]["effect"])
        ok &= abs(rate - want) <= tol
        rates_p.append(f"{cell} {want:.4f}")
        rates_g.append(f"{cell} {rate:.4f} ({len(_missed(below))}/{len(below)})")
    split = sum(GG.above_bar(float(r["contig_n50"])) != (r["spans_gene"] == "1")
                for r in rows)
    ok &= split == 0
    return agree(ok, "; ".join(pub) + "; below bar " + ", ".join(rates_p),
                 "; ".join(got) + "; below bar " + ", ".join(rates_g),
                 f"the bar at {_P.value('genomes.contiguity_bar_bp'):,.0f} bp splits "
                 f"{'every cell as' if not split else f'{split} cells unlike'} the "
                 "table's spans_gene column", **_grid_data(rows))


@register("P3.contiguity_tests", "retention",
          "The miss rate is 140/923 (0.1517) for the IP3R cells and 42/309 "
          "(0.1359) for the ryanodine control, indistinguishable (Fisher "
          "p = 0.578); the odds of finding the gene rise 8.10-fold per "
          "tenfold increase in contig N50 in the IP3R series and 19.99-fold "
          "in the ryanodine series.",
          "Every test in contiguity_tests.tsv recomputed from the per-cell "
          "table with this project's statistics: Wilson intervals, Fisher "
          "exact tests (the two series; chromosome-level against the rest), "
          "a logistic regression of found on log10 contig N50 (IRLS, Wald "
          "p), and the one-sided Mann-Whitney U of found against missed N50.",
          "rederived", (GG.CONTIG, TESTS))
def contiguity_tests():
    rows = G.read_tsv(GG.CONTIG)
    tests = _tests()
    tol, ptol, alpha = (_P.value("check.stat_tol"), _P.value("check.log_p_tol"),
                        _P.value("check.alpha"))

    def p_ok(a, b):
        return (a == b == 0) or (a > 0 and b > 0 and
                                 abs(math.log10(a) - math.log10(b)) <= ptol)

    ok, lines, fits = True, [], {}
    counts = {}
    for s in SERIES:
        cells = _series(rows, s)
        k, n = len(_missed(cells)), len(cells)
        counts[s] = (k, n)
        row = tests[("false_negative_rate", s)]
        lo, hi = wilson(k, n, alpha)
        ci = [float(x) for x in re.findall(r"\d\.\d+", row["detail"])]
        ok &= abs(k / n - float(row["effect"])) <= tol and \
            abs(lo - ci[0]) <= tol and abs(hi - ci[1]) <= tol
        x = np.log10(_n50(cells))
        y = np.array([r["false_negative"] == "0" for r in cells], float)
        f = logistic_fit(x, y)
        row = tests[("logit_found~log10(contig_n50)", s)]
        ok &= f["converged"] and _as_printed(f["b1"], _printed(row["effect"], "beta")) \
            and _as_printed(math.exp(f["b1"]), _printed(row["detail"], "OR per 10x")) \
            and p_ok(f["p"], float(row["p"]))
        fits[s] = {"b0": f["b0"], "b1": f["b1"], "x": x.tolist(), "y": y.tolist()}
        found = _n50([r for r in cells if r["false_negative"] == "0"])
        mw = mann_whitney_greater(found, _n50(_missed(cells)))
        row = tests[("mannwhitney_contig_n50_found_vs_missed", s)]
        ok &= mw["u"] == _number(row["detail"], "U") and p_ok(mw["p"], float(row["p"]))
        chrom = [r for r in cells if r["level"] == "Chromosome"]
        rest = [r for r in cells if r["level"] != "Chromosome"]
        fe = fisher_exact([[len(_missed(chrom)), len(chrom) - len(_missed(chrom))],
                           [len(_missed(rest)), len(rest) - len(_missed(rest))]])
        row = tests[("fisher_chromosome_vs_lower", s)]
        ok &= _as_printed(fe["odds"], _printed(row["effect"], "odds")) \
            and p_ok(fe["p"], float(row["p"]))
        lines.append(f"{s}: {k}/{n} missed ({k / n:.4f}, CI {lo:.4f}–{hi:.4f}); "
                     f"OR per 10× {math.exp(f['b1']):.2f} (β {f['b1']:.3f}, p {f['p']:.3g}); "
                     f"U {mw['u']:.0f} (one-sided p {mw['p']:.3g}); chromosome vs rest "
                     f"p {fe['p']:.3g}")
    (a, n1), (c, n2) = counts["itpr_present"], counts["ryr_sister"]
    fe = fisher_exact([[a, n1 - a], [c, n2 - c]])
    row = tests[("fisher_itpr_vs_ryr_false_negative", "the two control series")]
    ok &= p_ok(fe["p"], float(row["p"])) and fe["p"] > alpha
    lines.append(f"ITPR vs RyR Fisher odds {fe['odds']:.2f}, p {fe['p']:.3f}")
    pub = (f"{len(tests)} tests; e.g. p {row['p']} for ITPR vs RyR; OR per 10× "
           + " and ".join(_number(tests[("logit_found~log10(contig_n50)", s)]["detail"],
                                  "OR per 10x").__format__(".2f") for s in SERIES))
    return agree(ok, pub, "; ".join(lines), fits=fits,
                 bar=_P.value("genomes.contiguity_bar_bp"))


@register("P4.recovery_channels", "archive",
          "257 of 309 ITPR1, 260 of 307 ITPR2 and 227 of 307 ITPR3 genes have "
          "no record that resolves to them (ryanodine control 196 of 309). "
          "Of the 744 unreachable, 289 sit in species with no reference "
          "proteome, 186 with only fragmentary records, 254 with full-length "
          "records none of which resolves to the paralogue, and 15 with no "
          "family record; 179 are held both as DNA and as a protein record.",
          "Each cell's recovery channel rebuilt from gene_recovery.tsv's "
          "count columns (reference proteome? any family record? any "
          "full-length record resolving to a cell? to this cell?), compared "
          "with the table's own channel label in every cell, with "
          "gene_recovery_by_cell.tsv, and with the numbers Paper 4's claims "
          "ledger states (AR12–AR16).",
          "rederived", (GG.RECOVERY, BY_CELL, LEDGER4))
def recovery_channels():
    rows = G.read_tsv(GG.RECOVERY)
    ours = [GG.recovery_channel(r) for r in rows]
    label_diff = sum(not (r["recovery_channel"].startswith(c) if c != GG.NO_GENE
                          else r["gene_present"] != "1")
                     for r, c in zip(rows, ours))
    by_cell = {r["cell"]: r for r in G.read_tsv(BY_CELL)}
    ok, pub, got = label_diff == 0, [], []
    per_cell = {}
    for cell in GG.CELLS:
        ch = [c for r, c in zip(rows, ours) if r["cell"] == cell and c != GG.NO_GENE]
        inv = sum(c != GG.REACHABLE for c in ch)
        want = (int(by_cell[cell]["n_protein_db_invisible"]),
                int(by_cell[cell]["n_genes_present"]))
        ok &= (inv, len(ch)) == want
        pub.append(f"{cell} {want[0]}/{want[1]}")
        got.append(f"{cell} {inv}/{len(ch)}")
        per_cell[cell] = {k: sum(c == k for c in ch) for k in
                          (GG.REACHABLE, GG.NO_PROTEOME, GG.FRAGMENTS_ONLY,
                           GG.OTHER_PARALOG, GG.NO_RECORD)}
    itpr = [c for r, c in zip(rows, ours) if r["cell"] != "RYR" and c != GG.NO_GENE]
    reasons = {"AR12": GG.NO_PROTEOME, "AR13": GG.FRAGMENTS_ONLY,
               "AR14": GG.OTHER_PARALOG, "AR15": GG.NO_RECORD, "AR16": GG.REACHABLE}
    ledger = {r["id"]: r["expected"] for r in G.read_tsv(LEDGER4) if r["id"] in reasons}
    r_pub = [int(ledger[k]) for k in reasons]
    r_got = [sum(c == v for c in itpr) for v in reasons.values()]
    ok &= r_pub == r_got
    pub.append("reasons " + "/".join(map(str, r_pub)))
    got.append("reasons " + "/".join(map(str, r_got)))
    return agree(ok, "; ".join(pub), "; ".join(got),
                 f"the rebuilt channel differs from the table's label in {label_diff} "
                 f"of {len(rows)} cells", per_cell=per_cell)
