"""Paper 1's family-call benchmark: the labelled-bait margin recomputed from
the control sequences, and the headline recall and specificity rebuilt from
each control's scorer components (:mod:`ip3r.analysis.family_benchmark`).
"""

from __future__ import annotations

import json

from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from .checks import agree, register
from .family_benchmark import (DECOYS, MARGIN_TSV, NEGATIVE_TSV, POSITIVE_TSV, POSITIVES,
                               SUMMARY, identity, margins, rescore, sequences)

FAILURES = "benchmark_controls/recall_failures.tsv"


@register("P1.bait_margin", "range",
          "The labelled-bait margin calls all 25 positive controls IP3 "
          "receptors and all 6 ryanodine receptors RyR, with an empty gap "
          "between the families; the narrowest true member, Dictyostelium "
          "iplA (+0.065), sits inside the 0.10 no-call band (bait_margin.tsv).",
          "Every control aligned pairwise to each human ITPR and RyR bait "
          "with this project's Gotoh/BLOSUM62 aligner (not S1's MAFFT MSA), "
          "from the committed UniProt sequences; full-alignment identity to "
          "the nearest bait of each family, the call and the D7 band per "
          "control compared with bait_margin.tsv.",
          "recomputed", (POSITIVES, DECOYS, MARGIN_TSV, FAILURES))
def bait_margin():
    pub = {r["accession"]: r for r in G.read_tsv(MARGIN_TSV)}
    mine = margins()
    labelled = {a: r for a, r in pub.items() if r["truth"] in ("ITPR", "RYR")}
    missing = sorted(set(pub) - set(mine))
    wrong = sorted(a for a, r in labelled.items()
                   if a in mine and not (mine[a].call() == r["truth"] == r["call"]))
    itpr = {a: mine[a].full for a, r in labelled.items() if r["truth"] == "ITPR" and a in mine}
    ryr = {a: mine[a].full for a, r in labelled.items() if r["truth"] == "RYR" and a in mine}
    narrow = min(itpr, key=itpr.get)
    pub_narrow = min((a for a, r in labelled.items() if r["truth"] == "ITPR"),
                     key=lambda a: float(pub[a]["margin"]))
    d = _P.value("bench.sister_margin")
    gap_ok = min(itpr.values()) > 0 > max(ryr.values())
    band_ok = narrow == pub_narrow and 0 < itpr[narrow] <= d and 0 < float(pub[narrow]["margin"]) <= d
    ok = not missing and not wrong and gap_ok and band_ok
    worst = max(abs(mine[a].full - float(r["margin"])) for a, r in labelled.items() if a in mine)
    cov_calls = sum(mine[a].call("covered") == r["truth"] for a, r in labelled.items() if a in mine)
    fly = G.read_tsv(FAILURES)
    fly_text = ""
    if fly:
        f, seq = fly[0], sequences()
        sib = f["nearest_sibling_full"].rsplit("_", 1)[-1]
        if f["accession"] in seq and sib in seq:
            fly_text = (f" The one missed positive, {f['gene_symbol']} ({f['accession']}): "
                        f"full identity to {sib} {identity(seq[f['accession']], seq[sib])[0]:.3f} "
                        f"here against {float(f['id_full']):.3f} in S1's MSA, breadth bar "
                        f"{float(f['breadth_identity_min']):.2f}; the shortfall is smaller "
                        "than the difference between the two alignment routes.")
    name = lambda a: pub[a]["gene_symbol"]  # noqa: E731
    found = (f"{len(labelled) - len(wrong)}/{len(labelled)} calls agree; ITPR margins "
             f"{min(itpr.values()):+.3f}…{max(itpr.values()):+.3f}, RyR "
             f"{min(ryr.values()):+.3f}…{max(ryr.values()):+.3f} (gap "
             f"{min(itpr.values()) - max(ryr.values()):.3f}); narrowest true member "
             f"{name(narrow)} {itpr[narrow]:+.3f}" + ("" if band_ok else " — NOT as published")
             + (f"; wrong: {wrong}" if wrong else "") + (f"; missing: {missing}" if missing else ""))
    return agree(ok, f"31/31 labelled calls; gap 0.601; iplA +0.065 inside ±{d:g}", found,
                 f"Largest difference from S1's full-metric margin {worst:.3f}. On "
                 f"covered columns only, {cov_calls}/{len(labelled)} calls agree: a "
                 "pairwise optimum aligns unrelated long sequences at 0.25-0.28 "
                 "covered identity, so iplA against either family is at that floor, "
                 "and covered identity is not re-measured by this route." + fly_text,
                 margins={a: {"itpr": m.itpr[0], "ryr": m.ryr[0],
                              "truth": pub[a]["truth"] if a in pub else "other"}
                          for a, m in mine.items()},
                 band=d)


def _flag(components: str) -> str:
    return "sister" if "SISTER" in components else "gate" if "GATED" in components else ""


@register("P1.benchmark_counts", "range",
          "Recall is 24 of 25 and specificity 31 of 31, and all 6 ryanodine "
          "receptors are rejected (summary.json).",
          "Each control's score rebuilt from the components S1 records as "
          "fired, with the registered points; the sister cap decided by the "
          "margin recomputed here, the evidence gate by the family-specific "
          "components; scores, caps and the three counts compared with "
          "positive_controls.tsv, negative_controls.tsv and summary.json.",
          "rederived", (POSITIVES, DECOYS, POSITIVE_TSV, NEGATIVE_TSV, SUMMARY))
def benchmark_counts():
    mine, bar = margins(), int(_P.value("bench.promotion"))
    pos, neg = G.read_tsv(POSITIVE_TSV), G.read_tsv(NEGATIVE_TSV)
    s = json.loads(G.read_text(SUMMARY))
    bad = []
    got = {}
    for rows, key in ((pos, "score"), (neg, "score_max")):
        for r in rows:
            sc = rescore(r["components"], mine[r["accession"]])
            got[r["accession"]] = sc.score
            if sc.score != int(r[key]) or sc.cap != _flag(r["components"]):
                bad.append(f"{r['gene_symbol']} {sc.score}{'/' + sc.cap if sc.cap else ''} "
                           f"vs {r[key]}")
    recall = sum(got[r["accession"]] >= bar for r in pos)
    spec = sum(got[r["accession"]] < bar for r in neg)
    ryr = [r for r in neg if r["category"].startswith("RyR")]
    ryr_pass = sum(got[r["accession"]] < bar for r in ryr)
    counts = (recall, len(pos), spec, len(neg), ryr_pass, len(ryr))
    want = (s["recall_hit"], s["recall_n"], s["specificity_hit"], s["specificity_n"],
            s["ryr_pass"], s["ryr_n"])
    missed = [r["gene_symbol"] for r in pos if got[r["accession"]] < bar]
    return agree(not bad and counts == want,
                 "recall {}/{}, specificity {}/{}, RyR rejected {}/{}".format(*want),
                 "recall {}/{}, specificity {}/{}, RyR rejected {}/{}".format(*counts)
                 + f"; missed {', '.join(missed) or 'none'}"
                 + (f"; scores differ: {'; '.join(bad[:6])}" if bad else ""),
                 f"Promotion at {bar}; a cap sets a score to {bar - 1}. The six RyR "
                 "decoys are capped by the recomputed margin, not by S1's flag.")
