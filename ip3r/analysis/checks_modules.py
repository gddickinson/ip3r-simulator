"""Paper 6's module contrast, re-derived from the deep alignments.

The paper compares the IP3-binding core with the pore by pairing the two
inside each vertebrate orthologue. These checks rebuild both modules from
this project's annotation (:mod:`ip3r.core.modules`), map them onto S17's
deep alignments by walking the reference row, recompute every tip's two
identities and run the paired tests with this project's statistics
(:mod:`ip3r.analysis.module_contrast`). The published ``module_contrast.tsv``
is only read to compare against.
"""

from __future__ import annotations

import math

from ..config import PARALOGS
from ..core import genes_data as G
from ..core.modules import module
from ..parameters import PARAMETERS as _P
from .checks import agree, register
from .module_contrast import alignment_path, paired_contrast

MAP = "ligand_site/module_map.tsv"
CONTRAST = "ligand_site/module_contrast.tsv"
ALIGNMENTS = tuple(alignment_path(p) for p in PARALOGS)
CORE = "contact_span"
ALPHA = 0.05


def _published(pore_def: str) -> dict[str, dict]:
    return {r["paralog"]: r for r in G.read_tsv(CONTRAST)
            if r["core_definition"] == CORE and r["pore_definition"] == pore_def}


def _compare(pore_def: str) -> tuple[bool, list[str], dict]:
    """Re-derive one pore definition's three rows; compare every field."""
    tol, ptol = _P.value("check.stat_tol"), _P.value("check.log_p_tol")
    pub = _published(pore_def)
    ok, lines, data = True, [], {}
    for gene in PARALOGS:
        pc = paired_contrast(gene, CORE, pore_def)
        s, p = pc.stats(), pub[gene]
        counts = ("n_tips", "n_dropped", "n_core_more_conserved",
                  "n_pore_more_conserved", "n_ties")
        bad = [k for k in counts if s[k] != int(p[k])]
        bad += [k for k in ("mean_core_identity", "mean_pore_identity",
                            "mean_difference")
                if abs(s[k] - float(p[k])) > tol]
        dp = abs(math.log10(s["p_wilcoxon"]) - math.log10(float(p["p_wilcoxon"])))
        if dp > ptol:
            bad.append("p_wilcoxon")
        if s["direction"] != p["direction"]:
            bad.append("direction")
        ok &= not bad
        lines.append(f"{gene}: {s['mean_difference']:+.4f} on {s['n_tips']} tips "
                     f"({s['n_core_more_conserved']} core / "
                     f"{s['n_pore_more_conserved']} pore), p = {s['p_wilcoxon']:.3g}"
                     + (f" — differs in {', '.join(bad)}" if bad else ""))
        data[gene] = {"core": pc.core.tolist(), "pore": pc.pore.tolist(),
                      "p": s["p_wilcoxon"], "mean_difference": s["mean_difference"]}
    return ok, lines, data


@register("P6.module_map", "ligand",
          "The ligand core is the span of the ten IP3 contacts and the pore "
          "module is PF00520 less the luminal loop, with the spans in "
          "module_map.tsv.",
          "Both modules rebuilt from this project's imported sites and "
          "domain map (core.modules), each validated to hold all ten "
          "contacts / both filter and gate residues and nothing of the "
          "other; spans compared with module_map.tsv.",
          "rederived", (MAP,))
def module_map():
    pub = G.read_tsv(MAP)
    ok, lines = True, []
    for r in pub:
        if r["definition"] == "ibc_literature":      # a literature boundary,
            continue                                 # not a rule to rebuild
        m = module(r["paralog"], r["definition"])
        excl = ";".join(f"{a}-{b}" for a, b in m.excluded)
        mine = (m.start, m.end, excl, len(m.residues))
        theirs = (int(r["start"]), int(r["end"]), r["excluded"], int(r["n_residues"]))
        ok &= mine == theirs
        if mine != theirs:
            lines.append(f"{r['paralog']} {r['definition']}: {mine} vs {theirs}")
    core, pore = module("ITPR3", CORE), module("ITPR3", "channel_minus_luminal")
    found = "; ".join(lines) or (
        f"all nine spans identical (ITPR3: core {core.start}–{core.end}, "
        f"pore {pore.start}–{pore.end} less {pore.excluded[0][0]}–"
        f"{pore.excluded[0][1]})")
    return agree(ok, "spans as in module_map.tsv", found)


@register("P6.module_contrast", "ligand",
          "Paired per orthologue, the pore is more conserved than the ligand "
          "core in ITPR1 and ITPR3 (223 vs 32, 218 vs 42 tips) and level in "
          "ITPR2 (q = 0.13).",
          "Every tip's identity to the human reference in each module "
          "recomputed from S17's deep alignments (own column map, own "
          "modules, 50 % coverage floor); sign counts, means and a "
          "tie-corrected signed-rank test (own code) compared with "
          "module_contrast.tsv.",
          "rederived", (CONTRAST,) + ALIGNMENTS)
def module_contrast():
    ok, lines, data = _compare("channel_minus_luminal")
    sig = {g: d["p"] < ALPHA for g, d in data.items()}
    claim = (sig["ITPR1"] and sig["ITPR3"] and not sig["ITPR2"]
             and data["ITPR1"]["mean_difference"] < 0
             and data["ITPR3"]["mean_difference"] < 0)
    if not claim:
        lines.append("the pattern stated (pore ahead in ITPR1 and ITPR3, "
                     "ITPR2 level) does not hold")
    return agree(ok and claim, "pore > core in ITPR1, ITPR3; ITPR2 level",
                 "; ".join(lines), paralogs=data)


@register("P6.loop_reverses", "ligand",
          "Counting the luminal loop as pore reverses the answer in all three "
          "paralogues: the core leads by about five identity points.",
          "The same paired re-derivation with the pore module taken as the "
          "whole PF00520 span, compared with module_contrast.tsv.",
          "rederived", (CONTRAST,) + ALIGNMENTS)
def loop_reverses():
    ok, lines, data = _compare("channel_all")
    lead = [d["mean_difference"] for d in data.values()]
    claim = all(x > 0 and d["p"] < ALPHA for x, d in zip(lead, data.values()))
    if not claim:
        lines.append("the core does not lead in all three")
    return agree(ok and claim, "core > pore in all three with the loop in",
                 "; ".join(lines), paralogs=data)
