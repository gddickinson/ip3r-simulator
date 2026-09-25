"""Paper 3's lesion lead re-derived from the per-locus integrity table.

S15b §8.2 says ITPR3's identity-matched indel excess is a bird result, and
§8.3 that most of the bird pairs sit below the contiguity bar. Both are
rebuilt here by :mod:`.lesion_strata` from ``integrity_loci.tsv``: the
pairs, the per-class sign tests and their BH q-values, and the split at this
project's registered bar (not the table's ``contig_spans_gene``). The
outcome carries the strata for the exhibit; the Genomes tab's lesion layer
draws the same pairs genome by genome.
"""

from __future__ import annotations

import math

from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from . import genome_grid as GG
from . import lesion_strata as L
from .checks import agree, register


def _close(a: float, printed: str) -> bool:
    """A p or q against its printed value (six significant figures)."""
    if printed == "":
        return not math.isfinite(a)
    b = float(printed)
    return math.isfinite(a) and abs(a - b) <= 1e-5 * max(abs(b), 1e-300)


def _half(t: tuple, pos: str, neg: str, p: str) -> bool:
    return (t[0], t[1]) == (int(pos), int(neg)) and _close(t[2], p)


@register("P3.lesion_strata", "retention",
          "Stratified by class, ITPR3's identity-matched indel excess is 25 "
          "genomes to 2 in Aves (q = 4.52e-05) and 7 to 6 in Actinopteri "
          "(p = 1); 21 of the 27 bird pairs are in assemblies below the "
          "contiguity bar (19:2, p = 0.0002), and the 6 above all point the "
          "same way. The bird ITPR1 and ITPR2 deficits are the same "
          "comparison read from the other side.",
          "Pairs rebuilt from integrity_loci.tsv: each genome's best-covered "
          "locus per cell at lesion.coverage_bar, lesions per kilo-residue "
          "from frameshifts and stops, against the median of its siblings "
          "within lesion.identity_window. Exact sign test per cell × class, "
          "BH over strata with ≥ lesion.min_untied untied pairs, each "
          "significant stratum split at genomes.contiguity_bar_bp; compared "
          "with integrity_pairs.tsv, lesion_by_class.tsv and "
          "lesion_class_controls.tsv.",
          "rederived", (L.LOCI, L.PAIRS, L.BY_CLASS, L.CONTROLS, GG.CONTIG))
def lesion_strata():
    loci, contig = G.read_tsv(L.LOCI), G.read_tsv(GG.CONTIG)
    vclass = {r["accession"]: r["vclass"] for r in contig}
    above = {r["accession"]: GG.above_bar(float(r["contig_n50"])) for r in contig}
    prs = L.pairs(loci)
    ours = {(p.accession, p.cell): p.diff for p in prs}
    pub = {(r["accession"], r["cell"]): float(r["diff"])
           for r in G.read_tsv(L.PAIRS) if r["matched"] == "1"}
    pair_bad = (set(ours) ^ set(pub)) | {k for k in ours.keys() & pub.keys()
                                         if abs(ours[k] - pub[k]) > 1e-6}
    strata = L.stratify(prs, vclass)
    table = {(r["cell"], r["vclass"]): r for r in G.read_tsv(L.BY_CLASS)}
    strat_bad = set(table) ^ {(s.cell, s.vclass) for s in strata}
    for s in strata:
        r = table.get((s.cell, s.vclass))
        if r is None:
            continue
        counts = tuple(int(r[k]) for k in ("n_genomes", "n_pos", "n_neg", "n_ties"))
        if (counts != (s.n_genomes, s.n_pos, s.n_neg, s.n_ties)
                or r["direction"] != s.direction or not _close(s.p, r["p"])
                or not _close(s.q, r["q_bh"])):
            strat_bad.add((s.cell, s.vclass))
    splits = L.controls(strata, prs, vclass, above)
    ctab = {(r["cell"], r["vclass"]): r for r in G.read_tsv(L.CONTROLS)}
    split_bad = set(ctab) ^ {(c.stratum.cell, c.stratum.vclass) for c in splits}
    for c in splits:
        r = ctab.get((c.stratum.cell, c.stratum.vclass))
        if r is None:
            continue
        if not (_half(c.above, r["pos_above"], r["neg_above"], r["p_above"])
                and _half(c.below, r["pos_below"], r["neg_below"], r["p_below"])
                and ";".join(c.siblings) == r["sibling_of"]
                and round(c.median_identity, 4) == float(r["median_identity"])
                and round(c.median_others_identity, 4)
                == float(r["median_others_identity"])):
            split_bad.add((c.stratum.cell, c.stratum.vclass))
    by = {(s.cell, s.vclass): s for s in strata}
    bird, fish = by[("ITPR3", "Aves")], by[("ITPR3", "Actinopteri")]
    bsplit = next(c for c in splits if (c.stratum.cell, c.stratum.vclass)
                  == ("ITPR3", "Aves"))
    alpha = _P.value("check.alpha")
    pattern = (bird.direction == "excess" and bird.q <= alpha
               and not (fish.testable and fish.p <= alpha)
               and bsplit.below[2] <= alpha and bsplit.above[1] == 0
               and set(bsplit.siblings) == {"ITPR1", "ITPR2"})
    ok = pattern and not (pair_bad or strat_bad or split_bad)
    t = table[("ITPR3", "Aves")]
    ct = ctab.get(("ITPR3", "Aves"), {})
    published = (f"ITPR3 Aves {t['n_pos']}:{t['n_neg']} (q {float(t['q_bh']):.3g}); "
                 f"below the bar {ct.get('pos_below')}:{ct.get('neg_below')} "
                 f"(p {float(ct.get('p_below') or 'nan'):.2g}), above "
                 f"{ct.get('pos_above')}:{ct.get('neg_above')}; {len(pub)} pairs, "
                 f"{len(table)} strata")
    found = (f"ITPR3 Aves {bird.n_pos}:{bird.n_neg} (q {bird.q:.3g}), Actinopteri "
             f"{fish.n_pos}:{fish.n_neg} (p {fish.p:.2g}); below the bar "
             f"{bsplit.below[0]}:{bsplit.below[1]} (p {bsplit.below[2]:.2g}), above "
             f"{bsplit.above[0]}:{bsplit.above[1]}; {len(ours)} pairs, "
             f"{len(strata)} strata")
    detail = (f"{len(pair_bad)} pairs, {len(strat_bad)} strata and {len(split_bad)} "
              f"bar splits differ from the tables; siblings of the bird ITPR3 "
              f"stratum: {', '.join(bsplit.siblings) or 'none'}")
    return agree(ok, published, found, detail, strata=[
        {"cell": s.cell, "vclass": s.vclass, "pos": s.n_pos, "neg": s.n_neg,
         "q": s.q} for s in strata if s.testable],
        splits=[{"cell": c.stratum.cell, "vclass": c.stratum.vclass,
                 "above": c.above[:2], "below": c.below[:2]} for c in splits])
