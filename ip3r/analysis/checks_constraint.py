"""Constraint and ligand-site findings, re-derived from the per-residue tables.

Paper 5's inputs are the per-residue conservation tables
(``constraint_<gene>_<acc>.tsv``) and the harvested variants; its outputs are
the per-element means, the classifier AUCs and the rankings stated in prose.
These checks recompute the outputs from the inputs with this project's code,
using **this project's own** element assignment (``core.annotations``, from
the domain map) rather than the ``element`` column S17 wrote — so a
disagreement in how residues were assigned would surface here.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from ..config import PARALOG_ACC
from ..core import genes_data as G
from ..core.annotations import residue_elements
from ..parameters import PARAMETERS as _P
from .checks import agree, register
from .stats import auc

BY_ELEMENT = "constraint/constraint_by_element.tsv"
VARIANTS = "constraint/variants.tsv"
VARTEST = "constraint/variant_constraint_test.tsv"
CONTACTS = "ligand_site/contact_test.tsv"
MODULES = "ligand_site/module_map.tsv"
OMEGA = "selection/omega_table.tsv"
LAYERS = ("deep", "shallow", "vert", "family")


def _table(gene: str) -> str:
    return f"constraint/constraint_{gene}_{PARALOG_ACC[gene]}.tsv"


@lru_cache(maxsize=None)
def per_residue(gene: str) -> list[dict]:
    return G.read_tsv(_table(gene))


def element_means(gene: str, metric: str = "deep_jsd") -> dict[str, tuple[float, int]]:
    """Mean of a reliable deep-layer metric per element (this project's
    element assignment). ``metric`` is ``deep_jsd`` or ``deep_frac_modal``."""
    els = residue_elements(gene)
    acc: dict[str, list[float]] = {}
    for row in per_residue(gene):
        if not G.as_bool(row["deep_reliable"]):
            continue
        v = G.as_float(row[metric])
        if np.isfinite(v):
            acc.setdefault(els[int(row["resi"]) - 1], []).append(v)
    return {k: (float(np.mean(v)), len(v)) for k, v in acc.items()}


def _published_means() -> dict[tuple[str, str], dict]:
    return {(r["paralog"], r["element"]): r for r in G.read_tsv(BY_ELEMENT)}


@register("P5.element_means", "constraint",
          "Per-element mean deep-layer conservation (JSD) of each human "
          "paralog, as committed in constraint_by_element.tsv.",
          "Means recomputed from the per-residue tables over reliable "
          "residues, with residues assigned to elements by this project's "
          "rule from the domain map (not S17's element column).",
          "rederived", (BY_ELEMENT, "constraint/constraint_*.tsv"))
def element_means_check():
    pub = _published_means()
    tol = _P.value("check.stat_tol")
    worst, n, lines = 0.0, 0, []
    for gene in PARALOG_ACC:
        for el, (m, k) in element_means(gene).items():
            row = pub.get((gene, el))
            if row is None or el == "linker":
                continue
            d = abs(m - float(row["mean_jsd"]))
            worst = max(worst, d)
            n += 1
            if d > tol or k != int(row["n_sites"]):
                lines.append(f"{gene} {el}: {m:.4f} (n={k}) vs "
                             f"{row['mean_jsd']} (n={row['n_sites']})")
    return agree(not lines, f"{n} paralog × element means",
                 f"largest |Δ| {worst:.4f}" + ("; " + "; ".join(lines) if lines else
                                             ", every site count identical"),
                 "Linkers are excluded: S17 names each linker by its flanking "
                 "domains and this project pools them.")


def _ranking(gene: str, metric: str = "deep_jsd") -> list[tuple[str, float]]:
    pub = _published_means()
    means = element_means(gene, metric)
    ranked = [(el, m) for el, (m, _) in means.items()
              if el != "linker" and (gene, el) in pub
              and not G.as_bool(pub[(gene, el)]["is_control"])]
    return sorted(ranked, key=lambda x: -x[1])


def _rank_of(ranked, element: str, tie: float = 5e-5) -> tuple[int, int]:
    """(best, worst) 1-based rank of ``element``, ties sharing a range."""
    v = dict(ranked)[element]
    above = sum(m > v + tie for _, m in ranked)
    level = sum(abs(m - v) <= tie for _, m in ranked)
    return above + 1, above + level


@register("P5.gate_filter_most_conserved", "constraint",
          "On the modal-residue fraction the gate and filter are the two most "
          "conserved elements in all three paralogs (ITPR1's filter tied "
          "with RIH-associated); on the divergence metric the gate ranks "
          "first in ITPR2 and ITPR3 and second in ITPR1, and the filter "
          "second to fourth (paper 5, Results).",
          "Both metrics' element means recomputed from the per-residue "
          "tables and ranked; ties within 5e-5 share a rank range.",
          "rederived", (BY_ELEMENT, "constraint/constraint_*.tsv"))
def gate_filter_top():
    lines, ok = [], True
    want_gate_jsd = {"ITPR1": 2, "ITPR2": 1, "ITPR3": 1}
    for gene in PARALOG_ACC:
        fm = _ranking(gene, "deep_frac_modal")
        jsd = _ranking(gene, "deep_jsd")
        g_fm, f_fm = _rank_of(fm, "gate"), _rank_of(fm, "selectivity_filter")
        g_j, f_j = _rank_of(jsd, "gate"), _rank_of(jsd, "selectivity_filter")
        ok &= g_fm[0] <= 2 and f_fm[0] <= 2
        ok &= g_j[0] == want_gate_jsd[gene] and 2 <= f_j[0] <= 4
        lines.append(f"{gene}: modal-fraction ranks gate {g_fm[0]}, filter "
                     f"{f_fm[0]}{'-' + str(f_fm[1]) if f_fm[1] != f_fm[0] else ''}; "
                     f"JSD ranks gate {g_j[0]}, filter {f_j[0]}")
    return agree(ok, "as stated", "; ".join(lines))


@register("P5.report_both_metrics", "constraint",
          "S17 report §5: \"The gate and the selectivity filter are the most "
          "constrained elements of the protein, on both metrics and in all "
          "three paralogs.\"",
          "The same rankings, read for the stronger statement: gate and "
          "filter the top two on JSD as well as on the modal fraction, in "
          "every paralog.",
          "rederived", (BY_ELEMENT, "constraint/constraint_*.tsv"))
def report_both_metrics():
    lines, ok = [], True
    for gene in PARALOG_ACC:
        for metric in ("deep_frac_modal", "deep_jsd"):
            top = {el for el, _ in _ranking(gene, metric)[:2]}
            hit = top == {"gate", "selectivity_filter"}
            if metric == "deep_frac_modal" and not hit:
                r = _rank_of(_ranking(gene, metric), "selectivity_filter")
                hit = r[0] <= 2 and _rank_of(_ranking(gene, metric), "gate")[0] <= 2
            ok &= hit
            if not hit:
                lines.append(f"{gene} {metric.replace('deep_', '')}: top two "
                             f"{', '.join(el for el, _ in _ranking(gene, metric)[:2])}")
    return agree(ok, "top two on both metrics everywhere",
                 "; ".join(lines) or "holds",
                 "The paper's Results state the narrower version, which "
                 "P5.gate_filter_most_conserved checks.")


@register("P5.luminal_loop_least", "constraint",
          "The luminal loop inside the pore domain is the least conserved "
          "element of the receptor.",
          "Elements ranked by the recomputed mean deep JSD.",
          "rederived", (BY_ELEMENT, "constraint/constraint_*.tsv"))
def luminal_last():
    lines, ok = [], True
    for gene in PARALOG_ACC:
        ranked = _ranking(gene)
        ok &= ranked[-1][0] == "luminal_loop"
        lines.append(f"{gene}: last {ranked[-1][0]} ({ranked[-1][1]:.3f}), "
                     f"next {ranked[-2][0]} ({ranked[-2][1]:.3f})")
    return agree(ok, "luminal loop last in all three", "; ".join(lines))


@register("P5.gate_identical", "constraint",
          "The gate is identical between the three human paralogs.",
          "Gate residues of each paralog joined to the others through the "
          "shared msa_v2 alignment column in the per-residue tables.",
          "rederived", ("constraint/constraint_*.tsv",))
def gate_identical():
    by_col = {}
    for gene in PARALOG_ACC:
        els = residue_elements(gene)
        by_col[gene] = {int(r["msa_col"]): (int(r["resi"]), r["aa"])
                        for r in per_residue(gene)
                        if r["msa_col"] not in ("", "-1")
                        and els[int(r["resi"]) - 1] == "gate"}
    cols = sorted(set.intersection(*(set(v) for v in by_col.values())))
    strings = {g: "".join(by_col[g][c][1] for c in cols) for g in PARALOG_ACC}
    ok = bool(cols) and len(set(strings.values())) == 1
    return agree(ok, "identical", ", ".join(f"{g} {s}" for g, s in strings.items())
                 + f" over {len(cols)} shared columns")


def _score(site: dict | None, layer: str) -> float | None:
    """A layer's score at a residue, or None below the occupancy floor."""
    if not site:
        return None
    v = G.as_float(site.get(f"{layer}_jsd", ""))
    occ = G.as_float(site.get(f"{layer}_occupancy", ""))
    if not np.isfinite(v) or (np.isfinite(occ)
                              and occ < _P.value("constraint.min_occupancy")):
        return None
    return v


def _sites(gene: str) -> dict[int, dict]:
    return {int(r["resi"]): r for r in per_residue(gene)}


def position_scores(rows, gene: str | None, layer: str, bucket: str,
                    require_all: bool = False) -> list[float]:
    """One score per labelled **position** (not per variant), ClinVar only.

    Positions, because two alleles at one residue are one observation of
    that residue's constraint; ClinVar only, because the curated UniProt
    variants are the harvest's positive control, not part of the test set.
    """
    seen, out = set(), []
    sites = {g: _sites(g) for g in PARALOG_ACC}
    for r in rows:
        if r["source"] != "clinvar" or r["class_bucket"] != bucket:
            continue
        if gene and r["gene"] != gene:
            continue
        key = (r["gene"], int(r["resi"]))
        if key in seen:
            continue
        site = sites[r["gene"]].get(int(r["resi"]))
        if require_all and any(_score(site, la) is None for la in LAYERS):
            continue
        x = _score(site, layer)
        if x is not None:
            seen.add(key)
            out.append(x)
    return out


def _variant_auc(rows, gene, layer, require_all=False):
    pos = position_scores(rows, gene, layer, "P/LP", require_all)
    neg = position_scores(rows, gene, layer, "B/LB", require_all)
    return auc(pos, neg), len(pos), len(neg)


@register("P5.variant_auc", "constraint",
          "Conservation separates pathogenic from benign ClinVar positions: "
          "the AUCs of variant_constraint_test.tsv, all three contrasts "
          "that compare labelled sets.",
          "AUC recomputed per position from the per-residue tables with a "
          "Mann-Whitney implementation of this project's own (ties count "
          "one half), under S17's stated rules: ClinVar only, one score per "
          "position, occupancy at least 0.5.",
          "rederived", (VARIANTS, VARTEST, "constraint/constraint_*.tsv"))
def variant_auc():
    rows = G.read_tsv(VARIANTS)
    tol = _P.value("check.stat_tol")
    bad, n = [], 0
    for t in G.read_tsv(VARTEST):
        if t["contrast"] == "P/LP vs whole protein" or not t["auc"]:
            continue
        gene = None if t["gene"] == "POOLED" else t["gene"]
        a, npos, nneg = _variant_auc(rows, gene, t["layer"],
                                     require_all="all layers" in t["contrast"])
        n += 1
        if abs(a - float(t["auc"])) > tol or (npos, nneg) != (
                int(t["n_positive"]), int(t["n_negative"])):
            bad.append(f"{t['gene']} {t['layer']} [{t['contrast']}]: {a:.4f} "
                       f"({npos}/{nneg}) vs {t['auc']} ({t['n_positive']}/"
                       f"{t['n_negative']})")
    return agree(not bad, f"{n} AUCs", "all agree within tolerance, with "
                 "identical position counts" if not bad else "; ".join(bad))


@register("P5.deep_ranks_third", "constraint",
          "On the positions every layer scores, the family-wide layer is the "
          "best classifier and the deep within-paralog layer ranks third of "
          "four.",
          "Pooled AUC per layer recomputed on the all-layers-scorable "
          "positions, then ranked.",
          "rederived", (VARIANTS, VARTEST, "constraint/constraint_*.tsv"))
def deep_third():
    rows = G.read_tsv(VARIANTS)
    aucs = {la: _variant_auc(rows, None, la, require_all=True)[0] for la in LAYERS}
    order = sorted(aucs, key=lambda la: -aucs[la])
    ok = order[0] == "family" and order.index("deep") == 2
    return agree(ok, "family first, deep third",
                 ", ".join(f"{la} {aucs[la]:.4f}" for la in order), aucs=aucs)


@register("P5.vus_count", "constraint",
          "The per-residue tables stratify 1,546 variants of uncertain "
          "significance.",
          "Variants whose class bucket is VUS counted in variants.tsv.",
          "rederived", (VARIANTS,))
def vus_count():
    n = sum(r["class_bucket"] == "VUS" for r in G.read_tsv(VARIANTS))
    return agree(n == 1546, "1,546", f"{n:,}")


@register("P5.omega_range", "constraint",
          "All three paralogs are held far below neutrality: one-ratio ω "
          "0.02 to 0.04.",
          "The M0 ω of each paralog read from omega_table.tsv and rounded to "
          "two decimals.",
          "read", (OMEGA,))
def omega_range():
    om = {r["set"]: float(r["omega"]) for r in G.read_tsv(OMEGA)
          if r["job"] in {f"m0_{g}" for g in PARALOG_ACC}}
    ok = len(om) == 3 and all(0.02 <= round(v, 2) <= 0.04 for v in om.values())
    return agree(ok, "0.02-0.04", ", ".join(f"{g} {v:.4f}" for g, v in sorted(om.items())))


@register("P6.contacts_vs_core", "ligand",
          "The ten IP3 contacts are more conserved than the rest of the "
          "ligand core (deep layer), by the differences in contact_test.tsv.",
          "Mean deep JSD of the contact residues and of the rest of the "
          "primary ligand-core span recomputed from the per-residue tables.",
          "rederived", (CONTACTS, MODULES, "constraint/constraint_*.tsv"))
def contacts_vs_core():
    spans = {r["paralog"]: (int(r["start"]), int(r["end"]))
             for r in G.read_tsv(MODULES)
             if r["module"] == "ligand_core" and G.as_bool(r["is_primary"])}
    pub = {r["paralog"]: r for r in G.read_tsv(CONTACTS)
           if r["contact_set"] == "s0_contact" and r["background"] == "rest_of_core"
           and r["layer"] == "deep"}
    tol = _P.value("check.stat_tol")
    lines, ok = [], True
    for gene, (lo, hi) in spans.items():
        c, b = [], []
        for r in per_residue(gene):
            resi = int(r["resi"])
            v = G.as_float(r["deep_jsd"])
            if not (lo <= resi <= hi) or not np.isfinite(v):
                continue
            (c if G.as_bool(r["ip3_contact"]) else b).append(v)
        diff = float(np.mean(c) - np.mean(b))
        p = float(pub[gene]["difference"])
        ok &= abs(diff - p) <= tol and diff > 0
        lines.append(f"{gene}: {diff:+.4f} (published {p:+.4f})")
    return agree(ok, "contacts more conserved in all three", "; ".join(lines))
