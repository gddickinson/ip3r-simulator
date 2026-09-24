"""Every findings check must be able to change its mind.

For each registered check: copy the ip3r_genes tables it **declares** as
sources into a temporary project, run it (it must give the same verdict as on
the real project — which proves the source list is complete), then plant a
change in one table and require the verdict to flip. A check whose verdict no
planted input can flip asserts nothing (the PIEZO1 project's rule: calibrate a
checking instrument before believing it).

A new check without an entry in ``PLANTS`` fails
``test_every_check_has_a_calibration``.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from ip3r.analysis import checks as C
from ip3r.analysis import checks_constraint, module_contrast, shell_constraint
from ip3r.config import GENES_DIR
from ip3r.parameters import PARAMETERS
from conftest import needs_genes, needs_structure


# ------------------------------------------------------------ table editing

def _edit(path: Path, fn) -> None:
    """Rewrite a TSV row by row: ``fn(row) -> row``; rows are dicts."""
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        cols, rows = r.fieldnames, list(r)
    rows = [fn(dict(x)) for x in rows]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows([x for x in rows if x is not None])


def _set(where: dict, **new):
    def fn(row):
        if all(row.get(k) == v for k, v in where.items()):
            row.update({k: str(v) for k, v in new.items()})
        return row
    return fn


def _first(where: dict, **new):
    done = []

    def fn(row):
        if not done and all(row.get(k) == v for k, v in where.items()):
            row.update({k: str(v) for k, v in new.items()})
            done.append(1)
        return row
    return fn


def _json(path: Path, fn) -> None:
    d = json.loads(path.read_text())
    fn(d)
    path.write_text(json.dumps(d))


R = "results/"
CON3 = R + "constraint/constraint_ITPR3_Q14573.tsv"
CON1 = R + "constraint/constraint_ITPR1_Q14643.tsv"
META = R + "s0_baseline/review_figures/structure_meta.json"


def _swap_tips(d: Path, a: str, b: str) -> None:
    """Exchange two tip labels in the committed tree (an input plant)."""
    p = d / R / "phylogeny/rooted.nwk"
    t = p.read_text()
    assert a in t and b in t
    p.write_text(t.replace(a, "\0").replace(b, a).replace("\0", b))


def _gate_filter_boost(d: Path):
    for g, a in (("ITPR1", "Q14643"), ("ITPR2", "Q14571"), ("ITPR3", "Q14573")):
        _edit(d / f"results/constraint/constraint_{g}_{a}.tsv",
              lambda r: r if r["element"] not in ("gate", "selectivity_filter")
              else {**r, "deep_jsd": "0.99", "deep_frac_modal": "1.0"})


def _pore_to_reference(d: Path, gene: str) -> None:
    """Rewrite one deep alignment so every tip matches the reference over
    the whole channel domain (loop included)."""
    from ip3r.core.modules import module
    p = d / R / f"constraint/aln_{gene}.fasta"
    recs, lab = {}, None
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            lab = line[1:]
            recs[lab] = []
        else:
            recs[lab].append(line.strip())
    recs = {k: "".join(v) for k, v in recs.items()}
    ref = next(v for k, v in recs.items() if k.startswith("REF|"))
    cols = [i for i, c in enumerate(ref) if c != "-"]
    m = module(gene, "channel_all")
    span = set(cols[r - 1] for r in m.residues)
    out = []
    for k, v in recs.items():
        row = "".join(ref[i] if i in span else c for i, c in enumerate(v))
        out += [f">{k}", row]
    p.write_text("\n".join(out) + "\n")


def _contact_shell_conserved(d: Path) -> None:
    """Every contact-shell residue of ITPR3 made invariant: a step at 4.5 Å."""
    contact = {str(r.resi) for r in shell_constraint.measured_shells()
               if r.shell == "contact"}
    _edit(d / CON3, lambda r: {**r, "deep_jsd": "0.99"} if r["resi"] in contact else r)


def _drop_below_bar(d: Path) -> None:
    """The least contiguous genome above the bar re-assembled just under it:
    one genome and its cells move across the split."""
    p = d / R / "methods/contiguity_cells.tsv"
    with open(p, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    bar = PARAMETERS.value("genomes.contiguity_bar_bp")
    acc = min((r for r in rows if float(r["contig_n50"]) >= bar),
              key=lambda r: float(r["contig_n50"]))["accession"]
    _edit(p, _set({"accession": acc}, contig_n50=int(bar) - 1))


def _rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _call_in(d: Path, rank: str, name: str) -> None:
    """One proteome of a named phylum or class given an IP3R call."""
    taxa = {r["taxid"] for r in _rows(d / R / "s20_sweep/taxonomy.tsv") if r[rank] == name}
    assert taxa
    done = []

    def fn(row):
        if not done and row["taxid"] in taxa:
            row["n_itpr"] = "1"
            done.append(1)
        return row
    _edit(d / R / "s20_sweep/proteome_presence.tsv", fn)


def _gene_in(d: Path, phylum: str) -> None:
    """A complete gene model moved onto the first genome of a phylum."""
    acc = next(m["accession"] for m in _rows(d / R / "s23_scope/genome_manifest_s23.tsv")
               if m["phylum"] == phylum)
    _edit(d / R / "s23_scope/copies.tsv", _first({}, accession=acc))


def _drop_first(path: Path) -> None:
    done = []

    def fn(row):
        if not done:
            done.append(1)
            return None
        return row
    _edit(path, fn)


PLANTS = {
    "P5.element_means": lambda d: _edit(d / R / "constraint/constraint_by_element.tsv",
                                        _set({"paralog": "ITPR1", "element": "MIR"}, mean_jsd=0.9)),
    "P5.gate_filter_most_conserved": lambda d: _edit(
        d / CON3, lambda r: {**r, "deep_jsd": "0.1", "deep_frac_modal": "0.1"}
        if r["element"] == "gate" else r),
    "P5.report_both_metrics": _gate_filter_boost,
    "P5.luminal_loop_least": lambda d: _edit(
        d / CON1, lambda r: {**r, "deep_jsd": "0.99"} if r["element"] == "luminal_loop" else r),
    "P5.gate_identical": lambda d: _edit(
        d / CON3, _first({"element": "gate"}, aa="W")),
    "P5.variant_auc": lambda d: _edit(d / R / "constraint/variant_constraint_test.tsv",
                                      _first({"gene": "ITPR1", "layer": "deep"}, auc=0.1)),
    "P5.deep_ranks_third": lambda d: _edit(
        d / R / "constraint/variants.tsv",
        lambda r: {**r, "class_bucket": {"P/LP": "B/LB", "B/LB": "P/LP"}.get(
            r["class_bucket"], r["class_bucket"])}),
    "P5.vus_count": lambda d: _edit(d / R / "constraint/variants.tsv",
                                    _first({"class_bucket": "VUS"}, class_bucket="other")),
    # An input, not the published table: proves the rows are rebuilt.
    "P5.vus_stratification": lambda d: _edit(
        d / R / "constraint/variants.tsv",
        _first({"gene": "ITPR3", "class_bucket": "P/LP", "source": "clinvar"},
               class_bucket="B/LB")),
    "P5.omega_range": lambda d: _edit(d / R / "selection/omega_table.tsv",
                                      _set({"job": "m0_ITPR1"}, omega=0.2)),
    "P6.contacts_vs_core": lambda d: _edit(
        d / R / "ligand_site/contact_test.tsv",
        _set({"contact_set": "s0_contact", "background": "rest_of_core",
              "paralog": "ITPR2", "layer": "deep"}, difference=0.5)),
    "P2.sister_pair": lambda d: _edit(d / R / "phylogeny/membership_audit.tsv",
                                      _first({"assigned_to": "ITPR1", "rule": "core"},
                                             assigned_to="ITPR3")),
    "P2.au_test": lambda d: _edit(d / R / "phylogeny/au_test.tsv",
                                  _set({"tree": "H1_12"}, p_AU=0.3)),
    # Tree plants are made in the tree itself: move a tip, or raise one support.
    "P2.paralog_clades": lambda d: _swap_tips(d, "ITPR1_Bos_taurus_Bovine_Q9TU34",
                                              "ITPR2_Homo_sapiens_Human_Q14571"),
    "P2.cyclostome_lineages": lambda d: _swap_tips(
        d, "vertebrate_basal_Petromyzon_marinus_Sea_lamprey_A0ACM8C056",
        "ITPR3_Homo_sapiens_Human_Q14573"),
    "P2.support_bar": lambda d: (d / R / "phylogeny/rooted.nwk").write_text(
        (d / R / "phylogeny/rooted.nwk").read_text().replace("83.1/77", "83.1/97")),
    "P2.teleost_itpr1": lambda d: _edit(d / R / "duplication/teleost_copies.tsv",
                                        _set({"group": "teleost"}, ITPR1_copies=1)),
    "P3.no_absent_cells": lambda d: _edit(d / R / "loss_dynamics/character_matrix.tsv",
                                          _first({}, state="absent")),
    "P3.dollo_zero": lambda d: _edit(d / R / "loss_counts/dollo_counts.tsv",
                                     _set({"coding": "family", "is_base": "1"},
                                          dollo_losses_max=1)),
    "P3.false_negatives": lambda d: _edit(d / R / "methods/contiguity_cells.tsv",
                                          _first({"control": "itpr_present",
                                                  "false_negative": "0"}, false_negative=1)),
    # Papers 3/4 grid checks: input plants in the per-cell tables.
    "P3.miss_by_contiguity": _drop_below_bar,
    "P3.contiguity_tests": lambda d: _edit(d / R / "methods/contiguity_cells.tsv",
                                           _first({"control": "ryr_sister",
                                                   "false_negative": "0"}, false_negative=1)),
    "P4.recovery_channels": lambda d: _edit(
        d / R / "methods/gene_recovery.tsv",
        _first({"recovery_channel": "protein_database_and_genome"},
               n_records_resolving_to_cell=0)),
    "P4.unreachable": lambda d: _edit(d / R / "methods/gene_recovery.tsv",
                                      _first({"cell": "ITPR1", "gene_present": "1"},
                                             recovery_channel="protein_db")),
    # Paper 1: input plants in the per-proteome, per-record and per-genome rows.
    "P1.presence_range": lambda d: _edit(d / R / "s20_sweep/proteome_presence.tsv",
                                         _first({"group": "bacteria_genus"}, n_itpr=1)),
    "P1.kingdom_absences": lambda d: _call_in(d, "phylum", "Streptophyta"),
    "P1.relaxed_controls": lambda d: _edit(
        d / R / "s20_sweep/relaxed_hits.tsv",
        _set({"profile": "PF08709", "group": "viridiplantae"},
             full_evalue="1e-10", hmm_coverage="0.9")),
    "P1.absence_targets": lambda d: _call_in(d, "class", "Cestoda"),
    "P1.absences": lambda d: _gene_in(d, "Ascomycota"),
    "P1.copy_number": lambda d: _drop_first(d / R / "s23_scope/copies.tsv"),
    "P1.record_chase": lambda d: _edit(d / R / "s20_sweep/plant_fungal_verdicts.tsv",
                                       _first({"verdict": "real_gene"}, length=1500)),
    "LEDGER.claims": lambda d: _edit(d / "manuscript/claims_check.tsv",
                                     _first({}, verdict="fail")),
    "S0.c4_symmetry": lambda d: _json(d / META, lambda m: m.update(c4_residual_rmsd_A=0.5)),
    "S0.selectivity_filter": lambda d: _json(d / META, lambda m: m.update(filter_min_radius_A=6.0)),
    "S0.gate": lambda d: _json(d / META, lambda m: m.update(gate_lining_residues=["PHE2513"])),
    "S0.pore_profile": lambda d: _edit(
        d / R / "s0_baseline/review_figures/structure_pore.tsv",
        lambda r: {**r, "min_heavy_atom_radius": str(float(r["min_heavy_atom_radius"]) + 1)}),
    "S0.ip3_contacts": lambda d: _json(d / META, lambda m: m["ip3_contacts"][
        "same_subunit"].__setitem__(0, "LYS266")),
    "P6.shell_agreement": lambda d: _edit(
        d / R / "ligand_site/shell_agreement.tsv",
        lambda r: {**r, "n_contact_le_4.5A": "0"} if r["pdb_id"] == "8TKG" else r),
    "P6.module_map": lambda d: _edit(d / R / "ligand_site/module_map.tsv",
                                     _set({"paralog": "ITPR2", "definition": "channel_all"},
                                          end=2552)),
    # The pattern, not only a number: ITPR2 made significant in the table.
    "P6.module_contrast": lambda d: _edit(
        d / R / "ligand_site/module_contrast.tsv",
        _set({"paralog": "ITPR2", "core_definition": "contact_span",
              "pore_definition": "channel_minus_luminal"}, p_wilcoxon=0.001)),
    # An input plant: every ITPR3 tip given the reference's pore, so the
    # core no longer leads with the loop counted in.
    "P6.loop_reverses": lambda d: _pore_to_reference(d, "ITPR3"),
    # The real verdict is a discrepancy; the flip is a looser cutoff that
    # takes R503 in, run with the modified registry explicitly allowed.
    "P6.contacts_heavy_atom": ("param", "ligand.contact_cutoff", 5.0),
    "P6.shell_distances": lambda d: _edit(d / R / "ligand_site/ligand_shells.tsv",
                                          _first({"shell": "third"}, shell="second")),
    "P6.shell_constraint": lambda d: _edit(d / R / "ligand_site/shell_constraint.tsv",
                                           _set({"paralog": "ITPR2", "shell": "fourth"},
                                                mean_jsd=0.70)),
    # the pattern: ITPR1's trend made significant in the published table
    "P6.shell_trend": lambda d: _edit(d / R / "ligand_site/shell_trend.tsv",
                                      _set({"paralog": "ITPR1"}, p_distance_vs_jsd=0.001)),
    # an input plant: no published table states the step, so make one
    "P6.no_contact_step": _contact_shell_conserved,
}


def _clear_caches() -> None:
    checks_constraint.per_residue.cache_clear()
    module_contrast.clear_caches()
    shell_constraint.clear_caches()


def _copy_sources(check, dest: Path) -> None:
    for rel in check.sources:
        src_root = GENES_DIR / "results"
        if rel.startswith("../"):
            src_root, rel = GENES_DIR, rel[3:]
        for p in sorted(src_root.glob(rel)):
            target = dest / p.resolve().relative_to(GENES_DIR.resolve())
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(p, target)


def test_every_check_has_a_calibration():
    missing = [c.id for c in C.all_checks() if c.id not in PLANTS]
    assert not missing, f"checks with no planted flip: {missing}"


@needs_genes
@pytest.mark.parametrize("check_id", sorted(PLANTS))
def test_check_flips_on_planted_input(check_id, tmp_path, monkeypatch):
    check = C.REGISTRY.get(check_id) or {c.id: c for c in C.all_checks()}[check_id]
    if check.structures:
        needs = needs_structure(*check.structures)
        if needs.args[0]:
            pytest.skip(needs.kwargs["reason"])
    (tmp_path / "results").mkdir()
    _copy_sources(check, tmp_path)
    monkeypatch.setenv("IP3R_GENES_DIR", str(tmp_path))
    _clear_caches()
    before = C.run_check(check).status
    assert before in ("confirmed", "discrepancy"), \
        f"{check_id} did not run on its declared sources: {C.run_check(check).outcome.detail}"
    plant = PLANTS[check_id]
    try:
        if isinstance(plant, tuple):
            PARAMETERS.set_value(plant[1], plant[2])
            after = C.run_check(check, allow_modified=True).status
        else:
            plant(tmp_path)
            _clear_caches()
            after = C.run_check(check).status
    finally:
        PARAMETERS.reset()
        _clear_caches()
    assert after != before, f"{check_id}: planted change did not flip {before}"


def test_missing_project_is_not_run(tmp_path, monkeypatch):
    monkeypatch.setenv("IP3R_GENES_DIR", str(tmp_path / "nowhere"))
    _clear_caches()
    out = C.run_check({c.id: c for c in C.all_checks()}["P5.vus_count"])
    assert out.status == "not_run"


def test_modified_registry_refuses(monkeypatch):
    try:
        PARAMETERS.set_value("check.stat_tol", 0.05)
        out = C.run_check({c.id: c for c in C.all_checks()}["P5.vus_count"])
        assert out.status == "not_run" and "parameters differ" in out.outcome.detail
    finally:
        PARAMETERS.reset()
