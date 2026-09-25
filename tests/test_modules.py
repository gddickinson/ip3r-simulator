"""The Paper 6 modules and the alignment column map they are measured on."""

import numpy as np
import pytest

from ip3r.analysis import module_contrast as MC
from ip3r.core.annotations import functional_sites, reference_sequence
from ip3r.core.modules import Module, ModuleRefusal, _validate, module, modules
from conftest import needs_genes


@pytest.mark.parametrize("gene", ["ITPR1", "ITPR2", "ITPR3"])
def test_primary_modules_hold_their_sites_and_are_disjoint(gene):
    core, pore = modules(gene)
    fs = functional_sites(gene)
    assert set(fs["ip3_contact"]) <= set(core.residues)
    assert set(fs["filter_lining"]) | set(fs["gate_lining"]) <= set(pore.residues)
    assert not set(core.residues) & set(pore.residues)
    # the luminal loop is a hole, not a truncation
    assert pore.end == module(gene, "channel_all").end
    assert len(pore.residues) < len(module(gene, "channel_all").residues)


def test_validation_refuses_a_core_missing_a_contact():
    c = functional_sites("ITPR3")["ip3_contact"]
    short = Module("ITPR3", "ligand_core", "contact_span", min(c) + 1, max(c))
    with pytest.raises(ModuleRefusal):
        _validate(short, functional_sites("ITPR3"))


@needs_genes
@pytest.mark.parametrize("gene", ["ITPR1", "ITPR3"])
def test_column_map_lands_on_the_reference_residue(gene):
    MC.clear_caches()
    ref = MC._reference_row(gene)[1]
    seq = reference_sequence(gene)
    cmap = MC.reference_columns(gene)
    for r in (1, 266, 2400, len(seq)):
        assert ref[cmap[r]].upper() == seq[r - 1].upper()


@needs_genes
def test_tip_identity_matches_brute_force():
    MC.clear_caches()
    gene, d = "ITPR2", "contact_span"
    label, row = MC.read_alignment(gene)[5]
    ref = MC._reference_row(gene)[1]
    cols = [MC.reference_columns(gene)[r] for r in module(gene, d).residues]
    cov = [c for c in cols if row[c] != "-"]
    want = sum(row[c].upper() == ref[c].upper() for c in cov) / len(cov)
    ident, coverage = MC.tip_identities(gene, d)[label]
    assert abs(ident - want) < 1e-12 and abs(coverage - len(cov) / len(cols)) < 1e-12


@needs_genes
def test_wrong_reference_row_is_refused(monkeypatch):
    MC.clear_caches()
    real = MC.read_alignment("ITPR3")
    tampered = tuple((lab, row[:100] + row[101:] + "-") if lab.startswith("REF|")
                     else (lab, row) for lab, row in real)
    monkeypatch.setattr(MC, "read_alignment", lambda gene: tampered)
    MC.reference_columns.cache_clear()
    with pytest.raises(MC.AlignmentMismatch):
        MC.reference_columns("ITPR3")
    MC.reference_columns.cache_clear()


@needs_genes
def test_coverage_floor_drops_truncated_tips():
    MC.clear_caches()
    pc = MC.paired_contrast("ITPR1")
    assert len(pc.dropped) == 3 and len(pc.tips) == 260
    assert np.all(np.isfinite(pc.diff))


def test_literature_core_carried_and_follows_its_parameters():
    from ip3r.parameters import PARAMETERS
    # S22's module_map.tsv via S17's transfer; ours via our own alignment
    assert [(module(g, "ibc_literature").start, module(g, "ibc_literature").end)
            for g in ("ITPR1", "ITPR2", "ITPR3")] == [(224, 604), (224, 604), (225, 604)]
    try:
        PARAMETERS.set_value("ligand.ibc_end", 560)       # cuts contacts off
        with pytest.raises(ModuleRefusal):
            module("ITPR1", "ibc_literature")
    finally:
        PARAMETERS.reset()
    assert module("ITPR1", "ibc_literature").end == 604
