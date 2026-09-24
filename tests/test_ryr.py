"""RyR1: the curated resource, the curation rules, then the real panel.

The resource is checked against anchors the literature fixes independently
of the curation (the GGGIGD filter at 4894-4899, the I4937 gate), and the
rules are checked on hand-made entries, each built to fail one rule. On the
real panel the findings are pinned, including the one that refutes the
salt-bridge reading: D4899 is bridged, yet D4899Q cuts the conductance 5x.
"""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from ip3r.config import RESOURCE_DIR
from ip3r.core.annotations import (element_of, functional_sites, is_ryr,
                                   reference_sequence, residue_constraint)
from ip3r.io.registry import load_registry
from conftest import needs_structure

ROOT = Path(__file__).resolve().parent.parent
PANEL = ("7TDG", "7TDI", "8RRX", "9HEO", "9OL6", "9R8O")


def _curate():
    spec = importlib.util.spec_from_file_location("curate_ryr", ROOT / "scripts" / "curate_ryr.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_resource_has_provenance_and_the_literature_anchors():
    d = json.loads((RESOURCE_DIR / "ryr1.json").read_text())
    assert len(d["provenance"]["sources"]) == 4
    assert all(len(h) == 64 for h in d["provenance"]["sources"].values())
    seq = reference_sequence("RYR1")
    assert len(seq) == 5037
    assert seq[4893:4899] == "GGGIGD" and seq[4936] == "I"
    assert element_of("RYR1", 4899) == "channel"
    assert is_ryr("RYR1") and not is_ryr("ITPR3")


def test_ryr_has_no_publication_annotation_and_paints_grey():
    assert functional_sites("RYR1") == {}
    assert np.isnan(residue_constraint("RYR1")).all()
    assert len(residue_constraint("RYR1")) == 5037


def test_registry_holds_one_deposit_per_state_and_a_same_paper_morph_pair():
    ryr = [e for e in load_registry() if e.family == "RyR"]
    assert tuple(sorted(e.pdb_id for e in ryr)) == PANEL
    assert all(e.paralog == "RYR1" and e.uniprot == "P11716" for e in ryr)
    start = next(e for e in ryr if "morph_start" in e.roles)
    end = next(e for e in ryr if "morph_end" in e.roles)
    assert (start.pdb_id, end.pdb_id) == ("9R8O", "9HEO")
    raw = {s["pdb_id"]: s for s in json.loads((RESOURCE_DIR / "ryr1.json").read_text())["structures"]}
    assert raw["9R8O"]["doi"] == raw["9HEO"]["doi"]


def _entry(title="Open-state RyR1", ligands=("CA",), mutation=None, n=4,
           length=5037, modelled=17000, res=3.2, method="ELECTRON MICROSCOPY"):
    ref = {"database_accession": "P11716"}
    return {"struct": {"title": title}, "exptl": [{"method": method}],
            "rcsb_entry_info": {"resolution_combined": [res],
                                "deposited_modeled_polymer_monomer_count": modelled},
            "nonpolymer_entities": [{"nonpolymer_comp": {"chem_comp": {"id": x}}} for x in ligands],
            "polymer_entities": [{
                "rcsb_polymer_entity": {"pdbx_mutation": mutation, "pdbx_number_of_molecules": n},
                "rcsb_polymer_entity_container_identifiers": {"reference_sequence_identifiers": [ref]},
                "entity_poly": {"rcsb_sample_sequence_length": length}}]}


def test_each_curation_rule_rejects_its_own_violation():
    c = _curate()
    assert c.verdict(_entry()) is None
    assert c.verdict(_entry(method="X-RAY DIFFRACTION")).startswith("1")
    assert c.verdict(_entry(length=536)).startswith("1")             # a domain crystal
    assert c.verdict(_entry(modelled=4000)).startswith("1")
    assert c.verdict(_entry(title="RyR1 (Local Refinement of TMD)")).startswith("1/2")
    assert c.verdict(_entry(mutation="Y523S")).startswith("2")
    assert c.verdict(_entry(ligands=("CA", "U1C"))).startswith("3")  # dantrolene
    assert c.verdict(_entry(res=4.3)).startswith("4")


def test_state_words_do_not_read_closed_inactivated_as_inactivated():
    c = _curate()
    assert c.state_of("in closed-inactivated conformation") == "closed-inactivated"
    assert c.state_of("in inactivated conformation") == "inactivated"
    assert c.state_of("in detergent in close state") == "closed"
    assert c.state_of("Open-state RyR1") == "open"
    assert c.state_of("membrane embedded skeletal muscle ryanodine receptor") is None


@needs_structure(*PANEL)
def test_panel_is_in_rabbit_numbering_and_gates_at_i4937():
    from ip3r.structure.states import state_panel
    rows = state_panel("RYR1")
    assert len(rows) == len(PANEL)
    assert all(r.summary.numbering.paralog == "RYR1" and r.summary.numbering.identity > 0.99
               for r in rows)
    shut = [r for r in rows if r.state in ("closed", "primed")]
    assert all("ILE4937" in r.summary.constrictions["gate"].residues for r in shut)
    open_ = max(rows, key=lambda r: r.radius("gate"))
    assert open_.pdb_id == "9HEO" and open_.radius("gate") > 1.5 * max(
        r.radius("gate") for r in rows if r is not open_)


@needs_structure("9HEO")
def test_open_ryr1_falls_short_and_the_bridged_filter_aspartate_refutes_pairing():
    from ip3r.io import loader
    from ip3r.physics.ryr_mutants import mutant_panel
    from ip3r.physics.unitary import published
    wt, rows = mutant_panel(loader.load("9HEO"))
    assert wt.bath == pytest.approx(0.25)
    measured = published("RYR1")["Xu 2006 (bilayer)"]
    # The same shortfall as ITPR3: the neutral continuum is several-fold short.
    assert wt.neutral.conductance_pS < measured / 4
    by = {r.name: r for r in rows}
    assert not by["E4955Q"].lining and by["E4955Q"].charged_ratio == pytest.approx(1.0)
    for name in ("D4899Q", "E4900N", "D4938N", "D4945N"):
        assert by[name].lining and by[name].charged_ratio < 1.0     # direction right
    d = by["D4899Q"]
    assert d.bridged and d.measured_ratio < 0.3
    assert d.paired_ratio > 0.95          # pairing predicts no effect: refuted
