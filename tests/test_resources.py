"""The imported resources are internally consistent and match ip3r_genes."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from ip3r.config import PARALOGS, RESOURCE_DIR
from ip3r.core.annotations import (element_of, elements, functional_sites,
                                   reference_sequence, residue_constraint, variants)
from ip3r.io.registry import load_registry
from conftest import needs_genes

ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("gene", PARALOGS)
def test_sequences_and_sites(gene):
    seq = reference_sequence(gene)
    assert 2600 < len(seq) < 2800
    fs = functional_sites(gene)
    assert len(fs["ip3_contact"]) == 10
    # The filter motif sits where the domain map puts the filter.
    f = next(e for e in elements(gene) if e.name == "selectivity_filter")
    assert "GGGVGD" in seq[f.start - 3:f.end + 3]
    assert len(residue_constraint(gene)) == len(seq)


def test_specific_elements_win():
    g = next(e for e in elements("ITPR3") if e.name == "gate")
    assert all(element_of("ITPR3", r) == "gate" for r in range(g.start, g.end + 1))
    assert element_of("ITPR3", 1) in ("linker", "nterm_trefoil")


def test_registry():
    ids = {e.pdb_id for e in load_registry()}
    assert {"6DQN", "8TKG", "8TKF", "7LHF", "9YKK"} <= ids
    assert sum(e.ip3_bound for e in load_registry() if not e.is_control) == 6
    assert [e.pdb_id for e in load_registry() if e.is_control] == ["7T3T"]


def test_variants_resource():
    assert sum(v["class_bucket"] == "VUS" for v in variants()) == 1546


def test_every_resource_records_its_sources():
    for name in ("sequences", "domains", "sites", "structures", "constraint", "variants"):
        prov = json.loads((RESOURCE_DIR / f"{name}.json").read_text())["provenance"]
        assert prov["sources"]


@needs_genes
def test_resources_in_step_with_ip3r_genes():
    r = subprocess.run([sys.executable, str(ROOT / "scripts/sync_genes.py"), "--check"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout
