"""The extended IP3R deposits (``scripts/curate_ip3r.py``): each rule
rejects its own violation, the state words read the titles they are for,
and adding the deposits moves no panel the publication's numbers come from."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from ip3r.config import RESOURCE_DIR
from ip3r.io.registry import load_registry
from conftest import needs_structure

ROOT = Path(__file__).resolve().parent.parent
EXTENDED = tuple(e.pdb_id for e in load_registry() if e.is_extended)


def _curate():
    sys.path.insert(0, str(ROOT))
    spec = importlib.util.spec_from_file_location("curate_ip3r", ROOT / "scripts" / "curate_ip3r.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _entry(method="ELECTRON MICROSCOPY", acc="Q14573", n=4, length=2671,
           modelled=8800, title="Human IP3R3 in the resting state", mutation=None,
           res=3.3, second=None):
    ent = lambda a: {"rcsb_polymer_entity": {"pdbx_mutation": mutation,       # noqa: E731
                                             "pdbx_number_of_molecules": n},
                     "rcsb_polymer_entity_container_identifiers": {
                         "reference_sequence_identifiers": [{"database_accession": a}]},
                     "entity_poly": {"rcsb_sample_sequence_length": length}}
    return {"rcsb_id": "0XXX", "struct": {"title": title}, "exptl": [{"method": method}],
            "rcsb_entry_info": {"resolution_combined": [res],
                                "deposited_modeled_polymer_monomer_count": modelled},
            "polymer_entities": [ent(acc)] + ([ent(second)] if second else [])}


def test_each_rule_rejects_its_own_violation():
    c = _curate()
    assert c.verdict(_entry()) is None
    assert c.verdict(_entry(acc="P29994")) is None                    # rat ITPR1, as 7LHF
    assert c.verdict(_entry(method="X-RAY DIFFRACTION")).startswith("1")
    assert c.verdict(_entry(acc="P11881")).startswith("1")            # mouse: not a numbering
    assert c.verdict(_entry(second="Q14571")).startswith("1")         # two paralogs
    assert c.verdict(_entry(n=5)).startswith("1")                     # 9YNO's dimer
    assert c.verdict(_entry(length=2453)).startswith("1")             # 6UQK's deletion
    assert c.verdict(_entry(modelled=4700)).startswith("1")           # 9YMZ's partial model
    assert c.verdict(_entry(mutation="D2477E")).startswith("2")
    assert c.verdict(_entry(title="IP3R3 mutant K508A")).startswith("2")
    assert c.verdict(_entry(res=4.12)).startswith("3")


def test_state_words_read_the_titles_they_are_for():
    c = _curate()
    assert c.state_of("in the inactive state") == "inactive"          # not "active"
    assert c.state_of("in the active state") == "active"
    assert c.state_of("Preactivated+Ca2+ State (+IP3/ATP)") == "preactivated + Ca2+"
    assert c.state_of("Preactivated State (+IP3/ATP)") == "preactivated"
    assert c.state_of("in the pre-active B state") == "preactivated (B)"
    assert c.state_of("Labile Resting State 2 (+IP3/ATP)") == "labile resting (2)"
    assert c.state_of("Class 2 IP3-bound human type 3") == "IP3-bound (class 2)"
    assert c.state_of("reconstituted into lipid nanodisc in the apo-state") == "apo (nanodisc)"
    assert c.state_of("a receptor revealing a self-binding peptide") is None


def test_resource_is_disjoint_from_the_publication_and_carries_its_rules():
    d = json.loads((RESOURCE_DIR / "ip3r_extended.json").read_text())
    rules = d["provenance"]["rules"]
    assert all(len(h) == 64 for h in d["provenance"]["sources"].values())
    others = set()
    for name in ("structures.json", "ip3r_controls.json"):
        others |= {s["pdb_id"] for s in json.loads((RESOURCE_DIR / name).read_text())["structures"]}
    rows = d["structures"]
    assert rows and not {r["pdb_id"] for r in rows} & others
    assert all(r["roles"] == ["extended"] and r["resolution"] <= rules["max_resolution_A"]
               for r in rows)
    assert all(e.family == "IP3R" and e.state for e in load_registry() if e.is_extended)


def test_no_panel_moves():
    from ip3r.physics.shortfall import open_entries
    from ip3r.structure.states import panel_entries
    assert {e.pdb_id for e in panel_entries("ITPR3")} == {
        "6DQJ", "6DQN", "7T3P", "8TKF", "8TKG", "8TKH", "8TLA"}
    assert not {e.pdb_id for e in panel_entries("ITPR3")} & set(EXTENDED)
    assert set(EXTENDED) <= {e.pdb_id for p in ("ITPR1", "ITPR2", "ITPR3")
                             for e in load_registry() if e.paralog == p}
    assert {e.pdb_id for e in panel_entries("ITPR3", extended=True)} > {
        e.pdb_id for e in panel_entries("ITPR3")}
    assert {e.pdb_id for e in open_entries()} == {"8TKF", "7T3T", "9HEO"}


def test_human_extended_deposits_superpose_on_their_paralog():
    from ip3r.structure.superpose import candidates
    assert "8TK8" in dict(candidates("8TKG", "ITPR3"))
    assert set(dict(candidates("9YKK", "ITPR2"))) == {"9YKY", "9YLI"}


@needs_structure(*EXTENDED)
def test_every_extended_deposit_measures_in_its_numbering_and_none_is_open():
    from ip3r.io import loader
    from ip3r.structure.channel import measure_channel
    for e in load_registry():
        if not e.is_extended:
            continue
        s = measure_channel(loader.load(e.pdb_id))
        want = e.paralog if e.human else None           # rat fits no human numbering
        assert (s.numbering.paralog if s.numbering else None) == want, e.pdb_id
        # So the open panel stays 8TKF and 7T3T: no extended pore is open.
        assert s.constrictions["gate"].radius < 3.0, e.pdb_id


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
