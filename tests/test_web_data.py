"""The web site's data: precomputed files, fixtures and API answers are the Python code's.

* The JS fixtures (``tests/fixtures/anm_reduced.json``, ``gating_curves.json``)
  are still what the repo's Python gives (numbers compared with a tolerance:
  libm differs between macOS and Linux in the last digit).
* The precomputed files decode to what they claim, and — where the deposits
  are downloaded (``ref/``) — are still what ``ip3r.web.export_*`` produce.
* The API refuses what it does not have, and the export stays small.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from ip3r.io import loader  # noqa: E402
from ip3r.web.app import DEPOSITS, MODE_DEPOSITS, PRECOMPUTED, app  # noqa: E402
from ip3r.web.codec import decode_coords, decode_i16, encode_coords, encode_i16  # noqa: E402
from ip3r.web.reduced_anm import FIXTURE_CASES, reduced_modes  # noqa: E402

FIXTURES = ROOT / "tests" / "fixtures"
client = TestClient(app)


def _json(name):
    return json.loads((PRECOMPUTED / name).read_text())


def test_codec_round_trip_and_range():
    x = np.array([[0.0, -1.234, 327.0], [12.345, -200.0, 0.005]])
    assert np.allclose(decode_coords(encode_coords(x)), np.round(x / 0.01) * 0.01)
    assert np.allclose(decode_i16(encode_i16([0.5, -0.25], 1 / 30000), 1 / 30000), [0.5, -0.25], atol=1e-4)
    with pytest.raises(ValueError):
        encode_coords(np.array([[400.0, 0, 0]]))


def test_anm_fixture_is_current():
    fx = json.loads((FIXTURES / "anm_reduced.json").read_text())
    assert [(c["pdb_id"], c["stride"], c["cutoff"]) for c in fx["cases"]] == [c[:3] for c in FIXTURE_CASES]
    for case, (pdb, stride, cutoff, n) in zip(fx["cases"], FIXTURE_CASES):
        now = reduced_modes(_json(f"structure_{pdb}.json"), stride, cutoff, n)
        assert now["n_sites"] == case["n_sites"] and now["residues"] == case["residues"]
        assert now["symmetry"] == case["symmetry"]
        assert np.allclose(now["eigenvalues"], case["eigenvalues"], rtol=1e-6, atol=1e-12)
        assert np.allclose(now["collectivity"], case["collectivity"], atol=1e-6)


def test_anm_fixture_is_well_posed():
    # One connected network and no near-zero internal modes, so every mode is defined (E pairs as a pair).
    for case in json.loads((FIXTURES / "anm_reduced.json").read_text())["cases"]:
        assert case["n_components"] == 1
        assert min(case["eigenvalues"]) > 1e-6 * max(case["eigenvalues"])
        assert set(case["symmetry"]) <= {"A", "B", "E"}


def test_gating_fixture_is_current():
    import make_anm_fixture
    fx = json.loads((FIXTURES / "gating_curves.json").read_text())
    now = make_anm_fixture.gating_fixture()
    assert now["constants"] == fx["constants"]
    for m in fx["curves"]:
        for p in fx["curves"][m]:
            assert np.allclose(now["curves"][m][p], fx["curves"][m][p], rtol=1e-9, atol=1e-15)
    for m, rows in fx["bells"].items():
        for a, b in zip(now["bells"][m], rows):
            assert np.isclose(a["c_peak"], b["c_peak"], rtol=1e-6)
            assert np.isclose(a["po_peak"], b["po_peak"], rtol=1e-9)


def test_precomputed_files_decode():
    for pdb in DEPOSITS:
        rec = _json(f"structure_{pdb}.json")
        assert len(rec["subunits"]) == 4 and rec["chains"] == [s["chain"] for s in rec["subunits"]]
        for s in rec["subunits"]:
            xyz = decode_coords(s["xyz"])
            assert xyz.shape == (len(s["resi"]), 3)
            assert len(s["bfactor"]) == len(s["resi"])
            if rec["painted_in"]:
                assert len(s["element"]) == len(s["jsd"]) == len(s["resi"])
            else:
                assert "element" not in s
        assert {"gate", "filter"} <= set(rec["constrictions"])
    for pdb in MODE_DEPOSITS:
        m = _json(f"modes_{pdb}.json")
        assert decode_coords(m["sites"]).shape == (m["n_sites"], 3)
        assert all(len(decode_i16(v, m["vector_scale"])) == 3 * m["n_sites"] for v in m["vectors"])
        assert len(m["symmetry"]) == len(m["eigenvalues"]) == len(m["vectors"])
    t = _json("transition_8TKG_8TKF.json")
    n = 4 * t["n_residues"]
    assert all(decode_coords(f).shape == (n, 3) for f in t["morph"]["frames"])
    assert t["morph"]["kept"][0] == 0 and t["morph"]["kept"][-1] == t["morph"]["n_frames"] - 1


def test_precomputed_say_what_the_readme_says():
    """The README's state table, transition and overlap numbers, read from the precomputed files."""
    gate = {r["pdb_id"]: round(r["gate"]["radius"], 2) for r in _json("states_ITPR3.json")["rows"]}
    assert gate == {"8TKH": 1.95, "8TLA": 2.44, "7T3P": 2.53, "6DQN": 2.55, "6DQJ": 2.69, "8TKG": 2.73,
                    "8TKF": 5.85, "7T3T": 5.34}
    t = _json("transition_8TKG_8TKF.json")
    assert t["n_residues"] == 2194
    assert round(t["gate_path"]["half_open"], 2) == 0.42
    a = t["overlap"]["A_subspace"]
    assert round(a["cumulative"][-1], 2) == 0.67 and round(t["overlap"]["null_cumulative"][-1], 2) == 0.04


def test_structures_painted_only_in_their_numbering():
    for pdb in DEPOSITS:
        rec = _json(f"structure_{pdb}.json")
        num = rec["numbering"]
        assert rec["painted_in"] == (num["paralog"] if num and num["passed"] else None)


@pytest.mark.skipif(not all(loader.is_local(p) for p in ("6DQN", "8TKG")), reason="deposits not downloaded (ref/)")
def test_precomputed_are_current_where_the_deposits_are_here():
    from ip3r.web.export_dynamics import modes_record
    from ip3r.web.export_structures import structure_record
    now, saved = structure_record("6DQN"), _json("structure_6DQN.json")
    assert now["source"]["sha256"] == saved["source"]["sha256"]
    for k in ("constrictions", "ip3_contacts", "chains", "painted_in"):
        assert now[k] == saved[k], k
    assert all(np.allclose(decode_coords(a["xyz"]), decode_coords(b["xyz"]), atol=0.011)
               for a, b in zip(now["subunits"], saved["subunits"]))
    m, ms = modes_record("8TKG"), _json("modes_8TKG.json")
    assert m["symmetry"] == ms["symmetry"]
    assert np.allclose(m["eigenvalues"], ms["eigenvalues"], rtol=1e-5)


def test_api_refuses_unknown_deposits():
    assert client.get("/api/structure/1ABC").status_code == 404
    assert client.get("/api/modes/8TKF").status_code == 404
    assert client.get("/api/structure/6dqn").status_code == 200


def test_gating_api_cites_every_constant():
    g = client.get("/api/gating").json()
    for m in g["models"].values():
        assert set(m["values"]) == set(m["meta"])
        assert all(meta["citation"] for meta in m["meta"].values())


def test_precomputed_stay_small():
    total = sum(p.stat().st_size for p in PRECOMPUTED.glob("*.json"))
    assert total < 10e6, f"{total / 1e6:.1f} MB"
