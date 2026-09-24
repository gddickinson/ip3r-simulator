import json
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import needs_structure

from ip3r.parameters import PARAMETERS

ROOT = Path(__file__).resolve().parent.parent


def test_registry_builds_and_validates():
    r = subprocess.run([sys.executable, str(ROOT / "scripts/build_parameters.py"), "--check"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout


def test_json_matches_the_table():
    sys.path.insert(0, str(ROOT / "scripts"))
    from parameter_table import P
    assert {p["key"]: p["value"] for p in P} == {k: p.default for k, p in PARAMETERS.parameters.items()}


def test_overrides_are_tracked():
    try:
        PARAMETERS.set_value("gating.d1", 0.2)
        assert PARAMETERS.modified
    finally:
        PARAMETERS.reset()
    assert not PARAMETERS.modified


# ------------------------------------------------ editing (Round 5.1)


@pytest.fixture
def calls():
    """A listener that counts notifications; the registry is left at defaults."""
    seen = []
    cb = lambda: seen.append(1)  # noqa: E731
    PARAMETERS.subscribe(cb)
    yield seen
    PARAMETERS.unsubscribe(cb)
    PARAMETERS.reset()


def test_listeners_fire_once_per_effective_change(calls):
    d = PARAMETERS.default("gating.d1")
    PARAMETERS.set_value("gating.d1", d)            # no-op: already the default
    assert calls == []
    PARAMETERS.set_value("gating.d1", d * 2)
    PARAMETERS.set_value("gating.d1", d * 2)        # same value again: no-op
    assert calls == [1]
    PARAMETERS.reset("gating.d2")                   # not overridden: no-op
    assert calls == [1]
    PARAMETERS.reset()
    assert calls == [1, 1]
    PARAMETERS.apply({"gating.d1": d * 2, "gating.d2": PARAMETERS.default("gating.d2") * 2})
    assert calls == [1, 1, 1]                        # one notification for the set


def test_clamp_is_reported_by_the_applied_value(calls):
    p = PARAMETERS.get("gating.d1")
    if p.maximum is None:
        pytest.skip("gating.d1 has no upper bound")
    assert PARAMETERS.set_value("gating.d1", p.maximum * 10) == p.maximum


def test_overrides_round_trip_through_a_file(tmp_path, calls):
    PARAMETERS.set_value("gating.d1", PARAMETERS.default("gating.d1") * 3)
    path = tmp_path / "o.json"
    PARAMETERS.write_overrides(path)
    saved = PARAMETERS.overrides()
    PARAMETERS.reset()
    PARAMETERS.set_value("gating.d2", PARAMETERS.default("gating.d2") * 3)
    assert PARAMETERS.read_overrides(path) == []
    assert PARAMETERS.overrides() == saved          # replaced, not merged


def test_import_refuses_a_bad_file_without_changing_anything(tmp_path, calls):
    PARAMETERS.set_value("gating.d1", PARAMETERS.default("gating.d1") * 3)
    before = PARAMETERS.overrides()
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"gating.d2": "high"}))
    with pytest.raises(ValueError):
        PARAMETERS.read_overrides(bad)
    assert PARAMETERS.overrides() == before
    odd = tmp_path / "odd.json"
    odd.write_text(json.dumps({"no.such.key": 1.0}))
    assert PARAMETERS.read_overrides(odd) == ["no.such.key"]
    assert not PARAMETERS.modified                   # the file listed nothing known


def test_import_back_to_defaults_notifies(tmp_path, calls):
    PARAMETERS.set_value("gating.d1", PARAMETERS.default("gating.d1") * 3)
    empty = tmp_path / "e.json"
    empty.write_text("{}")
    n = len(calls)
    PARAMETERS.read_overrides(empty)
    assert not PARAMETERS.modified and len(calls) == n + 1


def test_filter():
    assert PARAMETERS.matches("gating.d1", "GATING")
    assert not PARAMETERS.matches("gating.d1", "no-such-text")
    assert not PARAMETERS.matches("gating.d1", "", only_modified=True)


@needs_structure("6DQN")
def test_a_result_computed_under_an_edit_is_not_served_after_reset():
    """The hazard the listeners exist for: a check's memoised measurement
    made under an edited value must be forgotten when the value changes."""
    from ip3r.analysis import checks_structure as cs
    try:
        PARAMETERS.set_value("pore.step", PARAMETERS.default("pore.step") * 2)
        edited = cs._summary("6DQN")
        assert cs._summary.cache_info().currsize == 1
        PARAMETERS.reset()
        assert cs._summary.cache_info().currsize == 0
        assert cs._summary("6DQN") is not edited
    finally:
        PARAMETERS.reset()
