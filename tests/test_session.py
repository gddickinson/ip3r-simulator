"""Sessions: the view round-trips; nothing that is a result can ride along;
a malformed file is refused by name before any restore starts."""

from __future__ import annotations

import json
from dataclasses import fields

import pytest

from ip3r.io.session import (SESSION_FORMAT, Session, load_session,
                             parameter_differences, save_session)
from ip3r.parameters import PARAMETERS


def _session() -> Session:
    return Session(structure="8TKG", n_atoms=12345, style="tube", color_by="chain",
                   visible_chains=["A", "C"], sites=["ip3_contact"], show_pore=True,
                   show_lumen=True,
                   tab="Channel", camera_rotation=[0.5, 0.5, 0.5, 0.5],
                   camera_pivot=[1.0, 2.0, 3.0], camera_distance=420.0,
                   camera_pan=[0.0, -4.0, 0.0],
                   transition={"end": "8TKF", "fit": "pore", "method": "restrained",
                               "frame": 7, "paint": True},
                   dynamics={"gating": "mak", "puff_n": 12.0, "puff_model": "park-drive",
                             "tab": "Puffs"},
                   modes={"index": 6, "amplitude": 12.0},
                   variants={"paralog": "ITPR3", "class": "VUS", "layer": None,
                             "draw": True},
                   parameters={"ligand.contact_cutoff": 4.0})


def test_round_trip(tmp_path):
    s = _session()
    path = save_session(s, tmp_path / "a" / "s.json")
    back = load_session(path)
    assert back.as_dict() == s.as_dict()
    assert back.saved_at and back.format_version == SESSION_FORMAT


def test_holds_the_view_and_its_inputs_only():
    """A new field is a decision: a result stored here would go stale."""
    assert {f.name for f in fields(Session)} == {
        "structure", "n_atoms", "style", "color_by", "layer", "show_ligands",
        "visible_chains", "sites", "show_pore", "show_lumen", "completeness", "tab", "camera_rotation",
        "camera_pivot", "camera_distance", "camera_pan", "camera_slab", "orthographic",
        "transition", "dynamics", "modes", "variants", "parameters", "notes", "format_version",
        "software_version", "saved_at"}


def test_unknown_keys_are_dropped_newer_format_refused():
    d = _session().as_dict() | {"from_the_future": 1}
    assert Session.from_dict(d).structure == "8TKG"
    with pytest.raises(ValueError, match="newer"):
        Session.from_dict(d | {"format_version": SESSION_FORMAT + 1})


@pytest.mark.parametrize("key, value", [
    ("camera_distance", "far"), ("camera_distance", float("nan")),
    ("show_pore", 1), ("show_lumen", "yes"), ("sites", ["ok", 3]), ("camera_pivot", [0, 0, True]),
    ("n_atoms", 1.5), ("structure", 6)])
def test_wrong_type_is_refused_by_name(key, value):
    with pytest.raises(ValueError, match=key):
        Session.from_dict(_session().as_dict() | {key: value})


@pytest.mark.parametrize("change, match", [
    ({"camera_rotation": [0, 0, 0, 0]}, "zero quaternion"),
    ({"camera_rotation": [1, 0, 0]}, "4 numbers"),
    ({"camera_distance": -1.0}, "positive"),
    ({"transition": {"end": "8TKF", "coords": []}}, "unknown keys"),
    ({"transition": {"frame": 3}}, "'end'"),
    ({"parameters": {"ligand.contact_cutoff": "4"}}, "number"),
    ({"dynamics": {"trace": [1.0, 2.0]}}, "dynamics holds"),     # a result, not a setting
    ({"variants": {"rows": {"a": 1}}}, "variants holds"),
    ({"dynamics": {"puff_n": float("inf")}}, "dynamics holds"),
    ({"modes": {"index": -1, "amplitude": 12}}, "modes needs"),
    ({"modes": {"index": True, "amplitude": 12}}, "modes needs"),
    ({"modes": {"index": 3}}, "modes needs")])
def test_malformed_is_refused(change, match):
    with pytest.raises(ValueError, match=match):
        Session.from_dict(_session().as_dict() | change)


def test_a_file_without_panel_views_opens_with_them_empty(tmp_path):
    """Round 7.5 added the panel views without a format bump."""
    d = {k: v for k, v in _session().as_dict().items()
         if k not in ("dynamics", "modes", "variants")}
    s = Session.from_dict(d)
    assert s.dynamics == s.modes == s.variants == {}
    assert "animating" not in s.describe()
    assert "mode #7 animating" in _session().describe()


def test_rotation_is_normalised():
    s = Session.from_dict(_session().as_dict() | {"camera_rotation": [2, 0, 0, 0]})
    assert s.camera_rotation == [1.0, 0.0, 0.0, 0.0]


def test_integers_are_accepted_where_floats_are_stored(tmp_path):
    d = _session().as_dict() | {"camera_distance": 400, "camera_pivot": [0, 1, 2]}
    (tmp_path / "s.json").write_text(json.dumps(d))
    s = load_session(tmp_path / "s.json")
    assert s.camera_distance == 400.0 and s.camera_pivot == [0.0, 1.0, 2.0]


def test_not_an_object_is_refused(tmp_path):
    (tmp_path / "s.json").write_text("[1, 2]")
    with pytest.raises(ValueError, match="object"):
        load_session(tmp_path / "s.json")
    with pytest.raises(FileNotFoundError):
        load_session(tmp_path / "missing.json")


def test_parameter_differences():
    assert parameter_differences({}, {}) == {}
    assert parameter_differences({"a": 1.0}, {"a": 1.0}) == {}
    assert parameter_differences({"a": 1.0, "b": 2.0}, {"b": 3.0, "c": 4.0}) == {
        "a": (1.0, None), "b": (2.0, 3.0), "c": (None, 4.0)}


def test_replace_is_the_whole_set_with_one_notification():
    key, other = "ligand.contact_cutoff", "ligand.shell_radius"
    calls = []
    listener = lambda: calls.append(1)          # noqa: E731
    PARAMETERS.subscribe(listener)
    try:
        PARAMETERS.set_value(other, PARAMETERS.default(other) + 1.0)
        calls.clear()
        target = PARAMETERS.default(key) + 0.25
        unknown = PARAMETERS.replace({key: target, "no.such.key": 1.0})
        assert unknown == ["no.such.key"]
        assert PARAMETERS.overrides() == {key: target}     # `other` went back
        assert len(calls) == 1
        PARAMETERS.replace({key: target})
        assert len(calls) == 1                             # no change, no call
    finally:
        PARAMETERS.unsubscribe(listener)
        PARAMETERS.reset()
