"""Round 7.18's GUI additions, headless: the HUD's scale bar, the selection
and its distances, the sequence window's model and its numbering rule, the
Analyses menu against the CLI, and the guide's shortcut table against the
menus."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from conftest import c4_tetramer, make_structure, needs_structure
from ip3r.cli import build_parser
from ip3r.ui.analyses import ANALYSES, GROUPS, command, stamp
from ip3r.ui.help_content import SHORTCUTS
from ip3r.ui.hud import NICE_LENGTHS, HudSettings, nice_scale_length
from ip3r.ui.selection import SelectionController, residue_mask
from ip3r.ui.sequence_model import TRACKS, chain_sequence, track

UI = Path(__file__).resolve().parent.parent / "ip3r" / "ui"


# ------------------------------------------------------------------- HUD

def test_the_scale_bar_takes_the_largest_round_length_that_fits():
    assert nice_scale_length(37.0) == 20.0
    assert nice_scale_length(50.0) == 50.0
    assert nice_scale_length(0.3) == NICE_LENGTHS[0]      # never zero
    assert nice_scale_length(1e6) == NICE_LENGTHS[-1]


def test_the_scale_label_names_the_pivot_only_in_perspective():
    from types import SimpleNamespace
    from ip3r.ui.hud import HudOverlay
    cam = SimpleNamespace(orthographic=False)
    hud = SimpleNamespace(viewport=SimpleNamespace(scene=SimpleNamespace(camera=cam)))
    assert HudOverlay.scale_label(hud, 50.0) == "50 Å at the pivot"
    cam.orthographic = True
    assert HudOverlay.scale_label(hud, 50.0) == "50 Å"


def test_hud_settings_round_trip_and_ignore_unknown_keys():
    s = HudSettings(gnomon=False)
    assert HudSettings.from_dict(s.as_dict()) == s
    assert HudSettings.from_dict({"gnomon": False, "bogus": 1}) == s


# ------------------------------------------------------------- selection

class _Batch:
    def __init__(self):
        self.data, self.alpha = None, 1.0

    def upload(self, *a):
        self.data = a


class _Scene:
    def __init__(self):
        self.batches = {}

    def remove(self, name):
        self.batches.pop(name, None)

    def spheres(self, name):
        return self.batches.setdefault(name, _Batch())

    cylinders = spheres


class _Viewport:
    labels = []

    def set_overlay_labels(self, labels):
        self.labels = list(labels)

    def update(self):
        pass


class _Hud:
    def __init__(self):
        self.readouts = {}

    def set_readout(self, k, v):
        if v:
            self.readouts[k] = v
        else:
            self.readouts.pop(k, None)


class _SceneCtl:
    def __init__(self, st):
        self.structure, self.view, self.fit_target = st, None, "all"
        self.scene, self.viewport = _Scene(), _Viewport()
        self.followers, self.on_clear = [], []

    def describe_atom(self, i):
        return f"atom {i}"


@pytest.fixture
def sel():
    st = make_structure(c4_tetramer(n_per=12))
    ctl = _SceneCtl(st)
    return SelectionController(ctl, _Hud(), lambda *_: None), ctl, st


def test_residue_mask_picks_each_pair_on_its_own_chain_only():
    st = make_structure(c4_tetramer(n_per=12))
    m = residue_mask(st, {("A", 3), ("C", 5)})
    assert m.sum() == 2
    assert set(zip(st.chain[m], st.res_seq[m])) == {("A", 3), ("C", 5)}


def test_a_click_selects_a_residue_and_shift_toggles(sel):
    s, ctl, st = sel
    i = int(np.flatnonzero((st.chain == "B") & (st.res_seq == 4))[0])
    s.pick(i)
    assert s.residues == {("B", 4)}
    assert "selection" in ctl.scene.batches
    j = int(np.flatnonzero((st.chain == "B") & (st.res_seq == 6))[0])
    s.pick(j, extend=True)
    s.pick(i, extend=True)
    assert s.residues == {("B", 6)}
    assert s.select_everywhere(6) == 4 and len(s.residues) == 4


def test_two_measured_atoms_give_their_distance_and_follow_a_frame(sel):
    s, ctl, st = sel
    s.arm(True)
    s.pick(0)
    assert not s.distances
    s.pick(5)
    d = float(np.linalg.norm(st.xyz[0] - st.xyz[5]))
    assert s.distances == [(0, 5)]
    assert f"{d:.2f} Å" in ctl.viewport.labels[0][1]
    assert f"{d:.2f} Å" in s.hud.readouts["distance"]
    # A frame moves atom 5 by 3 Å along x: the label follows the frame.
    xyz = st.xyz.astype(float).copy()
    xyz[5] += (3.0, 0.0, 0.0)
    for cb in ctl.followers:
        cb(xyz)
    assert f"{np.linalg.norm(xyz[0] - xyz[5]):.2f} Å" in ctl.viewport.labels[0][1]


def test_a_new_deposit_drops_the_selection_and_the_distances(sel):
    s, ctl, st = sel
    s.pick(0)
    s.add_point(0)
    s.add_point(3)
    for cb in ctl.on_clear:
        cb()
    assert not s.residues and not s.distances and not s.hud.readouts


def test_centring_stops_the_automatic_fit(sel):
    s, ctl, st = sel

    class _Cam:
        pivot, pan, distance, slab_front = np.zeros(3), np.ones(3), 500.0, 3.0
    ctl.scene.camera = _Cam()
    s.centre_on(7)
    assert np.allclose(ctl.scene.camera.pivot, st.xyz[7])
    assert ctl.fit_target is None and ctl.scene.camera.distance == 40.0


# ------------------------------------------------------ sequence window

def test_the_built_sequence_without_a_construct_is_all_resolved():
    st = make_structure(c4_tetramer(n_per=12))
    s = chain_sequence(st, "A")
    assert s.source == "built" and len(s) == 12 and s.n_resolved == 12
    assert s.letters == "A" * 12 and s.positions == tuple(range(1, 13))


def test_a_deposit_in_no_numbering_is_painted_grey():
    st = make_structure(c4_tetramer(n_per=12))
    s = chain_sequence(st, "A")
    for key in ("element", "constraint"):
        assert all(not d.background for d in track(s, key, None).values())


@needs_structure("6DQN", "7LHF")
def test_6dqn_keeps_its_unresolved_residues_and_paints_the_filter():
    from ip3r.core.annotations import functional_sites
    from ip3r.io import loader
    from ip3r.render.colormaps import MISSING
    st = loader.load("6DQN")
    s = chain_sequence(st, "A")
    built = set(st.res_seq[(st.chain == "A") & ~st.hetero].tolist())
    assert s.source == "construct" and s.n_resolved == len(built)
    assert len(s) > s.n_resolved                  # gaps kept, not closed
    assert list(s.positions) == sorted(s.positions)
    el = track(s, "element", "ITPR3")
    filt = functional_sites("ITPR3")["filter_lining"][0]
    assert el[filt].background and el[filt].underline
    grey = "#" + "".join(f"{int(round(255 * c)):02x}" for c in MISSING)
    js = track(s, "constraint", "ITPR3")
    assert all(d.background != grey for d in js.values())   # missing = no fill
    # Rat 7LHF: the residue-keyed tracks are empty when no numbering fits.
    rat = chain_sequence(loader.load("7LHF"), loader.load("7LHF").chains[0])
    assert all(not d.background for d in track(rat, "element", None).values())
    assert set(TRACKS) >= {"element", "constraint", "resolved"}


# -------------------------------------------------------- analyses menu

def test_every_analysis_is_a_command_the_cli_accepts():
    parser = build_parser()
    assert {a.group for a in ANALYSES} <= set(GROUPS)
    assert len({a.key for a in ANALYSES}) == len(ANALYSES)
    for a in ANALYSES:
        args = command(a, "8TKF", "ITPR3")
        ns = parser.parse_args(args)                 # SystemExit on a bad flag
        assert ns.cmd == args[0], a.key
        assert "{" not in " ".join(args)


def test_an_analysis_needing_a_deposit_refuses_without_one():
    info = next(a for a in ANALYSES if a.key == "info")
    with pytest.raises(ValueError):
        command(info, None, None)
    states = next(a for a in ANALYSES if a.key == "states")
    assert command(states, None, None)[-1] == "ITPR3"


def test_the_stamp_names_the_command_the_deposit_and_the_parameters():
    s = stamp(["lumen", "8TKF"], "8TKF", {})
    assert "python -m ip3r lumen 8TKF" in s and "registered defaults" in s
    s = stamp(["born"], None, {"born.reach": 8.0})
    assert "deposit on screen: none" in s and "born.reach = 8" in s


# ------------------------------------------------------------- the guide

def test_the_guide_lists_every_menu_shortcut_and_no_other():
    text = (UI / "menus.py").read_text()
    used = set(re.findall(r'"((?:Ctrl\+[\w+]+)|F1|Esc|Space)"', text))
    listed = {k for k, _ in SHORTCUTS if re.fullmatch(r"(Ctrl\+[\w+]+)|F1|Esc|Space", k)}
    assert used == listed
