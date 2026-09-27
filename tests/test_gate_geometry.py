"""Round 7.21: the gate widened. The move is calibrated on synthetic atoms
(exact radial shift, nothing else touched, C4 kept), then on 9HEO (the
gate opens by Δ, the filter stays), before the reading is trusted."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from ip3r.physics.gate_geometry import deltas, widen_gate, widening
from ip3r.structure.pore import Constriction
from ip3r.structure.symmetry import Frame
from conftest import c4_tetramer, make_structure, needs_structure

W = 8.0


def _summary(z_gate=0.0):
    frame = Frame(axis=np.array([0.0, 0.0, 1.0]), centre=np.zeros(3),
                  basis=np.eye(3), chains=list("ABCD"))
    return SimpleNamespace(frame=frame,
                           constrictions={"gate": Constriction("gate", z_gate, 3.0, [])})


def _cylindrical(xyz):
    return np.hypot(xyz[:, 0], xyz[:, 1]), np.arctan2(xyz[:, 1], xyz[:, 0]), xyz[:, 2]


def test_taper_is_full_at_the_gate_half_midway_and_zero_at_the_edge():
    z = np.array([0.0, W / 2, -W / 2, W, -W, 3 * W])
    np.testing.assert_allclose(widening(z, 0.0, 2.0, W), [2.0, 1.0, 1.0, 0, 0, 0], atol=1e-12)


def test_each_atom_moves_outward_by_the_taper_and_nothing_else():
    st = make_structure(c4_tetramer(n_per=200, seed=3))
    st.xyz[:, 2] = np.linspace(-20, 20, len(st.xyz))       # spread across the window
    wide = widen_gate(st, _summary(0.0), 2.5, W)
    r0, a0, z0 = _cylindrical(st.xyz.astype(float))
    r1, a1, z1 = _cylindrical(wide.xyz.astype(float))
    np.testing.assert_allclose(r1 - r0, widening(z0, 0.0, 2.5, W), atol=1e-4)
    np.testing.assert_allclose(z1, z0, atol=1e-5)
    np.testing.assert_allclose(np.cos(a1 - a0), 1.0, atol=1e-8)       # no turn
    outside = np.abs(z0) >= W
    assert outside.any() and np.array_equal(wide.xyz[outside], st.xyz[outside])


def test_c4_is_kept_and_a_ring_at_the_gate_opens_by_delta():
    blocks = c4_tetramer(n_per=40, seed=1)
    for b in blocks:
        b[:, 2] = np.linspace(-6, 6, len(b))
    wide = widen_gate(make_structure(blocks), _summary(0.0), 3.0, W)
    turn = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
    a, b = wide.xyz[:40].astype(float), wide.xyz[40:80].astype(float)
    np.testing.assert_allclose(a @ turn.T, b, atol=1e-4)
    ring = make_structure([np.array([[3.0 * np.cos(t), 3.0 * np.sin(t), 0.0]])
                           for t in np.pi / 2 * np.arange(4)])
    r, _, _ = _cylindrical(widen_gate(ring, _summary(0.0), 3.0, W).xyz.astype(float))
    np.testing.assert_allclose(r, 6.0, atol=1e-5)


def test_the_scan_runs_from_the_deposit():
    assert deltas()[0] == 0.0 and deltas() == sorted(deltas())


@needs_structure("9HEO")
def test_9heo_gate_opens_by_delta_and_the_filter_stays():
    from ip3r.io import loader
    from ip3r.structure.channel import measure_channel
    st = loader.load("9HEO")
    s = measure_channel(st)
    zg, zf = s.constrictions["gate"].z, s.constrictions["filter"].z
    wide = measure_channel(widen_gate(st, s, 2.0))
    assert wide.profile.at(zg) - s.profile.at(zg) == pytest.approx(2.0, abs=0.3)
    assert abs(wide.profile.at(zf) - s.profile.at(zf)) < 0.05
    assert np.degrees(np.arccos(abs(wide.frame.axis @ s.frame.axis))) < 0.1


@needs_structure("9HEO")
def test_9heo_without_its_gate_is_still_far_from_xu():
    """The finding, pinned at 1 Å: the gate stops holding Ca2+'s resistance,
    yet P_Ca:P_K moves little and the K+ conductance far less than 801 pS
    needs."""
    from ip3r.io import loader
    from ip3r.physics.gate_geometry import gate_scan
    scan = gate_scan(loader.load("9HEO"), widths=[0.0, 4.0],
                     readings=("neutral", "pb + csc"), spacing=1.0)
    shut, wide = scan.rows
    assert shut.gate_share["pb + csc"] > 0.4 and wide.gate_share["pb + csc"] < 0.1
    assert wide.ratios["neutral"] == pytest.approx(shut.ratios["neutral"], rel=0.05)
    assert shut.ratios["pb + csc"] < wide.ratios["pb + csc"] < 2.0 < scan.measured
    assert wide.conductance < 1.5 * shut.conductance
