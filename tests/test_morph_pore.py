"""The gate along a morph, with side chains interpolated.

The endpoint frames must be measured exactly as the deposits are (the same
rule, the same numbers), and the shortcut the viewer used before — side
chains riding their C-alpha — must be caught missing the end deposit's gate.
"""

import numpy as np
import pytest

from conftest import needs_structure
from ip3r.io.loader import load
from ip3r.structure.channel import measure_channel
from ip3r.structure.morph import morph
from ip3r.structure.morph_pore import AtomPath, GatePath, atom_path, gate_path
from ip3r.structure.transition import prepare_transition


def _path(offset_start, offset_end):
    k = len(offset_start)
    return AtomPath(np.zeros(k, int), np.asarray(offset_start, float),
                    np.asarray(offset_end, float), np.ones(k), np.array(["X1"] * k),
                    np.asarray(offset_start, float), np.arange(k))


def test_offsets_exact_at_both_ends():
    p = _path([[1.0, 0, 0]], [[0, 2.0, 0]])
    site = np.array([[5.0, 5, 5]])
    assert np.allclose(p.coords(site, 0.0), [[6, 5, 5]])
    assert np.allclose(p.coords(site, 1.0), [[5, 7, 5]])


def test_chord_of_a_flipped_side_chain():
    """A side chain swung 180° passes through its own C-alpha mid-path:
    the chord is its full length, and zero at both ends."""
    p = _path([[3.0, 0, 0]], [[-3.0, 0, 0]])
    assert p.offset_error(0.5) == pytest.approx(3.0)
    assert p.offset_error(0.0) == p.offset_error(1.0) == 0.0
    assert _path([[3.0, 0, 0]], [[3.0, 0, 0]]).offset_error(0.5) == 0.0


def _gate(values):
    n = len(values)
    z = np.zeros(n)
    return GatePath(np.linspace(0, 1, n), np.asarray(values, float), z, z, z,
                    [[]] * n, z, z)


def test_half_open_and_overshoot():
    assert _gate([2.0, 3.0, 4.0, 5.0, 6.0]).half_open() == pytest.approx(0.5)
    assert _gate([2.0, 2.0, 2.0, 2.0, 6.0]).half_open() == pytest.approx(0.875)
    assert _gate([6.0, 4.0, 2.0]).half_open() == pytest.approx(0.5)   # closing
    assert _gate([2.0, 7.0, 6.0]).overshoot() == pytest.approx(1.0)
    assert _gate([2.0, 4.0, 6.0]).overshoot() == 0.0


def _measured(a, b):
    sa, sb = load(a), load(b)
    tr = prepare_transition(sa, sb)
    path = atom_path(sa, sb, tr)
    return sa, sb, tr, path, gate_path(path, tr, morph(tr.start, tr.end))


@pytest.fixture(scope="module")
def itpr3():
    return _measured("8TKG", "8TKF")


@needs_structure("8TKG", "8TKF")
def test_last_frame_is_the_end_deposit_atom_for_atom(itpr3):
    sa, sb, tr, path, _ = itpr3
    assert path.meta["unmatched"] == 0 and path.meta["n_atoms"] > 60000
    r, t = tr.meta["end_transform"]
    last = path.coords(tr.end, 1.0)
    # the end deposit read by residue and name, independently of atom_path
    k = 0
    ch_a, res, name = sa.chain[path.index[k]], sa.res_seq[path.index[k]], sa.atom_name[path.index[k]]
    ch_b = tr.chains_end[tr.chains_start.index(ch_a)]
    j = np.flatnonzero((sb.chain == ch_b) & (sb.res_seq == res) & (sb.atom_name == name))[0]
    assert np.allclose(last[k], sb.xyz[j] @ r.T + t, atol=1e-4)
    assert np.allclose(path.coords(tr.start, 0.0), sa.xyz[path.index], atol=1e-4)


@pytest.mark.parametrize("pair", [
    pytest.param(("8TKG", "8TKF"), marks=needs_structure("8TKG", "8TKF")),
    pytest.param(("9R8O", "9HEO"), marks=needs_structure("9R8O", "9HEO"))])
def test_endpoint_frames_measure_as_the_deposits(pair):
    g = _measured(*pair)[4]
    for i, pdb in ((0, pair[0]), (-1, pair[1])):
        own = measure_channel(load(pdb), include_hetero=False).constrictions
        assert g.gate[i] == pytest.approx(own["gate"].radius, abs=0.05)
        assert g.filter[i] == pytest.approx(own["filter"].radius, abs=0.05)
        assert g.lining[i] == own["gate"].residues
    assert g.offset_error[0] == g.offset_error[-1] == 0.0


@needs_structure("8TKG", "8TKF")
def test_rigid_side_chains_miss_the_end_gate(itpr3):
    """The case that must fail: 8TKF's backbone lined by 8TKG's rotamers."""
    g = itpr3[4]
    assert g.rigid_gate[0] == pytest.approx(g.gate[0])
    assert abs(g.rigid_gate[-1] - g.gate[-1]) > 0.5
    assert g.gate[0] < 3.0 < 5.5 < g.gate[-1]
    assert g.overshoot() < 0.05 and 0.3 < g.half_open() < 0.6
