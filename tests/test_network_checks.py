"""Local modes located, and the network measured against its own choices.

Synthetic: a C4 tetramer of compact blocks with one planted flap (two sites
after an unresolved stretch, held by a few long springs) must produce local
modes on the flap and name the stretch; without the flap there must be none.
Real (skipped when not downloaded): 8TKG's stride-3 artefacts are the flap at
residue 86 after unresolved 77-85; the lowest-A-mode overlap moves with the
cutoff while the collective A modes together hold.
"""

import numpy as np
import pytest

from conftest import needs_structure
from ip3r.physics.anm import ANM
from ip3r.physics.network_checks import (_gap_beside, cutoff_grid, cutoff_scan,
                                         describe_local, local_modes, rmsip,
                                         stride_agreement)
from ip3r.structure.symmetry import rotation_matrix


def _tetramer(flap: bool):
    # Each subunit fills one quadrant, so the four make one solid slab.
    g = np.arange(5) * 3.8
    block = np.array([[x, y, z] for x in g for y in g for z in g[:4]]) + [3.0, 3.0, 0.0]
    residues = list(range(1, 101))
    if flap:
        # Two sites past the block's outer x face, numbered after a missing 101-110.
        block = np.vstack([block, [[3 + 15.2 + 9.0, 10.6, 5.7],
                                   [3 + 15.2 + 12.5, 10.6, 5.7]]])
        residues += [111, 112]
    z = np.array([0.0, 0.0, 1.0])
    coords = np.vstack([block @ rotation_matrix(z, k * np.pi / 2).T for k in range(4)])
    return coords, np.array(residues)


def _modes(coords):
    anm = ANM(coords, axis=np.array([0.0, 0.0, 1.0]), cutoff=15.0)
    return anm.label_symmetry(anm.calc_modes(20))


def test_planted_flap_is_located_with_its_gap():
    coords, res = _tetramer(flap=True)
    local = local_modes(_modes(coords), res, resolved=res)
    assert len(local) >= 4                          # one flap, four subunits
    assert {m.residue for m in local} <= {111, 112}
    # The gap is looked for only between the residue's sampled neighbours:
    # beside the flap's first site, not its tip.
    assert all(m.gap == ((101, 110) if m.residue == 111 else None) for m in local)
    assert _gap_beside(111, res, res) == (101, 110)
    assert {m.irrep for m in local[:4]} >= {"A", "B", "E"}
    assert all(line.split(" on residue ")[1].split(",")[0] in ("111", "112")
               for line in describe_local(local))


def test_no_flap_no_local_modes():
    coords, res = _tetramer(flap=False)
    assert local_modes(_modes(coords), res, resolved=res) == []


def test_gap_beside():
    sampled = np.array([70, 74, 86, 89])
    resolved = np.r_[70:77, 86:93]
    assert _gap_beside(86, sampled, resolved) == (77, 85)
    assert _gap_beside(89, sampled, resolved) is None
    assert _gap_beside(86, sampled, None) is None


def test_rmsip():
    e = np.eye(6)
    assert rmsip(e, e, 3) == pytest.approx(1.0)
    assert rmsip(e[:3], e[3:], 3) == pytest.approx(0.0)
    # The same subspace in another basis still scores 1.
    rot = np.linalg.qr(np.random.default_rng(0).normal(size=(3, 3)))[0]
    assert rmsip(e[:3], rot @ e[:3], 3) == pytest.approx(1.0)


def test_cutoff_grid_is_registered():
    grid = cutoff_grid()
    assert grid[0] == 12.0 and grid[-1] == 21.0 and len(grid) == 7


@pytest.fixture(scope="module")
def tr():
    from ip3r.io.loader import load
    from ip3r.structure.transition import prepare_transition
    return prepare_transition(load("8TKG"), load("8TKF"))


@needs_structure("8TKG", "8TKF")
def test_8tkg_local_modes_are_the_residue_86_flap(tr):
    from ip3r.physics.transition_modes import transition_overlap
    ov = transition_overlap(tr, stride=3)
    local = local_modes(ov.modes, ov.residues, tr.residues)
    assert [m.index + 1 for m in local] == [11, 12, 13, 14, 15]
    assert {(m.residue, m.gap) for m in local} == {(86, (77, 85))}
    assert local_modes(transition_overlap(tr, stride=2).modes,
                       tr.residues[::2], tr.residues) == []


@needs_structure("8TKG", "8TKF")
def test_single_mode_moves_with_cutoff_subspace_holds(tr):
    lo, hi = cutoff_scan(tr, cutoffs=(15.0, 21.0), stride=3)
    assert lo.lowest_a_overlap == pytest.approx(0.415, abs=0.005)   # the logged value
    assert lo.lowest_a_overlap - hi.lowest_a_overlap > 0.15
    assert abs(hi.collective_a - lo.collective_a) < 0.1
    assert hi.n_local == 0 and lo.n_local == 5


@needs_structure("8TKG", "8TKF")
def test_stride_two_agrees_with_full_network(tr):
    agree = stride_agreement(tr, strides=(2, 4))
    assert agree[2] > 0.95 and agree[4] < 0.85
