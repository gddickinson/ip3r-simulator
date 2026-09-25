"""Morph, transition basis and mode overlap — synthetic calibrations first.

Every measuring function here is shown to give the known answer on a field
built to have it (a pure A, B or E field; a displacement that *is* a mode;
a rigid rotation) before it is trusted on 8TKG -> 8TKF.
"""

import numpy as np
import pytest

from conftest import c4_tetramer
from ip3r.physics.anm import ANM, ModeSet, apply_generator
from ip3r.physics.transition_modes import (irrep_fractions, null_cumulative,
                                           remove_rigid_body, transition_overlap)
from ip3r.structure.morph import morph, peptide_pairs
from ip3r.structure.symmetry import Frame, rotation_matrix
from ip3r.structure.transition import Transition

Z = np.array([0.0, 0.0, 1.0])


def _helix(n=40):
    t = np.arange(n) * 100 * np.pi / 180
    return np.column_stack([12 + 2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)])


def _tetramer(block):
    return np.vstack([block @ rotation_matrix(Z, np.pi / 2 * k).T for k in range(4)])


# ------------------------------------------------------------------- morph

def test_restrained_removes_the_chord_artefact():
    start = _tetramer(_helix())
    end = start @ rotation_matrix(Z, np.radians(60)).T     # a large rigid swing
    lin, res = morph(start, end, "linear", 21), morph(start, end, "restrained", 21)
    assert lin.bond_error.max() > 0.3                      # chords shorten bonds (0.43 A)
    assert res.bond_error.max() < 0.05 * lin.bond_error.max()
    for tr in (lin, res):                                  # endpoints are experimental
        assert np.allclose(tr.frames[0], start) and np.allclose(tr.frames[-1], end)
    assert tr.nearest(1.0) == 20 and tr.nearest(0.49) == 10


def test_chain_breaks_are_not_restrained():
    block = _helix()
    block[20:] += [0, 0, 30.0]                             # a gap mid-chain
    x = _tetramer(block)
    pairs = peptide_pairs(x, x)
    assert len(pairs) == 4 * (len(block) - 2)
    with pytest.raises(ValueError):
        morph(x, x[:-1])


# ---------------------------------------------------------------- symmetry

def _fields():
    x = _tetramer(_helix())
    rng = np.random.default_rng(1)
    u = rng.normal(size=x.shape)
    powers = [apply_generator(u, Z, 4, k) for k in range(4)]
    a = sum(powers) / 4
    b = sum((-1) ** k * p for k, p in enumerate(powers)) / 4
    return x, a, b, u - a - b


def test_irrep_fractions_calibrated_on_pure_fields():
    _, a, b, e = _fields()
    for field, name in ((a, "A"), (b, "B"), (e, "E")):
        f = irrep_fractions(field, Z)
        assert f[name] == pytest.approx(1.0, abs=1e-9), (name, f)


def test_rigid_body_removed_exactly():
    x = _tetramer(_helix())
    rot = x @ rotation_matrix([1, 2, 3], 1e-3).T - x + [0.3, -0.2, 0.1]
    assert np.linalg.norm(remove_rigid_body(rot, x)) < 1e-3 * np.linalg.norm(rot)
    # a projector: applying it twice changes nothing
    _, a, _, _ = _fields()
    kept = remove_rigid_body(a, x)
    assert np.allclose(remove_rigid_body(kept, x), kept)
    # a purely internal move (radial breathing about the axis) is kept whole
    breathe = x * [1.0, 1.0, 0.0]
    kept = remove_rigid_body(breathe, x)
    assert np.linalg.norm(kept - breathe) < 1e-9 * np.linalg.norm(breathe)


def test_collectivity_bounds():
    n = 50
    even = np.ones((1, n, 3))
    one = np.zeros((1, n, 3))
    one[0, 7] = 1.0
    ms = ModeSet(np.ones(2), np.concatenate([even, one]))
    k = ms.collectivity()
    assert k[0] == pytest.approx(1.0) and k[1] == pytest.approx(1 / n)


def test_null_matches_irrep_dimension():
    sym = np.array(["A"] * 5)
    n_sites = 400
    got = null_cumulative(sym, {"A": 1.0, "B": 0.0, "E": 0.0}, n_sites)
    assert got[-1] == pytest.approx(np.sqrt(5 / (3 * n_sites / 4 - 2)))


# ------------------------------------------------------------ overlap

def _synthetic_transition(disp_fn):
    start = np.vstack(c4_tetramer(n_per=60, seed=3))
    per = len(start) // 4
    residues = np.arange(1, per + 1)
    frame = Frame(Z, np.zeros(3), np.eye(3), list("ABCD"))
    end = start + disp_fn(start)
    return Transition("S", "E", "ITPR3", residues, list("ABCD"), list("ABCD"),
                      start, end, "global", 0.0, 0.0, frame)


def test_overlap_is_one_when_the_move_is_a_mode():
    x = np.vstack(c4_tetramer(n_per=60, seed=3))
    anm = ANM(x, axis=Z)
    ms = anm.label_symmetry(anm.calc_modes(12))
    a = ms.first("A", collective=False)
    tr = _synthetic_transition(lambda s: 2.0 * ms.mode(a, 1.0))
    ov = transition_overlap(tr, stride=1, n_modes=12)
    assert ov.overlap[a] > 0.999
    assert ov.irrep_fraction["A"] > 0.999
    assert ov.cumulative[-1] > 0.999


def test_overlap_of_a_random_move_sits_at_the_null():
    rng = np.random.default_rng(5)
    tr = _synthetic_transition(lambda s: 0.5 * rng.normal(size=s.shape))
    ov = transition_overlap(tr, stride=1, n_modes=12)
    assert ov.cumulative[-1] < 4 * ov.null_cumulative[-1]


def test_a_subspace_holds_a_move_spread_over_a_modes():
    """Two A modes mixed: no single mode explains it, the A subspace does
    completely, and a B move leaves the A subspace empty (ceiling 0)."""
    x = np.vstack(c4_tetramer(n_per=60, seed=3))
    anm = ANM(x, axis=Z)
    ms = anm.label_symmetry(anm.calc_modes(12))
    a = np.flatnonzero(ms.symmetry == "A")[:2]
    b = ms.first("B", collective=False)
    unit = [ms.vectors[i] / np.linalg.norm(ms.vectors[i]) for i in a]
    mixed = _synthetic_transition(lambda s: 20.0 * (unit[0] + unit[1]))
    sub = transition_overlap(mixed, stride=1, n_modes=12).subspace("A", collective=False)
    ov = transition_overlap(mixed, stride=1, n_modes=12)
    assert ov.overlap.max() < 0.75            # each carries 1/sqrt(2)
    assert sub.total > 0.999 and sub.ceiling == pytest.approx(1.0, abs=1e-6)
    assert np.all(np.diff(sub.cumulative) >= 0) and sub.null[-1] < 0.2
    moved_b = _synthetic_transition(lambda s: 2.0 * ms.mode(b, 1.0))
    sub_b = transition_overlap(moved_b, stride=1, n_modes=12).subspace("A", collective=False)
    assert sub_b.total < 0.01 and sub_b.ceiling < 0.01
