"""The axis finders and the C4 residual, on tetramers built to be C4."""

import numpy as np

from ip3r.structure.symmetry import (axis_by_centroids, axis_by_superposition,
                                     c4_residual, kabsch, rotation_axis_angle,
                                     rotation_matrix, subunit_ca, tetramer_frame)
from conftest import c4_tetramer, make_structure


def test_kabsch_recovers_a_rotation():
    rng = np.random.default_rng(0)
    p = rng.normal(size=(30, 3))
    r0 = rotation_matrix([1, 2, 3], 0.7)
    r, t = kabsch(p, p @ r0.T + [1, 2, 3])
    assert np.allclose(r, r0, atol=1e-9) and np.allclose(t, [1, 2, 3])
    axis, angle = rotation_axis_angle(r)
    assert abs(angle - 0.7) < 1e-9
    assert abs(abs(axis @ np.array([1, 2, 3]) / np.sqrt(14)) - 1) < 1e-9


def test_exact_c4_tetramer():
    st = make_structure(c4_tetramer())
    ca = subunit_ca(st)
    axis, _, angle = axis_by_superposition(ca)
    assert abs(angle - 90.0) < 1e-6
    assert abs(abs(axis[2]) - 1.0) < 1e-9
    ax2, _, planar = axis_by_centroids(ca)
    assert abs(abs(ax2[2]) - 1.0) < 1e-9 and planar < 1e-9
    assert c4_residual(st, tetramer_frame(st)) < 1e-6


def test_residual_can_fail():
    """Break the symmetry and the residual must say so."""
    st = make_structure(c4_tetramer(noise=1.0))
    assert c4_residual(st, tetramer_frame(st)) > 0.5


def test_chains_ordered_right_handed():
    st = make_structure(c4_tetramer(), chains="ADBC")      # scrambled letters
    fr = tetramer_frame(st)
    cents = [fr.to_frame(np.mean(list(subunit_ca(st)[c].values()), 0)) for c in fr.chains]
    az = [np.degrees(np.arctan2(c[1], c[0])) % 360 for c in cents]
    steps = [(az[(i + 1) % 4] - az[i]) % 360 for i in range(4)]
    assert all(abs(s - 90) < 1e-6 for s in steps)
