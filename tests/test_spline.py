"""Rotation-minimising frames: the properties a cartoon depends on."""

import numpy as np

from ip3r.render.spline import parallel_transport_frames


def test_frames_are_orthonormal_on_a_random_walk_with_a_degenerate_stretch():
    rng = np.random.default_rng(3)
    p = np.cumsum(rng.normal(size=(300, 3)), axis=0)
    p[40:44] = p[40]                       # repeated points: zero tangents
    t, n, b = parallel_transport_frames(p)
    assert np.isfinite(n).all() and np.isfinite(b).all()
    ok = np.linalg.norm(t, axis=1) > 0.5
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0)
    assert np.allclose(np.einsum("ij,ij->i", t[ok], n[ok]), 0.0, atol=1e-9)
    assert np.allclose(np.cross(t, n), b)


def test_a_planar_curve_keeps_its_out_of_plane_normal():
    """Parallel transport does not spin about a plane curve: a normal that
    starts perpendicular to the plane stays there."""
    s = np.linspace(0, 3 * np.pi, 500)
    p = np.c_[np.cos(s) * (2 + np.sin(3 * s)), np.sin(s) * 2, np.zeros_like(s)]
    _, n, _ = parallel_transport_frames(p, initial_normal=np.array([0.0, 0.0, 1.0]))
    assert np.allclose(n, [0.0, 0.0, 1.0], atol=1e-9)


def test_a_helix_normal_turns_with_the_curve_not_about_it():
    s = np.linspace(0, 4 * np.pi, 800)
    p = np.c_[np.cos(s), np.sin(s), 0.3 * s]
    t, n, _ = parallel_transport_frames(p)
    # Rotation-minimising: consecutive normals differ only by the tangent's turn.
    turn = np.arccos(np.clip(np.einsum("ij,ij->i", t[1:], t[:-1]), -1, 1))
    step = np.arccos(np.clip(np.einsum("ij,ij->i", n[1:], n[:-1]), -1, 1))
    assert (step <= turn + 1e-9).all()
