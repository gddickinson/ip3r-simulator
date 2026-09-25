"""The voxelised pore and its 3-D ohmic solve (Round 7.6), calibrated on
shapes with a known answer before any deposit is read."""

import numpy as np
import pytest

from ip3r.physics.ohmic3d import (cylinder_volume, geometric_conductance,
                                  hall_cylinder)
from ip3r.structure.pore_volume import open_voxels, volume_from_mask
from conftest import needs_structure

H, BOX, BATH = 0.5, 15.0, 15.0


def _grid(length):
    n = int(round(BOX / H))
    xs = np.arange(-n, n + 1) * H
    zs = np.arange(-length / 2 - BATH, length / 2 + BATH + 1e-9, H)
    x, y, z = np.meshgrid(xs, xs, zs, indexing="ij")
    return xs, zs, np.hypot(x, y), z


def _pores(radius, length, centres):
    """Open everywhere except the membrane, pierced by cylinders at ``centres``."""
    xs, zs, _, z = _grid(length)
    x, y, _ = np.meshgrid(xs, xs, zs, indexing="ij")
    hole = np.zeros(z.shape, bool)
    for cx, cy in centres:
        hole |= np.hypot(x - cx, y - cy) <= radius
    acc = (np.abs(z) > length / 2) | hole
    return volume_from_mask(acc, xs, zs, (-length / 2, length / 2), 0.0,
                            seal=BOX * 2)


def test_cylinder_matches_hall():
    """Series resistance plus Hall's two access terms, to the voxel error
    (a smaller box adds a little access: the ratio sits just below 1)."""
    g = geometric_conductance(cylinder_volume(4.0, 20.0, H, BOX, BATH))
    assert g.converged
    assert 0.90 < g.g / hall_cylinder(4.0, 20.0) < 1.02


def test_grid_converges_on_cylinder():
    exact = hall_cylinder(4.0, 20.0)
    err = [abs(geometric_conductance(cylinder_volume(4.0, 20.0, h, 12.0, 12.0)).g
               / exact - 1) for h in (1.0, 0.5)]
    assert err[1] < err[0]


def test_two_pores_conduct_twice():
    one = geometric_conductance(_pores(3.0, 20.0, [(0.0, 0.0)])).g
    two = geometric_conductance(_pores(3.0, 20.0, [(-7.0, 0.0), (7.0, 0.0)])).g
    assert two / one == pytest.approx(2.0, rel=0.05)


def test_blind_hole_is_dropped():
    """A pocket reached from one bath only is a dead end: it stays joined to
    the bath but carries no current, so the conductance does not move."""
    length = 20.0
    xs, zs, rho, z = _grid(length)
    x = np.meshgrid(xs, xs, zs, indexing="ij")[0]
    pore = rho <= 3.0
    blind = (np.hypot(x - 8.0, np.meshgrid(xs, xs, zs, indexing="ij")[1]) <= 2.0) \
        & (z > 0)
    base = (np.abs(z) > length / 2) | pore
    plain = volume_from_mask(base, xs, zs, (-10, 10), 0.0, seal=BOX * 2)
    holed = volume_from_mask(base | blind, xs, zs, (-10, 10), 0.0, seal=BOX * 2)
    assert holed.accessible > plain.accessible
    assert geometric_conductance(holed).g == pytest.approx(
        geometric_conductance(plain).g, rel=1e-9)


def test_sideways_exit_conducts():
    """A cap on the axis just past the membrane: shut to a profile along the
    axis (radius 0 there), open in 3-D because ions go round it."""
    length = 20.0
    xs, zs, rho, z = _grid(length)
    base = (np.abs(z) > length / 2) | (rho <= 3.0)
    cap = (np.abs(z - (length / 2 + 3.0)) <= 1.5) & (rho <= 7.0)
    open_ = geometric_conductance(volume_from_mask(base, xs, zs, (-10, 10), 0.0,
                                                   seal=BOX * 2)).g
    capped = volume_from_mask(base & ~cap, xs, zs, (-10, 10), 0.0, seal=BOX * 2)
    g = geometric_conductance(capped).g
    assert np.all(rho[:, :, np.argmin(abs(zs - 13.0))][capped.mask[:, :, np.argmin(
        abs(zs - 13.0))]] > 7.0)
    assert 0.5 < g / open_ < 1.0


def test_membrane_seal_blocks_the_lipid_space():
    """With nothing but bath and membrane, the seal is all that stops
    current: a seal of 3 Å is a 3 Å pore, and one of 0 Å leaves the single
    axis column of voxels."""
    length = 20.0
    xs, zs, rho, z = _grid(length)
    acc = np.ones(z.shape, bool)
    column = geometric_conductance(volume_from_mask(acc, xs, zs, (-10, 10), 0.0,
                                                    seal=0.0)).g
    sealed = geometric_conductance(volume_from_mask(acc, xs, zs, (-10, 10), 0.0,
                                                    seal=3.0)).g
    drawn = geometric_conductance(_pores(3.0, length, [(0.0, 0.0)])).g
    assert sealed == pytest.approx(drawn, rel=1e-9)
    assert column < 0.05 * sealed


def test_open_voxels_hard_spheres():
    """One atom: a voxel is open exactly when it clears vdW + probe."""
    xs = np.arange(-4, 4.01, 0.5)
    zs = np.array([0.0])
    ok = open_voxels(np.zeros((1, 3)), np.array([1.5]), xs, zs, 1.0)
    x, y = np.meshgrid(xs, xs, indexing="ij")
    assert np.array_equal(ok[:, :, 0], np.hypot(x, y) >= 2.5)


def test_open_states_are_the_three():
    from ip3r.physics.shortfall import open_entries
    assert sorted(e.pdb_id for e in open_entries()) == ["7T3T", "8TKF", "9HEO"]


@needs_structure("8TKF", "7T3T")
def test_deposits_on_a_coarse_grid():
    """8TKF at 1 Å: the lumen conducts more than its inscribed circle, the
    seal is inside the protein (20 = 25 Å) and a seal past it leaks; the
    second open ITPR3 deposit reads within a third of the first."""
    from ip3r.io.loader import load
    from ip3r.physics.ohmic3d import conductance_3d
    from ip3r.physics.unitary import unitary
    from ip3r.structure.channel import measure_channel
    st = load("8TKF")
    s = measure_channel(st)
    g = conductance_3d(st, s, spacing=1.0)
    assert g.converged
    assert g.conductance_pS > 1.2 * unitary(st, s).neutral.conductance_pS
    assert conductance_3d(st, s, spacing=1.0, seal=20.0).conductance_pS == \
        pytest.approx(g.conductance_pS, rel=0.01)
    assert conductance_3d(st, s, spacing=1.0, seal=35.0).conductance_pS > \
        5 * g.conductance_pS
    other = load("7T3T")
    g2 = conductance_3d(other, measure_channel(other), spacing=1.0)
    assert g2.conductance_pS == pytest.approx(g.conductance_pS, rel=0.33)


@needs_structure("8TKF")
def test_window_does_not_matter():
    from ip3r.physics.shortfall import window_scan
    rows = window_scan("8TKF")
    g = [r[1] for r in rows]
    assert max(g) / min(g) < 1.05
    assert all(r[2] < 0.05 for r in rows)
