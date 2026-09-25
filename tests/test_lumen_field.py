"""Where the voltage falls (Round 7.10): the field and its drawn surface
calibrated on shapes with a known answer, then the real open deposits."""

import numpy as np
import pytest

from ip3r.physics.lumen_field import field_from_volume
from ip3r.physics.ohmic3d import cylinder_volume
from ip3r.structure.pore_volume import volume_from_mask

H, L, MARGIN, HALF = 0.5, 30.0, 15.0, 10.0


class _Frame:
    """Identity frame: frame coordinates are lab coordinates."""
    basis = np.eye(3)
    centre = np.zeros(3)

    def from_frame(self, xyz):
        return np.asarray(xyz, float)


def _necked(radius, neck_radius, neck_lo, neck_hi):
    """A cylinder through a slab, narrowed between ``neck_lo`` and ``neck_hi``."""
    n = int(round(HALF / H))
    xs = np.arange(-n, n + 1) * H
    zs = np.arange(-L / 2 - MARGIN, L / 2 + MARGIN + 1e-9, H)
    x, y, z = np.meshgrid(xs, xs, zs, indexing="ij")
    rho = np.hypot(x, y)
    r = np.where((z >= neck_lo) & (z <= neck_hi), neck_radius, radius)
    acc = (np.abs(z) > L / 2) | (rho <= r)
    return volume_from_mask(acc, xs, zs, (-L / 2, L / 2), 0.0, seal=radius)


def _profile(radius, neck_radius=None, lo=0.0, hi=0.0):
    z = np.arange(-L / 2 - MARGIN, L / 2 + MARGIN + 1e-9, H / 4)
    r = np.full(len(z), float(radius))
    if neck_radius is not None:
        r[(z >= lo) & (z <= hi)] = neck_radius
    return z, r


def test_uniform_cylinder_drops_linearly():
    vol = cylinder_volume(4.0, L, H, HALF, MARGIN)
    f = field_from_volume(vol, (-L / 2, L / 2), r_free=_profile(4.0))
    assert f.conducts and f.converged
    line = (f.z - f.z[0]) / (f.z[-1] - f.z[0])
    # Inside the slab the drop is linear; the ends bend within a radius of
    # the mouths (current converging from the bath), so check the middle.
    mid = np.abs(f.z) < L / 2 - 4.0
    assert np.max(np.abs(f.drop_3d - line)[mid]) < 0.03
    assert np.max(np.abs(f.drop_1d - line)) < 0.02
    assert abs(f.half_z("3d")) < 0.5 and abs(f.half_z("1d")) < 0.5
    # Every plane of the pore is one voxel disc of radius 4.
    disc = f.area_3d[mid]
    assert np.ptp(disc) == 0 and disc[0] == pytest.approx(np.pi * 16, rel=0.08)
    # The access resistance takes a share: 2 × 1/(4a) against L/(πa²).
    access = 2 / (4 * 4.0) / (L / (np.pi * 16) + 2 / (4 * 4.0))
    assert 1 - f.window_share == pytest.approx(access, abs=0.04)


def test_the_neck_holds_the_drop_on_both_routes():
    vol = _necked(5.0, 2.0, 2.0, 8.0)
    f = field_from_volume(vol, (-L / 2, L / 2), r_free=_profile(5.0, 2.0, 2.0, 8.0))
    d3, d1 = f.drop_across(5.0, 4.0)
    # 1-D by hand: the neck is 6 Å of radius 2 against 24 Å of radius 5,
    # read over the neck and 1 Å of wide pore either side.
    by_hand = (6 / 4 + 2 / 25) / (6 / 4 + 24 / 25)
    assert d1 == pytest.approx(by_hand, abs=0.02)
    assert d3 > 0.5 and abs(d3 - d1) < 0.1
    # An off-centre neck pulls the half-drop point towards it.
    assert 2.0 < f.half_z("3d") < 8.0 and 2.0 < f.half_z("1d") < 8.0


def test_a_blocked_pore_carries_no_field():
    vol = _necked(5.0, -1.0, -1.0, 1.0)      # radius 0 would keep the axis column
    f = field_from_volume(vol, (-L / 2, L / 2), r_free=_profile(5.0, 0.0, -1.0, 1.0))
    assert not f.conducts
    assert np.all(np.isnan(f.drop_3d)) and f.shut_1d
    assert "no" in f.summary() and "path" in f.summary()


def test_the_1d_route_needs_a_profile():
    vol = cylinder_volume(4.0, L, H, HALF, MARGIN)
    f = field_from_volume(vol, (-L / 2, L / 2))
    assert np.all(np.isnan(f.drop_1d)) and np.isfinite(f.half_z("3d"))


def test_the_surface_wraps_the_lumen():
    from ip3r.render.lumen_mesh import drawn_mask, lumen_mesh
    vol = cylinder_volume(4.0, L, H, HALF, MARGIN)
    f = field_from_volume(vol, (-L / 2, L / 2), r_free=_profile(4.0))
    mesh = lumen_mesh(f, _Frame(), radius=6.0, smoothing=0.0)
    p = mesh.positions
    rho = np.hypot(p[:, 0], p[:, 1])
    wall = np.abs(p[:, 2]) < L / 2 - 2
    # The unsmoothed surface sits half a voxel outside the last open voxel.
    assert np.all(np.abs(rho[wall] - 4.0) < 1.5 * H)
    assert np.all((p[:, 2] >= -L / 2 - H) & (p[:, 2] <= L / 2 + H))
    radial = p[wall, :2] / rho[wall, None]
    assert np.mean(np.sum(mesh.normals[wall, :2] * radial, axis=1)) > 0.8   # staircase facets are diagonal
    # φ follows the drop: the cytosolic end (+z) red, the luminal blue.
    assert np.corrcoef(p[wall, 2], mesh.phi[wall])[0, 1] > 0.99
    assert np.all(np.isfinite(mesh.phi))
    assert drawn_mask(f, radius=6.0).sum() < vol.mask.sum()


def test_nothing_to_draw_when_shut():
    from ip3r.render.lumen_mesh import lumen_mesh
    vol = _necked(5.0, -1.0, -1.0, 1.0)
    f = field_from_volume(vol, (-L / 2, L / 2))
    assert lumen_mesh(f, _Frame()) is None


@pytest.fixture(scope="module")
def real():
    from ip3r.io.loader import StructureUnavailable, load
    from ip3r.physics.lumen_field import lumen_field
    from ip3r.structure.channel import measure_channel
    out = {}
    for pdb in ("8TKF", "8TKG"):
        try:
            st = load(pdb)
        except StructureUnavailable:
            pytest.skip(f"{pdb} not fetched")
        s = measure_channel(st)
        out[pdb] = (s, lumen_field(st, s))
    return out


def test_8tkf_drops_where_the_1d_model_says(real):
    s, f = real["8TKF"]
    assert f.conducts and f.window_share > 0.95
    assert abs(f.half_z("3d") - f.half_z("1d")) < 2.0
    filt = s.constrictions["filter"]
    d3, d1 = f.drop_across(filt.z, 3.0)
    assert abs(d3 - d1) < 0.1 and d3 > 0.25


def test_the_resting_state_is_shut(real):
    _, f = real["8TKG"]
    assert not f.conducts
