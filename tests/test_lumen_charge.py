"""The charged lumen (Round 7.12): the readings the viewer colours by,
calibrated where the answer is known, then held to Round 7.11's numbers."""

import numpy as np
import pytest

from ip3r.physics.lumen_charge import ChargedLumen
from ip3r.physics.lumen_field import field_from_volume
from ip3r.physics.ohmic3d import cylinder_volume
from ip3r.render.colormaps import MISSING, ramp
from ip3r.render.lumen_mesh import lumen_mesh, wall_colors

H, L, MARGIN, HALF = 0.5, 30.0, 15.0, 10.0


class _Frame:
    basis = np.eye(3)

    def from_frame(self, xyz):
        return np.asarray(xyz, float)


def _cylinder():
    vol = cylinder_volume(4.0, L, H, HALF, MARGIN)
    z = np.arange(-L / 2 - MARGIN, L / 2 + MARGIN + 1e-9, H / 4)
    return field_from_volume(vol, (-L / 2, L / 2), r_free=(z, np.full(len(z), 4.0)))


def _charged(f, u_1d):
    nan = np.full(f.volume.mask.shape, np.nan)
    return ChargedLumen(neutral=f, closure="slice", pair_bridges=False,
                        wall=None, charge=None, u=nan, mu=nan, g=f.g,
                        u_3d=np.zeros(len(f.z)), u_1d=u_1d, mu_3d=f.phi_3d,
                        area_1d=f.area_1d)


def test_the_1d_charged_drop_reduces_to_the_neutral_one():
    """No charge, or a uniform one (a constant factor on every slice),
    leaves ∫dz/(A e^{-u}) normalised exactly as ∫dz/A."""
    f = _cylinder()
    for u in (np.zeros(len(f.z)), np.full(len(f.z), -2.0)):
        c = _charged(f, u)
        assert np.allclose(c.drop_1d, f.drop_1d)
        assert np.allclose(c.drop_3d, f.drop_3d)
    assert _charged(f, np.zeros(len(f.z))).ratio == 1.0


def test_a_cation_well_pulls_the_1d_drop_out_of_it():
    """A well (u < 0) raises K+ there, so less of the drop falls in it:
    by hand, a well of depth ln 4 over the luminal half puts 1/5 of the
    drop there instead of 1/2."""
    f = _cylinder()
    u = np.where(f.z < 0, -np.log(4.0), 0.0)
    c = _charged(f, u)
    assert np.interp(0.0, f.z, c.drop_1d) == pytest.approx(0.2, abs=0.02)
    assert c.half_z("1d") > 3.0 > f.half_z("1d")


def test_the_steepest_point_and_the_share_read_the_3d_drop():
    """A drop falling in one step is steepest there and holds it all."""
    f = _cylinder()
    c = _charged(f, np.zeros(len(f.z)))
    c.mu_3d = np.where(f.z < 2.0, 0.0, 1.0)
    assert abs(c.steepest_z() - 2.0) <= 0.5
    assert c.drop_across(2.0, 1.0) == pytest.approx(1.0, abs=0.01)
    assert c.drop_across(-8.0, 1.0) < 0.01


def test_the_wall_ramp_is_fixed_and_grey_where_missing():
    span = 5.0
    c = wall_colors(np.array([-span, 0.0, span, np.nan, -3 * span]), span=span)
    assert np.allclose(c[0], ramp(np.array([0.0]))[0])       # well: blue end
    assert np.allclose(c[1], ramp(np.array([0.5]))[0])       # zero: pale
    assert np.allclose(c[2], ramp(np.array([1.0]))[0])       # repulsive: red
    assert np.allclose(c[3], MISSING)
    assert np.allclose(c[4], c[0])                           # saturates, no rescale


def test_a_mesh_samples_its_own_voxels():
    """Recolouring reads the voxel each vertex was coloured from."""
    f = _cylinder()
    mesh = lumen_mesh(f, _Frame())
    assert np.array_equal(mesh.sample(f.phi), mesh.phi)
    assert np.all(np.isfinite(mesh.sample(f.phi)))


# ---------------------------------------------------------------- real data

@pytest.fixture(scope="module")
def tkf():
    from ip3r.io import loader
    from ip3r.structure.channel import measure_channel
    try:
        st = loader.load("8TKF")
    except loader.StructureUnavailable:
        pytest.skip("8TKF not fetched")
    return st, measure_channel(st)


@pytest.fixture(scope="module")
def tkf_field(tkf):
    from ip3r.physics.lumen_field import lumen_field
    return lumen_field(*tkf)


def test_8tkf_matches_round_7_11s_k_reading(tkf, tkf_field):
    """The drawn lumen is wall_3d's electrostatic volume: the same K+ g."""
    from ip3r.physics.charged3d import wall_3d
    from ip3r.physics.lumen_charge import charged_lumen
    st, s = tkf
    c = charged_lumen(st, tkf_field, "slice", s)
    w = wall_3d(st, s, closures=("slice",))
    assert c.converged
    assert c.g == pytest.approx(w.per_species["slice"]["K+"], rel=1e-6)
    assert tkf_field.g == pytest.approx(w.per_species["neutral"]["K+"], rel=1e-6)
    # Round 7.11's full-wall reading raises K+ (the co-ion Cl- pays).
    assert c.ratio > 1.0
    # The drawn potential is the solved one, on the lumen only.
    inside = tkf_field.volume.mask
    assert np.all(np.isfinite(c.u[inside])) and np.all(np.isnan(c.u[~inside]))
    assert np.all(np.isfinite(c.drop_3d))
    assert c.drop_3d[0] == 0.0 and c.drop_3d[-1] == pytest.approx(1.0)


def test_8tkf_pairing_turns_the_k_ratio_down(tkf, tkf_field):
    """Round 7.11: with D2478's salt bridges paired, PB lowers g."""
    from ip3r.physics.lumen_charge import charged_lumen
    st, s = tkf
    full = charged_lumen(st, tkf_field, "pb", s)
    paired = charged_lumen(st, tkf_field, "pb", s, pair_bridges=True)
    assert full.ratio > 2.0 and paired.ratio < 1.0
    u, _ = full.well()
    assert -20 < u < 0
    with pytest.raises(ValueError):
        charged_lumen(st, tkf_field, "donnan", s)



def test_8tkf_dielectric_is_round_7_13s_dipole_and_pair_omitted(tkf):
    """Round 7.14: the lumen box's dielectric closure is Round 7.13's field
    on the same volume (1 Å grid, as its tests): unpaired, the K+ g of the
    dipole reading; paired, of "pair omitted" (D2478 and R2471 left out).
    The box's charge is counted beyond the lining's."""
    from ip3r.physics.bridge_charge import READINGS, pair_readings
    from ip3r.physics.lumen_charge import charged_lumen
    from ip3r.physics.lumen_field import lumen_field
    st, s = tkf
    f = lumen_field(st, s, spacing=1.0)
    chosen = {k: READINGS[k] for k in ("dipole", "pair omitted")}
    ref = pair_readings(st, 2478, 2471, s, spacing=1.0, readings=chosen)
    for paired, label in ((False, "dipole"), (True, "pair omitted")):
        c = charged_lumen(st, f, "dielectric", s, pair_bridges=paired)
        assert c.converged
        assert c.g == pytest.approx(ref[label].per_species["dielectric"]["K+"],
                                    rel=1e-6), label
        assert c.ratio > 1.0
        assert c.wall.placed < c.charge.net_charge - 1.0   # more than the lining
        assert "in the box" in c.summary()
    inside = f.volume.mask
    assert np.all(np.isfinite(c.u[inside])) and np.all(np.isnan(c.u[~inside]))
