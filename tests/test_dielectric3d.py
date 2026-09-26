"""Round 7.13: the dielectric closure, calibrated before it is read.

Coulomb's law in a uniform medium, the flux condition at a dielectric
interface, no field without charge, charge conservation, the insulating
limit against Round 7.11's Poisson–Boltzmann, and a dipole behind a wall;
then the deposits: D2478 ionised by both pKa routes, and the pair's field
between its two limits.
"""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics.charge3d import local_density, poisson_boltzmann
from ip3r.physics.dielectric3d import dielectric_pb, point_density
from ip3r.physics.ohmic3d import cylinder_volume
from ip3r.physics.permeation import potassium_species
from ip3r.physics.pore_charge import AVOGADRO
from ip3r.physics.radial_pb import EPS0
from ip3r.structure.pore_volume import PoreVolume

_E, _KB = 1.602176634e-19, 1.380649e-23


def _box(n, nz, h, grounded="all"):
    """Every voxel solvent; Dirichlet on every face ("all") or on the two
    end planes ("ends")."""
    xs = (np.arange(n) - (n - 1) / 2) * h
    zs = (np.arange(nz) - (nz - 1) / 2) * h
    mask = np.ones((n, n, nz), bool)
    top = np.zeros_like(mask)
    top[:, :, -1] = True
    bottom = np.zeros_like(mask)
    bottom[:, :, 0] = True
    if grounded == "all":
        top[[0, -1], :, :] = True
        top[:, [0, -1], :] = True
    return PoreVolume(mask, xs, zs, h, top, bottom, (zs[0], zs[-1]), 0.0)


def _bjerrum_vacuum() -> float:
    t = _P.value("permeation.temperature")
    return _E ** 2 / (4 * np.pi * EPS0 * _KB * t) * 1e10          # Å


def _charge(vol, x):
    return float(x.sum() * vol.spacing ** 3 * 1e-30 * AVOGADRO)


# ------------------------------------------------------------------ calibration
def test_a_charge_in_a_uniform_medium_follows_coulomb():
    vol = _box(61, 61, 1.0)
    eps = 40.0
    x = point_density(vol, np.zeros((1, 3)), np.array([1.0]), width=1.0)
    u, ok, _ = dielectric_pb(vol, np.full(vol.mask.shape, eps), x, [])
    assert ok
    c = 30
    r1, r2 = 3, 6
    expected = _bjerrum_vacuum() / eps * (1 / r1 - 1 / r2)
    assert u[c + r1, c, c] - u[c + r2, c, c] == pytest.approx(expected, rel=0.03)
    # the same along the diagonal of a face: the grid is isotropic enough
    assert u[c + 3, c + 4, c] - u[c + 6, c + 8, c] == pytest.approx(
        _bjerrum_vacuum() / eps * (1 / 5 - 1 / 10), rel=0.03)


def test_the_displacement_is_continuous_across_an_interface():
    """A charged plane, ε 40 to one side of an interface and 4 beyond it:
    the field in the two uncharged layers scales inversely with ε, exactly."""
    vol = _box(3, 41, 0.5, grounded="ends")
    eps = np.full(vol.mask.shape, 40.0)
    eps[:, :, 25:] = 4.0
    x = np.zeros(vol.mask.shape)
    x[:, :, 10] = -50.0
    u, ok, _ = dielectric_pb(vol, eps, x, [])
    assert ok
    line = u[1, 1]
    s_water = line[20] - line[19]
    s_protein = line[31] - line[30]
    assert s_water * 40.0 == pytest.approx(s_protein * 4.0, rel=1e-6)
    assert np.allclose(np.diff(line[11:25]), s_water, rtol=1e-6)


def test_no_charge_no_field():
    vol = cylinder_volume(3.0, 20.0, 1.0, 10.0, 8.0)
    eps = np.where(vol.mask, 40.0, 4.0)
    u, ok, _ = dielectric_pb(vol, eps, np.zeros(vol.mask.shape),
                             potassium_species())
    assert ok and not np.any(np.abs(u) > 1e-12)


def test_point_density_conserves_charge_wherever_it_lies():
    vol = cylinder_volume(3.0, 20.0, 0.5, 10.0, 8.0)
    pos = np.array([[0.0, 0.0, 0.0],        # in the lumen
                    [7.0, 1.0, 2.0],        # in the protein
                    [0.0, 0.0, 99.0]])      # outside the box: dropped
    q = np.array([-1.0, 2.0, 5.0])
    for i, expected in ((0, -1.0), (1, 2.0)):
        x = point_density(vol, pos[i:i + 1], q[i:i + 1], width=1.0)
        assert _charge(vol, x) == pytest.approx(expected, rel=1e-12)
    assert not point_density(vol, pos[2:], q[2:], width=1.0).any()
    x = point_density(vol, pos[1:2], q[1:2], width=1.0)
    assert not x[vol.mask].any()            # 7 Å out, 1 Å wide: all in protein


def test_an_insulating_protein_is_round_7_11s_poisson_boltzmann():
    """Charge in the lumen and the protein's permittivity → 0: the whole-box
    solve becomes the lumen-only solve with no flux into the protein."""
    vol = cylinder_volume(3.0, 16.0, 0.5, 6.0, 6.0)
    sp = potassium_species()
    pos = np.array([[2.5, 0.0, -3.0], [-2.5, 0.0, 3.0]])
    x, unreached = local_density(vol, pos, np.array([-2.0, -1.0]), width=1.5)
    assert not unreached
    ref, ok_ref, _ = poisson_boltzmann(vol, x, sp, permittivity=40.0)
    eps = np.where(vol.mask, 40.0, 40.0 * 1e-4)
    u, ok, _ = dielectric_pb(vol, eps, x, sp)
    assert ok and ok_ref
    assert np.max(np.abs(u - ref)[vol.mask]) < 0.02
    assert np.min(ref[vol.mask]) < -1.0                  # a real well was tested


def test_a_partner_behind_the_wall_cancels_part_of_the_field():
    """An acid at the wall of a 3 Å pore and a base 3.6 Å further out. On
    the axis, the pair's potential in uniform ε is the Coulomb ratio
    1 − 4.2/7.8 = 0.46, raised a little because the salt screens the
    farther charge more. A low-ε protein cancels less: the buried base's
    field reaches the water less than the acid's. The voxel next to the acid
    is the acid's alone, whatever the partner."""
    vol = cylinder_volume(3.0, 20.0, 0.5, 10.0, 8.0)
    sp = potassium_species()
    c, k = int(np.argmin(np.abs(vol.xs))), int(np.argmin(np.abs(vol.zs)))
    acid, base = np.array([[4.2, 0.0, 0.0]]), np.array([[7.8, 0.0, 0.0]])
    ratio, nearest = {}, {}
    for protein in (40.0, 4.0):
        eps = np.where(vol.mask, 40.0, protein)
        alone, _, _ = dielectric_pb(vol, eps, point_density(vol, acid, [-1.0]), sp)
        both, _, _ = dielectric_pb(vol, eps, point_density(
            vol, np.vstack([acid, base]), [-1.0, 1.0]), sp)
        assert alone[c, c, k] < both[c, c, k] < 0.0
        ratio[protein] = both[c, c, k] / alone[c, c, k]
        nearest[protein] = both[vol.mask].min() / alone[vol.mask].min()
    assert 0.46 < ratio[40.0] < 0.56
    assert ratio[40.0] < ratio[4.0] < 0.7
    assert nearest[4.0] > 0.9


# -------------------------------------------------------------------- real data
def _deposit(pdb):
    from ip3r.io import loader
    try:
        return loader.load(pdb)
    except Exception:
        pytest.skip(f"{pdb} not fetched")


def test_real_8tkf_d2478_is_ionised_by_both_routes():
    """The network with R2471′ counted puts D2478's pKa below 3; without it,
    still below 5; PROPKA at 5.2. The base stays charged. At pH 7.3 the pair
    is two charges under every reading."""
    pytest.importorskip("propka")
    from ip3r.physics.bridge_charge import lining_bridges, titrate_pair
    st = _deposit("8TKF")
    assert (2478, 2471) in lining_bridges(st)
    rows = {t.reading: t for t in titrate_pair(st, 2478, 2471)}
    assert rows["network"].acid[1] < 3.0 < rows["network -base"].acid[0]
    assert all(t.acid[1] < 5.5 for t in rows.values())
    assert all(t.base[0] > 11.0 for t in rows.values() if t.base is not None)
    assert all(max(t.acid_charge()) < -0.98 for t in rows.values())


def test_real_8tkf_the_dipole_lies_between_its_two_limits():
    """The pair as two charges raises g less than D2478 alone and more than
    no pair at all; every reading converged, every charge in the box placed."""
    from ip3r.physics.bridge_charge import READINGS, pair_readings
    st = _deposit("8TKF")
    got = pair_readings(st, 2478, 2471, spacing=1.0)
    r = {k: w.ratio(READINGS[k][0]) for k, w in got.items()}
    assert all(w.converged for w in got.values())
    assert r["pair omitted"] < r["dipole"] < r["base omitted"]
    assert got["dipole"].fields["dielectric"].placed == pytest.approx(
        got["base omitted"].fields["dielectric"].placed + 4.0, abs=1e-6)
