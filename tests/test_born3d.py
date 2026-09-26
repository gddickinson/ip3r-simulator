"""Round 7.15: the image (Born) self-energy, calibrated before it is read.

The planar Debye–Hückel wall (the half-space Green's function's reflected
part), its unscreened limit by hand, zero in a uniform medium, the sign
when the wall is the higher permittivity, the box and reach cuts, parallel
= serial; the swept surface keeps every ion centre an ion radius from the
low-ε region; the image cost enters Poisson–Boltzmann and the Boltzmann
energy as z²W. Then 8TKF: the filter's axis costs more than the vestibule's.
"""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import born3d
from ip3r.physics.charge3d import WallField
from ip3r.physics.dielectric3d import dielectric_pb, point_density, swept_mask
from ip3r.physics.permeation import potassium_species
from ip3r.physics.radial_pb import EPS0
from ip3r.structure.pore_volume import PoreVolume

_E, _KB = 1.602176634e-19, 1.380649e-23
EW, EP = 40.0, 4.0


def _debye(eps, species) -> float:
    """Debye length, Å, by the textbook formula."""
    t = _P.value("permeation.temperature")
    ionic = sum(s.valence ** 2 * s.concentration * 1000.0 for s in species)
    return float(np.sqrt(eps * EPS0 * _KB * t / (_E ** 2 * 6.02214076e23 * ionic)) * 1e10)


def _wall(n, d, h=1.0, eps_wall=EP):
    """A cube of water with a planar wall: voxel planes below ``n - d`` are
    ``eps_wall``, so the centre voxel is (d + 0.5) h from the interface."""
    size = 2 * n + 1
    eps = np.full((size, size, size), EW)
    eps[:, :, :n - d] = eps_wall
    return eps


def test_bjerrum_length_by_hand():
    t = _P.value("permeation.temperature")
    assert born3d.bjerrum_vacuum() == pytest.approx(
        _E ** 2 / (4 * np.pi * EPS0 * _KB * t) * 1e10, rel=1e-9)


def test_unscreened_planar_image_is_the_textbook_one():
    d = 3.0
    lb = born3d.bjerrum_vacuum()
    assert born3d.planar_image(d, EW, EP) == pytest.approx(
        lb / (4 * EW * d) * (EW - EP) / (EW + EP), rel=1e-6)


def test_a_uniform_medium_costs_nothing():
    eps = np.full((25, 25, 25), EW)
    k2 = np.full(eps.shape, 0.01)
    w = born3d.self_energy(eps, k2, [(12, 12, 12)], 1.0, half=6.0,
                           eps_ref=EW, k2_ref=0.01)
    assert abs(w[0]) < 1e-6


@pytest.mark.parametrize("h,d,tol", [(1.0, 2, 0.05), (1.0, 3, 0.05),
                                     (0.5, 4, 0.02), (0.5, 6, 0.02)])
def test_the_screened_planar_wall_matches_debye_huckel(h, d, tol):
    """At the registered box and the 1 Å grid the deposits are read on,
    within 5 %; halving the grid in a 16 Å box, within 2 %."""
    species = potassium_species()
    lam = _debye(EW, species)
    n = int(round(16 / h))
    eps = _wall(n, d)
    t = _P.value("permeation.temperature")
    k2w = h * h * 1e-20 * (_E * 6.02214076e23) ** 2 / (EPS0 * 8.314462618 * t) * sum(
        s.valence ** 2 * s.concentration * 1000.0 for s in species)
    k2 = np.where(eps == EW, k2w, 0.0)
    half = _P.value("born.box_half_width") if h == 1.0 else 16.0
    w = born3d.self_energy(eps, k2, [(n, n, n)], h, half=half, eps_ref=EW,
                           k2_ref=k2w)[0]
    assert w == pytest.approx(born3d.planar_image((d + 0.5) * h, EW, EP, lam),
                              rel=tol)


def test_a_higher_permittivity_wall_attracts():
    n, d = 12, 2
    eps = _wall(n, d, eps_wall=400.0)
    k2 = np.where(eps == EW, 0.02, 0.0)
    w = born3d.self_energy(eps, k2, [(n, n, n)], 1.0, eps_ref=EW, k2_ref=0.02)[0]
    assert w < 0 and born3d.planar_image(d + 0.5, EW, 400.0, 7.0) < 0


def test_the_box_cut_is_converged_at_the_registered_half_width():
    """In the IP3R bath a 12 Å box and a 20 Å one agree within 3 %; the
    cut matters only because the salt screens the field."""
    n, d = 20, 2
    eps = _wall(n, d)
    k2w = born3d.screening(_slab(), np.ones((9, 9, 15), bool),
                           potassium_species()).max()
    k2 = np.where(eps == EW, k2w, 0.0)
    got = [born3d.self_energy(eps, k2, [(n, n, n)], 1.0, half=hw, eps_ref=EW,
                              k2_ref=k2w)[0]
           for hw in (_P.value("born.box_half_width"), 20.0)]
    assert got[0] == pytest.approx(got[1], rel=0.03)


def test_the_reach_cut_drops_only_what_a_wall_barely_touches():
    lam = _debye(EW, potassium_species())
    assert born3d.planar_image(_P.value("born.reach"), EW, EP, lam) < 0.02
    size = 31
    xs = np.arange(size) - 15.0
    eps = np.full((size, size, size), EW)
    eps[:, :, :3] = EP
    mask = eps == EW
    vol = PoreVolume(mask, xs, xs, 1.0, np.zeros_like(mask), np.zeros_like(mask),
                     (-15.0, 15.0), 0.0)
    pts = born3d.targets(vol, eps, reach=5.0, eps_water=EW)
    assert pts[:, 2].min() == 3 and pts[:, 2].max() == 3 + 4


def test_parallel_workers_give_the_serial_answer():
    n, d = 10, 1
    eps = _wall(n, d)
    k2 = np.where(eps == EW, 0.02, 0.0)
    pts = np.argwhere(eps == EW)
    pts = pts[(np.abs(pts - n) <= 3).all(1)][:210]
    one = born3d.self_energy(eps, k2, pts, 1.0, half=4.0, eps_ref=EW, k2_ref=0.02)
    two = born3d.self_energy(eps, k2, pts, 1.0, half=4.0, eps_ref=EW,
                             k2_ref=0.02, workers=2)
    np.testing.assert_allclose(one, two, atol=1e-12)


def test_the_swept_surface_keeps_centres_an_ion_radius_from_the_wall():
    reach = np.zeros((21, 21, 21), bool)
    reach[10, 10, 10] = True
    swept = swept_mask(reach, 2.0, 1.0)
    i = np.argwhere(swept) - 10
    assert np.linalg.norm(i, axis=1).max() == pytest.approx(2.0)
    assert swept.sum() == 33          # the lattice points within 2 of a point


def _slab(n=9, nz=15):
    xs = (np.arange(n) - (n - 1) / 2)
    zs = (np.arange(nz) - (nz - 1) / 2)
    mask = np.ones((n, n, nz), bool)
    top = np.zeros_like(mask)
    top[:, :, -1] = True
    bottom = np.zeros_like(mask)
    bottom[:, :, 0] = True
    return PoreVolume(mask, xs, zs, 1.0, top, bottom, (zs[0], zs[-1]), 0.0)


def test_an_image_cost_alone_leaves_a_symmetric_salt_neutral():
    vol = _slab()
    w = np.full(vol.mask.shape, 1.5)
    u, ok, _ = dielectric_pb(vol, np.full(vol.mask.shape, EW),
                             np.zeros(vol.mask.shape), potassium_species(),
                             self_energy=w)
    assert ok and np.abs(u).max() < 1e-9


def test_an_image_cost_weakens_the_screening_of_a_fixed_charge():
    vol = _slab()
    x = point_density(vol, np.zeros((1, 3)), np.array([-1.0]), width=1.0)
    species = potassium_species()
    eps = np.full(vol.mask.shape, EW)
    bare, _, _ = dielectric_pb(vol, eps, x, species)
    cost, _, _ = dielectric_pb(vol, eps, x, species,
                               self_energy=np.full(vol.mask.shape, 1.0))
    c = tuple(s // 2 for s in vol.mask.shape)
    assert cost[c] < bare[c] < 0


def test_the_boltzmann_energy_counts_z_squared():
    mask = np.ones((2, 2, 2), bool)
    w = np.full(mask.shape, 0.5)
    f = WallField("dielectric", np.zeros(mask.shape), np.full(mask.shape, -1.0),
                  mask, 0.0, self_energy=w)
    assert f.energy(1)[0, 0, 0] == pytest.approx(-0.5)
    assert f.energy(-1)[0, 0, 0] == pytest.approx(1.5)
    assert f.energy(2)[0, 0, 0] == pytest.approx(0.0)


# ------------------------------------------------------------------ the deposit
@pytest.fixture(scope="module")
def tkf():
    from ip3r.io import loader
    from ip3r.structure.channel import measure_channel
    from ip3r.structure.pore_volume import pore_volume
    try:
        st = loader.load("8TKF")
    except Exception as exc:                      # noqa: BLE001
        pytest.skip(f"8TKF not available: {exc}")
    s = measure_channel(st)
    species = potassium_species()
    vol = pore_volume(st, s.frame, s.span, species[0].radius, spacing=1.0)
    return st, s, vol, species


def test_8tkf_no_ion_centre_sits_within_its_radius_of_the_dielectric(tkf):
    from scipy import ndimage
    st, s, vol, _ = tkf
    eps, ions = born3d.ion_maps(st, s.frame, vol)
    d = ndimage.distance_transform_edt(eps == EW, sampling=vol.spacing)
    assert d[vol.mask].min() >= vol.probe - 1e-9


def test_8tkf_the_filter_axis_costs_more_than_the_vestibule(tkf):
    st, s, vol, species = tkf
    eps, ions = born3d.ion_maps(st, s.frame, vol)
    k2 = born3d.screening(vol, ions, species)
    r = vol.radius()
    pts = []
    for z in (s.constrictions["filter"].z, s.span[0] - 15.0):
        k = int(np.argmin(np.abs(vol.zs - z)))
        i = np.unravel_index(np.argmin(np.where(vol.mask[:, :, k], r, np.inf)), r.shape)
        pts.append((i[0], i[1], k))
    w = born3d.self_energy(eps, k2, pts, vol.spacing, k2_ref=float(k2.max()))
    assert 0.5 < w[0] < 3.0
    assert w[1] < 0.5 * w[0]


def test_closures_without_a_protein_refuse_the_image_cost():
    from ip3r.physics.charged3d import wall_3d
    with pytest.raises(ValueError, match="dielectric"):
        wall_3d(None, closures=("pb",), self_energy=np.zeros((2, 2, 2)))


def test_the_cube_template_builds_the_operator():
    rng = np.random.default_rng(3)
    eps = rng.choice([EP, EW], size=(7, 7, 7))
    k2 = rng.uniform(0, 0.1, size=eps.shape)
    a, _ = born3d._Template(7).build(eps.ravel(), k2.ravel())
    assert abs(a - born3d.operator(eps, k2)).max() < 1e-12


def test_the_registered_tolerance_is_tight_enough():
    n, d = 12, 2
    eps = _wall(n, d)
    k2 = np.where(eps == EW, 0.02, 0.0)
    loose = born3d.self_energy(eps, k2, [(n, n, n)], 1.0, eps_ref=EW, k2_ref=0.02)[0]
    before = _P.overrides().get("born.cg_tolerance")
    try:
        _P.set_value("born.cg_tolerance", 1e-10)
        tight = born3d.self_energy(eps, k2, [(n, n, n)], 1.0, eps_ref=EW,
                                   k2_ref=0.02)[0]
    finally:
        _P.reset("born.cg_tolerance") if before is None else _P.set_value(
            "born.cg_tolerance", before)
    assert abs(loose - tight) < 1e-4
