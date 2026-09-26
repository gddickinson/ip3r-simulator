"""Round 7.11: the charged pore in 3-D, calibrated before it is read.

The reduction (linear response = equilibrium + one weighted Laplace per
species) is held to the existing 1-D Gummel solver; the weighted Laplace to
its 1-D series by hand; the charge placements to conservation; the
Poisson–Boltzmann solve to its Donnan limit, its Debye decay and Gauss's law.
"""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics._pnp_kernels import F_FARADAY, R_GAS
from ip3r.physics.charge3d import (donnan_field, local_density,
                                   poisson_boltzmann, slice_density,
                                   wall_field)
from ip3r.physics.charged3d import linear_response_1d
from ip3r.physics.ohmic3d import (bernoulli_weight, cylinder_volume,
                                  geometric_conductance)
from ip3r.physics.permeation import potassium_species, solve_pnp
from ip3r.physics.pore_charge import AVOGADRO, ChargedGroup, PoreCharge
from ip3r.physics.radial_pb import EPS0
from ip3r.structure.pore_volume import volume_from_mask


def _tube(n=4, length=60.0, h=0.5):
    """A square tube filling the box: sides insulating, the end planes bath.
    Anything uniform in-plane is exactly 1-D on it."""
    xs = (np.arange(n) - (n - 1) / 2) * h
    zs = np.arange(-length / 2, length / 2 + 1e-9, h)
    acc = np.ones((n, n, len(zs)), bool)
    return volume_from_mask(acc, xs, zs, (zs[0] + h / 2, zs[-1] - h / 2), 0.0,
                            seal=1e3)


def _group(z, q, res=1):
    return ChargedGroup(res, "ASP" if q < 0 else "LYS", "A", z, 0.0, 0.0, q)


def _ring_profile():
    z = np.linspace(-30, 30, 241)
    r = 4.0 + 1.5 * np.cos(np.pi * z / 30)
    fixed = (-3000 * np.exp(-0.5 * ((z + 8) / 3) ** 2)
             + 2000 * np.exp(-0.5 * ((z - 8) / 3) ** 2))
    return z, r, fixed


# ------------------------------------------------ the linear-response reduction
def test_linear_response_is_the_gummel_solver_at_small_voltage():
    """Two rings of opposite sign: the Gummel loop at 1 mV and the series of
    Boltzmann-weighted resistors agree; the rings lower g, as in 1-D."""
    z, r, fixed = _ring_profile()
    sp = potassium_species()
    solved = solve_pnp(z, r, voltage=0.001, species=sp, fixed_charge=fixed)
    lr = linear_response_1d(z, r, sp, fixed)
    assert solved.converged
    assert lr == pytest.approx(solved.conductance, rel=0.01)
    assert linear_response_1d(z, r, sp) == pytest.approx(
        solve_pnp(z, r, voltage=0.001, species=sp).conductance, rel=0.01)
    assert lr < linear_response_1d(z, r, sp)


def test_bernoulli_weight_limits():
    assert bernoulli_weight(np.array([0.7]), np.array([0.7]))[0] == pytest.approx(np.exp(-0.7))
    # 1 / mean(e^E) with E linear 0 -> 2, by quadrature
    e = np.linspace(0, 2, 20001)
    assert bernoulli_weight(np.array([0.0]), np.array([2.0]))[0] == pytest.approx(
        1.0 / np.trapezoid(np.exp(e), e) * 2.0, rel=1e-6)
    assert bernoulli_weight(np.array([0.0]), np.array([1e-12]))[0] == pytest.approx(1.0)


def test_weighted_laplace_is_the_1d_series_on_a_tube():
    vol = _tube(n=3, length=20.0, h=0.5)
    rng = np.random.default_rng(1)
    e_plane = rng.normal(0, 1.5, len(vol.zs))
    energy = np.broadcast_to(e_plane, vol.mask.shape).copy()
    g = geometric_conductance(vol, tol=1e-13, energy=energy).g
    w = bernoulli_weight(e_plane[:-1], e_plane[1:])
    n_xy = 9
    expected = vol.spacing * 1e-10 / np.sum(1.0 / (n_xy * w))
    assert g == pytest.approx(expected, rel=1e-8)


def test_uniform_energy_scales_the_conductance():
    vol = cylinder_volume(4.0, 20.0, 1.0, 12.0, 8.0)
    g0 = geometric_conductance(vol).g
    g = geometric_conductance(vol, energy=np.full(vol.mask.shape, 1.3)).g
    assert g == pytest.approx(np.exp(-1.3) * g0, rel=1e-6)
    assert geometric_conductance(vol, energy=np.zeros(vol.mask.shape)).g == pytest.approx(g0, rel=1e-12)


# ------------------------------------------------------------ charge placement
def _placed(vol, x):
    return float(x.sum() * vol.spacing ** 3 * 1e-30 * AVOGADRO)


def test_placements_conserve_charge():
    vol = cylinder_volume(4.0, 30.0, 0.5, 10.0, 6.0)
    groups = [_group(-5.0, -1.0), _group(4.0, +1.0, 2), _group(9.0, -1.0, 3)]
    # the axial Gaussians are normalised by the trapezoid rule, summed by voxel
    assert _placed(vol, slice_density(vol, groups)) == pytest.approx(-1.0, rel=1e-4)
    pos = np.array([[4.5, 0.0, -5.0], [0.0, 4.5, 4.0], [-4.5, 0, 9.0]])
    x, unreached = local_density(vol, pos, np.array([-1.0, 1.0, -1.0]))
    assert unreached == []
    assert _placed(vol, x) == pytest.approx(-1.0, rel=1e-9)
    # the local charge sits by its centre, not across the plane
    k = int(np.argmin(np.abs(vol.zs + 5.0)))
    plane = x[:, :, k]
    i = int(np.argmin(np.abs(vol.xs - 3.5)))
    j = int(np.argmin(np.abs(vol.xs + 3.5)))
    zero = int(np.argmin(np.abs(vol.xs)))
    assert abs(plane[i, zero]) > 3 * abs(plane[j, zero])


def test_a_charge_that_reaches_no_voxel_is_reported():
    vol = cylinder_volume(3.0, 20.0, 1.0, 8.0, 4.0)
    x, unreached = local_density(vol, np.array([[80.0, 0.0, 0.0]]),
                                 np.array([-1.0]), ["ASP9/A"])
    assert unreached == ["ASP9/A"] and not x.any()


def test_no_charge_is_the_neutral_pore():
    vol = cylinder_volume(4.0, 20.0, 1.0, 10.0, 6.0)
    empty = PoreCharge(np.zeros(1), np.zeros(1))
    sp = potassium_species()
    for closure in ("slice", "local", "pb"):
        wf = wall_field(vol, closure, sp, empty, positions=np.zeros((0, 3)))
        assert not wf.potential.any() and wf.placed == 0.0


# ------------------------------------------------------------ Poisson–Boltzmann
def _kappa(sp, eps):
    ionic = sum(s.valence ** 2 * s.concentration * 1000.0 for s in sp)
    return np.sqrt(F_FARADAY ** 2 * ionic / (EPS0 * eps * R_GAS
                                             * _P.value("permeation.temperature")))


def test_pb_reaches_donnan_in_a_long_charged_tube():
    vol = _tube(n=3, length=120.0, h=0.5)
    x = np.zeros(vol.mask.shape)
    x[:, :, np.abs(vol.zs) <= 40.0] = -500.0            # mol/m³, 0.5 M
    sp = potassium_species()
    u, ok, _ = poisson_boltzmann(vol, x, sp, permittivity=40.0)
    donnan = donnan_field(vol, x, sp)
    mid = int(np.argmin(np.abs(vol.zs)))
    assert ok
    assert u[1, 1, mid] == pytest.approx(donnan[1, 1, mid], rel=1e-3)
    assert u[1, 1, mid] == pytest.approx(-np.arcsinh(0.5 / (2 * 0.14)), rel=1e-3)


def test_pb_decays_at_the_debye_length_and_obeys_gauss():
    vol = _tube(n=3, length=100.0, h=0.5)
    x = np.zeros(vol.mask.shape)
    x[:, :, np.abs(vol.zs) <= 0.5] = -5.0               # weak: linear regime
    sp = potassium_species()
    eps = 40.0
    u, ok, _ = poisson_boltzmann(vol, x, sp, permittivity=eps)
    assert ok
    kappa = _kappa(sp, eps) * 1e-10                      # 1/Å
    z1, z2 = 10.0, 25.0
    k1, k2 = (int(np.argmin(np.abs(vol.zs - z))) for z in (z1, z2))
    rate = np.log(u[1, 1, k1] / u[1, 1, k2]) / (vol.zs[k2] - vol.zs[k1])
    assert rate == pytest.approx(kappa, rel=0.02)
    ions = sum(s.valence * s.concentration * 1000.0 * np.exp(-s.valence * u[vol.mask])
               for s in sp)
    assert float(ions.sum()) == pytest.approx(-float(x[vol.mask].sum()), rel=1e-3)


# ---------------------------------------------------------------------- real data
def test_real_8tkf_linear_response_matches_the_solver():
    from ip3r.io import loader
    from ip3r.physics.unitary import unitary
    try:
        st = loader.load("8TKF")
    except Exception:
        pytest.skip("8TKF not fetched")
    u = unitary(st)
    sp = potassium_species(bath=u.bath)
    r = np.maximum(u.profile.r_free, 0.0)
    for fixed in (u.charge.density, u.paired_charge.density):
        solved = solve_pnp(u.profile.z, r, voltage=0.001, species=sp,
                           fixed_charge=fixed)
        assert linear_response_1d(u.profile.z, r, sp, fixed) == pytest.approx(
            solved.conductance, rel=0.01)


def test_a_narrow_cylinder_keeps_the_1d_junctions_under_every_closure():
    """Two opposite rings in a 3 Å cylinder (point ions): the slice closure
    is the 1-D reading, and the other two, the charge placed at the wall
    and screened, agree with it, because the lumen is narrow against both
    the Gaussian and the Debye length. The rings cut g ~37×."""
    from ip3r.physics.permeation import IonSpecies, _species_conductivity
    from ip3r.physics.pore_charge import map_charge
    temperature = _P.value("permeation.temperature")
    sp = [IonSpecies("K+", 1, 1.96e-9, 0.0, 0.14),
          IonSpecies("Cl-", -1, 2.03e-9, 0.0, 0.14)]
    radius, length = 3.0, 40.0
    groups = [ChargedGroup(1, "ASP", "A", -8.0, radius, radius, -4.0),
              ChargedGroup(2, "LYS", "A", 8.0, radius, radius, 4.0)]
    vol = cylinder_volume(radius, length, 0.5, 14.0, 10.0)
    g0 = geometric_conductance(vol).g
    sigma = {s.name: _species_conductivity(s, temperature) for s in sp}
    charge = PoreCharge(np.zeros(1), np.zeros(1), groups)
    pos = np.array([[radius, 0, -8.0], [radius, 0, 8.0]])
    ratio = {}
    for closure in ("slice", "local", "pb"):
        wf = wall_field(vol, closure, sp, charge, positions=pos)
        g = sum(sigma[s.name] * geometric_conductance(
            vol, energy=wf.energy(s.valence)).g for s in sp)
        ratio[closure] = g / (sum(sigma.values()) * g0)
    z = np.linspace(-length / 2, length / 2, 801)
    r = np.full_like(z, radius)
    one_d = (linear_response_1d(z, r, sp, map_charge(groups, z, r))
             / linear_response_1d(z, r, sp))
    assert ratio["slice"] == pytest.approx(one_d, rel=0.05)
    for closure in ("local", "pb"):
        assert ratio["slice"] <= ratio[closure] < 1.5 * ratio["slice"]


def _deposit(pdb):
    from ip3r.io import loader
    try:
        return loader.load(pdb)
    except Exception:
        pytest.skip(f"{pdb} not fetched")


def test_real_8tkf_paired_wall_lowers_g_under_every_closure():
    """With D2478's salt bridges paired, the ITPR3 wall lowers the 3-D
    conductance whatever the placement; the unpaired wall raises it once
    the charge sits at its own centres, and most when screened."""
    from ip3r.physics.charged3d import wall_3d
    st = _deposit("8TKF")
    paired = wall_3d(st, spacing=1.0, pair_bridges=True)
    assert paired.converged
    assert all(paired.ratio(c) < 0.7 for c in paired.charged)
    full = wall_3d(st, spacing=1.0)
    assert full.ratio("slice") < 1.0 < full.ratio("local") < full.ratio("pb")
    assert all(abs(f.placed - full.charge.net_charge) < 1e-3
               for f in full.fields.values())


def test_real_9heo_acidic_wall_raises_g_under_every_closure():
    from ip3r.physics.charged3d import wall_3d
    w = wall_3d(_deposit("9HEO"), spacing=1.0)
    assert all(w.ratio(c) > 1.2 for c in w.charged)
