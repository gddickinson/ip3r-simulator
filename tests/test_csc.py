"""Charge-space competition (Round 7.17).

Each term is held to a limit it must reach without the code under test:
hard spheres to Carnahan-Starling and to the derivative of their own free
energy; the MSA to Debye-Hueckel, to the restricted primitive model's
closed form and to the derivative of Blum & Hoye's free energy. The
partition with both terms off must be the Donnan partition, a bath must
come back as itself, and the wall's groups must carry exactly the charge
the fixed-charge map carries. Then the physics Gillespie 2008 reports: a
small monovalent crowds Ca2+ out of a charged slice, a large one lets it
in (his Fig. 8), and RyR1's filter reproduces his split of the binding
energy (Fig. 7).
"""

from __future__ import annotations

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import csc
from ip3r.physics._pnp_kernels import _donnan_potential
from ip3r.physics.permeation import IonSpecies, _accessible_area, solve_pnp
from ip3r.physics.pnp_closures import CscClosure, solve_offsets
from ip3r.physics.pore_charge import ChargedGroup, map_charge
from conftest import needs_structure

_TH = 0.0255
_SIG = np.array([2.76, 2.0, 3.62, 2.8])


def _fd(f, rho, i, h=1e-8):
    up, dn = rho.copy(), rho.copy()
    up[i] += h
    dn[i] -= h
    return float((f(up) - f(dn)) / (2 * h))


# ------------------------------------------------------------ hard spheres
@pytest.mark.parametrize("eta", [0.05, 0.3, 0.45])
def test_one_component_is_carnahan_starling(eta):
    rho = eta / (np.pi / 6 * 2.8 ** 3)
    mu = csc.hard_sphere(np.array([[rho]]), np.array([2.8]))[0, 0]
    assert mu == pytest.approx((8 * eta - 9 * eta ** 2 + 3 * eta ** 3)
                               / (1 - eta) ** 3, rel=1e-12)


def test_hard_sphere_mu_is_the_free_energy_derivative():
    rho = np.array([0.004, 0.002, 0.003, 0.02])
    mu = csc.hard_sphere(rho[:, None], _SIG)[:, 0]
    for i in range(4):
        d = _fd(lambda r: csc.hs_free_energy(r[:, None], _SIG)[0], rho, i, 1e-7)
        assert mu[i] == pytest.approx(d, rel=1e-6)


def test_no_density_no_excess():
    assert np.all(csc.hard_sphere(np.zeros((3, 2)), _SIG[:3]) == 0.0)


# --------------------------------------------------------------- the MSA
def test_msa_mu_is_the_free_energy_derivative():
    z = np.array([1.0, 2.0, -1.0, -0.5])            # a half-charged oxygen too
    rho = np.array([0.004, 0.002, 0.005, 0.006])    # neutral
    mu = csc.msa(rho[:, None], z, _SIG, 7.1)[:, 0]
    for i in range(4):
        d = _fd(lambda r: csc.msa_free_energy(r[:, None], z, _SIG, 7.1)[0],
                rho, i)
        assert mu[i] == pytest.approx(d, rel=1e-6)


def test_msa_dilute_limit_is_debye_hueckel():
    lb, c = 7.1, 1e-11
    mu = csc.msa(np.array([[c], [c]]), np.array([1.0, -1.0]),
                 np.array([3.0, 3.0]), lb)[:, 0]
    kappa = np.sqrt(4 * np.pi * lb * 2 * c)
    assert mu == pytest.approx(-lb * kappa / 2, rel=1e-3)


def test_msa_restricted_primitive_model_closed_form():
    lb, c, s = 7.1, 0.001, 3.0
    mu = csc.msa(np.array([[c], [c]]), np.array([1.0, -1.0]),
                 np.array([s, s]), lb)[:, 0]
    kappa = np.sqrt(4 * np.pi * lb * 2 * c)
    gamma = (np.sqrt(1 + 2 * kappa * s) - 1) / (2 * s)
    assert mu == pytest.approx(-lb * gamma / (1 + gamma * s), rel=1e-10)


def test_bjerrum_length_by_hand():
    # 7.0 A in water at 25 C is the textbook value.
    assert csc.bjerrum(78.4, 298.15) == pytest.approx(7.15, abs=0.01)


# -------------------------------------------------- the wall as a fluid
def _groups():
    return [ChargedGroup(1, "ASP", "A", 10.0, 4.0, 4.0, -1.0),
            ChargedGroup(2, "GLU", "B", 14.0, 4.0, 4.0, -0.6),
            ChargedGroup(3, "LYS", "C", 25.0, 5.0, 5.0, 1.0)]


def test_structural_charge_is_the_fixed_charge_map():
    z = np.linspace(0.0, 40.0, 161)
    r = 3.0 + 0.1 * z
    wall = csc.structural(_groups(), z, r)
    assert wall.names == ("O", "N")
    np.testing.assert_allclose(wall.charge(), map_charge(_groups(), z, r),
                               rtol=1e-12, atol=1e-9)
    # two oxygens per acid, whatever its charge; the partly protonated
    # Glu lowers the oxygens' valence, not their number
    per_len = np.trapezoid(wall.density[0] * np.pi * np.maximum(
        r, _P.value("permeation.radius_potassium")) ** 2, z) * 1e-30 * 6.02214076e23
    assert per_len == pytest.approx(4.0, rel=1e-6)
    assert wall.valence[0].min() >= -0.5 - 1e-12


def test_from_fixed_carries_the_net_charge():
    f = np.array([-3000.0, 0.0, 500.0])
    np.testing.assert_allclose(csc.Structural.from_fixed(f).charge(), f)


# ------------------------------------------------------- the partition
def _kcl_fluid(wall, **kw):
    sp = [IonSpecies("K+", 1, 1e-9, 1.38, 0.15),
          IonSpecies("Cl-", -1, 1e-9, 1.81, 0.15)]
    return sp, csc.Fluid.of(sp, wall, **kw)


def test_terms_off_is_the_donnan_partition():
    fixed = np.array([-20000.0, -3000.0, 0.0, 800.0])
    wall = csc.Structural.from_fixed(fixed)
    _, fl = _kcl_fluid(wall, hs=False, use_msa=False)
    bath = np.full((2, 4), 150.0)
    p = csc.partition(bath, fixed, fl, _TH, fl.water)
    assert p.converged
    np.testing.assert_allclose(p.psi, _donnan_potential([1, -1], bath, fixed,
                                                        _TH), atol=1e-9)
    np.testing.assert_allclose(p.mu, 0.0)


def test_a_bath_comes_back_as_itself():
    _, fl = _kcl_fluid(csc.Structural.empty(1))
    bath = np.array([150.0, 150.0])
    a, aw = csc.bath_activity(bath, fl)
    p = csc.partition(a[:, None], np.zeros(1), fl, _TH, aw)
    np.testing.assert_allclose(p.c[:, 0], bath, rtol=1e-8)
    assert p.water[0] == pytest.approx(fl.water, rel=1e-8)
    assert abs(p.psi[0]) < 1e-10


def _crowded(mono_diameter):
    """150 mM X+ and 1 mM Ca2+ against a -20 M slice: (X+, Ca2+) in M."""
    fixed = np.array([-20000.0])
    sp = [IonSpecies("X+", 1, 1e-9, mono_diameter / 2, 0.15),
          IonSpecies("Ca2+", 2, 1e-9, 1.0, 0.001),
          IonSpecies("Cl-", -1, 1e-9, 1.81, 0.152)]
    fl = csc.Fluid.of(sp, csc.Structural.from_fixed(fixed))
    bath = np.array([s.concentration * 1000.0 for s in sp])
    a, aw = csc.bath_activity(bath, fl)
    p = csc.partition(a[:, None], fixed, fl, _TH, aw,
                      initial=(bath[:, None], np.array([fl.water])))
    assert p.converged
    return p.c[0, 0] / 1000.0, p.c[1, 0] / 1000.0


def test_a_small_monovalent_crowds_calcium_out():
    # Gillespie 2008 Fig. 8: Li+ (1.33 A) against Cs+ (3.40 A).
    li, ca_li = _crowded(1.33)
    cs, ca_cs = _crowded(3.40)
    assert li > 2.0 * cs
    assert ca_cs > ca_li


def test_terms_off_through_the_gummel_loop_is_the_donnan_closure():
    z = np.linspace(0.0, 30.0, 61)
    r = np.full_like(z, 5.0)
    fixed = -8000.0 * np.exp(-0.5 * ((z - 15.0) / 3.0) ** 2)
    sp = [IonSpecies("K+", 1, 1e-9, 1.38, 0.25),
          IonSpecies("Cl-", -1, 1e-9, 1.81, 0.27, 0.25),
          IonSpecies("Ca2+", 2, 0.4e-9, 1.0, 0.01, 0.0)]
    ref = solve_pnp(z, r, voltage=0.005, species=sp, fixed_charge=fixed)
    zz, rr = z * 1e-10, r * 1e-10
    areas = {s.name: _accessible_area(rr, s.radius) for s in sp}
    th = 8.314462618 * _P.value("permeation.temperature") / 96485.33212
    rule = CscClosure(sp, fixed, th, hs=False, use_msa=False)
    got = solve_offsets(rule, zz, rr, 0.005, sp, areas, fixed, th,
                        _P.value("permeation.temperature"), 400, 1e-10, 0.4,
                        True, False)
    assert got.converged
    assert got.pore_current == pytest.approx(ref.pore_current, rel=1e-6)


def test_the_csc_closure_runs_and_reports():
    z = np.linspace(0.0, 30.0, 61)
    r = np.full_like(z, 5.0)
    fixed = -8000.0 * np.exp(-0.5 * ((z - 15.0) / 3.0) ** 2)
    res = solve_pnp(z, r, species=[IonSpecies("K+", 1, 1e-9, 1.38, 0.25),
                                   IonSpecies("Cl-", -1, 1e-9, 1.81, 0.25)],
                    fixed_charge=fixed, closure="csc")
    assert res.converged and res.meta["closure"] == "csc"
    assert res.meta["csc_terms"]["hs"] and res.meta["csc_terms"]["msa"]
    assert min(res.meta["water_M"]) < 55.5


# ------------------------------------------------------------ RyR1 (real)
@needs_structure("9HEO")
def test_ryr1_filter_binding_matches_gillespie():
    from ip3r.io import loader
    from ip3r.physics import csc_readings as cr
    e = cr.energetics(loader.load("9HEO"), 0.15, 0.001)
    # His Fig. 7: screening ~4 kT, excluded volume ~0.5-1 kT, and the
    # filter holds more Ca2+ than K+ at 1 mM Ca2+ against 150 mM K+.
    assert 3.0 < e.screening < 5.5
    assert 0.3 < e.excluded < 1.5
    assert e.ca > 3 * e.k
    assert e.packing < 0.6


@needs_structure("9HEO")
def test_ryr1_ratio_is_the_uncharged_gate_s():
    from ip3r.io import loader
    from ip3r.physics import csc_readings as cr
    from ip3r.physics.selectivity import ryr1_calcium_ratio
    w = cr.csc_wall(loader.load("9HEO"))
    _, donnan, _ = ryr1_calcium_ratio(w.z, w.radius, w.fixed)
    _, p, ok = ryr1_calcium_ratio(w.z, w.radius, w.fixed, closure="csc",
                                  structural=w.wall)
    assert ok
    assert donnan < p < 1.0                   # up, and nowhere near Xu's 7.0
    # The gate window, uncharged, caps it: its own ratio (by hand) is below
    # the whole pore's, and it holds most of Ca2+'s resistance.
    assert cr.gate_limit(w) < p < 1.5 * cr.gate_limit(w)
    s = cr.shares(w, "csc")
    assert s.gate["Ca2+"] > 0.5
    assert s.gate["Ca2+"] > 20 * s.filter["Ca2+"]
    assert s.filter["K+"] > cr.shares(w, "donnan").filter["K+"]


def test_moved_restores_the_parameter():
    before = _P.value("csc.structural_volume")
    from ip3r.physics.csc_readings import moved
    with moved("csc.structural_volume", 0.0):
        assert _P.value("csc.structural_volume") == 0.0
    assert _P.value("csc.structural_volume") == before
    assert _P.is_default("csc.structural_volume")
