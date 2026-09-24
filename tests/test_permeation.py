"""Unitary conductance: the solver against closed forms, then the real panel.

The calibrations are the ones a drift-diffusion solver can be held to without
trusting it: a cylinder has an exact ohmic conductance; a uniformly charged
cylinder has an exact Donnan partition; the solver must agree with the
series-resistor sum it was not used to derive; and it must refuse a pore
narrower than the ion.
"""

import numpy as np
import pytest

from ip3r.core.structure import Structure
from ip3r.parameters import PARAMETERS as _P
from ip3r.physics.permeation import (F_FARADAY, R_GAS, access_resistance,
                                     bulk_conductivity, potassium_species,
                                     series_conductance, solve_pnp)
from ip3r.physics.pore_charge import AVOGADRO, charged_groups, map_charge
from ip3r.structure.pore import PoreProfile
from ip3r.structure.symmetry import Frame
from conftest import needs_structure

_Z = np.arange(-20.0, 20.0 + 0.25, 0.5)          # a 40 A pore, A


def _analytic(radius, conc):
    """Ohmic cylinder + two Hall mouths, species concentrations in mol/m^3."""
    t = _P.value("permeation.temperature")
    species = potassium_species()
    per_len = sum(s.valence ** 2 * F_FARADAY ** 2 * s.diffusivity * c / (R_GAS * t)
                  * np.pi * ((radius - s.radius) * 1e-10) ** 2
                  for s, c in zip(species, conc))
    pore = (_Z[-1] - _Z[0]) * 1e-10 / per_len
    access = 2 * access_resistance(radius * 1e-10, bulk_conductivity(species, t))
    return 1e12 / (pore + access)


def test_cylinder_matches_the_closed_form():
    c = _P.value("permeation.bath_concentration") * 1000.0
    expected = _analytic(5.0, (c, c))
    r = np.full_like(_Z, 5.0)
    assert series_conductance(_Z, r)["conductance"] * 1e12 == pytest.approx(expected, rel=1e-3)
    assert solve_pnp(_Z, r).conductance_pS == pytest.approx(expected, rel=1e-3)


def test_neutral_solver_equals_the_series_sum_on_an_hourglass():
    r = 3.0 + 0.004 * _Z ** 2                       # 3 A waist, 4.6 A mouths
    s = series_conductance(_Z, r)["conductance"] * 1e12
    assert solve_pnp(_Z, r).conductance_pS == pytest.approx(s, rel=0.02)


def test_uniform_charge_gives_the_donnan_partition():
    c = _P.value("permeation.bath_concentration") * 1000.0
    x = -2.0 * c                                     # net negative wall charge
    k = -x / 2 + np.sqrt(x ** 2 / 4 + c ** 2)       # K+ enriched
    expected = _analytic(5.0, (k, c ** 2 / k))
    r = np.full_like(_Z, 5.0)
    got = solve_pnp(_Z, r, fixed_charge=np.full_like(_Z, x))
    assert got.converged
    assert got.conductance_pS == pytest.approx(expected, rel=0.01)
    assert got.conductance_pS > solve_pnp(_Z, r).conductance_pS     # enrichment


def test_opposite_rings_in_series_block_like_a_junction():
    """A K+-excluding zone then a Cl--excluding zone: each carrier must cross
    the zone where it is a co-ion, so the pair conducts less than no charge —
    the mechanism behind 8TKF's lysine rings."""
    c = _P.value("permeation.bath_concentration") * 1000.0
    x = np.where(_Z < 0, 3.0 * c, -3.0 * c)
    r = np.full_like(_Z, 5.0)
    assert solve_pnp(_Z, r, fixed_charge=x).conductance_pS < solve_pnp(_Z, r).conductance_pS


def test_too_narrow_a_pore_is_shut_with_its_reason():
    r = np.full_like(_Z, 5.0)
    r[40] = 1.0                                      # below the K+ radius
    res = solve_pnp(_Z, r)
    assert res.conductance == 0.0 and "sterically occluded" in res.blocked_by


def _charged_structure():
    """One ASP per subunit whose carboxylate lines a 5 A lumen, and one
    stubbed ASP (C-alpha only) behind the wall."""
    xyz, names, seq, chain = [], [], [], []
    for k, ch in enumerate("ABCD"):
        t = np.pi / 2 * k
        u = np.array([np.cos(t), np.sin(t), 0.0])
        for atom, rad in (("CA", 9.0), ("OD1", 5.5), ("OD2", 5.5)):
            xyz.append(u * rad + [0, 0, 2.0])
            names.append(atom)
            seq.append(10)
            chain.append(ch)
        xyz.append(u * 11.0 + [0, 0, -8.0])
        names.append("CA")
        seq.append(20)
        chain.append(ch)
    n = len(xyz)
    st = Structure(
        xyz=np.array(xyz, np.float32), element=np.array([a[0] for a in names], "U2"),
        atom_name=np.array(names, "U6"), res_name=np.full(n, "ASP", "U5"),
        res_seq=np.array(seq, np.int32), chain=np.array(chain, "U6"),
        hetero=np.zeros(n, bool), b_factor=np.zeros(n, np.float32),
        occupancy=np.ones(n, np.float32), alt_loc=np.full(n, ".", "U2"),
        entity=np.full(n, "1", "U6"), name="SYN")
    st._build_residue_index()
    return st


def test_charges_are_placed_from_side_chains_and_stubs_are_counted():
    st = _charged_structure()
    frame = Frame(np.array([0, 0, 1.0]), np.zeros(3), np.eye(3), list("ABCD"))
    prof = PoreProfile(_Z, np.full_like(_Z, 5.0), np.full_like(_Z, 3.5))
    groups, unplaced = charged_groups(st, frame, prof)
    assert [g.label() for g in groups] == ["ASP10"] * 4
    assert all(abs(g.z - 2.0) < 1e-6 and abs(g.margin - 0.5) < 1e-5 for g in groups)
    assert unplaced == [f"ASP20/{c}" for c in "ABCD"]
    density = map_charge(groups, _Z, prof.r_free)
    area = np.pi * prof.r_free ** 2
    total = np.trapezoid(density * area, _Z) * AVOGADRO / 1e30
    assert total == pytest.approx(-4.0, rel=1e-6)    # the kernel conserves charge


@needs_structure("8TKF", "8TKG", "8TKH", "8TLA", "7T3P", "6DQN", "6DQJ")
def test_only_the_activated_state_conducts_and_falls_short_of_the_measurement():
    from ip3r.physics.unitary import published, unitary_panel
    rows = unitary_panel("ITPR3", sweep=True)
    open_ = [u for u in rows if u.neutral.is_conducting]
    assert [u.name.upper() for u in open_] == ["8TKF"]
    u = open_[0]
    assert u.neutral.conductance_pS == pytest.approx(u.series_pS, rel=0.03)
    assert u.charged.converged
    assert {lab for lab, _, _ in u.charge.residues()} >= {"ASP2478", "LYS2482"}
    # The finding: no corner of the unmeasured-constant sweep reaches either
    # measured value, with or without the wall charge.
    ceiling = max(u.sweep["neutral"][1], u.sweep["charged"][1])
    assert ceiling < min(published().values())
