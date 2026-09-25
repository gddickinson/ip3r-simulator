"""Protonation readings of the pore (Round 7.4) on the real deposits.

The network and PROPKA are each held to their own truncations (a wider
network or context changes nothing), the four copies of a ring agree, the
wall accepts mean charges without losing any, Xu's ruler is the GHK
equation it says it is, and the findings the round reports are pinned:
the lysine rings stay charged at any permittivity, and no reading brings
8TKF within 20x of P_Ca:P_K 15.2.
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import pka
from ip3r.physics import protonation as pr
from ip3r.physics import selectivity as sel
from conftest import needs_structure

_HAS_PROPKA = importlib.util.find_spec("propka") is not None
needs_propka = pytest.mark.skipif(not _HAS_PROPKA, reason="propka not installed")


@pytest.fixture(scope="module")
def wall_8tkf():
    from ip3r.io import loader
    return pr.lining_wall(loader.load("8TKF"))


@pytest.fixture(scope="module")
def wall_9heo():
    from ip3r.io import loader
    return pr.lining_wall(loader.load("9HEO"))


def _lining_charge(wall, t):
    lining = {(g.chain, g.res_seq) for g in wall.groups}
    return {s.key: q for s, q in zip(t.sites, t.charge) if s.key in lining}


# ----------------------------------------------------------- the ruler
def test_xu_ruler_is_their_equation_1():
    """Xu's Eq. 1, E = (RT/F) ln(sqrt((K + 4 r Ca) / K)), is GHK with no Cl-
    term; their wild type's 9.5 mV gives 6.9 at 23 C (7.0 printed)."""
    k, ca = 0.25, 0.010
    cyt, lum = {"K+": k, "Cl-": k}, {"K+": k, "Cl-": k + 2 * ca, "Ca2+": ca}
    phi = sel.thermal_voltage()
    for r in (1.0, 3.3, 7.0):
        v = 0.5 * phi * np.log(1 + 4 * r * ca / k)
        assert sel.ghk_ratio(v, "Ca2+", lum, cyt, {"Cl-": 0.0}) == pytest.approx(r)
    assert sel.ghk_ratio(0.0095, "Ca2+", lum, cyt, {"Cl-": 0.0}) == pytest.approx(
        _P.value("selectivity.published_ryr1_pca_pk"), abs=0.15)


# ------------------------------------------------------------ the wall
@needs_structure("8TKF")
def test_wall_takes_mean_charges_without_losing_any(wall_8tkf):
    ch = {(g.chain, g.res_seq): 0.25 * g.charge for g in wall_8tkf.groups}
    fixed, net = wall_8tkf.density(ch)
    full, net_full = wall_8tkf.density()
    assert net == pytest.approx(0.25 * net_full)
    assert fixed == pytest.approx(0.25 * full)


# ------------------------------------------------------------ network
@needs_structure("8TKF")
def test_network_rings_agree_and_radius_is_wide_enough(wall_8tkf):
    """Sigmoidal permittivity: 25 A changes no lining charge by 0.01 e and the
    four copies agree to 0.01. The uniform eps-4 bound is long-ranged by
    construction (0.6 kT at 20 A): the radius moves it by < 0.04 e, and its
    copies differ by up to 0.15 e, because strong coupling amplifies the
    small differences between 8TKF's subunits (two seeds agree to 0.02)."""
    ph, ionic = pr.family_conditions(False)
    centres = wall_8tkf.centres()
    for eps, radius_tol, copies_tol in ((None, 0.01, 0.01), (4.0, 0.04, 0.2)):
        near = pka.titrate(pka.sites(wall_8tkf.structure, centres)[0], ph,
                           ionic, eps=eps)
        wide = pka.titrate(pka.sites(wall_8tkf.structure, centres, 25.0)[0],
                           ph, ionic, eps=eps)
        a, b = _lining_charge(wall_8tkf, near), _lining_charge(wall_8tkf, wide)
        assert max(abs(a[k] - b[k]) for k in a) < radius_tol, eps
        by: dict[int, list[float]] = {}
        for (_, seq), q in a.items():
            by.setdefault(seq, []).append(q)
        assert max(np.ptp(v) for v in by.values()) < copies_tol, eps


@needs_structure("8TKF")
def test_lysine_rings_stay_charged_at_any_permittivity(wall_8tkf):
    for eps in (None, 10.0, 4.0):
        charges, _ = pr.network_charges(wall_8tkf, eps)
        for g in wall_8tkf.groups:
            if g.res_name == "LYS":
                assert charges[(g.chain, g.res_seq)] > 0.99, (g.label(), eps)


# -------------------------------------------------------------- PROPKA
@needs_propka
@needs_structure("8TKF")
def test_propka_context_is_wide_enough(wall_8tkf):
    from ip3r.physics.pka_propka import propka_pka
    c = wall_8tkf.centres()
    a = propka_pka(wall_8tkf.structure, c, 25.0)
    b = propka_pka(wall_8tkf.structure, c, 35.0)
    for g in wall_8tkf.groups:
        k = (g.chain, g.res_seq)
        assert a[k][1] == pytest.approx(b[k][1], abs=0.01), k


@needs_propka
@needs_structure("9HEO")
def test_propka_misses_e4900_that_xu_shows_is_charged(wall_9heo):
    """E4900N moves RyR1's conductance (x0.63) and P_Ca:P_K (x0.64), so
    E4900 carries charge; PROPKA buries it (pKa ~8), the network does not."""
    ch, pk = pr.propka_charges(wall_9heo)
    net, _ = pr.network_charges(wall_9heo)
    for g in wall_9heo.groups:
        if g.res_seq == 4900:
            assert pk[(g.chain, 4900)][1] > 7.5
            assert net[(g.chain, 4900)] < -0.95


# ------------------------------------------------------------ findings
@needs_structure("8TKF")
def test_lysine_corners_bound_and_reproduce_the_acidic_reading(wall_8tkf):
    out = pr.corners(wall_8tkf, residues=[2482, 2529])
    best_off, best = out[0]
    assert best_off == frozenset({2482, 2529})
    assert best.pca_pk == pytest.approx(0.69, abs=0.02)   # selectivity's "acidic"
    assert best.pca_pk < _P.value("selectivity.published_pca_pk") / 20


@needs_propka
@needs_structure("8TKF")
def test_no_protonation_reading_approaches_vais(wall_8tkf):
    rows, _ = pr.readings(wall_8tkf)
    assert [r.label for r in rows] == ["formal", "network", "network eps 10",
                                       "network eps 4", "PROPKA"]
    assert all(r.converged for r in rows)
    assert max(r.pca_pk for r in rows) < 0.1
