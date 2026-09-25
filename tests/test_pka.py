"""The pKa network (Round 7.4), calibrated before it is read.

Closed forms: an isolated site is Henderson-Hasselbalch; a site beside a
charge that never titrates is shifted by exactly W / ln10; two coupled
sites follow their four-state partition function. The Monte Carlo is held
to exact enumeration on a strongly coupled 12-site double ring. The
electrostatic constants are held to textbook values.
"""

from __future__ import annotations

import numpy as np
import pytest

from ip3r.physics import pka

_T = 298.15


def _ring(n, radius, z, phase=0.0):
    a = phase + 2 * np.pi * np.arange(n) / n
    return np.c_[radius * np.cos(a), radius * np.sin(a), np.full(n, z)]


# ----------------------------------------------------------- constants
def test_bjerrum_and_debye_are_textbook():
    assert pka.bjerrum_vacuum(_T) == pytest.approx(560.2, abs=0.5)
    assert pka.bjerrum_vacuum(_T) / 78.4 == pytest.approx(7.15, abs=0.02)
    # 0.304 / sqrt(I) nm in water at 25 C
    assert pka.debye_angstrom(0.1, _T) == pytest.approx(9.61, rel=0.01)
    assert pka.debye_angstrom(0.0, _T) == np.inf


def test_sigmoidal_permittivity_runs_from_contact_to_water():
    r = np.linspace(0.0, 60.0, 601)
    e = pka.eps_sigmoidal(r)
    assert e[0] == pytest.approx(1.35, abs=0.01)
    assert e[-1] == pytest.approx(78.4, abs=0.1)
    assert np.all(np.diff(e) > 0)


# ------------------------------------------------------- closed forms
@pytest.mark.parametrize("q0,pk,ph", [(-1.0, 3.67, 4.5), (0.0, 10.4, 9.8),
                                      (-1.0, 4.25, 7.3)])
def test_isolated_site_is_henderson_hasselbalch(q0, pk, ph):
    hh = 1.0 / (1.0 + 10.0 ** (ph - pk))
    w = np.zeros((1, 1))
    assert pka.exact(w, [pk], [q0], ph)[0] == pytest.approx(hh, rel=1e-12)
    assert pka.monte_carlo(w, [pk], [q0], ph, sweeps=200)[0] == pytest.approx(
        hh, rel=1e-9)                   # the conditional estimator is exact here


def test_fixed_neighbour_shifts_pka_by_w_over_ln10():
    """An acid beside a base that never deprotonates (pKa 30) sees a fixed
    +1: its pKa falls by exactly W / ln10."""
    w12 = 2.0
    w = np.array([[0.0, w12], [w12, 0.0]])
    ph = 4.0
    theta = pka.exact(w, [4.0, 30.0], [-1.0, 0.0], ph)
    t = pka.Titration([pka.Site("A", 1, "ASP", (0, 0, 0)),
                       pka.Site("A", 2, "LYS", (5, 0, 0))], ph, theta, "exact")
    assert t.apparent_pka()[0] == pytest.approx(4.0 - w12 / np.log(10), abs=1e-9)


def test_two_coupled_acids_match_their_partition_function():
    w12, pk, ph = 1.5, 4.0, 4.3
    a = np.log(10) * (ph - pk)
    # states (th1, th2): charges (th - 1); G = a (th1 + th2) + W q1 q2
    g = {(0, 0): w12, (1, 0): a, (0, 1): a, (1, 1): 2 * a}
    z = sum(np.exp(-v) for v in g.values())
    p1 = (np.exp(-g[(1, 0)]) + np.exp(-g[(1, 1)])) / z
    w = np.array([[0.0, w12], [w12, 0.0]])
    assert pka.exact(w, [pk, pk], [-1, -1], ph) == pytest.approx([p1, p1], rel=1e-12)
    assert pka.monte_carlo(w, [pk, pk], [-1, -1], ph, sweeps=20000) == \
        pytest.approx([p1, p1], abs=0.005)


def test_monte_carlo_matches_enumeration_on_a_coupled_double_ring():
    """Two rings of four acids 6 A apart around the axis, and a ring of four
    bases between them, at a uniform permittivity of 4: every pair couples by
    kT or more, the regime where mean-field titration fails."""
    xyz = np.vstack([_ring(4, 5.0, 0.0), _ring(4, 6.0, 6.0, 0.4),
                     _ring(4, 8.0, 3.0, 0.8)])
    w = pka.interactions(xyz, 0.14, temperature=_T, eps=4.0)
    assert np.median(np.abs(w[np.triu_indices(12, 1)])) > 1.0
    pk = [3.67] * 4 + [4.25] * 4 + [10.4] * 4
    q0 = [-1.0] * 8 + [0.0] * 4
    rings = [np.arange(0, 4), np.arange(4, 8), np.arange(8, 12)]
    for ph in (4.0, 7.3):
        ex = pka.exact(w, pk, q0, ph)
        mc = pka.monte_carlo(w, pk, q0, ph, sweeps=3000, seed=3, blocks=rings)
        assert np.max(np.abs(mc - ex)) < 0.01, ph
    # without the ring move the sampler sticks in one of the ring's two
    # degenerate states at pH 7.3 (0.008-0.07 over seeds 1-5, against
    # <= 0.0005 with it): this is what the block move is for
    ex = pka.exact(w, pk, q0, 7.3)
    blocked = pka.monte_carlo(w, pk, q0, 7.3, sweeps=3000, seed=3, blocks=rings)
    single = pka.monte_carlo(w, pk, q0, 7.3, sweeps=3000, seed=3)
    assert np.max(np.abs(single - ex)) > 10 * np.max(np.abs(blocked - ex))
    # coupling matters here: the acids are far from Henderson-Hasselbalch
    hh = 1.0 / (1.0 + 10.0 ** (4.0 - np.array(pk)))
    assert np.max(np.abs(pka.exact(w, pk, q0, 4.0) - hh)[:8]) > 0.2


def test_like_charges_raise_an_acid_ring_pka():
    """Four acids on a 4 A ring protonate earlier than one alone."""
    xyz = _ring(4, 4.0, 0.0)
    w = pka.interactions(xyz, 0.14, temperature=_T)
    theta = pka.exact(w, [3.67] * 4, [-1.0] * 4, 4.0)
    assert np.all(theta > 1.0 / (1.0 + 10.0 ** (4.0 - 3.67)))


def test_charge_centres_too_close_are_refused():
    with pytest.raises(ValueError):
        pka.interactions(np.zeros((2, 3)), 0.14)


def test_exact_refuses_a_large_network():
    with pytest.raises(ValueError):
        pka.exact(np.zeros((21, 21)), [4.0] * 21, [-1.0] * 21, 7.0)
