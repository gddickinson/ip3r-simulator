"""The puff microdomain (Round 7.2): the pools balance by hand, the steady
single-channel rise is the closed form, the dye saturates where it must,
the cluster reproduces the receptor's stationary P_open when every receptor
sees the same Ca2+, and the IPI fit recovers known constants."""

import math
from dataclasses import replace

import numpy as np
import pytest

from ip3r.physics import microdomain as md
from ip3r.physics import park_drive as pdm
from ip3r.physics import puff_stats as ps
from ip3r.physics.puffs_domain import DomainPuffParams, simulate_cluster_domain

P = 0.2


def _settle(y, n_open, d, t=2.0, dt=1e-4):
    for _ in range(int(t / dt)):
        y = md.rk4_step(y, n_open, P, d, dt)
    return y


def test_start_store_is_the_codes():
    d = md.DomainParams()
    assert md.start_store(d) == pytest.approx(10 * (45 - 0.1 / 100 - 0.1))


@pytest.mark.parametrize("clamp", md.CLAMPS)
def test_rest_is_steady(clamp):
    d = md.DomainParams()
    rest = md.rest_state(P, d)
    dd = d if clamp == "none" else md.clamped(d, rest, clamp)
    dy = md.rhs(rest.y, 0, P, dd)
    # c, cb, b are steady; a free store's total drifts at Jin - Jpm only.
    assert max(abs(dy[0]), abs(dy[1]), abs(dy[3])) < 1e-10
    if clamp != "none":
        assert abs(dy[2]) < 1e-10


def test_store_balance_by_hand():
    """d(cs)/dt from the bookkeeping equals gamma2 (Jserca - Jleak - Jipr):
    the bound dye must be counted in the total, or this fails."""
    d = md.DomainParams()
    y = (0.3, 4.0, 46.0, 5.0)
    n_open, h = 3, 1e-7
    dy = md.rhs(y, n_open, P, d)
    y2 = tuple(a + h * b for a, b in zip(y, dy))
    got = (md.store(y2, d) - md.store(y, d)) / h
    c, cb, cs = y[0], y[1], md.store(y, d)
    serca = d.v_serca * c ** d.n_serca / (c ** d.n_serca + d.k_serca ** d.n_serca)
    want = d.gamma2 * (serca - d.k_leak * (cs - c) - d.k_ipr * n_open * (cs - cb))
    assert got == pytest.approx(want, rel=1e-5)


def test_bath_holds_the_cytosol():
    d = md.DomainParams()
    rest = md.rest_state(P, d)
    y = _settle(tuple(rest.y), 5, md.clamped(d, rest, "bath"), t=0.2)
    assert y[0] == pytest.approx(rest.c, abs=1e-12)


def test_one_open_receptor_steady_rise_and_fluorescence():
    d = md.DomainParams()
    rest = md.rest_state(P, d)
    dd = md.clamped(d, rest, "bath")
    y = _settle(tuple(rest.y), 1, dd, t=0.5)
    cb = (d.k_diff * rest.c + d.k_ipr * rest.cs) / (d.k_diff + d.k_ipr)
    assert y[1] == pytest.approx(cb, rel=1e-6)
    kd = d.fluo_koff / d.fluo_kon
    assert y[3] == pytest.approx(d.fluo_total * cb / (cb + kd), rel=1e-6)
    assert md.fluorescence(y[3], rest) - 1 == pytest.approx(ps.blip_df(rest, d),
                                                            rel=1e-6)
    # The mean-field model's coupling, to first order in k_ipr / k_diff.
    assert cb - rest.c == pytest.approx(md.coupling_per_open(rest, d), rel=1e-3)


def _steady_rise(n, scale):
    d = md.DomainParams().scaled(scale)
    rest = md.rest_state(P, d)
    y = _settle(tuple(rest.y), n, md.clamped(d, rest, "bath"), t=0.3)
    return y[1] - rest.c, float(md.fluorescence(y[3], rest)) - 1


def test_dye_bends_only_near_its_kd():
    """Ca2+ rises linearly with the number open at any release rate; dF/F0
    is compressed more the further the microdomain passes K_d = 2 uM
    (Cao 2013 Fig. 8 against Fig. S9)."""
    squeeze = {}
    for scale in (1.0, 4.0):
        ca1, f1 = _steady_rise(1, scale)
        ca20, f20 = _steady_rise(20, scale)
        assert ca20 / ca1 == pytest.approx(20, rel=0.05)
        squeeze[scale] = f20 / (20 * f1)
    assert 0 < squeeze[4.0] < squeeze[1.0] < 1.0
    assert squeeze[4.0] < 0.5 * squeeze[1.0]


def test_clamp_is_checked():
    d = md.DomainParams()
    with pytest.raises(ValueError):
        md.clamped(d, md.rest_state(P, d), "cytosol")


def test_cluster_samples_the_stationary_open_probability():
    """No release and the mouth set to the microdomain: every receptor sees
    rest Ca2+, open or shut, and the cluster must sample stationary P_open."""
    dp = DomainPuffParams(n_channels=400, clamp="bath")
    dp.domain = replace(dp.domain, k_ipr=0.0)
    rest = md.rest_state(P, dp.domain)
    dp.domain = replace(dp.domain, mouth_per_store=rest.c / rest.cs)
    tr = simulate_cluster_domain(P, 8.0, 1, dp)
    burn = len(tr.n_open) // 4
    got = tr.n_open[burn:].mean() / 400
    want = float(pdm.open_probability(rest.c, P))
    assert got == pytest.approx(want, rel=0.15)
    assert np.allclose(tr.ca, rest.c)


def test_seeded_and_recorded():
    dp = DomainPuffParams(n_channels=10)
    a = simulate_cluster_domain(0.1, 2.0, 3, dp)
    b = simulate_cluster_domain(0.1, 2.0, 3, dp)
    assert np.array_equal(a.n_open, b.n_open) and np.array_equal(a.f_ratio, b.f_ratio)
    assert np.all(a.peaks >= a.n_open)
    assert np.all(a.f_peak >= a.f_ratio - 1e-12)
    assert a.f_ratio[0] == pytest.approx(1.0)
    assert np.all(a.store == a.store[0])            # clamped by default


def test_free_store_depletes_under_release():
    dp = DomainPuffParams(n_channels=20, clamp="none")
    dp.domain = dp.domain.scaled(20.0)
    tr = simulate_cluster_domain(P, 5.0, 0, dp)
    assert tr.store.min() < 0.9 * tr.store[0]


# ------------------------------------------------------------------ IPIs

def test_thurley_density_is_normalised():
    t = np.linspace(0, 200, 400001)
    for lam, xi in ((0.3, 1.0), (0.2, 50.0)):
        assert np.trapezoid(ps.thurley_pdf(t, lam, xi), t) == pytest.approx(1, abs=1e-4)


def test_fit_recovers_a_refractory_process():
    rng = np.random.default_rng(0)
    t = ps.thurley_sample(0.3, 0.8, 2000, rng)
    lam, xi = ps.fit_thurley(t)
    assert lam == pytest.approx(0.3, rel=0.15)
    assert xi == pytest.approx(0.8, rel=0.2)
    assert ps.refractory_lr(t, lam, xi) > 3.84
    assert ps.cv(t) < 0.9


def test_an_exponential_process_is_not_called_refractory():
    rng = np.random.default_rng(1)
    t = rng.exponential(1 / 0.3, 2000)
    lam, xi = ps.fit_thurley(t)
    assert xi / lam > 20
    assert ps.refractory_lr(t, lam, xi) < 3.84
    assert ps.cv(t) == pytest.approx(1, abs=0.07)


def test_events_read_from_fluorescence():
    dp = DomainPuffParams(n_channels=10)
    tr = simulate_cluster_domain(0.1, 30.0, 0, dp)
    ev = ps.fluorescence_events(tr)
    assert ev and all(e["end"] >= e["start"] for e in ev)
    level = 0.5 * ps.blip_df(tr.rest, dp.domain)
    assert all(e["peak_df"] > level for e in ev)
    s = ps.puff_stats(tr)
    assert s.n_puffs <= s.n_events and len(s.ipis) == max(s.n_puffs - 1, 0)
    assert math.isfinite(s.blip_mean) and s.blip_mean > 0
