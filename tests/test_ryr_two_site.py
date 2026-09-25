"""The two-site inactivation gate: its states, its steady state, its fit."""

import numpy as np
import pytest

from ip3r.physics import ryr_gating as rg
from ip3r.physics import ryr_two_site as ts


@pytest.fixture(scope="module")
def fit():
    return ts.fit_two_site()


def test_exits_follow_the_state_code():
    # s = a + 2n: the activation gate flips a, binding adds 2, unbinding
    # takes 2 away; an impossible move is a self-loop (its rate is 0).
    for s in range(6):
        n = s // 2
        act, bind, unbind = ts.TwoSiteParams.dest[s]
        assert act == s ^ 1
        assert bind == (s + 2 if n < 2 else s)
        assert unbind == (s - 2 if n > 0 else s)
    rates = ts.two_site(5.0, 300.0, 100.0).exit_rates(np.arange(6), np.full(6, 20.0))
    bind, unbind = rates[6:12], rates[12:]
    assert np.all(bind[4:] == 0) and np.all(unbind[:2] == 0)


@pytest.mark.parametrize("mg", [0.0, 1000.0])
@pytest.mark.parametrize("c", [0.1, 5.0, 40.0, 300.0, 3000.0])
def test_null_space_is_the_closed_form(c, mg):
    sp = rg.with_mg(ts.two_site(5.0, 300.0, 100.0), mg, k_mg_a=700.0)
    pi = rg.stationary(c, sp)
    assert pi[sp.open_mask].sum() == pytest.approx(rg.open_probability(c, sp),
                                                   rel=1e-9)
    assert pi[sp.inact_mask].sum() == pytest.approx(
        1.0 - sp.uninactivated(c), rel=1e-9, abs=1e-15)


def test_slope_ruler_on_a_hill_function():
    k, n = 30.0, 1.7
    f = lambda c: c ** n / (c ** n + k ** n)
    assert ts.bell_slope(f, k) == pytest.approx(n / 2, rel=1e-6)


@pytest.mark.parametrize("ratio,lo,hi", [(1e3, 1.0, 1.05), (1e-3, 1.95, 2.0)])
def test_k2_over_k1_sets_the_gate_slope_between_one_and_two(ratio, lo, hi):
    # Slope at the half point is 2 - c_h/(K1 + c_h): K2 >> K1 puts c_h ~ K2
    # far above K1 (slope 1); K2 << K1 puts it at sqrt(K1 K2) far below.
    assert lo <= ts.gate_hill_slope(ts.two_site(5.0, 100.0, 100.0 * ratio)) <= hi


def test_fit_meets_all_three_targets(fit):
    m, b = rg.murayama_bell(), rg.bell_at(fit)
    assert b.c_half_act == pytest.approx(m.c_half_act, rel=1e-6)
    assert b.c_half_inh == pytest.approx(m.c_half_inh, rel=1e-6)
    want = ts.bell_slope(rg.murayama_activity, m.c_half_inh)
    got = ts.bell_slope(lambda c: rg.open_probability(c, fit), b.c_half_inh)
    assert got == pytest.approx(want, rel=1e-5)
    # The one-site fit meets the points but not the slope.
    one = rg.fit_to_bell()
    one_slope = ts.bell_slope(lambda c: rg.open_probability(c, one),
                              rg.bell_at(one).c_half_inh)
    assert abs(one_slope / want - 1) > 0.3
    assert fit.mg == 0.0 and 1.0 < ts.gate_hill_slope(fit) < 2.0


def test_steeper_gate_inactivates_less_at_cleft_calcium(fit):
    # Same half point, steeper flank: below it, less inactivation. Where a
    # cleft channel lives (tens of µM) the two-site gate is the weaker one.
    one = rg.fit_to_bell()
    for c in (20.0, 50.0):
        assert 1 - fit.uninactivated(c) < 0.5 * (1 - one.k_i / (c + one.k_i))


def test_trigger_opens_every_uninactivated_closed_state():
    t = ts.TwoSiteParams.trigger_map
    sp = ts.TwoSiteParams
    assert np.all(sp.open_mask[t[~sp.inact_mask]])
    assert np.array_equal(t[sp.inact_mask], np.flatnonzero(sp.inact_mask))
    assert np.array_equal(t[sp.open_mask], np.flatnonzero(sp.open_mask))


def test_uncoupled_array_samples_its_stationary_state(monkeypatch, fit):
    # The simulators read the scheme's own tables: six states, exact
    # Gillespie, independent channels at a clamped 150 µM (where the gate
    # inactivates ~20 %).
    from ip3r.physics import sparks_cleft as sc
    monkeypatch.setattr(sc, "couplings_for", lambda pp: np.zeros((30, 30)))
    tr = sc.simulate_sparks_cleft(0.0, 20.0, 3,
                                  sc.CleftSparkParams(ca_rest=150.0, gating=fit))
    pi = rg.stationary(150.0, fit)
    assert tr.n_open.mean() / 30 == pytest.approx(pi[fit.open_mask].sum(), rel=0.05)
    assert tr.n_inactivated.mean() / 30 == pytest.approx(
        pi[fit.inact_mask].sum(), rel=0.1)


def test_two_site_couplon_loses_control_further(fit):
    from ip3r.physics import ec_release as er
    one = er.configurations()["fitted, no Mg2+"]
    after = [er.summarise(er.ensemble(-30.0, er.with_gating(sp), 10)).c_after
             for sp in (one, fit)]
    assert after[1] > after[0] > 0.5
