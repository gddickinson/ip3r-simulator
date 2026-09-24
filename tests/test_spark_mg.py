"""Mg2+ in the RyR1 scheme, and what it does to a cleft spark.

* Mg2+ = 0 leaves the scheme exactly as it was, and the bell fit ignores
  any Mg2+ it is handed (the measured bell has none).
* With Mg2+, the closed form still equals the generator's null space. The
  activation gate's half point moves to Ka (1 + Mg/K_Mg,A), and at rest the
  inactivation gate holds Mg/(Mg + Ki) of the channels (the Ca2+-equal rule).
* The trigger opens exactly the channels that are not inactivated.
* The finding: fitted to the 25 C bell, a triggered array never shuts
  without Mg2+. Fibre Mg2+ at the activation site alone shuts every one, with
  no channel inactivated, so the spark ends by losing its feedback.
"""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as P
from ip3r.physics import ryr_gating as rg
from ip3r.physics import spark_mg as sm
from ip3r.physics.sparks_cleft import CleftSparkParams, simulate_sparks_cleft


@pytest.fixture(scope="module")
def fit():
    return rg.fit_to_bell()


def test_no_mg_is_the_scheme_unchanged():
    sp = rg.SternParams()
    for c in (0.1, 10.0, 300.0):
        assert np.array_equal(rg.generator(c, rg.with_mg(sp, 0.0)),
                              rg.generator(c, sp))


def test_the_fit_is_made_without_mg():
    sp = rg.SternParams()
    a, b = rg.fit_to_bell(sp=sp), rg.fit_to_bell(sp=rg.with_mg(sp, 1000.0))
    assert (b.mg, b.k_a, b.k_i) == (0.0, pytest.approx(a.k_a), pytest.approx(a.k_i))


@pytest.mark.parametrize("mg", [30.0, 1000.0])
def test_closed_form_matches_the_null_space_under_mg(fit, mg):
    sp = rg.with_mg(fit, mg)
    for c in (0.1, 5.0, 80.0, 900.0):
        assert rg.stationary(c, sp)[rg.OPEN] == pytest.approx(
            float(rg.open_probability(c, sp)), rel=1e-9)


def test_mg_moves_half_activation_by_the_competitive_factor(fit):
    sp = rg.with_mg(rg.with_constants(fit.k_a, 1e12, fit), 1000.0, mg_i=0.0)
    assert sp.k_a_eff == pytest.approx(fit.k_a * (1 + 1000.0 / sp.k_mg_a))
    assert float(rg.open_probability(sp.k_a_eff, sp)) == pytest.approx(0.5, rel=1e-6)


def test_mg_inactivates_as_ca_would_at_rest(fit):
    sp = rg.with_mg(fit, 1000.0)
    st = rg.stationary(1e-6, sp)
    assert st[2] + st[3] == pytest.approx(1000.0 / (1000.0 + fit.k_i), rel=1e-6)
    assert P.value("ryr.mg_i_relative") == 1.0


def test_the_trigger_opens_every_channel_not_inactivated(fit):
    sp = rg.with_mg(fit, 1000.0)
    for seed in range(3):
        tr = simulate_sparks_cleft(0.0, 0.002, seed, CleftSparkParams(gating=sp),
                                   trigger=True)
        assert tr.n_open[0] + tr.n_inactivated[0] == 30
        assert tr.n_inactivated[0] > 15          # ~80 % held by Mg2+ at rest


def test_without_mg_a_triggered_array_never_shuts(fit):
    r = sm.triggered(fit, "no Mg", trials=4, window=0.5)
    assert (r.opened, r.ended) == (30, 0)


def test_mg_at_the_activation_site_alone_ends_every_spark_uninactivated(fit):
    r = sm.triggered(rg.with_mg(fit, None, mg_i=0.0), "A only", trials=8)
    assert (r.opened, r.ended) == (30, 8)
    assert r.inactivated_end == 0.0
    assert 5.0 < r.duration_ms < 60.0


def test_the_ratio_reading_is_a_weaker_competitor(fit):
    k = sm.k_mg_a_by_ratio(fit)
    assert k == pytest.approx(fit.k_a * 54.0 / 0.51)
    assert rg.with_mg(fit, 1000.0, k).k_a_eff < rg.with_mg(fit, 1000.0).k_a_eff
