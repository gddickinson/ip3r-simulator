"""What ends a cleft spark: the refit, its bookkeeping and the finding.

* The fit lands on both measured flanks at 25 and 37 C, keeps Stern's on
  rates, and recovers Stern's own constants from Stern's own bell.
* Scaling the inactivation rates leaves the steady state untouched, which
  is why a bell cannot fix them and the rate has to be scanned.
* A spark still running at the trace's end is flagged, and one that ended
  is not (a hand-made trace).
* The finding: Stern's scheme ends its sparks in tens of ms; fitted to the
  25 C bell, sparks never end; fitted to 37 C, none start; and duration
  grows with Ki.
"""

import numpy as np
import pytest

from ip3r.physics import ryr_gating as rg
from ip3r.physics import spark_termination as st
from ip3r.physics.puff_compare import SPARK_FIT, params_for, spark_ends
from ip3r.physics.puffs import PuffTrace


@pytest.mark.parametrize("mp", [rg.MurayamaParams(), rg.MurayamaParams.at_37()])
def test_the_fit_lands_on_both_measured_flanks_and_keeps_the_on_rates(mp):
    target = rg.murayama_bell(mp)
    fit = rg.fit_to_bell(target)
    b = rg.bell_at(fit)
    assert b.c_half_act == pytest.approx(target.c_half_act, rel=1e-5)
    assert b.c_half_inh == pytest.approx(target.c_half_inh, rel=1e-5)
    sp = rg.SternParams()
    assert (fit.k_act_on, fit.k_inact_on) == (sp.k_act_on, sp.k_inact_on)


def test_the_fit_recovers_sterns_constants_from_sterns_own_bell():
    sp = rg.SternParams()
    far = rg.with_constants(20.0, 300.0, sp)             # start somewhere else
    fit = rg.fit_to_bell(rg.bell_at(sp), far)
    assert fit.k_a == pytest.approx(sp.k_a, rel=1e-5)
    assert fit.k_i == pytest.approx(sp.k_i, rel=1e-5)


def test_the_measured_bell_asks_for_a_25_fold_weaker_inactivation():
    fit = rg.fit_to_bell()
    assert 4.0 < fit.k_a < 6.0                           # activation barely moves
    assert fit.k_i / rg.SternParams().k_i > 20           # inactivation 25x weaker


def test_rate_scaling_leaves_the_steady_state_alone():
    sp = rg.SternParams()
    fast = rg.with_constants(sp.k_a, sp.k_i, sp, rate=30.0)
    assert fast.k_inact_on == pytest.approx(30 * sp.k_inact_on)
    for c in (1.0, 30.0, 300.0):
        assert rg.stationary(c, fast) == pytest.approx(rg.stationary(c, sp), rel=1e-9)


def _trace(peaks):
    peaks = np.asarray(peaks)
    t = np.arange(len(peaks)) * 1e-3
    return PuffTrace(t, peaks, np.zeros_like(t), params_for(SPARK_FIT), 0.0,
                     peaks, np.zeros_like(peaks))


def test_a_spark_running_at_the_end_is_flagged_and_an_ended_one_is_not():
    ended = spark_ends(_trace([0, 20, 25, 10, 0, 0]))
    running = spark_ends(_trace([0, 0, 20, 25, 28, 29]))
    assert [e["unterminated"] for e in ended] == [False]
    assert [e["unterminated"] for e in running] == [True]


def test_stern_ends_its_sparks_and_the_bell_fitted_scheme_does_not():
    stern, fit25, fit37 = st.refit(duration=4.0, seeds=2)
    assert stern.n_sparks > 0 and stern.unterminated == 0
    assert stern.duration_ms < 40
    assert fit25.n_sparks > 0 and fit25.unterminated == fit25.n_sparks
    assert fit25.open_fraction > 0.3                     # release never stops
    assert fit37.n_sparks == 0                           # and at 37 C never starts


def test_spark_duration_grows_with_ki():
    rows = st.ki_scan(duration=4.0, seeds=2, kis=[10.0, 60.0, 250.0])
    d = [r.duration_ms for r in rows]
    assert d[0] < d[1] < d[2]
    assert rows[0].unterminated == 0 and rows[2].unterminated == rows[2].n_sparks


def test_no_inactivation_rate_ends_the_fitted_sparks():
    rows = st.rate_scan(rg.fit_to_bell(), duration=3.0, seeds=1, rates=[0.1, 10.0])
    assert all(r.n_sparks and r.unterminated == r.n_sparks for r in rows)
