"""Round 6.13: the use gate at Rios & Pizarro 2026's millisecond speed.

Each claim is held against its bilayer-speed control, so a test that passes
for the fast gate must fail for the slow one."""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics.bell import measure_bell
from ip3r.physics.ryr_gating import fit_to_bell
from ip3r.physics.ryr_use import fit_with_use, no_ca_gate, with_use
from ip3r.physics.spark_mg import triggered
from ip3r.physics.use_speed import rios_rates, rios_ratio, speed_values

TRIALS = 6


@pytest.fixture(scope="module")
def fits():
    on, off = rios_rates()
    rho = off / on
    slow = _P.value("ryr.k_use_on")
    return {"slow": fit_with_use(slow, k_use_off=rho * slow),
            "fast": fit_with_use(on, k_use_off=off)}


def test_rates_are_the_tables_half_times():
    """Table 1's I1_ht and R1_ht are half times: k = ln2 / t."""
    on, off = rios_rates()
    assert on == pytest.approx(np.log(2) / 3.5e-3)
    assert off == pytest.approx(np.log(2) / 20e-3)
    assert rios_ratio() == pytest.approx(0.175)
    assert on / _P.value("ryr.k_use_on") > 3000


def test_speed_grid_ends_on_both_rates():
    k = speed_values()
    assert k[0] == pytest.approx(_P.value("ryr.k_use_on"))
    assert k[-1] == pytest.approx(rios_rates()[0])
    assert np.all(np.diff(k) > 0)


def test_the_bell_still_barely_sees_the_speed(fits):
    """Round 6.10's reduction at the bell survives to ~10 %: a gate 3500x
    faster at the same ratio moves the fitted Ki by more than the 1e-3 the
    bilayer-speed test allows, but by less than 20 %."""
    slow, fast = fits["slow"], fits["fast"]
    for k in ("k_a", "k_i"):
        rel = abs(getattr(fast, k) / getattr(slow, k) - 1.0)
        assert 0.005 < rel < 0.2, k


def test_rios_speed_ends_a_spark_without_the_ca_gate(fits):
    """The round's result. With the Ca2+ gate removed from each fit, the
    bilayer-speed gate leaves a triggered array open for hundreds of ms (or
    for good); at Rios's speed the use gate alone ends it on the measured
    6.3 ms scale."""
    slow = triggered(no_ca_gate(fits["slow"]), "slow", TRIALS)
    fast = triggered(no_ca_gate(fits["fast"]), "fast", TRIALS)
    assert fast.ended == TRIALS
    assert fast.duration_ms < 3.0 * _P.value("spark.published_release_duration")
    assert slow.ended < TRIALS or slow.duration_ms > 100.0


def test_rios_alone_leaves_the_bell_no_descending_limb():
    """Rios's claim that no Ca2+ role is needed, on Murayama's measurement:
    his gate beside the fitted activation gate caps P_open at the available
    share and never falls to half, so the bell's inhibition is not his
    gate's to explain."""
    on, off = rios_rates()
    sp = no_ca_gate(with_use(fit_to_bell(), on, off))
    b = measure_bell(sp.open_probability)
    assert np.isnan(b.c_half_inh)
    assert b.po_peak <= sp.available + 1e-9
    assert sp.available == pytest.approx(off / (on + off))


def test_speed_rescues_the_mixed_cluster():
    """Round 6.11: with the measured 0.8 of channels carrying the gate the
    population bell forces a weak Ca2+ gate and sparks lasted ~200 ms. At
    Rios's speed the carriers end the spark themselves."""
    on, off = rios_rates()
    rho = off / on
    slow_k = _P.value("ryr.k_use_on")
    slow = fit_with_use(slow_k, k_use_off=rho * slow_k, fraction=0.8)
    fast = fit_with_use(on, k_use_off=off, fraction=0.8)
    ts, tf = (triggered(s, n, TRIALS) for s, n in ((slow, "s"), (fast, "f")))
    assert tf.ended == TRIALS and tf.duration_ms < 40.0
    assert ts.ended < TRIALS or ts.duration_ms > 100.0


def test_a_fast_gate_suppresses_ignition(fits):
    """The speed acts on the way in as well as the way out: channels that
    inactivate within ms of opening cannot recruit, so spontaneous sparks
    become rarer as the gate speeds up at a fixed ratio."""
    from ip3r.physics.spark_termination import measure
    slow = measure(fits["slow"], "slow", 10.0, 2)
    fast = measure(fits["fast"], "fast", 10.0, 2)
    assert slow.n_sparks >= 4
    assert fast.n_sparks <= slow.n_sparks // 3


def test_a_row_without_a_use_gate_prints():
    from ip3r.physics.spark_mg import Triggered
    from ip3r.physics.spark_termination import Termination
    from ip3r.physics.use_speed import SpeedRow
    t = Termination("x", 1.0, 1.0, 1.0, 0, 0.0, float("nan"), 0,
                    float("nan"), 0.0)
    g = Triggered("x", 0.0, 1.0, 1, 1.0, 0, float("nan"), 0.0, float("nan"))
    assert "no gate" in SpeedRow("x", 0.0, 0.0, t, g).row()
