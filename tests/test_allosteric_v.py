"""The V channel (Rios et al. 1993) against its closed forms and against the
figures of the two papers that used it."""

from dataclasses import replace

import numpy as np
import pytest

from ip3r.physics import allosteric_v as av

# Stern et al. 1997 Fig. 11 (their stand-alone V check, step to 0 mV),
# digitised: ms after the step -> open probability. The printed axis labels
# repeat "0.0" where 0.1 belongs; the ticks are read by position.
FIG11_T = np.array([1, 2, 3, 5, 7.5, 10, 15, 20, 30]) * 1e-3
FIG11_PO = np.array([0.01, 0.026, 0.049, 0.115, 0.22, 0.292, 0.385, 0.424, 0.45])
# Stern Fig. 12 (the couplon runs), V at 0 mV, where the C curve is below it.
FIG12_T = np.array([10.1, 11.8, 13.4]) * 1e-3
FIG12_PO = np.array([0.426, 0.438, 0.447])
# Stern Fig. 12 V plateaus at 0, -30, -50 mV.
FIG12_PLATEAU = {0.0: 0.457, -30.0: 0.052, -50.0: 0.003}


@pytest.mark.parametrize("v", [-90.0, -50.0, -20.0, 0.0, 30.0])
def test_null_space_is_mwc(v):
    p = av.stationary(v)
    assert p[av.OPEN_STATES].sum() == pytest.approx(av.open_probability(v), rel=1e-9)
    moved = (np.arange(10) % 5) @ p / av.N_SENSORS
    assert moved == pytest.approx(av.charge(v), rel=1e-9)


def test_detailed_balance():
    q, p = av.generator(-35.0), av.stationary(-35.0)
    flow = p[:, None] * q
    np.fill_diagonal(flow, 0.0)
    assert np.allclose(flow, flow.T, atol=1e-12)


def test_saturating_po_is_the_papers():
    """'the simulated open probability reaches a maximum of 0.72 when all
    voltage sensors are in the activating position' (Rios 1993, fiber 827)."""
    rp = av.RiosParams()
    assert 1.0 / (1.0 + rp.big_l * rp.f ** 8) == pytest.approx(0.72, abs=0.005)
    assert av.open_probability(200.0) == pytest.approx(0.716, abs=0.002)


def test_rate_scale_leaves_equilibria():
    fast = av.RiosParams(rate_scale=7.0)
    assert np.allclose(av.stationary(-25.0, fast), av.stationary(-25.0))


def test_fiber_827_gives_sterns_plateaus():
    for v, po in FIG12_PLATEAU.items():
        assert av.open_probability(v) == pytest.approx(po, rel=0.1)


def _rms(rp, t, ref):
    return float(np.sqrt(np.mean((av.step_response(0.0, t, rp) - ref) ** 2)))


def test_fig11_follows_the_printed_rates_fig12_twice_them():
    """Stern's text doubles every rate. Their couplon figure follows that;
    their stand-alone check of the model does not."""
    printed, doubled = av.RiosParams(rate_scale=1.0), av.RiosParams(rate_scale=2.0)
    assert _rms(printed, FIG11_T, FIG11_PO) < 0.03
    assert _rms(doubled, FIG11_T, FIG11_PO) > 0.06
    assert _rms(doubled, FIG12_T, FIG12_PO) < 0.02
    assert _rms(printed, FIG12_T, FIG12_PO) > 0.1
    assert av.RiosParams().rate_scale == 2.0


def test_repolarisation_closes_the_channel():
    po = av.step_response(-90.0, [0.0, 0.01], holding=0.0)
    assert po[0] == pytest.approx(av.open_probability(0.0))
    assert po[1] < 0.01


def test_exits_table_matches_generator():
    dest, rate = av.exits(-40.0)
    q = av.generator(-40.0)
    for i in range(10):
        assert rate[i].sum() == pytest.approx(-q[i, i])
        for d, r in zip(dest[i], rate[i]):
            if r > 0:
                assert q[i, d] == pytest.approx(r)


def test_each_sensor_favours_opening_by_f_squared():
    rp = replace(av.RiosParams(), rate_scale=1.0)
    q = av.generator(-60.0, rp)
    for j in range(5):
        eq = q[j, 5 + j] / q[5 + j, j]
        assert eq == pytest.approx(rp.f ** (-2 * j) / rp.big_l)
