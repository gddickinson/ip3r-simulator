"""The couplon (V + C channels) and the release it gives under voltage clamp."""

import numpy as np
import pytest

from ip3r.physics import allosteric_v as av
from ip3r.physics import cleft, ec_release as er
from ip3r.physics.couplon import (CouplonParams, CouplonTrace, simulate_couplon,
                                  step_protocol, v_trajectory)
from ip3r.physics.ryr_gating import SternParams

# Stern et al. 1997 Fig. 12, C channels, read on the ticks' positions (the
# printed labels repeat "0.0" where 0.1 belongs): mV -> (peak, plateau).
FIG12_C = {0.0: (0.64, 0.105), -30.0: (0.39, 0.078)}


def _v_ensemble(v, n, t, rp=None):
    """Open fraction of ``n`` independent V channels at times ``t``."""
    rng = np.random.default_rng(1)
    rp = rp or av.RiosParams()
    prot = step_protocol(v, pulse=0.05)
    out = np.zeros(len(t))
    for _ in range(n):
        ts, ss = v_trajectory(prot, float(t[-1]) + 1e-9, rp, rng)
        ts, ss = np.array(ts), np.array(ss)
        out += np.array([ss[ts <= ti].sum() for ti in t])
    return out / n


def test_monte_carlo_is_the_master_equation():
    """Stern's own check (their Fig. 11): the V channels' Monte Carlo mean
    equals the integrated master equation, here through the step and after
    repolarisation."""
    t = np.array([0.003, 0.006, 0.012, 0.03, 0.049, 0.0515, 0.055])
    mc = _v_ensemble(0.0, 1500, t)
    on = av.step_response(0.0, t[t < 0.05])
    off_start = av.stationary(0.0)
    q = av.generator(-90.0)
    from scipy.linalg import expm
    off = [(off_start @ expm(q * (ti - 0.05)))[av.OPEN_STATES].sum()
           for ti in t[t >= 0.05]]
    assert np.allclose(mc, np.concatenate([on, off]), atol=0.035)


def test_v_to_c_coupling_is_reciprocal_and_scaled():
    """The Green's function is symmetric: Ca2+ at C_i from V_k equals Ca2+ at
    V_k from C_i, at equal current."""
    g = cleft.CleftGeometry.from_parameters()
    h = cleft.v_coupling_matrix(g.current, g)
    u, _ = cleft._solved(g)
    at_v = cleft._sample(g, u, cleft.couplon(g).v_sites)      # (V, C)
    assert np.allclose(h, at_v.T, rtol=0.05)
    assert np.allclose(cleft.v_coupling_matrix(0.1, g), h * 0.1 / g.current)


def test_nearest_v_gives_a_few_micromolar():
    h = cleft.v_coupling_matrix(0.1)
    top = np.sort(h, axis=1)[:, ::-1]
    assert 5.0 < np.median(top[:, 0]) < 8.0          # across the row
    assert np.median(top[:, 3]) < 0.3 * np.median(top[:, 0])


def test_no_v_current_no_trigger():
    """With V channels passing nothing, a resting array stays silent."""
    pp = CouplonParams(v_current=0.0)
    tr = simulate_couplon(step_protocol(0.0), 0.1, 0, pp)
    assert tr.v_open.max() > 5 and tr.c_open.max() <= 1


def test_stern_couplon_reproduces_fig12():
    """With Stern's constants the couplon is theirs: peak, plateau, and
    release stops at repolarisation (their 'control')."""
    for v, (peak, plateau) in FIG12_C.items():
        s = er.summarise(er.ensemble(v, er.with_gating(SternParams()), 40))
        assert s.c_peak == pytest.approx(peak, abs=0.08)
        assert s.c_plateau == pytest.approx(plateau, abs=0.04)
        assert s.c_after < 0.02
        assert s.peak_to_plateau > 2.0


def test_fitted_without_mg_loses_control():
    fit = er.configurations()["fitted, no Mg2+"]
    s = er.summarise(er.ensemble(-30.0, er.with_gating(fit), 10))
    assert s.c_after > 0.5 and s.peak_to_plateau < 1.2


def test_both_site_mg_keeps_control_but_no_peak():
    both = er.configurations()["fitted, Mg2+ both sites"]
    s = er.summarise(er.ensemble(-30.0, er.with_gating(both), 20))
    assert s.c_after < 0.005 and s.c_plateau < 0.05
    assert s.inactivated_end > 0.7


def test_events_read_by_hand():
    t = np.arange(12) * 1e-3
    c = np.array([0, 2, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0])
    tr = CouplonTrace(t, c, np.zeros(12, int), c, np.zeros(12, int),
                      CouplonParams(), ())
    ev = er.events(tr)
    assert [(round(e["start"], 3), round(e["duration"], 3), e["peak_open"])
            for e in ev] == [(0.001, 0.002, 3), (0.005, 0.001, 1)]
    assert tr.flux()[2] == pytest.approx(3 * 0.3)
