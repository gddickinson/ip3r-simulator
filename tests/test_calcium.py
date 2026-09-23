import numpy as np

from ip3r.physics.calcium import (CellParams, Trace, fluxes, oscillation_metrics,
                                  simulate)


def test_total_calcium_is_conserved():
    tr = simulate(0.5, t_end=50)
    cp = CellParams()
    assert np.allclose(tr.c + cp.c1 * tr.c_er, cp.c0, atol=1e-9)


def test_oscillation_window_edges():
    assert oscillation_metrics(simulate(0.5, t_end=200))["oscillates"]
    assert not oscillation_metrics(simulate(0.2, t_end=200))["oscillates"]
    assert not oscillation_metrics(simulate(1.0, t_end=200))["oscillates"]


def test_period_is_seconds_scale():
    m = oscillation_metrics(simulate(0.5, t_end=200))
    assert 5 < m["period"] < 30


def test_damped_spiral_is_not_an_oscillation():
    """The criterion must be able to say no to a decaying ring-down."""
    t = np.linspace(0, 100, 2001)
    c = 0.3 + 0.2 * np.exp(-t / 15) * np.sin(2 * np.pi * t / 10)
    tr = Trace(t, c, np.zeros_like(t), 0.5, np.zeros_like(t), np.zeros_like(t))
    assert not oscillation_metrics(tr)["oscillates"]
    sustained = Trace(t, 0.3 + 0.2 * np.sin(2 * np.pi * t / 10), *(np.zeros_like(t),),
                      0.5, np.zeros_like(t), np.zeros_like(t))
    m = oscillation_metrics(sustained)
    assert m["oscillates"] and abs(m["period"] - 10) < 0.2


def test_fluxes_balance_at_rest():
    f = fluxes(0.1, 0.9, 0.0)
    assert f["channel"] == 0.0                     # no IP3, no channel flux
