"""The park/drive cluster simulator: the step splitting reproduces the
stationary distribution, and the same test fails at a step too coarse."""

import numpy as np
import pytest

from ip3r.physics import park_drive as pd
from ip3r.physics.puffs_pd import ParkDrivePuffParams, simulate_cluster_pd


def _clamped_open_fraction(c, p, dt=None, n=400, duration=4.0, seed=1):
    pp = ParkDrivePuffParams()
    pp.n_channels = n
    if dt is not None:
        pp.dt = dt
    tr = simulate_cluster_pd(p, duration, seed, pp, clamp_ca=c)
    burn = len(tr.n_open) // 4            # start is all-parked, not stationary
    return tr.n_open[burn:].mean() / n


@pytest.mark.parametrize("c,p", [(0.5, 0.2), (1.0, 1.0)])
def test_clamped_cluster_reproduces_stationary_open_probability(c, p):
    expected = float(pd.open_probability(c, p))
    assert _clamped_open_fraction(c, p) == pytest.approx(expected, rel=0.03)


def test_a_coarse_step_is_caught():
    # The error grows with the step (measured at 0.5 uM Ca2+, 0.2 uM IP3:
    # 0.2 % at 0.1 ms, 1.6 % at 0.5 ms, 4.5 % at 2 ms, 10.5 % at 5 ms). At
    # 5 ms the same comparison must fail, or the test above proves nothing.
    c, p = 0.5, 0.2
    expected = float(pd.open_probability(c, p))
    got = _clamped_open_fraction(c, p, dt=5e-3, n=100, duration=100.0)
    assert abs(got / expected - 1) > 0.03


def test_no_ip3_only_park_flickers():
    pp = ParkDrivePuffParams()
    pp.n_channels = 100
    tr = simulate_cluster_pd(0.0, 2.0, 0, pp)
    g = pd.ParkDriveParams()
    park_po = g.q45 / (g.q45 + g.q54)
    assert tr.n_open.mean() / 100 == pytest.approx(park_po, rel=0.25)


def test_seeded_and_reproducible():
    a = simulate_cluster_pd(0.2, 1.0, seed=4)
    b = simulate_cluster_pd(0.2, 1.0, seed=4)
    assert np.array_equal(a.n_open, b.n_open) and np.array_equal(a.n_peak, b.n_peak)


def test_peak_bounds_snapshot():
    tr = simulate_cluster_pd(0.2, 2.0, seed=2)
    assert np.all(tr.n_peak >= tr.n_open)
