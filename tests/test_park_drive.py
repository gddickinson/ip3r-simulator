"""The park/drive receptor against its source (Cao et al. 2013/2014 code)
and against an independent route to its stationary distribution."""

import numpy as np
import pytest
from scipy.linalg import null_space

from ip3r.physics import park_drive as pd


def _full_generator(c, p):
    g = pd.ParkDriveParams()
    f = pd.ip3_functions(p, g)
    q24, q42 = pd.mode_rates(pd.gate_inf(c, f, g), f)
    q = pd.constant_generator(g)
    q[pd.C2, pd.C4] += q24
    q[pd.C2, pd.C2] -= q24
    q[pd.C4, pd.C2] += q42
    q[pd.C4, pd.C4] -= q42
    return q


@pytest.mark.parametrize("c,p", [(0.1, 0.2), (0.5, 0.2), (2.0, 1.0), (1.0, 10.0)])
def test_stationary_is_the_generator_null_space(c, p):
    pi = null_space(_full_generator(c, p).T)[:, 0]
    pi /= pi.sum()
    assert np.allclose(pd.stationary(c, p), pi, atol=1e-12)


def test_printed_ip3_functions_are_reproduced():
    # The code prints V24 = 62 + 880/(p^2+4) and a24 = 1 + 5/(p^2+0.5^2);
    # the registry holds them as base + amp x falling Hill.
    for p in (0.01, 0.2, 1.0, 10.0):
        f = pd.ip3_functions(p)
        assert f["V24"] == pytest.approx(62 + 880 / (p ** 2 + 4))
        assert f["a24"] == pytest.approx(1 + 5 / (p ** 2 + 0.5 ** 2))
        assert f["V42"] == pytest.approx(110 * p ** 2 / (p ** 2 + 0.1 ** 2))
        assert f["kn42"] == pytest.approx(0.41 + 25 * p ** 3 / (p ** 3 + 6.5 ** 3))


def test_modes_as_described():
    # "in drive mode, the steady-state open probability is around 70%";
    # "q54 is ~300 times greater than q45" (park almost never open).
    assert pd.drive_open_probability() == pytest.approx(0.70, abs=0.01)
    g = pd.ParkDriveParams()
    assert g.q54 / g.q45 == pytest.approx(300, rel=0.02)


def test_quiet_at_rest_and_drive_saturated_at_high_ip3():
    assert pd.park_fraction(0.1, 0.2) > 0.95
    assert pd.open_probability(0.1, 0.2) < 0.02
    # "For 10 uM IP3 ... the receptor is almost always in the drive mode"
    # for Ca2+ between 1 and ~50 uM: P_open sits at the drive-mode value.
    po = pd.open_probability(np.array([1.0, 5.0, 10.0]), 10.0)
    assert np.all(np.abs(po - pd.drive_open_probability()) < 0.02)


def test_bell_is_bell_shaped():
    b = pd.bell_at(0.2)
    assert b.c_half_act < b.c_peak < b.c_half_inh
    assert 0.2 < b.po_peak < 0.5


def test_no_ip3_never_drives():
    assert pd.park_fraction(0.5, 0.0) == pytest.approx(1.0)
