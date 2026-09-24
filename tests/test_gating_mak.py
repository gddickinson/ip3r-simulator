"""The Mak et al. (1998) Hill-type model: the registered constants are the
paper's, the curve has the shapes the paper reports, and the flank test that
separates it from De Young-Keizer is an instrument that can fail."""

import numpy as np
import pytest

from ip3r.physics import gating, gating_mak as mk
from ip3r.physics.bell import flank_shifts, measure_bell


def test_constants_are_the_papers():
    g = mk.MakParams()
    assert (g.p_max, g.k_act, g.h_act, g.h_inh) == (0.81, 0.21, 1.9, 3.9)
    assert (g.k_inf, g.k_ip3, g.h_ip3) == (52.0, 0.050, 4.0)


def test_k_inh_saturates_and_falls():
    # Eq. 2: ceiling K_inf at saturating IP3, half of it at K_IP3.
    assert mk.k_inh(10.0) == pytest.approx(52.0, rel=1e-6)
    assert mk.k_inh(0.050) == pytest.approx(26.0)
    # 33 nM: the Hill fit gives 8.3 uM; the paper's point was 9.5 uM.
    assert 7.5 < mk.k_inh(0.033) < 9.5
    assert mk.k_inh(0.0) == 0.0


def test_plateau_at_saturating_ip3():
    # "the open probability remained elevated (~0.8)" over 1-20 uM Ca2+.
    po = mk.open_probability(np.array([1.0, 3.0, 10.0, 20.0]), 10.0)
    assert np.all((po > 0.75) & (po <= 0.81))


def test_no_ip3_no_opening():
    assert float(mk.open_probability(0.3, 0.0)) == 0.0


def test_ip3_tunes_inhibition_alone():
    """33 nM -> 10 uM: inhibitory flank 6.2x, activating 1.016x (Mak);
    DYK moves both, 2.8x and 2.0x."""
    r = mk.compare_flanks()
    act, inh = r["Mak"]
    assert act < 1.05 and inh > 5.0
    act_d, inh_d = r["DYK"]
    assert act_d > 1.5                     # the DYK scheme fails the test


def test_flank_test_can_fail():
    """Plant IP3 dependence into K_act: the same ruler must now see the
    activating flank move."""
    g = mk.MakParams()

    def planted(p):
        h = mk.MakParams(**{**g.__dict__, "k_act": g.k_act * (1 + 0.1 / p)})
        return mk.bell_at(p, h)

    act, _ = flank_shifts(planted, 0.033, 10.0)
    assert act < 0.5


def test_bell_collapses_below_k_ip3():
    # "at 10 and 20 nM ... reductions of both the maximum Po and the range"
    b10, b20, b100 = (mk.bell_at(p) for p in (0.010, 0.020, 0.100))
    assert b10.po_peak < b20.po_peak < b100.po_peak
    assert b10.width_decades < b20.width_decades < b100.width_decades


def test_bell_ruler_agrees_with_closed_form():
    # A pure product of Hill terms has half points at K_act and K_inh when
    # they are far apart.
    f = lambda c: 1 / (1 + (1e-2 / c) ** 4) / (1 + (c / 1e2) ** 4)
    b = measure_bell(f)
    assert b.c_half_act == pytest.approx(1e-2, rel=1e-3)
    assert b.c_half_inh == pytest.approx(1e2, rel=1e-3)
    assert gating.bell_at(1.0).c_half_act < gating.bell_at(1.0).c_peak
