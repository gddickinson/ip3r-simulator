"""The mixed cluster (Round 6.11): only some RyR1s carry the use gate."""

import numpy as np
import pytest

from ip3r.physics.bell import measure_bell
from ip3r.physics.ryr_gating import SternParams, murayama_bell, stationary
from ip3r.physics.ryr_mixed import bind, fit_mixed, fraction_values, mixed
from ip3r.physics.ryr_use import OPEN, fit_with_use, with_use


@pytest.fixture(scope="module")
def use():
    return with_use(SternParams())


def test_mixture_is_the_weighted_mean_by_hand(use):
    """The population bell is f Po_use + (1 - f) Po_plain, computed here
    from the two schemes' own stationary states."""
    m = mixed(use, 0.7)
    for c in (0.5, 5.0, 80.0):
        po_use = stationary(c, use)[OPEN]
        po_plain = SternParams().open_probability(c)
        assert m.open_probability(c) == pytest.approx(
            0.7 * po_use + 0.3 * po_plain, rel=1e-9)


def test_limits_are_the_known_schemes(use):
    c = np.geomspace(0.1, 1e3, 7)
    assert np.allclose(mixed(use, 1.0).open_probability(c),
                       use.open_probability(c), rtol=1e-10)
    assert np.allclose(mixed(use, 0.0).open_probability(c),
                       SternParams().open_probability(c), rtol=1e-9)


def test_fraction_bounds_refused(use):
    with pytest.raises(ValueError):
        mixed(use, 1.2)


def test_bind_counts_and_unbound_is_refused(use):
    m = mixed(use, 0.8)
    with pytest.raises(ValueError, match="bound"):
        m.exit_rates(np.zeros(5, int), np.ones(5))
    b = m.for_cluster(25, np.random.default_rng(3))
    assert b.mask.sum() == 20
    assert bind(SternParams(), 25, np.random.default_rng(0)) == SternParams()
    with pytest.raises(ValueError, match="mask for 25"):
        b.exit_rates(np.zeros(4, int), np.ones(4))


def test_only_masked_channels_can_enter_the_use_gate(use):
    """Every channel open: the use-gate block is k_use on the masked ones
    and zero on the rest; a used-up channel still recovers."""
    b = mixed(use, 0.5).for_cluster(10, np.random.default_rng(1))
    state = np.full(10, OPEN)
    r = b.exit_rates(state, np.full(10, 5.0))[20:]
    assert np.allclose(r[b.mask], use.k_use_on)
    assert np.all(r[~b.mask] == 0.0)
    state[:] = 5                                   # OU: used up, open
    assert np.allclose(b.exit_rates(state, np.full(10, 5.0))[20:],
                       use.k_use_off)


def test_draw_never_starts_an_unmasked_channel_used_up(use):
    b = mixed(use, 0.5).for_cluster(200, np.random.default_rng(2))
    s = b.draw(20.0, np.random.default_rng(4))
    assert np.all(s[~b.mask] < 4)
    assert np.any(s[b.mask] >= 4)          # at 20 µM the gate is populated


def test_fraction_grid_spans_zero_to_one():
    f = fraction_values()
    assert f[0] == 0.0 and f[-1] == 1.0


def test_population_fit_lands_on_the_measured_bell():
    target = murayama_bell()
    m = fit_mixed(0.8)
    b = measure_bell(m.open_probability)
    assert b.c_half_act == pytest.approx(target.c_half_act, rel=1e-4)
    assert b.c_half_inh == pytest.approx(target.c_half_inh, rel=1e-4)


def test_fewer_carriers_weaken_the_shared_ca_gate():
    """The bell is a population measurement, so the fewer channels carry the
    use gate, the more of the descending limb the Ca2+ gate must take back:
    Ki climbs from the all-use fit toward the one-site fit's."""
    kis = [fit_mixed(f).k_i for f in (1.0, 0.8, 0.5)]
    assert kis[0] < kis[1] < kis[2]
    assert fit_mixed(1.0).k_i == pytest.approx(fit_with_use().k_i, rel=1e-6)


def test_the_bell_not_the_carriers_lengthens_sparks():
    """Round 6.11's decomposition. The non-inactivating channels alone
    (all-use gate, 0.8 carry) leave sparks as short as all carrying; the
    Ca2+ gate the population bell forces (0.8 fit, all carry) makes them
    several times longer."""
    from ip3r.physics.spark_termination import fraction_controls, measure
    alone, gate = fraction_controls(0.8, None, 10.0, 2)
    full = measure(mixed(fit_mixed(1.0), 1.0), "all", 10.0, 2)
    assert alone.unterminated == 0 and full.unterminated == 0
    assert alone.duration_ms < 2.0 * full.duration_ms
    assert gate.duration_ms > 5.0 * full.duration_ms
