"""The use-dependent inactivation gate: its tables, its null space against a
hand solution, and the non-equilibrium claim the module rests on."""

from dataclasses import replace

import numpy as np
import pytest

from ip3r.physics.bell import measure_bell
from ip3r.physics.ryr_gating import (SternParams, fit_to_bell, murayama_bell,
                                     stationary)
from ip3r.physics.ryr_two_site import fit_two_site
from ip3r.parameters import PARAMETERS as _P
from ip3r.physics.ryr_use import (OPEN, STATES, UseParams, cycle_flux,
                                  fit_with_use,
                                  no_ca_gate, open_probability_no_ca_gate,
                                  residual_bound, use_rate_at, with_use)

#: The solves below are held to each other relatively, not absolutely: the
#: rates span the activation on rate's c^2, so at 1 mM the generator's
#: condition number costs about eight digits.
CAS = (0.1, 1.0, 5.0, 30.0, 100.0)


def test_state_code_and_tables():
    """``dest`` flips exactly the gate's own bit, and the masks follow."""
    sp = UseParams()
    for s in range(8):
        a, i, u = s % 2, (s // 2) % 2, s // 4
        assert sp.dest[s, 0] == (1 - a) + 2 * i + 4 * u
        assert sp.dest[s, 1] == a + 2 * (1 - i) + 4 * u
        assert sp.dest[s, 2] == a + 2 * i + 4 * (1 - u)
        assert sp.open_mask[s] == (a == 1 and i == 0 and u == 0)
        assert sp.inact_mask[s] == (i == 1 or u == 1)
    assert STATES[OPEN] == "O"
    assert sp.open_mask.sum() == 1


def test_trigger_opens_only_available_channels():
    sp = UseParams()
    for s in range(8):
        want = OPEN if not sp.inact_mask[s] else s
        assert sp.trigger_map[s] == want


def test_use_gate_is_entered_only_from_the_open_state():
    """The whole mechanism: the entry rate is zero unless the channel
    conducts, and it does not depend on Ca2+."""
    sp = UseParams()
    states = np.arange(8)
    for c in CAS:
        r_use = sp.exit_rates(states, np.full(8, c)).reshape(3, 8)[2]
        for s in range(8):
            if s // 4 == 1:
                assert r_use[s] == pytest.approx(sp.k_use_off)
            elif s == OPEN:
                assert r_use[s] == pytest.approx(sp.k_use_on)
            else:
                assert r_use[s] == 0.0


def test_null_space_equals_the_hand_solution_without_the_ca_gate():
    """Removing the Ca2+ gate leaves C, O, OU and CU, whose stationary open
    probability was solved by hand in the module."""
    sp = no_ca_gate(with_use(fit_to_bell()))
    for c in CAS:
        pi = stationary(c, sp)
        assert pi[2] + pi[3] < 1e-12                      # CI, I unreachable
        assert pi[OPEN] == pytest.approx(open_probability_no_ca_gate(c, sp),
                                         rel=1e-8)


def test_removing_the_use_gate_recovers_sterns_product():
    """With ``k_use_on`` 0 the scheme must be Stern's exactly, so the third
    gate cannot have changed the other two."""
    base = SternParams()
    sp = with_use(base, k_use_on=0.0, k_use_off=0.1)
    for c in CAS:
        assert sp.open_probability(c) == pytest.approx(base.open_probability(c),
                                                       rel=1e-8)
        assert stationary(c, sp)[4:].sum() < 1e-8


def test_removing_the_ca_gate_leaves_the_activation_gate_alone():
    """``no_ca_gate`` on a plain Stern scheme is the activation gate only:
    the Ki-infinite limit that ``SternParams.open_probability`` must give
    without a 0/0."""
    base = SternParams()
    sp = no_ca_gate(base)
    assert sp.k_i == float("inf")
    for c in CAS:
        want = c ** 2 / (c ** 2 + base.k_a ** 2)
        assert sp.open_probability(c) == pytest.approx(want, rel=1e-12)
        assert stationary(c, sp)[OPEN] == pytest.approx(want, rel=1e-8)


def test_the_cycle_carries_net_flux_and_is_edge_consistent():
    """The non-equilibrium claim: a positive flux around
    C -> O -> OU -> CU -> C, equal across every edge of the isolated cycle.
    """
    sp = no_ca_gate(with_use(SternParams()))
    for c in CAS:
        pi, q = stationary(c, sp), sp.generator(c)
        j = cycle_flux(c, sp)
        assert j > 0
        for src, dst in ((0, 1), (1, 5), (5, 4)):
            net = pi[src] * q[src, dst] - pi[dst] * q[dst, src]
            assert net == pytest.approx(j, rel=1e-6)


def test_detailed_balance_holds_once_the_gate_can_be_entered_when_shut():
    """Calibration of the test above: let the use gate be entered from C as
    well, and the same cycle falls into detailed balance. The net flux is
    therefore diagnostic of the mechanism, not of the topology."""

    class Reversible(UseParams):
        def exit_rates(self, state, c):
            r = super().exit_rates(state, c).reshape(3, -1)
            r[2] = np.where(state // 4 == 1, self.k_use_off, self.k_use_on)
            return r.ravel()

    sp = no_ca_gate(Reversible())
    for c in CAS:
        pi, q = stationary(c, sp), sp.generator(c)
        for src, dst in ((0, 1), (1, 5), (5, 4), (4, 0)):
            gross = pi[src] * q[src, dst]
            net = gross - pi[dst] * q[dst, src]
            assert abs(net) < 1e-8 * gross


def test_the_use_gate_caps_the_open_probability_at_its_available_share():
    """Held open, the channel is available only ``k_use-/(k_use + k_use-)``
    of the time, so no Ca2+ can push the peak above it."""
    sp = no_ca_gate(with_use(SternParams()))
    assert sp.open_probability(1e3) == pytest.approx(sp.available, rel=1e-3)
    assert all(sp.open_probability(c) < sp.available for c in CAS)
    fast = no_ca_gate(with_use(SternParams(), k_use_on=1e4, k_use_off=0.1))
    assert fast.open_probability(1e3) < 1e-3


def test_ca_inactivation_relieves_the_use_gate():
    """The mechanism behind everything below. The use gate is entered only
    while conducting, so the share of channels it holds down tracks how much
    the channel conducts: it *rises* with Ca2+ up to the bell's peak and then
    falls away, because shutting the channel by the Ca2+ route protects it
    from the use gate. The two inactivations partly cancel."""
    sp = with_use(fit_to_bell())
    held = {c: stationary(c, sp)[4:].sum() for c in (5.0, 30.0, 300.0, 3000.0)}
    assert held[30.0] > held[5.0]                              # rising flank
    assert held[30.0] > held[300.0] > held[3000.0]             # then relieved
    assert held[3000.0] < 0.4 * held[30.0]


def test_the_use_gate_widens_the_bell_rather_than_hiding_in_it():
    """Measured here, and it is the round's finding. At the registered
    recovery ratio (0.2) the use gate crushes the peak and pushes the inhibitory flank
    out, because Ca2+ inactivation relieves it: the bell widens from 1.86 to
    2.9 decades. An equilibrium bell is therefore not blind to this gate."""
    fit = fit_to_bell()
    plain = measure_bell(fit.open_probability)
    used = measure_bell(with_use(fit).open_probability)
    assert plain.width_decades == pytest.approx(1.86, abs=0.02)
    assert used.po_peak < 0.2 * plain.po_peak
    assert used.c_half_act < 0.5 * plain.c_half_act
    assert used.c_half_inh > 4 * plain.c_half_inh
    assert used.width_decades == pytest.approx(2.90, abs=0.02)


def test_the_bell_alone_does_not_bound_the_ratio_from_below():
    """Round 6.9 reported that no Ca2+ gate fits the bell beside a use gate
    faster than tau ~0.3 s (at recovery 0.1 s^-1: ratio ~0.03), and read the
    edge as agreement with Laver & Lamb's 1-3 s. It was a stalled solver.
    Beside a gate that barely recovers, the Ca2+ gate *fitted without it*
    does lose its falling flank -- but refitted beside it, a gate exists all
    the way down (Ka up, Ki down by decades) and lands on both points."""
    fit = fit_to_bell()
    assert not np.isfinite(
        measure_bell(with_use(fit, ratio=0.02).open_probability).c_half_inh)
    target = murayama_bell()
    kis = []
    for rho in (0.02, 0.005, 0.002):
        sp = fit_with_use(ratio=rho)
        got = measure_bell(sp.open_probability)
        assert got.c_half_act == pytest.approx(target.c_half_act, rel=1e-5)
        assert got.c_half_inh == pytest.approx(target.c_half_inh, rel=1e-5)
        kis.append(sp.k_i)
    assert kis == sorted(kis, reverse=True) and kis[-1] < 1.0


def test_the_measured_bell_is_not_counted_twice():
    """``fit_with_use`` puts the Ca2+ gate and the use gate together against
    the measured half-points, so where a fit exists it really lands on them
    -- unlike fitting the Ca2+ gate alone and adding the gate afterwards."""
    target = murayama_bell()
    slow = fit_with_use(1.0 / 30.0)                # a gate slow enough to fit
    got = measure_bell(slow.open_probability)
    assert got.c_half_act == pytest.approx(target.c_half_act, rel=1e-5)
    assert got.c_half_inh == pytest.approx(target.c_half_inh, rel=1e-5)
    # and it is a different Ca2+ gate from the one fitted without it
    assert abs(slow.k_i / fit_to_bell().k_i - 1.0) > 1e-3


def test_with_use_refuses_the_two_site_gate():
    """Its second inactivation site is not carried, so it must not be
    silently dropped."""
    with pytest.raises(TypeError, match="one-site"):
        with_use(fit_two_site())


def test_registered_rates_are_the_defaults():
    sp = UseParams()
    assert sp.tau_use == pytest.approx(1.0 / sp.k_use_on)
    assert replace(sp, k_use_on=0.0).tau_use == float("inf")


def test_the_fitted_ki_tracks_the_recovery_ratio():
    """Round 6.9's load-bearing coupling, and why ``ryr.use_recovery_ratio``
    is not a free knob: the faster the use gate recovers relative to its
    entry, the less of the measured bell's descending limb it accounts for,
    so the more the refitted Ca2+ gate has to -- and Ki climbs back towards
    the 249 µM of the fit that ignores the use gate."""
    kis = [fit_with_use(ratio=r).k_i for r in (0.1, 0.32, 1.0, 3.2)]
    assert kis == sorted(kis)                       # rises with recovery
    assert kis[0] < 30.0                            # slow recovery: near Stern
    assert kis[-1] > 150.0                          # fast: back towards 249


def test_only_the_ratio_enters_the_fit():
    """Round 6.10's reduction. Both use-gate rates are seconds and the Ca2+
    gate's are sub-millisecond, so the stationary state -- and so the fitted
    Ka and Ki -- depends on the ratio alone. A 7x slower gate at the same
    ratio gives the same fit; the same gate at another ratio does not."""
    fast, slow = (fit_with_use(on, ratio=0.32) for on in (0.316, 0.045))
    assert slow.k_a == pytest.approx(fast.k_a, rel=1e-3)
    assert slow.k_i == pytest.approx(fast.k_i, rel=1e-3)
    other = fit_with_use(0.045, ratio=1.0)
    assert abs(other.k_i / fast.k_i - 1.0) > 0.5


def test_with_use_keeps_the_ratio_when_only_the_speed_is_given():
    sp = with_use(SternParams(), k_use_on=0.3)
    assert sp.k_use_off == pytest.approx(0.3 * _P.value("ryr.use_recovery_ratio"))
    assert with_use(SternParams(), 0.3, 0.01).k_use_off == 0.01


def test_the_reading_of_fig4_agrees_with_the_papers_own_numbers():
    """The registered rate is Fig. 4's intercept at 0 mV. Its line must
    give the abstract's tau 1-3 s at +40 mV, and its slope the text's
    z delta 1.14 +- 0.25 -- which it does only if the axis is log10 (a
    natural-log reading gives tau 8 s and z delta 0.51)."""
    assert 1.0 < 1.0 / use_rate_at(0.040) < 3.0
    kT_e = 1.380649e-23 * 295.0 / 1.602176634e-19        # V, room temperature
    z_delta = _P.value("ryr.use_rate_slope") * np.log(10.0) * kT_e
    assert abs(z_delta - 1.14) < 0.25
    assert use_rate_at(0.0) == _P.value("ryr.k_use_on")


def test_the_residual_bound_and_where_the_registered_ratio_sits():
    """R/(1-R) by hand, and the unverified ratio inside the span Fig. 8
    allows at +40 mV (the only potential at which the paper bounds it)."""
    assert residual_bound(0.5) == pytest.approx(1.0)
    lo, hi = (residual_bound(_P.value(f"ryr.use_residual_40mv_{k}"))
              for k in ("min", "max"))
    assert lo == pytest.approx(0.031, abs=0.001)
    assert hi == pytest.approx(0.515, abs=0.001)
    assert lo < _P.value("ryr.use_recovery_ratio") < hi


def test_a_stalled_start_is_not_reported_as_no_fit():
    """At ratio 0.18 a single fsolve start stalled between two ratios that
    fitted, and the stall was reported as "no Ca2+ gate reproduces the
    bell". The fit must be found and must land on the measured points."""
    target = murayama_bell()
    got = measure_bell(fit_with_use(0.056, ratio=0.18).open_probability)
    assert got.c_half_act == pytest.approx(target.c_half_act, rel=1e-5)
    assert got.c_half_inh == pytest.approx(target.c_half_inh, rel=1e-5)
