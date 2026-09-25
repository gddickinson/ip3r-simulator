"""The use-dependent inactivation gate: its tables, its null space against a
hand solution, and the non-equilibrium claim the module rests on."""

from dataclasses import replace

import numpy as np
import pytest

from ip3r.physics.bell import measure_bell
from ip3r.physics.ryr_gating import (SternParams, fit_to_bell, murayama_bell,
                                     stationary)
from ip3r.physics.ryr_two_site import fit_two_site
from ip3r.physics.ryr_use import (OPEN, STATES, UseParams, cycle_flux,
                                  fit_with_use,
                                  no_ca_gate, open_probability_no_ca_gate,
                                  with_use)

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
    sp = with_use(base, k_use_on=0.0)
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
    fast = no_ca_gate(with_use(SternParams(), k_use_on=1e4))
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
    """Measured here, and it is the round's finding. At the *measured* rate
    (tau 2 s) the use gate crushes the peak and pushes the inhibitory flank
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


def test_a_use_gate_faster_than_a_second_has_no_measurable_falling_flank():
    """Faster still and the descending limb leaves the searched range
    (10 mM), 37x above the measured KI: no Ca2+ gate can put it back."""
    fit = fit_to_bell()
    assert not np.isfinite(
        measure_bell(with_use(fit, 1.0 / 0.2).open_probability).c_half_inh)
    with pytest.raises(RuntimeError, match="no Ca2\\+ gate"):
        fit_with_use(1.0 / 0.2)


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


def test_the_fitted_ki_tracks_the_recovery_rate():
    """Round 6.9's load-bearing coupling, and why ``ryr.k_use_off`` is not a
    free knob: the faster the use gate recovers, the less of the measured
    bell's descending limb it accounts for, so the more the refitted Ca2+
    gate has to -- and Ki climbs back towards the 249 µM of the fit that
    ignores the use gate."""
    kis = [fit_with_use(1.0 / 3.162, k_use_off=off).k_i
           for off in (0.03, 0.1, 0.3, 1.0)]
    assert kis == sorted(kis)                       # rises with recovery
    assert kis[0] < 30.0                            # slow recovery: near Stern
    assert kis[-1] > 150.0                          # fast: back towards 249
