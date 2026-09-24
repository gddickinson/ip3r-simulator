"""SR content under the couplon: the pool, the hook, and Stern's Fig. 20."""

from dataclasses import replace

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS
from ip3r.physics import lumen
from ip3r.physics.couplon import CouplonParams, simulate_couplon, step_protocol
from ip3r.physics.ec_release import with_gating
from ip3r.physics.ryr_gating import SternParams

# Stern et al. 1997 Fig. 20, digitised on the tick positions (28-channel
# couplon, Stern's constants): mV -> SR content at 100 ms (C), and the
# release corrected for depletion over the plateau (D), pA per couplon.
FIG20_CONTENT = {0.0: 0.655, -30.0: 1.34}
FIG20_CORRECTED = {0.0: 1.2, -30.0: 0.4}


def test_conversion_by_hand():
    """4.8 couplons per um^3 and 1 pA for 1 ms: 1e-15 C / 2F = 5.18e-21 mol
    per couplon, times 4.8e15 couplons per litre = 24.9 uM."""
    lp = lumen.LumenParams()
    assert lumen.mm_per_pa_ms(lp) == pytest.approx(0.02487, rel=1e-3)


def test_content_path_closed_forms():
    lp = lumen.LumenParams(refill_tau=1e9)
    t = np.arange(0, 0.1, 1e-3)
    s = lumen.content_path(t, np.full(len(t), 1.0), lp)
    # no uptake, constant 1 pA: linear, 24.9 uM per ms
    assert s[50] == pytest.approx(lp.content - 50 * 0.02487, rel=1e-3)
    # uptake only: the deficit relaxes with tau exactly
    lp = lumen.LumenParams(refill_tau=0.05)
    s = lumen.content_path(t, np.r_[20.0, np.zeros(len(t) - 1)], lp)
    deficit = lp.content - s
    assert deficit[80] / deficit[30] == pytest.approx(np.exp(-0.05 / 0.05),
                                                      rel=1e-9)
    # constant current with uptake: S0 - k I tau (1 - exp(-t/tau)) to within
    # the one-bin splitting
    s = lumen.content_path(t, np.full(len(t), 0.5), lp)
    k = lumen.mm_per_pa_ms(lp) * 1e3
    want = lp.content - k * 0.5 * 0.05 * (1 - np.exp(-t / 0.05))
    assert np.allclose(s, want, rtol=0.02)
    # it never goes negative
    assert lumen.content_path(t, np.full(len(t), 1e3), lp).min() == 0.0


def test_full_scale_is_the_unscaled_couplon():
    """A path of ones changes nothing: the same couplon, event for event."""
    prot, pp = step_protocol(0.0), with_gating(SternParams())
    t = np.arange(0, 0.2, 1e-3)
    a = simulate_couplon(prot, 0.2, 3, pp)
    b = simulate_couplon(prot, 0.2, 3, pp, scale=(t, np.ones(len(t))))
    assert np.array_equal(a.c_open, b.c_open)
    assert np.array_equal(a.v_open, b.v_open)
    assert np.allclose(a.flux(), b.flux())


def test_empty_store_releases_nothing_and_triggers_nothing():
    """At zero content no current flows, so no V opening can raise the cleft
    Ca2+ and the C channels stay at their resting behaviour."""
    prot, pp = step_protocol(0.0), with_gating(SternParams())
    t = np.arange(0, 0.2, 1e-3)
    runs = [simulate_couplon(prot, 0.2, s, pp, scale=(t, np.zeros(len(t))))
            for s in range(10)]
    assert all(r.flux().max() == 0.0 for r in runs)
    assert max(r.c_open.max() for r in runs) <= 2
    assert max(r.v_open.max() for r in runs) > 5       # the V channels open


@pytest.fixture(scope="module")
def fig20():
    pp = lumen.fig20_couplon(with_gating(SternParams()))
    return {v: lumen.depleted_ensemble(v, pp, trials=80) for v in FIG20_CONTENT}


def test_path_is_consistent(fig20):
    lp = lumen.LumenParams()
    for d in fig20.values():
        assert d.change < lp.tolerance
        assert np.allclose(d.content, lumen.content_path(d.e.t, d.e.flux, lp))


def test_stern_constants_reproduce_fig20(fig20):
    """Content at 100 ms within 25 % of Stern's, and lower than his (first-
    order uptake under-refills during the pulse); the corrected plateau
    within 20 %; release stops at repolarisation."""
    for v, want in FIG20_CONTENT.items():
        d = fig20[v]
        assert d.at_ms(100) == pytest.approx(want, rel=0.25)
        assert d.at_ms(100) < want
        pulse, win = PARAMETERS.value("ec.pulse"), PARAMETERS.value("ec.plateau_window")
        plateau = (d.e.t >= pulse - win) & (d.e.t < pulse)
        corrected = float(np.nanmean(d.corrected()[plateau]))
        assert corrected == pytest.approx(FIG20_CORRECTED[v], rel=0.2)
        assert d.c_after < 0.01
    # depletion makes the raw plateau sag: the corrected record keeps
    # Stern's peak, the raw one exaggerates it
    d = fig20[0.0]
    assert d.peak_to_plateau > d.corrected_peak_to_plateau > 2.0


def test_record_grid_is_fine_enough():
    """The path is piecewise constant on the record grid; halving the bin
    moves the content at 100 ms by under 3 %."""
    pp = lumen.fig20_couplon(with_gating(SternParams()))
    coarse = lumen.depleted_ensemble(0.0, pp, trials=60)
    fine_cleft = replace(pp.cleft, record_dt=0.5 * pp.cleft.record_dt)
    fine = lumen.depleted_ensemble(0.0, replace(pp, cleft=fine_cleft), trials=60)
    assert fine.at_ms(100) == pytest.approx(coarse.at_ms(100), rel=0.03)


def test_default_couplon_is_sixty_channels():
    assert CouplonParams().cleft.n_channels == 30
    assert lumen.fig20_couplon().cleft.n_channels == 14


def test_depletion_peaks_only_by_overemptying():
    """The finding (Round 6.7). At the pool where Stern's constants release
    Rios's measured 50-60 % in a 100-ms pulse to +20 mV (twice Stern's
    2 mM), the scheme fitted to Murayama's bell releases far more, and its
    C channels keep almost none of the peak Stern's have."""
    from ip3r.physics.ec_release import configurations
    cfg = configurations()
    v = PARAMETERS.value("lumen.rios_voltage")
    lp = replace(lumen.LumenParams(), content=2 * PARAMETERS.value("lumen.content"))
    run = lambda sp: lumen.depleted_ensemble(
        v, lumen.fig20_couplon(with_gating(sp)), lp, trials=60)
    stern, fit = run(cfg["Stern 1997"]), run(cfg["fitted, no Mg2+"])
    lo = PARAMETERS.value("lumen.rios_released_low")
    hi = PARAMETERS.value("lumen.rios_released_high")
    assert lo - 0.05 < stern.released < hi + 0.05
    assert fit.released > hi + 0.15
    assert stern.po_peak_to_plateau > 3.5
    assert fit.po_peak_to_plateau < 2.0
