"""Round 7.24: each ion on the lumen at reversal.

Calibrated on tubes, where everything is 1-D: an uncharged junction's
concentrations are Planck's straight lines and electroneutral, and each
ion's drop is its own series resistance ∫e^{E}dz accumulated from the lumen
(the claim the viewer makes: the drop rises where the ion's resistance
lies). Then 8TKF under pb + csc.
"""

from __future__ import annotations

import numpy as np
import pytest

from ip3r.physics import reversal3d
from ip3r.physics.lumen_reversal import ion_grids
from ip3r.physics.selectivity import thermal_voltage

from conftest import needs_structure
from test_charged3d import _tube
from test_reversal3d import _band, _dom, _salt


def _at_reversal(vol, sp, fixed=None):
    fixed = np.zeros(vol.mask.shape) if fixed is None else fixed
    dom = _dom(vol, sp)
    v, st, _ = reversal3d.reversal(dom, sp, fixed)
    assert st.converged
    return dom, st


def _profile(grid):
    """The in-plane mean on every plane (a tube is uniform in-plane)."""
    return np.nanmean(grid, axis=(0, 1))


def test_uncharged_junction_is_plancks_straight_lines():
    """0.3 / 1.4 M KCl over 160 A (Round 7.23's Planck case): both ions'
    concentrations are one straight line from bath to bath, to 1 % of the
    range, and equal to each other (electroneutral)."""
    vol = _tube(n=4, length=160.0)
    sp = _salt(0.3, 1.4)
    dom, st = _at_reversal(vol, sp)
    conc, drop, baths = ion_grids(dom, st, sp)
    assert baths == {"K+": (0.3, 1.4), "Cl-": (0.3, 1.4)}
    z = vol.zs
    k, cl = _profile(conc["K+"]), _profile(conc["Cl-"])
    line = np.polyval(np.polyfit(z, k, 1), z)
    assert np.abs(k - line).max() < 0.01 * (1.4 - 0.3)
    assert np.abs(k - cl).max() < 0.01 * (1.4 - 0.3)
    assert k[0] == pytest.approx(0.3, rel=0.02)
    assert k[-1] == pytest.approx(1.4, rel=0.02)


def _resistance(dom, st, s, excess=None):
    """∫ e^{E} dz from the lumen, normalised: the ion's series resistance."""
    psi = st.v / thermal_voltage() * dom.phi0 + st.u
    e = _profile(np.where(dom.elec.mask, s.valence * psi, np.nan))
    r = np.cumsum(np.exp(e)) - np.exp(e) / 2 - np.exp(e[0]) / 2
    return r / r[-1]


@pytest.mark.parametrize("name", ["K+", "Cl-"])
def test_drop_is_each_ions_series_resistance(name):
    """A charged band (Round 7.23's) between unequal baths, at reversal:
    each ion's drop is its own accumulated resistance, so a counter-ion
    (K+ in the acidic band) drops little across the band and a co-ion (Cl-)
    most of its drop there."""
    vol = _tube(n=6, length=60.0)
    sp = _salt(0.03, 0.14)
    dom, st = _at_reversal(vol, sp, _band(vol))
    _, drop, _ = ion_grids(dom, st, sp)
    s = next(x for x in sp if x.name == name)
    d = _profile(drop[name])
    assert d[0] == pytest.approx(0.0, abs=1e-9)
    assert d[-1] == pytest.approx(1.0, abs=1e-9)
    assert np.abs(d - _resistance(dom, st, s)).max() < 0.005
    band = np.abs(vol.zs) < 8.0
    across = d[band][-1] - d[band][0]
    assert across < 0.15 if name == "K+" else across > 0.6


def test_an_ion_absent_from_one_bath_still_has_a_drop():
    """Ca2+ in the lumen only (Xu's and Vais's experiment): n_cyt = 0, the
    drop still runs 0 -> 1 and its concentration falls to nothing."""
    vol = _tube(n=4, length=60.0)
    sp = next(e for e in reversal3d.experiments(False) if e.name == "Ca2+").species()
    dom, st = _at_reversal(vol, sp)
    conc, drop, baths = ion_grids(dom, st, sp)
    assert baths["Ca2+"] == (0.01, 0.0)
    d, c = _profile(drop["Ca2+"]), _profile(conc["Ca2+"])
    assert d[0] == pytest.approx(0.0, abs=1e-9) and d[-1] == pytest.approx(1.0, abs=1e-9)
    assert c[0] == pytest.approx(0.01, rel=0.05) and c[-1] < 1e-4


# ------------------------------------------------------------------ 8TKF
@needs_structure("8TKF")
def test_real_8tkf_at_reversal():
    """8TKF, Vais's Ca2+ experiment under pb + csc: the reversal is Round
    7.23's (same code path, so the same V); Ca2+ gathers in the filter far
    above its bath, yet under a tenth of its drop falls there and more at
    the gate (7 % / 31 %); every ion's drop runs 0 -> 1 across the grid."""
    from ip3r.io import loader
    from ip3r.physics.lumen_reversal import reversal_lumen
    from ip3r.structure.channel import measure_channel
    st = loader.load("8TKF")
    s = measure_channel(st)
    r = reversal_lumen(st, "pb + csc", "Ca2+", s)
    assert r.converged
    assert r.species == ["K+", "Cl-", "Ca2+"]
    assert abs(sum(r.currents.values())) < 1e-3 * max(
        abs(c) for c in r.currents.values())
    peak, z = r.peak("Ca2+")
    assert peak > 100 * r.baths["Ca2+"][0]
    assert abs(z - s.constrictions["filter"].z) < 10.0
    for ion in r.species:
        d = r.drop_3d(ion)
        assert np.all(np.isfinite(d)) and d.min() > -0.05 and d.max() < 1.05
    # the finding: the well carries little of Ca2+'s drop, the gate more
    w = 3.0
    at = {k: r.drop_across("Ca2+", c.z, w) for k, c in s.constrictions.items()}
    assert at["filter"] < 0.1 < at["gate"]
    assert f"{r.v * 1e3:+.2f} mV" in r.summary()
    with pytest.raises(ValueError, match="reading must be"):
        reversal_lumen(st, "slice", "Ca2+", s)


@needs_structure("6DQN")
def test_a_shut_deposit_is_refused_by_name():
    """6DQN admits no K+ or Cl- path (the viewer's default deposit; this
    was a bare numpy error). 8TKG admits a 1 A Ca2+ path too narrow to
    reverse within the bracket, which reversal3d refuses by name itself."""
    from ip3r.io import loader
    from ip3r.physics.lumen_reversal import reversal_lumen
    with pytest.raises(ValueError, match="6DQN is shut: no K"):
        reversal_lumen(loader.load("6DQN"), "neutral", "Ca2+")
