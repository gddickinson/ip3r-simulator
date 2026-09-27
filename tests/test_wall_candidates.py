"""Round 7.26: Round 7.25's candidate walls built for the lumen box, and
read at reversal as the search read them."""

import numpy as np
import pytest

from conftest import needs_structure
from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import wall_candidates as wc
from ip3r.physics.wall_search import RingResult, regions


@pytest.fixture(scope="module")
def pore_8tkf():
    from ip3r.io import loader
    from ip3r.physics.reversal3d import wall
    from ip3r.physics.selectivity3d import prepare
    try:
        st = loader.load("8TKF")
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"8TKF unavailable: {e}")
    pore = prepare(st, spacing=_P.value("reversal3d.spacing"))
    return pore, wall(pore).fixed


def test_span_wells_are_the_searchs(pore_8tkf):
    """The well is −d on exactly the span's lumen voxels and nothing else;
    uncharged it carries no charge, with the charge it is the deposit's
    map unchanged."""
    pore, fixed = pore_8tkf
    d = _P.value("wallsearch.gui_well_depth")
    lo, hi = regions(pore.summary)["span"]
    w = wc.candidate_wall(pore, "span well", fixed)
    zs = pore.elec.zs[None, None, :]
    inside = pore.elec.mask & (zs >= lo) & (zs <= hi)
    assert np.all(w.well[inside] == -d) and np.all(w.well[~inside] == 0.0)
    assert not np.any(w.fixed)
    assert set(w.excess) == {"Ca2+"}
    c = wc.candidate_wall(pore, "span well + charge", fixed)
    assert c.fixed is fixed and np.array_equal(c.well, w.well)
    assert "uncharged" in w.description and "deposit's charge" in c.description


def test_ring_pair_is_neutral_and_memoised(pore_8tkf, monkeypatch):
    """The ring pair takes the highest-B converged ring of the search (a
    stand-in search here), its charge sums to zero, and the choice is
    remembered until a parameter changes."""
    pore, fixed = pore_8tkf
    calls = []

    def fake(st, spacing=None, progress=None):
        calls.append(spacing)
        return [RingResult(-90.0, 1.0, 1.0, 1.0, 0.9, True),
                RingResult(-70.0, 2.0, 1.1, 0.9, 1.3, True),
                RingResult(-60.0, 2.0, 1.2, 0.9, 5.0, False)]
    monkeypatch.setattr("ip3r.physics.wall_search.ring_search", fake)
    wc._RINGS.clear()
    w = wc.candidate_wall(pore, "ring pair", fixed)
    assert "z -70.0" in w.description and "±2 e" in w.description
    assert w.well is None and w.excess is None
    assert abs(w.fixed.sum()) < 1e-6 * np.abs(w.fixed).sum()
    assert np.any(w.fixed > 0) and np.any(w.fixed < 0)
    wc.candidate_wall(pore, "ring pair", fixed)
    assert len(calls) == 1
    key = "wallsearch.gui_well_depth"
    _P.set_value(key, 5.0)
    try:
        assert not wc._RINGS
    finally:
        _P.reset(key)
    wc._RINGS.clear()


def test_refusals(pore_8tkf):
    pore, fixed = pore_8tkf
    with pytest.raises(ValueError, match="candidate must be"):
        wc.candidate_wall(pore, "gate well", fixed)


@needs_structure("9HEO")
def test_a_ryr1_deposit_is_refused():
    from ip3r.io import loader
    from ip3r.physics.lumen_reversal import reversal_lumen
    with pytest.raises(ValueError, match="IP3R deposit"):
        reversal_lumen(loader.load("9HEO"), "span well", "Ca2+")


@needs_structure("8TKF")
def test_8tkf_span_well_at_reversal_is_round_725s():
    """The lumen box's span well reverses where Round 7.25's search put it
    (+8.09 mV, peak Ca2+ 0.326 M), and its Ca2+ drop is steepest at the
    well's cytosolic edge, not in the filter: leaving the well is the
    resistance the well adds."""
    from ip3r.io import loader
    from ip3r.physics.lumen_reversal import reversal_lumen
    from ip3r.structure.channel import measure_channel
    st = loader.load("8TKF")
    s = measure_channel(st)
    r = reversal_lumen(st, "span well", "Ca2+", s)
    assert r.candidate and r.converged
    assert r.v * 1e3 == pytest.approx(8.09, abs=0.01)
    assert r.peak("Ca2+")[0] == pytest.approx(0.326, abs=0.002)
    hi = regions(s)["span"][1]
    assert abs(r.steepest_z("Ca2+") - hi) < 3.0
    assert "span well (Ca2+-only well of" in r.summary()
