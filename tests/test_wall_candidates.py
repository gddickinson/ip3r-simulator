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


def test_ca_site_is_round_727s_crossing_site(pore_8tkf):
    """Round 7.28: the site candidate is uncharged, carries no fixed excess,
    and is 7.27's compensated blocking site over exactly the span's lumen
    at the drawn depth; its coupling spreads casite.sites over that band."""
    from ip3r.physics.ca_site import SiteCoupling, site_density
    from ip3r.physics.pnp3d import domain
    pore, fixed = pore_8tkf
    w = wc.candidate_wall(pore, "Ca2+ site", fixed)
    assert not np.any(w.fixed) and w.excess is None
    assert w.site.compensated and w.site.region == "span"
    assert w.site.depth == _P.value("casite.gui_depth")
    assert w.site.block == _P.value("casite.block")
    lo, hi = regions(pore.summary)["span"]
    zs = pore.elec.zs[None, None, :]
    assert np.array_equal(w.mask, pore.elec.mask & (zs >= lo) & (zs <= hi))
    c = w.coupling(domain(pore.elec, pore.vols))
    assert isinstance(c, SiteCoupling)
    assert c.s == site_density(w.mask, pore.spacing, _P.value("casite.sites"))
    assert wc.candidate_wall(pore, "span well", fixed).coupling(None) is None
    assert "compensated site" in w.description


def test_block_colours_are_fixed_and_grey_where_missing():
    from ip3r.render.colormaps import ramp
    from ip3r.render.lumen_mesh import block_colors
    top = _P.value("display.lumen_block_range")
    e = np.array([0.0, top / 2, top, 5 * top, np.nan])
    assert np.allclose(block_colors(e), ramp(np.array([0.0, 0.5, 1.0, 1.0,
                                                       np.nan])))


@needs_structure("8TKF")
def test_8tkf_site_at_reversal_is_round_727s():
    """The drawn site reverses where Round 7.27's root put it (+18.17 mV at
    4.4128 kT; +18.15 at the drawn 4.41), 2.01 Ca2+ held; θ is zero outside
    the span and at most 1; K+'s block is −ln(1 − fθ) voxel by voxel; the
    Cl− experiment carries no site reading (it holds no Ca2+). Six minutes."""
    from ip3r.io import loader
    from ip3r.physics.lumen_reversal import reversal_lumen
    from ip3r.structure.channel import measure_channel
    st = loader.load("8TKF")
    s = measure_channel(st)
    r = reversal_lumen(st, "Ca2+ site", "Ca2+", s)
    assert r.candidate and r.converged
    assert r.v * 1e3 == pytest.approx(18.15, abs=0.02)
    assert r.held == pytest.approx(2.01, abs=0.02)
    lo, hi = regions(s)["span"]
    zs = r.lumen.volume.zs[None, None, :]
    th = r.occupancy
    outside = np.isfinite(th) & ((zs < lo) | (zs > hi))
    assert np.all(th[outside] == 0.0) and np.nanmax(th) <= 1.0
    f = _P.value("casite.block")
    on = np.isfinite(th)
    assert np.allclose(r.k_block[on], -np.log1p(-f * th[on]))
    assert "the site holds 2.0" in r.summary()
    cl = reversal_lumen(st, "Ca2+ site", "Cl-", s)
    assert cl.occupancy is None and np.isnan(cl.held)
