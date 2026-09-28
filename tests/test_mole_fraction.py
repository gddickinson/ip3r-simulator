"""Round 7.29: the mole-fraction sweep, calibrated on a tube.

The baths are Vais's with the CaCl₂ replaced. The grid holds Vais's
10 mM, and a sweep point at 10 mM is Round 7.27's reversal exactly.
Without a site the dilute limit is the Ca²⁺-free pore and i_Ca at 0 mV is
linear in c. Only the block stops K⁺: the same site without it leaves K⁺
alone.
"""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import mole_fraction as mf
from ip3r.physics import pnp3d

from conftest import needs_structure
from test_charged3d import _tube


def _setup():
    vol = _tube(n=4, length=40.0)
    dom = pnp3d.domain(vol, {k: vol for k in ("K+", "Cl-", "Ca2+")})
    band = vol.mask & (np.abs(vol.zs) <= 10.0)[None, None, :]
    return vol, dom, band


def _curve(reading, cas, pcl=0.0, depth=6.0):
    vol, dom, band = _setup()
    return mf.curve(dom, band, 50.0, 0.5, reading, depth, "span",
                    np.asarray(cas), -0.04, pcl)


# ------------------------------------------------------------- bookkeeping
def test_baths_are_vais_with_the_cacl2_replaced():
    kcl = _P.value("permeation.bath_concentration")
    free = mf.experiment(0.0)
    assert free.lumen == free.cytosol == {"K+": kcl, "Cl-": kcl}
    e = mf.experiment(0.003)
    assert e.lumen == {"K+": kcl, "Cl-": kcl + 0.006, "Ca2+": 0.003}
    assert e.cytosol == {"K+": kcl, "Cl-": kcl}


def test_grid_holds_vais_10_mm():
    c = mf.concentrations()
    assert len(c) == 7
    assert c[0] == pytest.approx(1e-4) and c[-1] == pytest.approx(0.1)
    assert np.isclose(c, _P.value("selectivity.cacl2_lumen")).any()
    assert np.isclose(c, 1e-3).any()


def test_half_point_by_hand():
    ca = np.array([1e-3, 1e-2, 1e-1])
    assert mf.half_point(ca, np.array([1.0, 0.75, 0.25])) == pytest.approx(
        np.sqrt(1e-2 * 1e-1))
    assert np.isnan(mf.half_point(ca, np.array([1.0, 0.9, 0.8])))


def test_readings_are_the_site_with_and_without_its_block():
    s = mf.reading_site("site", 4.0, "span")
    assert s.compensated and s.block == _P.value("casite.block")
    n = mf.reading_site("site, no block", 4.0, "span")
    assert n.compensated and n.block == 0.0 and n.depth == 4.0
    assert mf.reading_site("uncharged", 4.0, "span") is None


# ----------------------------------------------------------------- on a tube
def test_point_at_10_mm_is_round_727s_reversal():
    """The sweep's reversal at 10 mM is the reversal the 7.27 test reads
    (same baths, same coupling)."""
    from test_ca_site import _coupling, _pca
    from ip3r.physics.csc_readings import _mixed
    vol, _, _ = _setup()
    cp = _coupling(vol, _mixed(0.14, 0.01), 6.0, comp=True, block=1.0)
    p = _curve("site", [0.01]).points[0]
    assert p.converged
    assert p.pca_pk == pytest.approx(_pca(vol, cp), rel=1e-6)


def test_uncharged_dilute_limit_and_linear_i_ca():
    c = _curve("uncharged", [1e-5, 1e-4, 3e-4, 1e-3])
    assert c.k_fraction[0] == pytest.approx(1.0, abs=1e-3)
    assert c.total_fraction[0] == pytest.approx(1.0, abs=1e-3)
    per = np.array([p.i_ca / p.ca for p in c.points])
    assert np.ptp(per) / per.mean() < 0.01
    assert all(p.occupancy == 0.0 and p.held == 0.0 for p in c.points)


def test_only_the_block_stops_k():
    cas = [1e-4, 1e-3, 1e-2, 1e-1]
    block, plain = _curve("site", cas), _curve("site, no block", cas)
    assert np.all(np.diff(block.k_fraction) < 0)
    assert block.k_fraction[-1] < 0.01
    assert plain.k_fraction.min() > 0.9
    assert np.all(block.k_fraction < plain.k_fraction)
    # the site fills as c passes K_d either way
    assert np.all(np.diff([p.occupancy for p in block.points]) > 0)
    assert block.points[-1].occupancy > 0.99


def test_i_ca_slope_through_zero():
    c = _curve("uncharged", [1e-4, 1e-3])
    per = c.points[0].i_ca / c.points[0].ca
    assert c.ica_slope(1e-4, 1e-3) == pytest.approx(per, rel=0.01)
    assert np.isnan(c.ica_slope(1.0, 2.0))


# ------------------------------------------------------------- 8TKF, pinned
@needs_structure("8TKF")
def test_8tkf_site_at_10_mm_is_the_crossing():
    """At Vais's 10 mM the sweep's site is Round 7.27's crossing: P_Ca:P_K
    15.2 at casite.gui_depth, the K⁺ current partly blocked."""
    from ip3r.io import loader
    keys = ("molefrac.ca_min", "molefrac.ca_max")
    try:
        for k in keys:
            _P.set_value(k, _P.value("selectivity.cacl2_lumen"))
        m = mf.sweep(loader.load("8TKF"), readings=("site",))
    finally:
        for k in keys:
            _P.reset(k)
    c = m.curves["site"]
    (p,) = c.points
    assert p.converged
    assert p.pca_pk == pytest.approx(
        _P.value("selectivity.published_pca_pk"), rel=0.02)
    assert 0.0 < c.k_fraction[0] < 0.9
