"""Round 7.27: the saturable Ca²⁺ site (compensated or not) and the K⁺
block its occupancy sets, calibrated on a tube before it is read.

A hidden energy moves transport and leaves Poisson exactly as it was; at
equilibrium the site holds c(1 + s/(c + K_d)) by hand; a dilute site is
7.25's well; the solve's energies are the coupling's own fixed point; the
Cl⁻ experiment cannot see the site; and at reversal the block and the
compensation each raise P_Ca:P_K in the direction the mechanism says.
"""

import numpy as np
import pytest

from ip3r.physics import ca_site, pnp3d, reversal3d
from ip3r.physics.csc_readings import _mixed
from ip3r.physics.selectivity import ghk_ratio, ions

from conftest import needs_structure
from test_charged3d import _tube


def _dom(vol, species):
    return pnp3d.domain(vol, {s.name: vol for s in species})


def _band(vol, half=10.0):
    return vol.mask & (np.abs(vol.zs) <= half)[None, None, :]


def _coupling(vol, sp, depth, comp=False, block=0.0, sites=50.0):
    """``sites``: the site density, mol/m³ (the tube's band is ~80 Å³, so
    four sites there would be 80 M)."""
    band = _band(vol)
    s = sites
    return ca_site.SiteCoupling(_dom(vol, sp), ca_site.Site(depth, comp, block),
                                band, s)


# ------------------------------------------------------------ the hidden term
def test_hidden_energy_moves_transport_and_not_poisson():
    vol = _tube(n=4, length=40.0)
    sp = _mixed(0.14, 0.01)
    dom, zero = _dom(vol, sp), np.zeros(vol.mask.shape)
    well = {"Ca2+": -4.0 * _band(vol)}
    plain = pnp3d.steady_state(dom, sp, 0.0, zero)
    hid = pnp3d.steady_state(dom, sp, 0.0, zero, hidden=well)
    seen = pnp3d.steady_state(dom, sp, 0.0, zero, excess=well)
    assert hid.converged and seen.converged
    assert np.abs(hid.u - plain.u).max() < 1e-12       # Poisson unmoved
    assert np.abs(seen.u).max() > 0.05                 # a visible well is not
    c = reversal3d.concentration(dom, hid, "Ca2+", 2, hidden=well)
    c0 = reversal3d.concentration(dom, plain, "Ca2+", 2)
    inside = _band(vol)
    assert np.allclose(c[inside] / c0[inside], np.exp(4.0), rtol=1e-10)


# ------------------------------------------------------------- the site alone
def test_compensated_site_holds_its_langmuir_total_by_hand():
    """Uncharged tube, v = 0, identical baths: ψ = 0, so c_f is the bath's and
    the total is c (1 + s / (c + K_d)), θ = c / (c + K_d)."""
    vol = _tube(n=4, length=40.0)
    sp = _mixed(0.14, 0.01)
    cp = _coupling(vol, sp, 6.0, comp=True)
    st = pnp3d.steady_state(_dom(vol, sp), sp, 0.0, np.zeros(vol.mask.shape),
                            coupling=cp)
    assert st.converged
    assert np.abs(st.u).max() < 1e-10
    c = 10.0                                           # mol/m3
    kd = cp.s / np.expm1(6.0)
    total = reversal3d.concentration(cp.dom, st, "Ca2+", 2)
    inside = _band(vol)
    assert np.allclose(total[inside], c * (1 + cp.s / (c + kd)), rtol=1e-4)
    assert np.allclose(total[~inside & vol.mask], c, rtol=1e-8)
    theta = cp.occupancy(cp.free(st.u, st.n["Ca2+"]))
    assert np.allclose(theta[inside], c / (c + kd), rtol=1e-4)


def test_dilute_site_is_the_round_725_well():
    """A site far from full (many sites, trace Ca²⁺) is a fixed −d well."""
    vol = _tube(n=4, length=40.0)
    sp = _mixed(0.14, 1e-6)
    dom, zero = _dom(vol, sp), np.zeros(vol.mask.shape)
    cp = _coupling(vol, sp, 4.0, sites=1e5)
    site = pnp3d.steady_state(dom, sp, 0.0, zero, coupling=cp)
    well = pnp3d.steady_state(dom, sp, 0.0, zero,
                              excess={"Ca2+": -4.0 * _band(vol)})
    a = reversal3d.concentration(dom, site, "Ca2+", 2)
    b = reversal3d.concentration(dom, well, "Ca2+", 2, {"Ca2+": -4.0 * _band(vol)})
    assert np.allclose(a[vol.mask], b[vol.mask], rtol=1e-5)


def test_solved_energies_are_the_couplings_fixed_point():
    """Uncompensated, blocking, away from equilibrium: the state's energies
    are what the coupling returns for that state, within the tolerance."""
    vol = _tube(n=4, length=40.0)
    sp = ions({"K+": 0.14, "Cl-": 0.16, "Ca2+": 0.01}, {"K+": 0.14, "Cl-": 0.14})
    cp = _coupling(vol, sp, 6.0, block=1.0)
    dom = cp.dom
    st = pnp3d.steady_state(dom, sp, 0.005, np.zeros(vol.mask.shape),
                            coupling=cp)
    assert st.converged
    psi = 0.005 / reversal3d.thermal_voltage() * dom.phi0 + st.u
    vis, _ = cp(psi, st.n)
    for k in ("Ca2+", "K+"):
        assert np.abs(vis[k] - st.coupled[0][k]).max() < 2e-4
    assert st.coupled[0]["K+"].max() > 0.1             # the block is on


def test_chloride_experiment_cannot_see_the_site():
    vol = _tube(n=4, length=40.0)
    sp = ions({"K+": 0.03, "Cl-": 0.03}, {"K+": 0.14, "Cl-": 0.14})
    dom, zero = _dom(vol, sp), np.zeros(vol.mask.shape)
    cp = _coupling(vol, _mixed(0.14, 0.01), 8.0, comp=True, block=1.0)
    a = pnp3d.steady_state(dom, sp, 0.01, zero)
    b = pnp3d.steady_state(dom, sp, 0.01, zero, coupling=cp)
    assert b.flux == a.flux
    assert np.array_equal(a.u, b.u)


# ------------------------------------------------- the mechanism at reversal
def _pca(vol, cp):
    lum = {"K+": 0.14, "Cl-": 0.16, "Ca2+": 0.01}
    cyt = {"K+": 0.14, "Cl-": 0.14}
    sp = ions(lum, cyt)
    v, st, _ = reversal3d.reversal(_dom(vol, sp), sp,
                                   np.zeros(vol.mask.shape), coupling=cp)
    assert st.converged
    return ghk_ratio(v, "Ca2+", lum, cyt, {"Cl-": 0.0})


def test_block_and_compensation_each_raise_the_ratio():
    vol = _tube(n=4, length=40.0)
    sp = _mixed(0.14, 0.01)
    plain = _pca(vol, _coupling(vol, sp, 8.0))
    comp = _pca(vol, _coupling(vol, sp, 8.0, comp=True))
    block = _pca(vol, _coupling(vol, sp, 8.0, block=1.0))
    both = _pca(vol, _coupling(vol, sp, 8.0, comp=True, block=1.0))
    assert comp > plain
    assert block > plain
    assert both > max(comp, block)


def test_site_density_by_hand():
    mask = np.ones((10, 10, 10), bool)                 # 1000 Å3 at 1 Å
    s = ca_site.site_density(mask, 1.0, 4.0)
    assert s == pytest.approx(4 / (ca_site.N_AVOGADRO * 1e-27))


def test_labels_and_grid():
    assert ca_site.Site(4.0, True, 1.0).label == (
        "compensated site 4 kT at span, K+ block 1")
    d = ca_site.depths()
    assert d[0] == pytest.approx(2.0) and d[-1] == pytest.approx(12.0)


# ------------------------------------------------------------- 8TKF, pinned
@needs_structure("8TKF")
def test_8tkf_compensated_blocking_site_passes_vais_with_cl_untouched():
    """Four compensated sites over 8TKF's span at 8 kT, blocking K⁺, nearly
    full at reversal: P_Ca:P_K past 15.2 while P_Cl:P_K stays the uncharged
    pore's exactly (the Cl⁻ experiment holds no Ca²⁺)."""
    from ip3r.io import loader
    from ip3r.parameters import PARAMETERS as _P
    from ip3r.physics.selectivity3d import prepare
    from ip3r.physics.wall_search import Candidate, _Solver, regions
    pore = prepare(loader.load("8TKF"), spacing=_P.value("reversal3d.spacing"))
    solver, bands = _Solver(pore), regions(pore.summary)
    neutral = solver.read(Candidate("charge", 0.0), bands)
    r = ca_site.read(solver, ca_site.Site(8.0, True, 1.0), bands)
    assert r.converged
    assert r.pcl_pk == neutral.pcl_pk
    assert r.pca_pk > _P.value("selectivity.published_pca_pk")
    assert r.occupancy > 0.95 and 3.8 < r.bound <= 4.0
