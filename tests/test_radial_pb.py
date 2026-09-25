"""The radial Poisson-Boltzmann closure (``physics/radial_pb.py``).

Calibrated against what it must reduce to without the code under test:
local Donnan where the slice is narrow against the Debye length, the
linearised (Debye-Hueckel) cylinder with its Bessel-function closed form,
and Gauss's law, which the discretisation must hold exactly. Then the whole
solver: the Donnan limit is the Donnan solve, a neutral pore does not see
the closure, and the grid is converged. Last, the finding on 8TKF.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import i0, i1

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import radial_pb
from ip3r.physics import selectivity as sel
from ip3r.physics._pnp_kernels import F_FARADAY, R_GAS, _donnan_potential
from ip3r.physics.permeation import solve_pnp
from conftest import needs_structure

_TH = R_GAS * 296.0 / F_FARADAY
_KCL = np.array([[140.0] * 4, [140.0] * 4])          # mol/m^3, K+ and Cl-
_X = np.array([-10000.0, -500.0, 300.0, 3000.0])      # mol/m^3
_R = np.array([4.0, 6.0, 10.0, 14.0]) * 1e-10


def _partition(**kw):
    return radial_pb.radial_partition([1, -1], _KCL, _X, _R, _TH, **kw)


def test_narrow_slice_limit_is_local_donnan():
    """R << lambda_D (here by a very polarisable medium): the radial
    potential is flat and every species' offset is the Donnan potential."""
    donnan = _donnan_potential([1, -1], _KCL, _X, _TH)
    p = _partition(permittivity=1e5)
    assert np.allclose(p.offset, donnan[None, :], atol=1e-4 * _TH)
    assert np.all(np.abs(p.wall - p.axis) < 0.01 * _TH)


def test_linear_regime_matches_the_bessel_closed_form():
    """A weak wall on a 1:1 reservoir: psi = A I0(kappa r) with
    A kappa I1(kappa R) = F X R / (2 eps)."""
    eps = _P.value("permeation.permittivity_pore") * radial_pb.EPS0
    x = np.array([-1.0, 0.5])                          # mol/m^3, weak
    radius = np.array([6.0, 12.0]) * 1e-10
    res = np.array([[140.0, 140.0], [140.0, 140.0]])
    p = radial_pb.radial_partition([1, -1], res, x, radius, _TH, cells=256)
    kappa = np.sqrt(2 * 140.0 * F_FARADAY / (eps * _TH))
    amp = F_FARADAY * x * radius / (2 * eps * kappa * i1(kappa * radius))
    # psi is held at cell centres: the outermost is half a cell in.
    assert p.wall == pytest.approx(amp * i0(kappa * radius * p.rho[-1]),
                                   rel=1e-3)
    assert p.axis == pytest.approx(amp * i0(kappa * radius * p.rho[0]),
                                   rel=1e-3)


def test_gauss_law_holds_exactly_in_every_slice():
    p = _partition()
    gamma = np.exp(-np.array([1, -1])[:, None] * p.offset / _TH)
    net = (np.array([1, -1])[:, None] * _KCL * gamma).sum(axis=0) + _X
    assert np.all(np.abs(net) <= 1e-8 * np.abs(_X))


def test_divalents_gather_at_a_negative_wall_more_than_donnan_says():
    """What only a radial solution can show: at a wide negative wall Ca2+
    sits deeper in the wall's potential than K+ does."""
    res = np.array([[140.0], [10.0], [160.0]])        # K+, Ca2+, Cl-
    p = radial_pb.radial_partition([1, 2, -1], res, np.array([-2000.0]),
                                   np.array([10e-10]), _TH)
    k, ca, cl = p.offset[:, 0]
    assert ca < k < 0 and cl > k


def test_grid_is_converged():
    coarse = _partition()
    fine = _partition(cells=4 * int(_P.value("permeation.radial_cells")))
    assert np.max(np.abs(coarse.offset - fine.offset)) < 1e-4


# -------------------------------------------------------- through the solver
_ZS = np.linspace(0.0, 40.0, 81)
_RS = np.full_like(_ZS, 4.0)
_XS = np.where(np.abs(_ZS - 20.0) < 5.0, -2000.0, 0.0)


def test_solver_in_the_donnan_limit_is_the_donnan_solve():
    _P.set_value("permeation.permittivity_pore", _P.get(
        "permeation.permittivity_pore").maximum)
    try:
        a = solve_pnp(_ZS, _RS, fixed_charge=_XS, closure="radial")
    finally:
        _P.reset("permeation.permittivity_pore")
    b = solve_pnp(_ZS, _RS, fixed_charge=_XS)
    # eps = 80 is not infinite: R/lambda_D is still ~0.5 here.
    assert a.converged and a.meta["closure"] == "radial"
    assert a.conductance == pytest.approx(b.conductance, rel=0.02)


def test_neutral_pore_does_not_see_the_closure():
    a = solve_pnp(_ZS, _RS, closure="radial")
    b = solve_pnp(_ZS, _RS)
    assert a.conductance == b.conductance
    with pytest.raises(ValueError):
        solve_pnp(_ZS, _RS, fixed_charge=_XS, closure="poisson")


def test_radial_screening_raises_a_wide_charged_pore_toward_neutral():
    """Screening can only weaken the wall: the radial answer lies between
    Donnan's and the uncharged pore's."""
    radius = np.full_like(_ZS, 10.0)
    fixed = np.where(np.abs(_ZS - 20.0) < 3.0, 3000.0, 0.0)
    g0 = solve_pnp(_ZS, radius).conductance
    gd = solve_pnp(_ZS, radius, fixed_charge=fixed).conductance
    gr = solve_pnp(_ZS, radius, fixed_charge=fixed, closure="radial").conductance
    assert gd < gr < g0


# ------------------------------------------------------------ the finding
@needs_structure("8TKF")
def test_8tkf_vestibule_screening_is_not_the_calcium_barrier():
    from ip3r.io import loader
    st = loader.load("8TKF")
    donnan = {r.label: r for r in sel.selectivity_panel(st)}
    radial = {r.label: r for r in sel.selectivity_panel(st, closure="radial")}
    pub = sel.published()
    assert all(r.selectivity.converged for r in radial.values())
    # The closure lifts the charged reading, but not within 100x of 15.2.
    assert radial["charged"].selectivity.pca_pk > donnan["charged"].selectivity.pca_pk
    assert radial["charged"].selectivity.pca_pk < pub["pca_pk"] / 100
    # The acidic rings sit in narrow slices: the closure leaves them alone.
    assert radial["acidic"].selectivity.pca_pk == pytest.approx(
        donnan["acidic"].selectivity.pca_pk, abs=0.02)
    # Still anion-tight.
    assert radial["charged"].selectivity.pcl_pk < pub["pcl_pk"] / 5
