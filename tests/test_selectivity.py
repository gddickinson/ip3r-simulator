"""Selectivity and i_Ca (Vais 2010's protocols through the pore model).

The ruler (GHK Eq. 1) is calibrated by a round trip and by Vais's own
arithmetic. The solver's asymmetric-bath path, including the per-mouth
Donnan jumps and an impermeant bath ion, is calibrated against closed
forms: Planck's liquid junction, Teorell-Meyer-Sievers for a uniformly
charged pore, and the Donnan jump an excluded NMDG+ sets.
"""

from __future__ import annotations

import numpy as np
import pytest

from ip3r.physics import selectivity as sel
from ip3r.physics.permeation import F_FARADAY, IonSpecies, solve_pnp
from conftest import needs_structure

_Z = np.linspace(0.0, 40.0, 81)             # A, a uniform cylinder
_R = np.full_like(_Z, 6.0)
_CYT = {"K+": 0.14, "Cl-": 0.14}
_LUM = {"K+": 0.03, "Cl-": 0.03}


def _pair(d_plus, d_minus, lumen, cytosol, radius=1.4):
    """K+ and Cl- with chosen diffusivities and one radius (equal areas)."""
    return [IonSpecies("K+", 1, d_plus, radius, lumen, cytosol),
            IonSpecies("Cl-", -1, d_minus, radius, lumen, cytosol)]


# ------------------------------------------------------------ the ruler
def test_ghk_ratio_inverts_ghk_reversal():
    kcl_lum = {"K+": 0.03, "Cl-": 0.14}
    v = sel.ghk_reversal({"K+": 1.0, "Cl-": 0.27}, kcl_lum, _CYT)
    assert sel.ghk_ratio(v, "Cl-", kcl_lum, _CYT, {}) == pytest.approx(0.27, rel=1e-6)
    ca_lum = {"K+": 0.14, "Cl-": 0.16, "Ca2+": 0.01}
    v = sel.ghk_reversal({"K+": 1.0, "Cl-": 0.27, "Ca2+": 15.2}, ca_lum, _CYT)
    assert v > 0                         # Ca2+ into the cytosol: V_cyt must rise
    assert sel.ghk_ratio(v, "Ca2+", ca_lum, _CYT, {"Cl-": 0.27}) == pytest.approx(
        15.2, rel=1e-6)


def test_vais_arithmetic_reproduced():
    """P_Ca from the i_Ca slope is 1.5e-18 m^3/s; GHK from 545 pS and the
    ratios (at the ~104 mM activity they used) is ~1.5e-17, ten times more."""
    slope = 0.30e-12 / 1e-3                      # A/M
    assert slope / (2 * F_FARADAY * 1000.0) == pytest.approx(1.5e-18, rel=0.05)
    p_ca, _ = sel.ghk_calcium_permeability(545e-12, 0.27, 15.2, 0.104)
    assert p_ca == pytest.approx(1.5e-17, rel=0.15)


# ------------------------------------------------------------ the solver
def _planck(d_plus, d_minus, c_lum, c_cyt):
    u = (d_plus - d_minus) / (d_plus + d_minus)
    return u * sel.thermal_voltage() * np.log(c_lum / c_cyt)


def test_uncharged_pore_gives_plancks_junction():
    d_plus, d_minus = 2.0e-9, 0.5e-9
    v, ok = sel.reversal_potential(_Z, _R, _pair(d_plus, d_minus, 0.03, 0.14))
    assert ok
    assert v == pytest.approx(_planck(d_plus, d_minus, 0.03, 0.14), abs=2e-5)
    assert v < 0                         # the faster K+ leaves the lumen positive


def _tms(x, d_plus, d_minus, c_lum, c_cyt):
    """Teorell-Meyer-Sievers: Donnan at each mouth, electroneutral Planck
    inside, for a 1:1 salt and a uniform fixed charge x (mol/m^3, signed)."""
    phi = sel.thermal_voltage()

    def inside(c):                      # counter/co-ion in the pore mouth
        cp = (-x + np.sqrt(x * x + 4 * c * c)) / 2
        return cp, -phi * np.log(cp / c)
    cl_, dl = inside(c_lum * 1000)
    cr_, dr = inside(c_cyt * 1000)
    a = d_minus * x / (d_plus + d_minus)
    u = (d_plus - d_minus) / (d_plus + d_minus)
    interior = -u * phi * np.log((cr_ + a) / (cl_ + a))
    return interior + dl - dr


@pytest.mark.parametrize("x", [-20.0, -300.0, -5000.0, 400.0])
def test_uniformly_charged_pore_gives_tms(x):
    d_plus, d_minus = 2.0e-9, 0.5e-9
    fixed = np.full_like(_Z, x)
    v, ok = sel.reversal_potential(_Z, _R, _pair(d_plus, d_minus, 0.03, 0.14),
                                   fixed)
    assert ok
    assert v == pytest.approx(_tms(x, d_plus, d_minus, 0.03, 0.14), abs=5e-5)


def test_strong_negative_wall_approaches_potassium_nernst():
    v, _ = sel.reversal_potential(_Z, _R, _pair(2e-9, 2e-9, 0.03, 0.14),
                                  np.full_like(_Z, -20000.0))
    nernst = sel.thermal_voltage() * np.log(0.03 / 0.14)
    assert v == pytest.approx(nernst, abs=2e-4)


def test_impermeant_bath_ion_sets_a_mouth_donnan_jump():
    """Lumen 30 KCl + 110 NMDG-Cl with NMDG+ left out: the luminal mouth
    holds K+ = Cl- = sqrt(30 x 140), and Planck runs from there."""
    d_plus, d_minus = 2.0e-9, 0.5e-9
    sp = [IonSpecies("K+", 1, d_plus, 1.4, 0.03, 0.14),
          IonSpecies("Cl-", -1, d_minus, 1.4, 0.14, 0.14)]
    v, ok = sel.reversal_potential(_Z, _R, sp)
    assert ok
    mouth = np.sqrt(0.03 * 0.14)
    jump = -sel.thermal_voltage() * np.log(mouth / 0.03)
    assert v == pytest.approx(_planck(d_plus, d_minus, mouth, 0.14) + jump,
                              abs=5e-5)


def test_symmetric_baths_unchanged_by_the_per_mouth_donnan():
    """The fix only applies to different baths: symmetric charged current
    equals the mean-bath arithmetic it replaced (zero at 0 V, and ohmic)."""
    fixed = np.where((_Z > 10) & (_Z < 30), -2000.0, 0.0)
    sp = _pair(1.96e-9, 2.03e-9, 0.14, 0.14)
    assert abs(solve_pnp(_Z, _R, 0.0, sp, fixed).pore_current) < 1e-20
    a = solve_pnp(_Z, _R, 0.01, sp, fixed).pore_current
    b = solve_pnp(_Z, _R, -0.01, sp, fixed).pore_current
    assert a == pytest.approx(-b, rel=1e-3)


def test_reversal_is_independent_of_the_diffusivity_scale():
    """Why selectivity tests the wall and not the unmeasured diffusivity."""
    z = np.linspace(0, 60, 121)
    r = 4.0 + 3.0 * np.abs(np.sin(z / 12))
    fixed = -3000.0 * np.exp(-((z - 30) / 5) ** 2)
    lum = {"K+": 0.14, "Cl-": 0.16, "Ca2+": 0.01}
    out = []
    for scale in (1.0, 0.1):
        sp = [IonSpecies(s.name, s.valence, s.diffusivity * scale, s.radius,
                         s.concentration, s.concentration_right)
              for s in sel.ions(lum, _CYT)]
        out.append(sel.reversal_potential(z, r, sp, fixed)[0])
    assert out[0] == pytest.approx(out[1], abs=1e-7)


def test_calcium_current_flows_lumen_to_cytosol_in_an_acidic_pore():
    z = np.linspace(0, 60, 121)
    fixed = -3000.0 * np.exp(-((z - 30) / 5) ** 2)
    c = sel.calcium_current(z, np.full_like(z, 5.0), fixed)
    assert c.converged
    assert np.all(c.current > 0)
    assert np.all(np.diff(c.current) > 0)          # grows with the gradient


def test_the_instrument_can_report_calcium_selectivity():
    """The 8TKF shortfall must not be a ceiling of the method: a pore charged
    along its length is Ca2+-selective past Vais's 15.2, and the same wall
    concentrated in one ring is not (Ca2+ must cross the neutral stretches
    without Donnan enrichment)."""
    z = np.linspace(0, 60, 121)
    r = np.full_like(z, 3.5)
    tract = sel.selectivity(z, r, np.full_like(z, -30000.0))
    assert tract.pca_pk > 15.2 and tract.pcl_pk < 0.01
    ring = sel.selectivity(z, r, -30000.0 * np.exp(-((z - 30) / 5) ** 2))
    assert ring.pca_pk < 1.0


# ------------------------------------------------------------ the finding
@needs_structure("8TKF")
def test_8tkf_cannot_reach_the_measured_selectivity():
    from ip3r.io import loader
    rows = {r.label: r for r in sel.selectivity_panel(loader.load("8TKF"))}
    pub = sel.published()
    assert set(rows) == {"neutral", "charged", "paired", "acidic"}
    assert all(r.selectivity.converged and r.calcium.converged
               for r in rows.values())
    # No reading comes within a factor of 10 of P_Ca:P_K = 15.2.
    assert max(r.selectivity.pca_pk for r in rows.values()) < pub["pca_pk"] / 10
    # The lining bases are the Ca2+ barriers: with them, Ca2+ barely passes.
    assert rows["charged"].calcium.slope_pA_per_mM < 1e-3
    assert rows["acidic"].selectivity.pca_pk > 10 * max(
        rows["charged"].selectivity.pca_pk, 0.01)
    # The wall makes the pore anion-tight, where the channel passes Cl-.
    assert rows["charged"].selectivity.pcl_pk < pub["pcl_pk"] / 5
    # No reading reproduces Vais's GHK excess: the model obeys GHK
    # to within 10 %, where the channel falls ~8-10x short of it.
    a = rows["acidic"]
    assert a.calcium.slope == pytest.approx(a.ghk_slope, rel=0.10)
