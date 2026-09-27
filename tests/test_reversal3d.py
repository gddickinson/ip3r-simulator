"""Round 7.23: the steady 3-D Poisson–Nernst–Planck solve and the bi-ionic
reversal read from it, calibrated before it is read.

At v = 0 between identical baths the solve is Round 7.11's Poisson–
Boltzmann equilibrium exactly, and a small v gives Round 7.19's linear-
response conductances. On a tube (anything uniform in-plane is 1-D there)
the reversal reaches Planck's junction when uncharged, Teorell–Meyer–
Sievers when charged, and the Donnan jump an excluded NMDG⁺ sets, each
closer as the Debye length shrinks against the tube.
"""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import pnp3d, reversal3d
from ip3r.physics.charge3d import poisson_boltzmann, slice_density
from ip3r.physics.csc_readings import _mixed
from ip3r.physics.ohmic3d import geometric_conductance
from ip3r.physics.permeation import IonSpecies
from ip3r.physics.pore_charge import ChargedGroup
from ip3r.physics.selectivity import _conditions, thermal_voltage

from conftest import needs_structure
from test_charged3d import _tube
from test_selectivity import _planck, _tms

_DP, _DM = 2.0e-9, 0.5e-9


def _salt(lumen, cytosol):
    return [IonSpecies("K+", 1, _DP, 1.4, lumen, cytosol),
            IonSpecies("Cl-", -1, _DM, 1.4, lumen, cytosol)]


def _dom(vol, species):
    return pnp3d.domain(vol, {s.name: vol for s in species})


def _band(vol):
    groups = [ChargedGroup(i, "ASP", "A", z, 0.0, 0.0, -1.0)
              for i, z in enumerate(np.linspace(-8, 8, 5))]
    return slice_density(vol, groups)


# ------------------------------------------------------------ the two limits
def test_equilibrium_is_poisson_boltzmann_exactly():
    vol = _tube(n=6, length=60.0)
    x, sp = _band(vol), _mixed(0.25, 0.01)
    s = pnp3d.steady_state(_dom(vol, sp), sp, 0.0, x)
    u, _, _ = poisson_boltzmann(vol, x, sp)
    assert s.converged
    assert np.abs(s.u - u).max() < 1e-10
    assert all(f == 0.0 for f in s.flux.values())


def test_small_voltage_is_linear_response():
    """Each species' flux at 0.1 mV = D g c z v / (kT/e), g Round 7.19's
    Boltzmann-weighted Laplace conductance in the equilibrium field."""
    vol = _tube(n=6, length=60.0)
    x, sp = _band(vol), _mixed(0.25, 0.01)
    dom = _dom(vol, sp)
    eq = pnp3d.steady_state(dom, sp, 0.0, x)
    v = 1e-4
    s = pnp3d.steady_state(dom, sp, v, x, initial=eq)
    for sp_i in sp:
        g = geometric_conductance(vol, energy=sp_i.valence * eq.u).g
        linear = (sp_i.diffusivity * g * sp_i.concentration * 1000.0
                  * sp_i.valence * v / thermal_voltage())
        assert s.flux[sp_i.name] == pytest.approx(linear, rel=2e-5)


# ------------------------------------------------------ the reversal on a tube
def _reversal(vol, sp, fixed=None, impermeant=()):
    fixed = np.zeros(vol.mask.shape) if fixed is None else fixed
    v, s, _ = reversal3d.reversal(_dom(vol, sp), sp, fixed,
                                  impermeant=impermeant)
    assert s.converged
    assert abs(s.current) < 1e-3 * max(abs(c) for c in s.currents.values())
    return v


def test_uncharged_tube_reaches_plancks_junction():
    """Planck's junction needs the gradient long against the Debye length:
    0.03/0.14 M over 60 A misses by 0.6 mV, 0.3/1.4 M over 160 A by 0.02."""
    want = _planck(_DP, _DM, 0.03, 0.14)
    short = _reversal(_tube(n=4, length=60.0), _salt(0.03, 0.14))
    long_ = _reversal(_tube(n=4, length=160.0), _salt(0.3, 1.4))
    assert abs(long_ - want) < 5e-5
    assert abs(long_ - want) < abs(short - want) / 10


def _charged(x, length):
    vol = _tube(n=4, length=length)
    fixed = np.zeros(vol.mask.shape)
    fixed[:, :, np.abs(vol.zs) < length / 2 - 3.0] = x
    return _reversal(vol, _salt(0.3, 1.4), fixed) - _tms(x, _DP, _DM, 0.3, 1.4)


@pytest.mark.parametrize("x", [-3000.0, 4000.0])
def test_charged_tube_reaches_tms(x):
    """A tube charged all along but 3 A at each end: Donnan at each end of
    the charge and Planck between, as the 1-D model's TMS test (the same
    ratios at ten times the salt, so the Debye length is small). The ends'
    double layers are a 1/length correction: 300 A misses by 0.1 mV (a
    cation-selective wall) and 1.0 mV (an anion-selective one, whose
    counter-ion is the slow Cl-), and doubling the length halves it."""
    miss = _charged(x, 300.0)
    assert abs(miss) < 1.2e-3
    if x > 0:
        assert abs(_charged(x, 600.0)) < 0.6 * abs(miss)


def test_excluded_nmdg_sets_the_donnan_jump():
    """Vais's substitute: lumen K+ 0.3, Cl- 1.4, NMDG+ 1.1 kept below the
    tube's luminal end. Past it K+ = Cl- = sqrt(0.3 x 1.4) (Donnan), and
    Planck runs from there, as the 1-D model's test at ten times the salt.
    The NMDG region's own polarisation is the correction: kept to 20 A of
    a 200 A tube it misses by 1.5 mV, to 5 A of 400 A by 0.2 mV."""
    vol = _tube(n=4, length=400.0)
    sp = [IonSpecies("K+", 1, _DP, 1.4, 0.3, 1.4),
          IonSpecies("Cl-", -1, _DM, 1.4, 1.4, 1.4)]
    region = vol.mask & (vol.zs[None, None, :] < -195.0)
    v = _reversal(vol, sp, impermeant=[pnp3d.Impermeant("NMDG+", 1, 1.1, region)])
    mouth = np.sqrt(0.3 * 1.4)
    jump = -thermal_voltage() * np.log(mouth / 0.3)
    assert v == pytest.approx(_planck(_DP, _DM, mouth, 1.4) + jump, abs=3e-4)


def test_no_root_is_refused():
    vol = _tube(n=4, length=60.0)
    sp = _salt(0.03, 0.14)
    with pytest.raises(ValueError, match="no reversal"):
        reversal3d.reversal(_dom(vol, sp), sp, np.zeros(vol.mask.shape),
                            bracket=0.002)


# ------------------------------------------------------------- the protocols
def test_protocols_are_the_1d_ones():
    xu, = reversal3d.experiments(True)
    kcl = _P.value("permeation.ryr1_bath_concentration")
    ca = _P.value("selectivity.ryr1_cacl2_lumen")
    assert xu.lumen == {"K+": kcl, "Cl-": kcl + 2 * ca, "Ca2+": ca}
    assert xu.cytosol == {"K+": kcl, "Cl-": kcl}
    cl, ca_exp = reversal3d.experiments(False)
    cyt, kcl_lum, ca_lum = _conditions()
    assert (cl.lumen, cl.cytosol, ca_exp.lumen, ca_exp.cytosol) == (
        kcl_lum, cyt, ca_lum, cyt)
    assert cl.impermeant == {"NMDG+": _P.value("selectivity.nmdg_cl")}
    # every reference bath is neutral among the permeant species
    for e in (xu, cl, ca_exp):
        bath = e.lumen if e.reference == "lumen" else e.cytosol
        assert sum(c * {"K+": 1, "Cl-": -1, "Ca2+": 2}[k]
                   for k, c in bath.items()) == pytest.approx(0.0, abs=1e-12)


def test_one_reference_bath_leaves_out_little():
    for ryr in (True, False):
        assert max(reversal3d.bath_gamma_difference(ryr).values()) < 0.02


# ------------------------------------------------------------- real deposits
@needs_structure("9HEO")
def test_real_9heo_reversal_does_not_rescue_p_ca():
    """RyR1 at Xu's reversal (1 Å grid). The uncharged pore obeys GHK: read
    with its own P_Cl:P_K, its reversal gives its linear-response ratio
    back (Xu's ruler, which takes P_Cl as 0, does not). The charged wall
    with the charge-space fluid reads within a fifth of linear response,
    and far below Xu's 7.0."""
    from ip3r.io import loader
    r = reversal3d.reversal_3d(loader.load("9HEO"),
                               readings=("neutral", "pb + csc"))
    assert r.converged and r.ryr
    neutral, full = r.readings["neutral"], r.readings["pb + csc"]
    assert neutral.pca_pk_model_cl == pytest.approx(r.linear["neutral"], rel=0.03)
    assert neutral.pca_pk < 0.2
    assert full.pca_pk == pytest.approx(r.linear["pb + csc"], rel=0.2)
    assert full.pca_pk < r.measured / 5
    assert full.peak_ca > 5.0                    # the filter still binds Ca2+
