"""Ion selectivity and the unitary Ca2+ current, read with Vais 2010's ruler.

The conductance model (:mod:`ip3r.physics.unitary`) is 2.4x short of the
measured K+ conductance at best. That shortfall is an absolute scale:
diffusivity, ion radius, where the profile ends. **Selectivity is a ratio.**
The in-pore diffusivity fraction cancels from it exactly (every flux scales
by the same factor, so the reversal potential does not move). A permeability
ratio therefore tests the wall charge, which is the part of the model the
conductance could not test cleanly.

The experiments are Vais et al. 2010's lum-out nuclear patches, run through
the same drift-diffusion solver, in their own solutions:

* **P_Cl : P_K**: cytosol 140 mM KCl, lumen 30 mM KCl + 110 mM NMDG-Cl.
  NMDG+ is the impermeant substitute, so it is not a species in the pore;
  it only sets the lumen's Cl- (and, through the Donnan jump at the luminal
  mouth, the partition there).
* **P_Ca : P_K**: 140 mM KCl both sides, with 10 mM CaCl2 added in the lumen.
* **i_Ca**: 140 mM KCl both sides, 3 uM Ca2+ in the cytosol, and 0.16 /
  0.55 / 1.1 mM free Ca2+ in the lumen, at 0 mV. At 0 mV only the Ca2+
  gradient drives current, so the pore current is i_Ca. It is read before
  the access correction, because ``current`` is V / R_total, and that is
  0 at 0 V by construction. The error is a 0.1 pA current through about
  1 GOhm of access, about 0.1 mV, which is negligible.

Each reversal potential is found by root-finding on the solver's pore
current. It is converted to a ratio with :func:`ghk_ratio`, Vais's Eq. 1
(the general GHK equation), with P_Cl : P_K taken from the first experiment
exactly as they did. The model has concentrations, not activities, so the
ruler is applied to the concentrations the model saw.

**Frame.** ``z[0]`` is the luminal mouth (the channel frame has the cytosol
at +z). The applied voltage is the cytosolic bath relative to the luminal
one, which is Vais's V_app (pipette = cytoplasmic side, measured against
the bath = luminal side in a lum-out patch). Current is positive from lumen
to cytosol for a cation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from ..parameters import PARAMETERS as _P
from ._pnp_kernels import F_FARADAY, R_GAS
from .permeation import IonSpecies, solve_pnp

__all__ = ["ions", "ghk_ratio", "ghk_reversal", "reversal_potential",
           "Selectivity", "selectivity", "calcium_current", "CalciumCurrent",
           "ghk_calcium_permeability", "published", "thermal_voltage",
           "Reading", "selectivity_panel", "slow_nmdg_bound",
           "ryr1_calcium_ratio"]

_VALENCE = {"K+": 1, "Cl-": -1, "Ca2+": 2}


def thermal_voltage() -> float:
    return R_GAS * _P.value("permeation.temperature") / F_FARADAY


def ions(lumen: dict[str, float], cytosol: dict[str, float],
         cation_radius: float | None = None) -> list[IonSpecies]:
    """Species for two baths (M, keyed ``K+``, ``Cl-``, ``Ca2+``); the lumen
    is ``z[0]``. ``cation_radius`` overrides K+ and Ca2+ together (the sweep).
    Every species sees the same registered in-pore diffusivity fraction."""
    scale = _P.value("permeation.diffusion_scale")
    diff = {"K+": _P.value("permeation.diffusion_potassium"),
            "Cl-": _P.value("permeation.diffusion_chloride"),
            "Ca2+": _P.value("permeation.diffusion_calcium")}
    radius = {"K+": _P.value("permeation.radius_potassium"),
              "Cl-": _P.value("permeation.radius_chloride"),
              "Ca2+": _P.value("permeation.radius_calcium")}
    if cation_radius is not None:
        radius["K+"] = radius["Ca2+"] = cation_radius
    names = [n for n in _VALENCE if lumen.get(n, 0.0) > 0 or cytosol.get(n, 0.0) > 0]
    return [IonSpecies(n, _VALENCE[n], diff[n] * scale, radius[n],
                       lumen.get(n, 0.0), cytosol.get(n, 0.0)) for n in names]


# ---------------------------------------------------------------- the ruler
def _ghk_term(v: float, z: int, lumen: float, cytosol: float,
              thermal: float) -> float:
    """Vais's bracket, ``([Y]i - [Y]o e^{-zV/phi}) / (1 - e^{-zV/phi})``,
    with i = cytosol (pipette) and o = lumen (bath); V = V_cyt - V_lumen."""
    a = -z * v / thermal
    if abs(a) < 1e-12:
        if lumen != cytosol:
            raise ValueError("GHK term is singular at 0 V with a gradient")
        return cytosol
    return (cytosol - lumen * np.exp(a)) / (-np.expm1(a))


def ghk_ratio(v_rev: float, target: str, lumen: dict, cytosol: dict,
              known: dict[str, float]) -> float:
    """Vais 2010 Eq. 1: P_target / P_K from a reversal potential, given the
    other permeant species' ratios to K+ (``known``; K+ itself is 1)."""
    thermal = thermal_voltage()
    ratios = {"K+": 1.0, **known}
    rest = sum(ratios[x] * _VALENCE[x] ** 2
               * _ghk_term(v_rev, _VALENCE[x], lumen.get(x, 0.0),
                           cytosol.get(x, 0.0), thermal)
               for x in _VALENCE if x != target
               and (lumen.get(x, 0.0) > 0 or cytosol.get(x, 0.0) > 0))
    own = _VALENCE[target] ** 2 * _ghk_term(
        v_rev, _VALENCE[target], lumen.get(target, 0.0),
        cytosol.get(target, 0.0), thermal)
    return -rest / own


def ghk_reversal(ratios: dict[str, float], lumen: dict, cytosol: dict) -> float:
    """The inverse: where GHK's zero-current condition sits for given
    permeability ratios (used to calibrate :func:`ghk_ratio`), V."""
    thermal = thermal_voltage()

    def total(v):                       # the GHK current, up to F^2/RT
        return v * sum(p * _VALENCE[x] ** 2
                   * _ghk_term(v, _VALENCE[x], lumen.get(x, 0.0),
                               cytosol.get(x, 0.0), thermal)
                   for x, p in ratios.items())
    return brentq(total, -0.2 + 1e-9, 0.2, xtol=1e-12)


# ------------------------------------------------------------- the model
def reversal_potential(z_A, radius_A, species, fixed_charge=None,
                       bracket: float = 0.1,
                       closure: str = "donnan") -> tuple[float, bool]:
    """Voltage at which the model's pore current is zero, V; and whether
    every solve on the way converged."""
    ok = [True]

    def current(v):
        r = solve_pnp(z_A, radius_A, voltage=v, species=species,
                      fixed_charge=fixed_charge, closure=closure)
        ok[0] &= r.converged
        return r.pore_current
    lo, hi = -bracket, bracket
    if current(lo) * current(hi) > 0:
        raise ValueError("no reversal within +/- %.0f mV" % (bracket * 1e3))
    return brentq(current, lo, hi, xtol=1e-6), ok[0]


def _conditions():
    kcl = _P.value("permeation.bath_concentration")
    dilute = _P.value("selectivity.kcl_dilute")
    cacl2 = _P.value("selectivity.cacl2_lumen")
    cyt = {"K+": kcl, "Cl-": kcl}
    kcl_lum = {"K+": dilute, "Cl-": dilute + _P.value("selectivity.nmdg_cl")}
    ca_lum = {"K+": kcl, "Cl-": kcl + 2.0 * cacl2, "Ca2+": cacl2}
    return cyt, kcl_lum, ca_lum


@dataclass
class Selectivity:
    """Reversal potentials (V) and the ratios Vais's Eq. 1 reads from them."""

    label: str
    v_kcl: float
    v_ca: float
    pcl_pk: float
    pca_pk: float
    converged: bool

    def row(self) -> str:
        flag = "" if self.converged else "  (n.c.)"
        return (f"{self.label:8s} V_rev {self.v_kcl * 1e3:+6.1f} / "
                f"{self.v_ca * 1e3:+6.1f} mV   P_Cl:P_K {self.pcl_pk:5.2f}   "
                f"P_Ca:P_K {self.pca_pk:6.2f}{flag}")


def selectivity(z_A, radius_A, fixed_charge=None, label: str = "",
                cation_radius: float | None = None,
                closure: str = "donnan") -> Selectivity:
    """Vais's two bi-ionic experiments on one profile and wall charge."""
    cyt, kcl_lum, ca_lum = _conditions()
    v_kcl, ok1 = reversal_potential(
        z_A, radius_A, ions(kcl_lum, cyt, cation_radius), fixed_charge,
        closure=closure)
    pcl = ghk_ratio(v_kcl, "Cl-", kcl_lum, cyt, {})
    v_ca, ok2 = reversal_potential(
        z_A, radius_A, ions(ca_lum, cyt, cation_radius), fixed_charge,
        closure=closure)
    pca = ghk_ratio(v_ca, "Ca2+", ca_lum, cyt, {"Cl-": pcl})
    return Selectivity(label, v_kcl, v_ca, pcl, pca, ok1 and ok2)


def ryr1_calcium_ratio(z_A, radius_A, fixed_charge=None,
                       closure: str = "donnan") -> tuple[float, float, bool]:
    """Xu 2006's protocol on one wall: symmetric KCl (the RyR1 bath) with
    CaCl2 added on the luminal side, read with their Eq. 1, which is GHK
    with no Cl- term (``known`` P_Cl = 0). Returns (V_rev, P_Ca:P_K,
    converged)."""
    kcl = _P.value("permeation.ryr1_bath_concentration")
    ca = _P.value("selectivity.ryr1_cacl2_lumen")
    cyt = {"K+": kcl, "Cl-": kcl}
    lum = {"K+": kcl, "Cl-": kcl + 2.0 * ca, "Ca2+": ca}
    v, ok = reversal_potential(z_A, radius_A, ions(lum, cyt), fixed_charge,
                               closure=closure)
    return v, ghk_ratio(v, "Ca2+", lum, cyt, {"Cl-": 0.0}), ok


@dataclass
class CalciumCurrent:
    """i_Ca at 0 mV against the luminal Ca2+ step (A per M of gradient)."""

    label: str
    delta: np.ndarray                  # M, lumen - cytosol
    current: np.ndarray                # A, lumen -> cytosol positive
    slope: float                       # A/M, least squares through the origin
    conductance: float                 # S, symmetric KCl, same wall (for GHK)
    converged: bool
    meta: dict = field(default_factory=dict)

    @property
    def slope_pA_per_mM(self) -> float:
        return self.slope * 1e12 / 1e3


def calcium_current(z_A, radius_A, fixed_charge=None, label: str = "",
                    conductance: float = float("nan"),
                    closure: str = "donnan") -> CalciumCurrent:
    """Vais's i_Ca protocol: symmetric KCl, a luminal Ca2+ step, 0 mV."""
    kcl = _P.value("permeation.bath_concentration")
    ca_cyt = _P.value("selectivity.ca_cytosol")
    cyt = {"K+": kcl, "Cl-": kcl + 2.0 * ca_cyt, "Ca2+": ca_cyt}
    deltas, currents, ok = [], [], True
    for key in ("selectivity.ca_lumen_low", "selectivity.ca_lumen_mid",
                "selectivity.ca_lumen_high"):
        ca = _P.value(key)
        lum = {"K+": kcl, "Cl-": kcl + 2.0 * ca, "Ca2+": ca}
        r = solve_pnp(z_A, radius_A, voltage=0.0, species=ions(lum, cyt),
                      fixed_charge=fixed_charge, closure=closure)
        ok &= r.converged
        deltas.append(ca - ca_cyt)
        currents.append(r.pore_current)
    d, i = np.array(deltas), np.array(currents)
    slope = float(d @ i / (d @ d))
    return CalciumCurrent(label, d, i, slope, conductance, ok)


def ghk_calcium_permeability(conductance: float, pcl_pk: float, pca_pk: float,
                             kcl: float) -> tuple[float, float]:
    """Vais's estimate: P_K from the symmetric-KCl conductance by GHK at 0 V,
    ``g = (F^2/RT) c (P_K + P_Cl)``, then P_Ca = ratio x P_K. Returns
    (P_Ca m^3/s, the i_Ca slope it predicts, A/M = 2 F P_Ca x 1000)."""
    thermal = thermal_voltage()
    p_k = conductance * thermal / (F_FARADAY * kcl * 1000.0 * (1.0 + pcl_pk))
    p_ca = pca_pk * p_k
    return p_ca, 2.0 * F_FARADAY * p_ca * 1000.0


def published() -> dict[str, float]:
    """Vais 2010's measured values (ratios, and the i_Ca slope in pA/mM)."""
    return {"pcl_pk": _P.value("selectivity.published_pcl_pk"),
            "pca_pk": _P.value("selectivity.published_pca_pk"),
            "ica_slope": _P.value("selectivity.published_ica_slope"),
            "pca_measured": _P.value("selectivity.published_pca_slope"),
            "conductance": _P.value("permeation.published_itpr3_dt40")}


# ------------------------------------------------------------- the panel
@dataclass
class Reading:
    """One wall-charge reading of one deposit, through all three protocols."""

    label: str
    net_charge: float
    conductance: float                 # S, symmetric KCl
    selectivity: Selectivity
    calcium: CalciumCurrent
    ghk_slope: float                   # A/M: i_Ca GHK predicts from g and ratios

    def row(self) -> str:
        return (f"{self.selectivity.row()}   g {self.conductance * 1e12:6.1f} pS"
                f"   i_Ca {self.calcium.slope_pA_per_mM:+.4f} pA/mM "
                f"(GHK from g: {self.ghk_slope * 1e9:+.4f})"
                f"   wall {self.net_charge:+.0f} e")


def _reading(label, z, radius, fixed, net, closure="donnan") -> Reading:
    g = solve_pnp(z, radius, fixed_charge=fixed, closure=closure).conductance
    sel = selectivity(z, radius, fixed, label, closure=closure)
    ca = calcium_current(z, radius, fixed, label, conductance=g,
                         closure=closure)
    _, ghk = ghk_calcium_permeability(
        g, max(sel.pcl_pk, 0.0), max(sel.pca_pk, 0.0),
        _P.value("permeation.bath_concentration"))
    return Reading(label, net, g, sel, ca, ghk)


def selectivity_panel(st, summary=None,
                      closure: str = "donnan") -> list[Reading]:
    """A deposit under every wall-charge reading: none (neutral), its lining
    charges (charged), less salt bridges (paired), and acidic only (every
    lining base neutralised, which is where the charged reading's Ca2+
    barriers turned out to be). ``closure`` is the charged slice's closure
    (:data:`ip3r.physics.permeation.CLOSURES`)."""
    from ..structure.channel import measure_channel
    from .pore_charge import pore_charge
    from .unitary import permeation_profile
    summary = summary or measure_channel(st)
    prof = permeation_profile(st, summary)
    radius = np.maximum(prof.r_free, 0.0)
    charged = pore_charge(st, summary.frame, prof)
    paired = pore_charge(st, summary.frame, prof, pair_bridges=True)
    bases = frozenset(g.res_seq for g in charged.groups if g.charge > 0)
    acidic = pore_charge(st, summary.frame, prof, neutralise=bases)
    rows = [_reading("neutral", prof.z, radius, None, 0.0, closure)]
    for label, ch in (("charged", charged), ("paired", paired),
                      ("acidic", acidic)):
        rows.append(_reading(label, prof.z, radius, ch.density, ch.net_charge,
                             closure))
    return rows


def slow_nmdg_bound(z_A, radius_A, fraction: float | None = None) -> float:
    """P_Cl : P_K of an uncharged pore with NMDG+ *inside* it at a small
    mobility, rather than excluded at the mouth: the other bound on how the
    impermeant substitute is treated (the truth, exclusion at the
    constriction, lies between)."""
    fraction = (_P.value("selectivity.nmdg_slow_fraction") if fraction is None
                else fraction)
    cyt, kcl_lum, _ = _conditions()
    sp = ions(kcl_lum, cyt)
    k = next(s for s in sp if s.name == "K+")
    nmdg = IonSpecies("NMDG+", 1, k.diffusivity * fraction, k.radius,
                      _P.value("selectivity.nmdg_cl"), 0.0)
    v, _ = reversal_potential(z_A, radius_A, sp + [nmdg])
    return ghk_ratio(v, "Cl-", kcl_lum, cyt, {})
