"""P_Ca:P_K at bi-ionic reversal in 3-D (Round 7.23).

Rounds 7.19 and 7.21 read P_Ca:P_K from linear response between identical
baths and found it at or below 1.8 on every open wall, gate or no gate.
Xu 2006 (RyR1, 7.0) and Vais 2010 (ITPR3, 15.2) measured it differently:
as the voltage at which no current flows when the two baths differ,
converted to a ratio by GHK. In 1-D the two readings already differ
(9HEO csc 0.87 against 0.64). Binding in the filter can set a reversal
potential without setting the conductance ratio, so the protocol itself
is the last difference between model and measurement that Rounds 7.19
and 7.21 left untested in 3-D.

Here each protocol runs through :mod:`.pnp3d`'s steady state on the same
voxels, and the reversal is the root of its net current:

* ``xu`` (RyR1): the RyR1 KCl on both sides, CaCl₂ added in the lumen,
  read with Xu's Eq. 1 (GHK with no Cl⁻ term), as
  :func:`.selectivity.ryr1_calcium_ratio`;
* ``vais`` (IP3R): P_Cl:P_K from dilute luminal KCl with impermeant
  NMDG-Cl, then P_Ca:P_K from luminal CaCl₂ read with that P_Cl:P_K, as
  :func:`.selectivity.selectivity`. NMDG⁺ is confined to the luminal bath
  below the membrane span (the 1-D model excludes it at the luminal mouth).

Readings (:data:`READINGS`): ``neutral`` (no wall charge, Poisson still
on: the junction needs it), ``pb`` (Round 7.11's local charge placement,
point ions) and ``pb + csc`` (Round 7.19's: each species' excess from the
locally neutral charge–space fluid, held fixed). The csc reference bath is
the experiment's side that holds every permeant species and is neutral
without the impermeant one (Xu and Vais's Ca²⁺: the lumen; Vais's KCl: the
cytosol); its excess is zero outside the wall groups' reach on both sides,
so the two baths' small difference in activity coefficient is left out
(:func:`bath_gamma_difference` measures it).

Local-neutrality closures (``slice``, ``local``, their csc forms) have no
non-equilibrium form here: they are Round 7.19's linear-response readings
only.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .charge3d import group_positions, local_density
from .csc3d import local_csc, wall_fluid
from .pnp3d import Domain, Impermeant, Steady, domain, steady_state
from .pore_charge import pore_charge
from .selectivity import _conditions, ghk_ratio, ions, thermal_voltage
from .selectivity3d import Pore, measured_ratio, prepare

__all__ = ["READINGS", "Experiment", "experiments", "Wall", "wall",
           "reversal", "RevReading", "Reversal3D", "reversal_3d",
           "mutant_panel", "charge_scan", "bath_gamma_difference"]

READINGS = ("neutral", "pb", "pb + csc")


@dataclass
class Experiment:
    name: str
    lumen: dict[str, float]            # M
    cytosol: dict[str, float]
    impermeant: dict[str, float] = field(default_factory=dict)  # lumen only
    reference: str = "lumen"           # the csc fluid's bath

    def species(self):
        return ions(self.lumen, self.cytosol)

    def reference_species(self):
        bath = self.lumen if self.reference == "lumen" else self.cytosol
        return ions(bath, bath)


def experiments(ryr: bool) -> list[Experiment]:
    """The family's protocol: Xu's one experiment, or Vais's two."""
    if ryr:
        kcl = _P.value("permeation.ryr1_bath_concentration")
        ca = _P.value("selectivity.ryr1_cacl2_lumen")
        return [Experiment("Ca2+", {"K+": kcl, "Cl-": kcl + 2.0 * ca,
                                    "Ca2+": ca}, {"K+": kcl, "Cl-": kcl})]
    cyt, kcl_lum, ca_lum = _conditions()
    return [Experiment("Cl-", kcl_lum, cyt,
                       {"NMDG+": _P.value("selectivity.nmdg_cl")}, "cytosol"),
            Experiment("Ca2+", ca_lum, cyt)]


@dataclass
class Wall:
    """One deposit's wall on the voxels: its charge map and the charge–space
    fluid of its groups (both Round 7.19's ``local`` placement)."""

    fixed: np.ndarray
    fluid: object
    net: float
    unreached: list[str]


def wall(pore: Pore, neutralise: frozenset[int] = frozenset(),
         pair_bridges: bool = False) -> Wall:
    from .charged3d import _profile
    st, summary = pore.st, pore.summary
    charge = pore_charge(st, summary.frame, _profile(st, summary),
                         pair_bridges=pair_bridges, neutralise=neutralise)
    pos = group_positions(st, summary.frame, charge.groups)
    q = np.array([g.charge for g in charge.groups])
    fixed, unreached = local_density(
        pore.elec, pos, q, [f"{g.label()}/{g.chain}" for g in charge.groups])
    return Wall(fixed, wall_fluid(pore.elec, charge.groups, pos, "local"),
                charge.net_charge, unreached)


def _impermeant(pore: Pore, exp: Experiment) -> list[Impermeant]:
    lo = pore.summary.span[0]
    below = pore.elec.mask & (pore.elec.zs[None, None, :] < lo)
    return [Impermeant(name, 1, c, below) for name, c in exp.impermeant.items()]


def reversal(dom: Domain, species, fixed, excess=None, impermeant=(),
             bracket: float | None = None) -> tuple[float, Steady, int]:
    """The voltage (V) at which the steady net current is zero; the steady
    state there; and how many steady states it took. Each solve starts
    from the nearest one already found."""
    bracket = _P.value("reversal3d.bracket") if bracket is None else bracket
    tol = _P.value("reversal3d.voltage_tolerance")
    done: dict[float, Steady] = {}

    def current(v: float) -> float:
        if v not in done:
            near = min(done, key=lambda w: abs(w - v)) if done else None
            done[v] = steady_state(dom, species, v, fixed, excess, impermeant,
                                   initial=done.get(near))
        return done[v].current

    probe = 0.01
    f0, f1 = current(0.0), current(probe)
    guess = float(np.clip(-f0 * probe / (f1 - f0) if f1 != f0 else 0.0,
                          -bracket, bracket))
    half = 0.005
    while True:
        a, b = max(guess - half, -bracket), min(guess + half, bracket)
        if current(a) * current(b) <= 0:
            break
        if a <= -bracket and b >= bracket:
            raise ValueError(f"no reversal within +/- {bracket * 1e3:.0f} mV")
        half *= 2
    v = brentq(current, a, b, xtol=tol)
    current(v)
    return v, done[v], len(done)


@dataclass
class RevReading:
    label: str
    v: dict[str, float]                # experiment -> V_rev, V
    pca_pk: float
    pcl_pk: float                      # nan under Xu's ruler (P_Cl taken 0)
    currents: dict[str, dict[str, float]]   # experiment -> species -> A at V_rev
    peak_ca: float                     # M, highest Ca2+ in the lumen at V_rev
    converged: bool
    solves: int
    # Xu's V_rev read with the model's own P_Cl:P_K (linear response)
    # instead of Xu's 0: what the ruler's assumption costs (RyR1 only)
    pca_pk_model_cl: float = float("nan")


@dataclass
class Reversal3D:
    name: str
    ryr: bool
    spacing: float
    readings: dict[str, RevReading]
    measured: float
    linear: dict[str, float] = field(default_factory=dict)   # Round 7.19's
    one_d: dict[str, float] = field(default_factory=dict)
    net: float = 0.0
    unreached: list[str] = field(default_factory=list)

    @property
    def converged(self) -> bool:
        return all(r.converged for r in self.readings.values())


def _excess(pore: Pore, w: Wall, exp: Experiment, **fluid_kw):
    ref = local_csc(pore.elec, w.fixed, exp.reference_species(), w.fluid,
                    **fluid_kw)
    return ref.excess, ref.converged


def _reading(pore, dom, w, label, exps, **fluid_kw) -> RevReading:
    fixed = np.zeros_like(w.fixed) if label == "neutral" else w.fixed
    v, cur, ok, solves = {}, {}, True, 0
    peak = 0.0
    for exp in exps:
        excess = None
        if label.endswith("csc"):
            excess, c_ok = _excess(pore, w, exp, **fluid_kw)
            ok &= c_ok
        sp = exp.species()
        v[exp.name], st, n = reversal(dom, sp, fixed, excess,
                                      _impermeant(pore, exp))
        ok &= st.converged
        solves += n
        cur[exp.name] = st.currents
        if exp.name == "Ca2+":
            psi = v[exp.name] / thermal_voltage() * dom.phi0 + st.u
            e = 2.0 * psi + (excess["Ca2+"] if excess else 0.0)
            c = st.n["Ca2+"] * np.exp(-np.clip(e, -40, 40))
            peak = float(c[pore.elec.mask].max() / 1000.0)
    by = {e.name: e for e in exps}
    if pore.ryr:
        e = by["Ca2+"]
        pcl, pca = float("nan"), ghk_ratio(v["Ca2+"], "Ca2+", e.lumen,
                                           e.cytosol, {"Cl-": 0.0})
    else:
        ec, ea = by["Cl-"], by["Ca2+"]
        pcl = ghk_ratio(v["Cl-"], "Cl-", ec.lumen, ec.cytosol, {})
        pca = ghk_ratio(v["Ca2+"], "Ca2+", ea.lumen, ea.cytosol, {"Cl-": pcl})
    return RevReading(label, v, float(pca), float(pcl), cur, peak, bool(ok),
                      solves)


def reversal_3d(st_or_pore, readings=READINGS, spacing: float | None = None,
                neutralise: frozenset[int] = frozenset(),
                linear: bool = True, **fluid_kw) -> Reversal3D:
    """Every reading of one deposit under its family's protocol, with
    Round 7.19's linear-response ratios on the same grid beside them
    (``linear``) and the 1-D readings."""
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = (st_or_pore if isinstance(st_or_pore, Pore)
            else prepare(st_or_pore, spacing=h))
    dom = domain(pore.elec, pore.vols)
    w = wall(pore, neutralise)
    exps = experiments(pore.ryr)
    out = {label: _reading(pore, dom, w, label, exps, **fluid_kw)
           for label in readings}
    lin, one_d = {}, {}
    if linear:
        from .selectivity3d import selectivity_3d
        s = selectivity_3d(pore, readings, neutralise=neutralise, **fluid_kw)
        lin = {k: r.ratio for k, r in s.readings.items()}
        one_d = s.one_d
        if pore.ryr:
            e = exps[0]
            for k, r in out.items():
                r.pca_pk_model_cl = float(ghk_ratio(
                    r.v["Ca2+"], "Ca2+", e.lumen, e.cytosol,
                    {"Cl-": s.readings[k].pcl_pk}))
    return Reversal3D(pore.st.name, pore.ryr, pore.spacing, out,
                      measured_ratio(pore.ryr), lin, one_d, w.net,
                      w.unreached)


def charge_scan(st: Structure, scales, spacing: float | None = None,
                progress=None) -> list[tuple[float, RevReading]]:
    """The ``pb`` reading with every wall charge multiplied by each scale:
    the path (P_Cl:P_K, P_Ca:P_K) a mean-field wall of this shape can take,
    from uncharged (0) through the deposit's (1)."""
    from dataclasses import replace
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = prepare(st, spacing=h)
    dom = domain(pore.elec, pore.vols)
    w = wall(pore)
    exps = experiments(pore.ryr)
    out = []
    for i, f in enumerate(scales):
        if progress:
            progress(i, len(scales), f)
        out.append((float(f), _reading(pore, dom, replace(w, fixed=w.fixed * f),
                                       "pb", exps)))
    return out


def mutant_panel(st: Structure | None = None, readings=("pb", "pb + csc"),
                 spacing: float | None = None, progress=None, **kw
                 ) -> list[tuple[str, float, Reversal3D]]:
    """Xu 2006's wild type and five mutants at Xu's reversal on the open
    RyR1 deposit: (name, measured P_Ca:P_K, reading); volumes shared."""
    from ..io import loader
    from .ryr_mutants import mutants, open_deposit
    st = st or loader.load(open_deposit())
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = prepare(st, spacing=h)
    key = "selectivity.published_ryr1_pca_pk"
    todo = [("wild type", frozenset(), _P.value(key))]
    for name, (res, _) in sorted(mutants().items(), key=lambda kv: kv[1][0]):
        todo.append((name, frozenset({res}), _P.value(f"{key}_{name.lower()}")))
    out = []
    for i, (name, off, pca) in enumerate(todo):
        if progress:
            progress(i, len(todo), name)
        out.append((name, pca, reversal_3d(pore, readings, neutralise=off,
                                           linear=False, **kw)))
    return out


def bath_gamma_difference(ryr: bool, **fluid_kw) -> dict[str, float]:
    """|ln γ_lumen − ln γ_cytosol| per species (kT) in the Ca²⁺ experiment's
    two baths, under the csc fluid with no wall: the size of what holding
    one reference bath leaves out."""
    from . import csc
    from .csc3d import fluid_for
    exp = next(e for e in experiments(ryr) if e.name == "Ca2+")
    out = {}
    lum, cyt = ions(exp.lumen, exp.lumen), ions(exp.cytosol, exp.cytosol)
    fl_l, a_l, _, b_l = fluid_for(lum, csc.Structural.empty(1), **fluid_kw)
    fl_c, a_c, _, b_c = fluid_for(cyt, csc.Structural.empty(1), **fluid_kw)
    gl = {s.name: np.log(a / b) for s, a, b in zip(lum, a_l, b_l)}
    gc = {s.name: np.log(a / b) for s, a, b in zip(cyt, a_c, b_c)}
    for name in gc:
        out[name] = float(abs(gl[name] - gc[name]))
    return out
