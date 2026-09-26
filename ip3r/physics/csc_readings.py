"""Round 7.17's readings: charge-space competition against the measured
selectivities.

* :func:`ryr1_panel`: Xu 2006's wild type and five mutants, P_Ca:P_K by
  their protocol and Eq. 1 and the K+ conductance in 250 mM KCl, under
  the Donnan closure and the csc closure with and without the wall
  groups' volume (``READINGS``).
* :func:`vais_reading`: an IP3R deposit under Vais 2010's protocols, the
  same readings.
* :func:`energetics`: the binding selectivity of the filter's most charged
  slice, split as Gillespie 2008 splits it (Fig. 7): the Ca2+ advantage
  over K+ from the mean potential, from screening (MSA) and from excluded
  volume (hard spheres), each relative to the bath.
* :func:`shares`: where each ion's resistance lies, in linear response
  about an equilibrium in a symmetric mixed bath (the integral
  1 / P_i = int exp(z_i w_i) / (D_i A_i) dz, w_i the species' offset),
  at the filter and at the gate, and the P_Ca:P_K that integral gives.
* :func:`gate_limit`: the same ratio over the gate window alone with no
  charge there, by hand: (D_Ca / D_K) x the areas' ratio.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from . import csc
from ._pnp_kernels import _donnan_potential
from .permeation import _accessible_area, potassium_species, solve_pnp
from .selectivity import ions, ryr1_calcium_ratio, selectivity, thermal_voltage

__all__ = ["READINGS", "CscWall", "csc_wall", "moved", "PanelRow",
           "ryr1_panel", "vais_reading", "Energetics", "energetics",
           "filter_slice",
           "Shares", "shares", "gate_limit", "SCAN_KEYS", "ryr1_scan"]

#: label -> (closure, csc.structural_volume or None).
READINGS = {"donnan": ("donnan", None), "csc": ("csc", 1.0),
            "csc, wall volume 0": ("csc", 0.0)}

#: The constants ``csc --scan`` moves.
SCAN_KEYS = ("csc.structural_volume", "csc.permittivity",
             "csc.water_diameter", "csc.oxygen_diameter")


@contextmanager
def moved(key: str, value):
    """One registered constant at ``value`` (None = leave it), restored after."""
    if value is None:
        yield
        return
    before = _P.overrides().get(key)
    _P.set_value(key, value)
    try:
        yield
    finally:
        if before is None:
            _P.reset(key)
        else:
            _P.set_value(key, before)


@dataclass
class CscWall:
    """A deposit's pore with its wall charge and the groups as fluid."""

    name: str
    z: np.ndarray
    radius: np.ndarray
    fixed: np.ndarray
    wall: csc.Structural
    ryr: bool
    constrictions: dict
    net: float = 0.0


def csc_wall(st: Structure, neutralise: frozenset[int] = frozenset(),
             base=None) -> CscWall:
    """The formal lining wall of ``st`` (less ``neutralise``)."""
    from .pore_charge import pore_charge
    from .protonation import lining_wall
    base = base or lining_wall(st)
    ch = pore_charge(st, base.summary.frame, base.profile,
                     neutralise=neutralise)
    z = np.asarray(base.profile.z, float)
    return CscWall(st.name, z, base.radius, ch.density,
                   csc.structural(ch.groups, z, base.radius), base.ryr,
                   base.summary.constrictions, ch.net_charge)


def _read(w: CscWall, reading: str, fn):
    closure, volume = READINGS[reading]
    with moved("csc.structural_volume", volume):
        return fn(closure, None if closure == "donnan" else w.wall)


# ------------------------------------------------------------ the panels
@dataclass
class PanelRow:
    name: str
    measured_pca: float
    measured_g: float                  # pS (NaN when not measured)
    pca: dict = field(default_factory=dict)    # reading -> P_Ca:P_K
    g: dict = field(default_factory=dict)      # reading -> pS
    converged: bool = True


def _ryr1_row(w: CscWall, name, pca, g_meas, readings) -> PanelRow:
    bath = _P.value("permeation.ryr1_bath_concentration")
    row = PanelRow(name, pca, g_meas)
    for r in readings:
        _, p, ok = _read(w, r, lambda cl, s: ryr1_calcium_ratio(
            w.z, w.radius, w.fixed, closure=cl, structural=s))
        g = _read(w, r, lambda cl, s: solve_pnp(
            w.z, w.radius, species=potassium_species(bath=bath),
            fixed_charge=w.fixed, closure=cl, structural=s))
        row.pca[r], row.g[r] = float(p), g.conductance_pS
        row.converged &= bool(ok and g.converged)
    return row


def ryr1_panel(st: Structure, readings=tuple(READINGS), progress=None
               ) -> list[PanelRow]:
    """Wild type and Xu's mutants, measured and under each reading."""
    from .protonation import lining_wall
    from .ryr_mutants import mutants
    base = lining_wall(st)
    key = "selectivity.published_ryr1_pca_pk"
    todo = [("wild type", frozenset(), _P.value(key),
             _P.value("permeation.published_ryr1"))]
    for name, (res, pS) in sorted(mutants().items(), key=lambda kv: kv[1][0]):
        todo.append((name, frozenset({res}),
                     _P.value(f"{key}_{name.lower()}"), pS))
    rows = []
    for i, (name, off, pca, pS) in enumerate(todo):
        rows.append(_ryr1_row(csc_wall(st, off, base), name, pca, pS,
                              readings))
        if progress:
            progress(i + 1, len(todo))
    return rows


def vais_reading(st: Structure, readings=tuple(READINGS)) -> dict:
    """reading -> (P_Ca:P_K, P_Cl:P_K, K+ g pS, converged), Vais's protocols."""
    w = csc_wall(st)
    out = {}
    for r in readings:
        sel = _read(w, r, lambda cl, s, label=r: selectivity(
            w.z, w.radius, w.fixed, label, closure=cl, structural=s))
        g = _read(w, r, lambda cl, s: solve_pnp(
            w.z, w.radius, fixed_charge=w.fixed, closure=cl, structural=s))
        out[r] = (sel.pca_pk, sel.pcl_pk, g.conductance_pS,
                  bool(sel.converged and g.converged))
    return out


# ------------------------------------------------- binding at the filter
@dataclass
class Energetics:
    """Gillespie 2008's split of the Ca2+ advantage over K+ at one slice,
    kT (positive favours Ca2+), and the slice's contents (M)."""

    z: float
    fixed: float
    k: float
    ca: float
    water: float
    packing: float
    mean: float
    screening: float
    excluded: float
    gamma: dict                        # bath ln(gamma): {"K+": .., "Ca2+": ..}


def _mixed(kcl: float, ca: float):
    bath = {"K+": kcl, "Cl-": kcl + 2.0 * ca, "Ca2+": ca}
    return ions(bath, bath)


def _split(c, water, fluid: csc.Fluid, wall: csc.Structural):
    """(hard-sphere mu, MSA mu) of the mobile species, kT."""
    from dataclasses import replace
    hs, _ = csc.excess(c, water, replace(fluid, msa=False), wall)
    es, _ = csc.excess(c, water, replace(fluid, hs=False), wall)
    return hs, es


def filter_slice(w: CscWall) -> int:
    """The most charged slice within ``csc.filter_reach`` of the filter."""
    zf = w.constrictions["filter"].z
    near = np.abs(w.z - zf) <= _P.value("csc.filter_reach")
    return int(np.flatnonzero(near)[np.argmax(np.abs(w.fixed[near]))])


def energetics(st_or_wall, kcl: float, ca: float) -> Energetics:
    """The filter slice in equilibrium with KCl + CaCl2 on both sides."""
    w = (st_or_wall if isinstance(st_or_wall, CscWall)
         else csc_wall(st_or_wall))
    k = filter_slice(w)
    sp = _mixed(kcl, ca)
    names = [s.name for s in sp]
    bath = np.array([s.concentration * 1000.0 for s in sp])
    wall = w.wall.take([k])
    fl = csc.Fluid.of(sp, wall)
    a, aw = csc.bath_activity(bath, fl)
    p = csc.partition(a[:, None], w.fixed[[k]], fl, thermal_voltage(), aw,
                      initial=(bath[:, None], np.array([fl.water])))
    empty = csc.Structural.empty(1)
    # The bath's water is csc.water_concentration by definition (its
    # activity was taken there).
    hs_b, es_b = _split(bath[:, None], np.array([fl.water]), fl, empty)
    hs_f, es_f = _split(p.c, p.water, fl, wall)
    i_k, i_ca = names.index("K+"), names.index("Ca2+")

    def adv(f, b):
        return float(-((f[i_ca, 0] - b[i_ca, 0]) - (f[i_k, 0] - b[i_k, 0])))
    return Energetics(float(w.z[k]), float(w.fixed[k] / 1000.0),
                      float(p.c[i_k, 0] / 1000.0),
                      float(p.c[i_ca, 0] / 1000.0),
                      float(p.water[0] / 1000.0),
                      float(csc.packing(p.c, p.water, fl, wall)[0]),
                      float(-p.psi[0] / thermal_voltage()),
                      adv(es_f, es_b), adv(hs_f, hs_b),
                      {n: float(hs_b[i, 0] + es_b[i, 0])
                       for i, n in enumerate(names) if n != "Cl-"})


# ------------------------------------------------ where the resistance is
@dataclass
class Shares:
    """Each ion's resistance at the filter and the gate (fractions), and
    the P_Ca:P_K the linear-response integral gives."""

    closure: str
    filter: dict
    gate: dict
    ratio: float


def _offsets(w: CscWall, closure: str, sp):
    bath = np.array([s.concentration * 1000.0 for s in sp])
    n, th = len(w.z), thermal_voltage()
    if closure == "donnan":
        psi = _donnan_potential([s.valence for s in sp],
                                np.repeat(bath[:, None], n, axis=1),
                                w.fixed, th)
        return np.repeat(psi[None, :], len(sp), axis=0), np.zeros(len(sp))
    fl = csc.Fluid.of(sp, w.wall)
    a, aw = csc.bath_activity(bath, fl)
    p = csc.partition(np.repeat(a[:, None], n, axis=1), w.fixed, fl, th, aw,
                      initial=(np.repeat(bath[:, None], n, axis=1),
                               np.full(n, fl.water)))
    # The bath's own excess, so exp(z w) is the partition relative to it.
    return p.offset(fl.valence, th), np.log(a / bath)


def shares(w: CscWall, closure: str = "csc", kcl: float | None = None,
           ca: float | None = None) -> Shares:
    """Linear-response resistance of K+ and Ca2+ along the pore."""
    kcl = (_P.value("permeation.ryr1_bath_concentration") if w.ryr
           else _P.value("permeation.bath_concentration")) if kcl is None else kcl
    ca = _P.value("selectivity.ryr1_cacl2_lumen") if ca is None else ca
    sp = _mixed(kcl, ca)
    off, ref = _offsets(w, closure, sp)
    th, half = thermal_voltage(), _P.value("lumen.constriction_half_width")
    dens, frac = {}, {"filter": {}, "gate": {}}
    for i, s in enumerate(sp):
        if s.name == "Cl-":
            continue
        area = np.maximum(_accessible_area(w.radius * 1e-10, s.radius), 1e-40)
        dens[s.name] = np.exp(s.valence * off[i] / th - ref[i]) / (
            s.diffusivity * area)
        total = np.trapezoid(dens[s.name], w.z)
        for name in frac:
            c = w.constrictions[name].z
            m = np.abs(w.z - c) <= half
            frac[name][s.name] = float(np.trapezoid(dens[s.name][m], w.z[m])
                                       / total)
    ratio = float(np.trapezoid(dens["K+"], w.z)
                  / np.trapezoid(dens["Ca2+"], w.z))
    return Shares(closure, frac["filter"], frac["gate"], ratio)


def gate_limit(w: CscWall) -> float:
    """P_Ca:P_K of the gate window alone, uncharged: the integral of
    1/(D A) for K+ over that for Ca2+, by hand from the profile."""
    half = _P.value("lumen.constriction_half_width")
    m = np.abs(w.z - w.constrictions["gate"].z) <= half
    r = w.radius[m]
    zz = w.z[m]
    rk = _P.value("permeation.radius_potassium")
    rc = _P.value("permeation.radius_calcium")
    ik = np.trapezoid(1.0 / (_P.value("permeation.diffusion_potassium")
                             * (r - rk) ** 2), zz)
    ic = np.trapezoid(1.0 / (_P.value("permeation.diffusion_calcium")
                             * (r - rc) ** 2), zz)
    return float(ik / ic)


def ryr1_scan(st: Structure, key: str, values) -> list[tuple[float, float,
                                                             float]]:
    """(value, wild-type P_Ca:P_K, D4899Q / wild type) under csc, with one
    registered constant moved (restored afterwards)."""
    if key not in SCAN_KEYS:
        raise ValueError(f"key must be one of {SCAN_KEYS}, not {key!r}")
    from .protonation import lining_wall
    base = lining_wall(st)
    out = []
    for v in values:
        with moved(key, v):
            p = []
            for off in (frozenset(), frozenset({4899})):
                w = csc_wall(st, off, base)
                p.append(ryr1_calcium_ratio(w.z, w.radius, w.fixed,
                                            closure="csc",
                                            structural=w.wall)[1])
        out.append((float(v), float(p[0]), float(p[1] / p[0])))
    return out
