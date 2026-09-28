"""The mole-fraction prediction (Round 7.29): P_Ca:P_K, the K⁺ current
and i_Ca against luminal Ca²⁺ through Round 7.27's crossing site.

Round 7.27 found one wall that gives Vais 2010's pair: a compensated
Ca²⁺ site over the span whose occupancy blocks K⁺. It crosses 15.2 at
d = 4.41 kT on 8TKF (K_d 3.3 mM) and is half full at Vais's 10 mM. A site
that works through occupancy predicts that what is measured depends on
how full the site is, so it depends on the luminal Ca²⁺. The uncharged
pore, and a site that does not block, are the controls that show which
part of the curve comes from the block.

For each luminal CaCl₂ c (``molefrac.ca_min``–``molefrac.ca_max``, added
to Vais's 140 mM KCl, cytosol 140 mM KCl, as in their P_Ca:P_K
experiment):

* **V_rev → P_Ca:P_K** by Vais's Eq. 1 with the pore's own P_Cl:P_K (the
  Cl⁻ experiment holds no Ca²⁺, so the site leaves it alone);
* **i_Ca at 0 mV** (Vais's i_Ca protocol: 0.30 pA/mM over 0.16–1.1 mM);
* **i_K at ``molefrac.voltage``**, against the same pore with no Ca²⁺ at
  all. Its fall is the block, and the c at which it halves is the
  prediction an experiment can test.

``v`` is cytosol − lumen throughout (:mod:`.pnp3d`). At a negative ``v``
both cations flow lumen → cytosol, so their currents are positive.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .ca_site import Site, SiteCoupling, load, site_density, site_regions
from .pnp3d import steady_state
from .reversal3d import Experiment, reversal
from .selectivity import _conditions, ghk_ratio
from .selectivity_bound import band_mask

__all__ = ["READINGS", "concentrations", "experiment", "Point", "Curve",
           "MoleFraction", "reading_site", "curve", "sweep", "half_point"]

#: The site as 7.27 found it, the same site not blocking K⁺, and no site.
READINGS = ("site", "site, no block", "uncharged")


def concentrations() -> np.ndarray:
    lo, hi = _P.value("molefrac.ca_min"), _P.value("molefrac.ca_max")
    n = int(round(np.log10(hi / lo) * _P.value("molefrac.points_per_decade")))
    return np.geomspace(lo, hi, n + 1)


def experiment(ca: float) -> Experiment:
    """Vais's P_Ca:P_K experiment with ``ca`` M CaCl₂ in the lumen (none
    when 0)."""
    cyt, _, lum = _conditions()
    lumen = {"K+": lum["K+"], "Cl-": lum["K+"] + 2.0 * ca}
    if ca > 0:
        lumen["Ca2+"] = ca
    return Experiment("Ca2+" if ca > 0 else "KCl", lumen, dict(cyt))


@dataclass
class Point:
    ca: float                          # M, luminal CaCl2
    v_rev: float                       # V
    pca_pk: float
    i_zero: dict[str, float]           # A per species at 0 mV
    i_hold: dict[str, float]           # A per species at molefrac.voltage
    occupancy: float                   # the band's mean theta at the hold
    held: float                        # Ca2+ ions held there
    converged: bool

    @property
    def i_ca(self) -> float:
        return self.i_zero.get("Ca2+", 0.0)


@dataclass
class Curve:
    reading: str
    site: Site | None
    kd: float                          # mM (nan without a site)
    free: dict[str, float]             # A per species at the hold, no Ca2+
    points: list[Point] = field(default_factory=list)

    @property
    def ca(self) -> np.ndarray:
        return np.array([p.ca for p in self.points])

    @property
    def k_fraction(self) -> np.ndarray:
        """i_K at the hold relative to the Ca²⁺-free pore's."""
        return np.array([p.i_hold["K+"] / self.i_k0 for p in self.points])

    @property
    def i_k0(self) -> float:
        return self.free["K+"]

    @property
    def total_fraction(self) -> np.ndarray:
        """The net current at the hold relative to the Ca²⁺-free pore's."""
        total = sum(self.free.values())
        return np.array([sum(p.i_hold.values()) / total for p in self.points])

    def ica_slope(self, lo: float, hi: float) -> float:
        """i_Ca at 0 mV per M over [lo, hi] (least squares through 0),
        A/M; Vais's protocol."""
        sel = [p for p in self.points if lo * (1 - 1e-9) <= p.ca <= hi * (1 + 1e-9)]
        if not sel:
            return float("nan")
        c = np.array([p.ca for p in sel])
        i = np.array([p.i_ca for p in sel])
        return float(c @ i / (c @ c))


@dataclass
class MoleFraction:
    name: str
    voltage: float                     # V, the hold
    pcl_pk: float
    region: str
    curves: dict[str, Curve] = field(default_factory=dict)


def reading_site(reading: str, depth: float, region: str) -> Site | None:
    if reading == "uncharged":
        return None
    block = _P.value("casite.block") if reading == "site" else 0.0
    return Site(depth, True, block, region)


def _hold(dom, exp, v, coupling, initial=None):
    fixed = np.zeros(dom.elec.mask.shape)
    return steady_state(dom, exp.species(), v, fixed, initial=initial,
                        coupling=coupling)


def curve(dom, mask: np.ndarray, density: float, spacing: float,
          reading: str, depth: float, region: str, cas: np.ndarray,
          voltage: float, pcl: float, progress=None) -> Curve:
    """One reading's points: the site (if any) of ``density`` mol/m³ over
    ``mask`` on the :mod:`.pnp3d` domain ``dom`` (grid ``spacing`` Å)."""
    site = reading_site(reading, depth, region)

    def coupling_for():
        return None if site is None else SiteCoupling(dom, site, mask, density)

    free = _hold(dom, experiment(0.0), voltage, None)
    out = Curve(reading, site, density / np.expm1(depth) if site else np.nan,
                free.currents)
    fixed = np.zeros(dom.elec.mask.shape)
    for i, ca in enumerate(cas):
        if progress:
            progress(i, len(cas), f"{reading}, {ca * 1e3:.3g} mM")
        exp = experiment(float(ca))
        cp = coupling_for()
        v, st, _ = reversal(dom, exp.species(), fixed, None, (), coupling=cp)
        pca = ghk_ratio(v, "Ca2+", exp.lumen, exp.cytosol, {"Cl-": pcl})
        zero = _hold(dom, exp, 0.0, cp, initial=st)
        hold = _hold(dom, exp, voltage, cp, initial=zero)
        theta, held = (load(cp, hold, spacing) if cp is not None
                       else (0.0, 0.0))
        out.points.append(Point(float(ca), v, float(pca), zero.currents,
                                hold.currents, theta, held,
                                bool(st.converged and zero.converged
                                     and hold.converged)))
    return out


def sweep(st: Structure, readings=READINGS, depth: float | None = None,
          region: str = "span", spacing: float | None = None,
          progress=None) -> MoleFraction:
    """Every reading over :func:`concentrations` on one IP3R deposit.
    ``depth`` defaults to ``casite.gui_depth`` (8TKF's crossing)."""
    from .selectivity3d import prepare
    from .wall_search import _Solver
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = prepare(st, spacing=h)
    if pore.ryr:
        raise ValueError("the sweep runs Vais 2010's protocol: give an IP3R "
                         "deposit")
    depth = _P.value("casite.gui_depth") if depth is None else depth
    solver = _Solver(pore)
    pcl, _, _ = solver.chloride(0.0)
    voltage = _P.value("molefrac.voltage")
    out = MoleFraction(st.name, voltage, pcl, region)
    mask = band_mask(pore, *site_regions(pore.summary)[region])
    density = site_density(mask, h, _P.value("casite.sites"))
    for r in readings:
        out.curves[r] = curve(solver.dom, mask, density, h, r, depth, region,
                              concentrations(), voltage, pcl, progress)
    return out


def half_point(ca: np.ndarray, frac: np.ndarray, level: float = 0.5
               ) -> float:
    """The c at which ``frac`` first falls through ``level``, interpolated
    in log c; NaN when it does not."""
    for (c0, f0), (c1, f1) in zip(zip(ca, frac), zip(ca[1:], frac[1:])):
        if f0 >= level > f1:
            t = (f0 - level) / (f0 - f1)
            return float(np.exp(np.log(c0) + t * (np.log(c1) - np.log(c0))))
    return float("nan")
