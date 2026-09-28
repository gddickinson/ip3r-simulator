"""A saturable Ca²⁺ site, compensated or not, and the K⁺ block its
occupancy sets, at Vais's reversal (Round 7.27).

Round 7.25 found that no mean-field wall reaches Vais 2010's pair
(P_Cl:P_K 0.27, P_Ca:P_K 15.2): for point ions in series,
(P_Cl:P_K)(P_Ca:P_K)² over the uncharged pore's is at most 1. A
Ca²⁺-only well leaves Cl⁻ alone, but it peaked at 4.48 on 8TKF because the
Ca²⁺ it gathers is uncompensated charge. Its emergent item names two
interactions the well lacks. This module adds both.

**The site.** ``casite.sites`` sites spread evenly over one band of the
lumen (density s over the band's voxels), each binding one Ca²⁺. Its
affinity is given as the dilute-limit depth d (kT): an empty site lowers
Ca²⁺'s energy by d, so K_d = s / (e^d − 1). With fast binding and bound
Ca²⁺ moving with the free (the well picture, as in 7.25), the total is
c_f (1 + s / (c_f + K_d)). This is a Ca²⁺ energy
μ = −ln(1 + s / (c_f + K_d)) that saturates as the site fills. The free
Ca²⁺ is c_f = n e^{−2ψ} in Slotboom form, whatever μ is, so μ is a
:mod:`.pnp3d` coupling. The occupancy is θ = c_f / (c_f + K_d).

* **Compensated**: each bound Ca²⁺ comes with −2e fixed, so filling the
  site costs no field. μ is then ``hidden``: transport sees it, Poisson
  counts only c_f.
* **Uncompensated**: Poisson counts the bound Ca²⁺ too (7.25's well, now
  saturable).

**The block.** A K⁺ passing through the band finds a fraction θ of it
held by Ca²⁺ and passes an occupied site with probability 1 − f
(f = ``casite.block``). In the mean field this is the energy
−ln(1 − f θ) on K⁺ over the band: the anomalous-mole-fraction mechanism
in its simplest continuum form. Poisson counts the K⁺ it excludes.
Neither term touches Cl⁻, and Vais's Cl⁻ experiment holds no Ca²⁺, so it
reads the uncharged pore's P_Cl:P_K exactly (0.29 on 8TKF, measured 0.27).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .pnp3d import Domain
from .reversal3d import _impermeant, concentration, reversal
from .selectivity import ghk_ratio, thermal_voltage
from .selectivity_bound import band_mask, bound_product
from .wall_search import _Solver, regions, score

__all__ = ["N_AVOGADRO", "Site", "SiteCoupling", "SiteResult",
           "SiteSearch", "site_density", "read", "load", "depths", "search",
           "required_depth", "site_regions"]

N_AVOGADRO = 6.02214076e23
#: The largest block energy (kT) on K⁺ where the site is full and f = 1.
_BLOCK_CLIP = 30.0


@dataclass(frozen=True)
class Site:
    depth: float                       # kT, the empty site's pull on Ca2+
    compensated: bool = False
    block: float = 0.0                 # f: K+ stopped by an occupied site
    region: str = "span"

    @property
    def label(self) -> str:
        kind = "compensated" if self.compensated else "uncompensated"
        blk = f", K+ block {self.block:g}" if self.block else ""
        return f"{kind} site {self.depth:g} kT at {self.region}{blk}"


def site_density(mask: np.ndarray, spacing: float, sites: float) -> float:
    """``sites`` spread evenly over the mask's voxels (Å grid): mol/m³."""
    volume = float(mask.sum()) * (spacing * 1e-10) ** 3
    return sites / (N_AVOGADRO * volume)


class SiteCoupling:
    """The site's energies as a :mod:`.pnp3d` coupling: Ca²⁺'s saturable
    μ (hidden when compensated) and K⁺'s block, from each iterate's free
    Ca²⁺."""

    def __init__(self, dom: Domain, site: Site, mask: np.ndarray,
                 density: float):
        self.dom, self.site, self.mask = dom, site, mask.astype(bool)
        self.s = density
        self.kd = density / np.expm1(site.depth)

    def free(self, psi: np.ndarray, n_ca: np.ndarray) -> np.ndarray:
        return n_ca * np.exp(-np.clip(2.0 * psi, -40.0, 40.0))

    def occupancy(self, c_free: np.ndarray) -> np.ndarray:
        return np.where(self.mask, c_free / (c_free + self.kd), 0.0)

    def energies(self, c_free: np.ndarray | None) -> tuple[dict, dict]:
        if c_free is None:                        # the dilute limit
            mu = -self.site.depth * self.mask
            theta = np.zeros(self.mask.shape)
        else:
            mu = -np.log1p(self.s / (c_free + self.kd)) * self.mask
            theta = self.occupancy(c_free)
        ca = {"Ca2+": mu}
        vis, hid = ({}, ca) if self.site.compensated else (ca, {})
        if self.site.block:
            k = -np.log1p(-self.site.block * np.minimum(theta, 1.0))
            vis["K+"] = np.minimum(k, _BLOCK_CLIP) * self.mask
        return vis, hid

    def start(self) -> tuple[dict, dict]:
        return self.energies(None)

    def __call__(self, psi, n) -> tuple[dict, dict]:
        if "Ca2+" not in n:                       # the Cl- experiment
            return {}, {}
        return self.energies(self.free(psi, n["Ca2+"]))


@dataclass
class SiteResult:
    site: Site
    pca_pk: float
    pcl_pk: float
    v: dict[str, float]
    peak_ca: float                     # M, total Ca2+ (free + bound)
    occupancy: float                   # the band's mean θ at the Ca2+ reversal
    bound: float                       # Ca2+ ions held by the site there
    converged: bool
    score: float = float("nan")
    product: float = float("nan")      # B relative to the uncharged pore


@dataclass
class SiteSearch:
    name: str
    measured: tuple[float, float]
    neutral: object                    # wall_search.Result, uncharged
    density: float                     # mol/m3 over the span
    results: list[SiteResult] = field(default_factory=list)
    required: dict[str, float] = field(default_factory=dict)

    @property
    def best(self) -> SiteResult:
        return min(self.results, key=lambda r: r.score)


def site_regions(summary) -> dict[str, tuple[float, float]]:
    """Round 7.25's bands, plus the luminal vestibule: from the span's
    luminal end to the filter band (Round 7.29)."""
    bands = regions(summary)
    bands["vestibule"] = (summary.span[0], bands["filter"][0])
    return bands


def depths() -> np.ndarray:
    step = _P.value("casite.depth_step")
    return np.arange(step, _P.value("casite.depth_max") + step / 2, step)


def read(solver: _Solver, site: Site, bands: dict) -> SiteResult:
    """One site at Vais's reversal on an uncharged wall (the Cl⁻ experiment
    is the uncharged pore's)."""
    pore, dom = solver.pore, solver.dom
    pcl, v_cl, ok_cl = solver.chloride(0.0)
    mask = band_mask(pore, *bands[site.region])
    density = site_density(mask, pore.spacing, _P.value("casite.sites"))
    coupling = SiteCoupling(dom, site, mask, density)
    fixed = np.zeros_like(solver.fixed)
    v, st, _ = reversal(dom, solver.ca.species(), fixed, None,
                        _impermeant(pore, solver.ca), coupling=coupling)
    pca = ghk_ratio(v, "Ca2+", solver.ca.lumen, solver.ca.cytosol,
                    {"Cl-": pcl})
    theta, held = load(coupling, st, pore.spacing)
    conc = concentration(dom, st, "Ca2+", 2)
    return SiteResult(site, float(pca), pcl, {"Cl-": v_cl, "Ca2+": v},
                      float(conc[pore.elec.mask].max() / 1000.0),
                      theta, held, bool(ok_cl and st.converged))


def load(coupling: SiteCoupling, st, spacing: float) -> tuple[float, float]:
    """The band's mean occupancy θ in steady state ``st``, and the Ca²⁺
    ions the site holds (0 with no Ca²⁺ in either bath)."""
    if "Ca2+" not in st.n:
        return 0.0, 0.0
    psi = st.v / thermal_voltage() * coupling.dom.phi0 + st.u
    theta = coupling.occupancy(coupling.free(psi, st.n["Ca2+"]))
    on = coupling.mask
    held = float(theta[on].sum() * coupling.s * (spacing * 1e-10) ** 3
                 * N_AVOGADRO)
    return float(theta[on].mean()), held


def _finish(r: SiteResult, neutral, measured) -> SiteResult:
    r.score = score(r.pca_pk, r.pcl_pk, measured)
    r.product = bound_product(r.pcl_pk, r.pca_pk, neutral.pcl_pk,
                              neutral.pca_pk)
    return r


def search(st: Structure, spacing: float | None = None, region: str = "span",
           kinds=((False, 0.0), (True, 0.0), (False, None), (True, None)),
           solve_required: bool = True, progress=None) -> SiteSearch:
    """Every depth of each kind (compensated, block; None = ``casite.block``)
    on one IP3R deposit at Vais's reversal; then, for each kind that
    crosses 15.2, the depth at which it does."""
    from .selectivity3d import prepare
    from .wall_search import Candidate
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = prepare(st, spacing=h)
    if pore.ryr:
        raise ValueError("the site is read against Vais 2010's pair: give "
                         "an IP3R deposit")
    measured = (_P.value("selectivity.published_pcl_pk"),
                _P.value("selectivity.published_pca_pk"))
    solver, bands = _Solver(pore), site_regions(pore.summary)
    neutral = solver.read(Candidate("charge", 0.0), bands)
    neutral.score = score(neutral.pca_pk, neutral.pcl_pk, measured)
    mask = band_mask(pore, *bands[region])
    out = SiteSearch(st.name, measured, neutral,
                     site_density(mask, h, _P.value("casite.sites")))
    f_block = _P.value("casite.block")
    todo = [Site(float(d), comp, f_block if f is None else f, region)
            for comp, f in kinds for d in depths()]
    for i, s in enumerate(todo):
        if progress:
            progress(i, len(todo), s.label)
        out.results.append(_finish(read(solver, s, bands), neutral, measured))
    if solve_required:
        for comp, f in kinds:
            f = f_block if f is None else f
            rows = [r for r in out.results if r.site.compensated == comp
                    and r.site.block == f and r.converged]
            d = required_depth(solver, bands, rows, measured[1])
            if d is not None:
                out.required[Site(d, comp, f, region).label] = d
    return out


def required_depth(solver, bands, rows: list[SiteResult], target: float
                   ) -> float | None:
    """The depth at which P_Ca:P_K crosses ``target``, by root-finding
    between the first grid pair that brackets it; None when none does."""
    rows = sorted(rows, key=lambda r: r.site.depth)
    for a, b in zip(rows, rows[1:]):
        if (a.pca_pk - target) * (b.pca_pk - target) <= 0:
            def f(d, base=a.site):
                s = Site(d, base.compensated, base.block, base.region)
                return read(solver, s, bands).pca_pk - target

            return float(brentq(f, a.site.depth, b.site.depth,
                                xtol=_P.value("casite.depth_tolerance")))
    return None
