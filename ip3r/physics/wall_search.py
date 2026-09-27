"""A search over wall models scored on P_Cl:P_K and P_Ca:P_K at once
(Round 7.25).

Each candidate is a wall on one deposit's voxels, read at the family's
bi-ionic reversal (Round 7.23's steady Poisson–Nernst–Planck, point ions)
and scored against both measured ratios:
S = |ln(P_Ca:P_K / 15.2)| + |ln(P_Cl:P_K / 0.27)| (Vais 2010). The families:

* ``charge``: the deposit's own lining charge scaled (a geometric grid
  from ``wallsearch.scale_min`` to 1): the mean-field wall, where
  :mod:`.selectivity_bound` says B ≤ 1 in series;
* ``well``: no charge, and a Ca²⁺-only energy −d over one band of the path
  (the span, the filter, the gate, or filter to gate). This is binding that
  is not the mean potential. Vais's Cl⁻ experiment holds no Ca²⁺, so it
  reads the uncharged pore's P_Cl:P_K whatever the well;
* ``well + charge``: the best well with the deposit's charge scaled in
  (a binding site is itself charged);
* ``rings`` (linear response only; :func:`ring_search`): the one route
  round the series bound for point ions, parallel paths. C4 rings of
  opposite charge at one height, offset by 45°, so that cations and anions
  could each have their own side of the lumen. It reports B.

The Cl⁻ experiment depends only on the charge map, so it is solved once
per charge scale.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .charge3d import local_density
from .pnp3d import domain
from .reversal3d import _impermeant, concentration, experiments, reversal, wall
from .selectivity import ghk_ratio
from .selectivity3d import prepare
from .selectivity_bound import band_mask, bound_product, linear_ratios

__all__ = ["FAMILIES", "Candidate", "Result", "Search", "regions",
           "scales", "depths", "score", "search", "RingResult", "ring_search"]

FAMILIES = ("charge", "well", "well + charge")


@dataclass(frozen=True)
class Candidate:
    family: str
    scale: float = 0.0                 # the deposit's charge multiplied by
    region: str = ""                   # where the Ca2+ well is
    depth: float = 0.0                 # kT

    @property
    def label(self) -> str:
        parts = []
        if self.region:
            parts.append(f"Ca well {self.depth:g} kT at {self.region}")
        if self.scale or not self.region:
            parts.append(f"charge x{self.scale:g}")
        return ", ".join(parts)


@dataclass
class Result:
    candidate: Candidate
    pca_pk: float
    pcl_pk: float
    v: dict[str, float]                # experiment -> V_rev, V
    peak_ca: float                     # M in the lumen at the Ca2+ reversal
    converged: bool
    score: float = float("nan")
    bound: float = float("nan")        # B relative to the uncharged pore


@dataclass
class Search:
    name: str
    measured: tuple[float, float]      # (P_Cl:P_K, P_Ca:P_K)
    neutral: Result
    results: list[Result] = field(default_factory=list)

    @property
    def best(self) -> Result:
        return min(self.results, key=lambda r: r.score)

    def family(self, name: str) -> list[Result]:
        return [r for r in self.results if r.candidate.family == name]


def regions(summary) -> dict[str, tuple[float, float]]:
    """The bands a well may cover: the membrane span, each constriction
    ± ``lumen.constriction_half_width``, and filter to gate."""
    half = _P.value("lumen.constriction_half_width")
    zf = summary.constrictions["filter"].z
    zg = summary.constrictions["gate"].z
    lo, hi = sorted((zf, zg))
    return {"span": tuple(summary.span), "filter": (zf - half, zf + half),
            "gate": (zg - half, zg + half),
            "filter to gate": (lo - half, hi + half)}


def scales() -> np.ndarray:
    n = int(_P.value("wallsearch.scale_points"))
    return np.geomspace(_P.value("wallsearch.scale_min"), 1.0, n)


def depths() -> np.ndarray:
    step = _P.value("wallsearch.depth_step")
    return np.arange(step, _P.value("wallsearch.depth_max") + step / 2, step)


def score(pca: float, pcl: float, measured: tuple[float, float]) -> float:
    """|ln| distance to the measured pair (P_Cl:P_K, P_Ca:P_K)."""
    if not (pca > 0 and pcl > 0):
        return float("inf")
    return float(abs(np.log(pca / measured[1])) + abs(np.log(pcl / measured[0])))


class _Solver:
    """One deposit's reversal machinery, the Cl⁻ experiment memoised by
    charge scale."""

    def __init__(self, pore):
        self.pore = pore
        self.dom = domain(pore.elec, pore.vols)
        self.fixed = wall(pore).fixed
        self.cl, self.ca = experiments(False)
        self._cl: dict[float, tuple[float, float, bool]] = {}

    def _v(self, exp, fixed, excess):
        v, st, _ = reversal(self.dom, exp.species(), fixed, excess,
                            _impermeant(self.pore, exp))
        return v, st

    def chloride(self, scale: float) -> tuple[float, float, bool]:
        if scale not in self._cl:
            v, st = self._v(self.cl, self.fixed * scale, None)
            pcl = ghk_ratio(v, "Cl-", self.cl.lumen, self.cl.cytosol, {})
            self._cl[scale] = (float(pcl), v, bool(st.converged))
        return self._cl[scale]

    def read(self, c: Candidate, bands: dict) -> Result:
        pcl, v_cl, ok_cl = self.chloride(c.scale)
        excess = None
        if c.region:
            excess = {"Ca2+": -c.depth * band_mask(self.pore, *bands[c.region])}
        v, st = self._v(self.ca, self.fixed * c.scale, excess)
        pca = ghk_ratio(v, "Ca2+", self.ca.lumen, self.ca.cytosol, {"Cl-": pcl})
        conc = concentration(self.dom, st, "Ca2+", 2, excess)
        peak = float(conc[self.pore.elec.mask].max() / 1000.0)
        return Result(c, float(pca), pcl, {"Cl-": v_cl, "Ca2+": v}, peak,
                      bool(ok_cl and st.converged))


def search(st: Structure, families=FAMILIES, spacing: float | None = None,
           progress=None) -> Search:
    """Every candidate of ``families`` on one IP3R deposit at Vais's
    reversal, scored against the measured pair."""
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = prepare(st, spacing=h)
    if pore.ryr:
        raise ValueError("the search reads Vais 2010's Cl- and Ca2+ pair: "
                         "give an IP3R deposit")
    measured = (_P.value("selectivity.published_pcl_pk"),
                _P.value("selectivity.published_pca_pk"))
    solver, bands = _Solver(pore), regions(pore.summary)
    todo = []
    if "charge" in families:
        todo += [Candidate("charge", float(s)) for s in scales()]
    if "well" in families:
        todo += [Candidate("well", 0.0, r, float(d)) for r in bands
                 for d in depths()]
    neutral = solver.read(Candidate("charge", 0.0), bands)
    out = Search(st.name, measured, _finish(neutral, neutral, measured))
    for i, c in enumerate(todo):
        if progress:
            progress(i, len(todo), c.label)
        out.results.append(_finish(solver.read(c, bands), neutral, measured))
    if "well + charge" in families:
        wells = [r for r in out.family("well") if r.converged] or [
            _finish(solver.read(Candidate("well", 0.0, "span",
                                          float(depths()[0])), bands),
                    neutral, measured)]
        best = min(wells, key=lambda r: r.score).candidate
        for s in scales():
            c = Candidate("well + charge", float(s), best.region, best.depth)
            if progress:
                progress(len(todo), len(todo), c.label)
            out.results.append(_finish(solver.read(c, bands), neutral, measured))
    return out


def _finish(r: Result, neutral: Result, measured) -> Result:
    r.score = score(r.pca_pk, r.pcl_pk, measured)
    r.bound = bound_product(r.pcl_pk, r.pca_pk, neutral.pcl_pk, neutral.pca_pk)
    return r


@dataclass
class RingResult:
    z: float
    charge: float                      # e per site, each sign
    pca_pk: float                      # relative to the uncharged pore
    pcl_pk: float
    bound: float
    converged: bool


def ring_positions(summary, z: float, angle0: float = 0.0) -> np.ndarray:
    """Eight sites at the wall's radius at ``z``: the four C4 copies of one
    site and, 45° on, of its opposite partner (positive first)."""
    r = summary.profile.at(z)
    angles = angle0 + np.deg2rad(np.r_[0, 90, 180, 270, 45, 135, 225, 315])
    return np.column_stack([r * np.cos(angles), r * np.sin(angles),
                            np.full(8, z)])


def ring_search(st: Structure, spacing: float | None = None,
                progress=None) -> list[RingResult]:
    """Opposite-charge C4 rings at each height through the span, each
    charge magnitude, in linear response: B for the parallel-path route."""
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = prepare(st, spacing=h)
    base = linear_ratios(pore)
    lo, hi = pore.summary.span
    step = _P.value("wallsearch.ring_step")
    heights = np.arange(lo + step / 2, hi, step)
    qs = (1.0, _P.value("wallsearch.ring_charge_max"))
    out = []
    for i, z in enumerate(heights):
        for q in qs:
            if progress:
                progress(i, len(heights), f"z {z:+.1f} A, {q:g} e")
            charges = np.r_[np.full(4, q), np.full(4, -q)]
            fixed, _ = local_density(pore.elec, ring_positions(pore.summary, z),
                                     charges)
            r = linear_ratios(pore, fixed)
            a, b = r.pca_pk / base.pca_pk, r.pcl_pk / base.pcl_pk
            out.append(RingResult(float(z), float(q), a, b,
                                  float(b * a * a), r.converged))
    return out
