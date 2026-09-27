"""What a mean field can do to P_Cl:P_K and P_Ca:P_K at once (Round 7.25).

Round 7.23 found that the deposits' charged walls shut Cl⁻ out (P_Cl:P_K
0.003 against Vais 2010's 0.27) long before they select Ca²⁺, and that
scaling the charge moves the pair along a path that never nears the
measured one. This module shows that the path is not an accident of the
walls' shape.

**The series bound.** In linear response a species of valence z crossing a
single-file path (cross-section A(z), energy zψ(z) in kT, point ions) has
permeability P_z ∝ D_z / ∫ e^{zψ}/A dz. Relative to the same path
uncharged, ρ_z = 1 / E_w[X^z] with X = e^{ψ} and w the normalised 1/A
weight. Hölder's inequality (log E[X^t] is convex in t) gives

    ρ_Cl ρ_Ca² ≤ ρ_K³,   i.e.   B = b a² ≤ 1,

with a = ρ_Ca/ρ_K and b = ρ_Cl/ρ_K the two measured ratios each divided by
the uncharged pore's own. Equality holds only for a uniform potential
(a tract of constant charge density). No profile of point-ion potential in
series, however charged or placed, raises P_Ca:P_K without lowering
P_Cl:P_K at least as the square: Vais's pair needs B ≈ 3 × 10³.

Three things lie outside the bound, and the round measures each:

* **parallel paths** in 3-D. Two paths of opposite sign side by side give
  log-convex conductance (:func:`parallel_ratios`), so B > 1. A lumen a
  few Debye lengths wide cannot hold them apart (:mod:`.wall_search`'s
  ring family);
* **non-linearity** at reversal. The wall's screening differs between
  the Cl⁻ and the Ca²⁺ experiment (7.23's scan reaches B 1.2);
* **species-specific energies**, i.e. binding that is not the mean
  potential. A Ca²⁺-only well leaves the Cl⁻ experiment, which holds no
  Ca²⁺, untouched. In linear response its ceiling is the share of Ca²⁺'s
  resistance it covers: :func:`well_ceiling`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["profile_ratios", "bound_product", "parallel_ratios",
           "LinearRatios", "linear_ratios", "band_mask", "well_ceiling"]


def profile_ratios(psi: np.ndarray, area: np.ndarray, dz: float = 1.0
                   ) -> dict[int, float]:
    """1-D linear response: each valence's permeability relative to the
    same path uncharged (ψ in kT/e, point ions, one diffusivity)."""
    psi, w = np.asarray(psi, float), dz / np.asarray(area, float)
    return {z: float(w.sum() / (w * np.exp(z * psi)).sum()) for z in (1, 2, -1)}


def bound_product(pcl: float, pca: float, pcl0: float, pca0: float) -> float:
    """B = b a²: the Cl⁻ and Ca²⁺ ratios each relative to the uncharged
    pore's (same protocol). At most 1 for any series point-ion profile."""
    return float((pcl / pcl0) * (pca / pca0) ** 2)


def parallel_ratios(paths: list[tuple[np.ndarray, np.ndarray]], dz: float = 1.0
                    ) -> dict[int, float]:
    """Relative permeabilities of 1-D paths ``(psi, area)`` in parallel."""
    def g(z, zero=False):
        return sum(1.0 / (dz / np.asarray(a, float)
                          * np.exp(0.0 if zero else z * np.asarray(p, float))).sum()
                   for p, a in paths)
    return {z: float(g(z) / g(z, zero=True)) for z in (1, 2, -1)}


@dataclass
class LinearRatios:
    pca_pk: float
    pcl_pk: float
    converged: bool


def linear_ratios(pore, fixed: np.ndarray | None = None,
                  extra: dict[str, np.ndarray] | None = None) -> LinearRatios:
    """Round 7.19's 3-D linear response (``pb``, point ions) for any fixed
    charge map on ``pore.elec`` (mol/m³, None = uncharged), plus
    ``extra`` species-specific energies (kT on the grid)."""
    from .charge3d import poisson_boltzmann
    from .selectivity3d import _read
    shape = pore.elec.mask.shape
    energies = {s.name: np.zeros(shape) for s in pore.species}
    ok = True
    if fixed is not None and np.any(fixed):
        u, ok, _ = poisson_boltzmann(pore.elec, fixed, pore.species)
        energies = {s.name: s.valence * u for s in pore.species}
    for name, e in (extra or {}).items():
        energies[name] = energies[name] + e
    r = _read(pore, "search", energies)
    return LinearRatios(r.ratio, r.pcl_pk, bool(ok and r.converged))


def band_mask(pore, z_lo: float, z_hi: float) -> np.ndarray:
    """The electrostatic lumen voxels with z_lo ≤ z ≤ z_hi (float 0/1)."""
    zs = pore.elec.zs
    band = (zs >= z_lo) & (zs <= z_hi)
    return (pore.elec.mask & band[None, None, :]).astype(float)


def well_ceiling(pore, z_lo: float, z_hi: float,
                 depth: float | None = None) -> float:
    """How far a Ca²⁺-only well over [z_lo, z_hi] can raise P_Ca:P_K in
    linear response, relative to the uncharged pore: at ``depth`` kT
    (``wallsearch.ceiling_depth``) the well's own resistance is gone, and
    1 / ceiling is the share of Ca²⁺'s resistance left outside it. NaN
    when the weighted Laplace solve does not converge."""
    from ..parameters import PARAMETERS as _P
    depth = _P.value("wallsearch.ceiling_depth") if depth is None else depth
    base = linear_ratios(pore)
    deep = linear_ratios(pore, extra={"Ca2+": -depth * band_mask(pore, z_lo, z_hi)})
    if not (base.converged and deep.converged):
        return float("nan")
    return deep.pca_pk / base.pca_pk
