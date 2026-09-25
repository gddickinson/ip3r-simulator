"""Radial Poisson-Boltzmann partition: the closure resolved across the pore.

The 1-D solver's charged closure is **local Donnan**: at each slice the wall
charge is neutralised by ions spread uniformly over the cross-section, at one
potential. That is exact when the slice is narrow against the Debye length,
and it overstates the charge wherever the slice is wide: a ring of lysines in
an 11.6 A vestibule, with a Debye length of 6-8 A, is screened near the wall
and leaves a less perturbed core the ions can pass through.

This module solves, in each slice independently, the cylindrical
Poisson-Boltzmann equation with the slice's wall charge on its wall
(Gauss's law as the boundary condition)::

    (1/r) d/dr (r dpsi/dr) = -(F / eps) sum_i z_i cbar_i exp(-z_i psi / phi_T)
    dpsi/dr (0) = 0,     dpsi/dr (R) = F X R / (2 eps)

with ``X`` the slice's fixed charge as a molar density over pi R^2 (exactly
the quantity :func:`ip3r.physics.pore_charge.map_charge` returns and the
Donnan closure uses), and ``cbar_i`` the concentration the ions would have at
zero potential (the reservoir the slice is in equilibrium with, radially).
Integrating the equation over the disc gives ``sum_i z_i <c_i> + X = 0``
for any ``cbar``: the slice is neutral as a whole, which is the Donnan
condition, but the potential is no longer uniform across it.

What each species sees is its cross-section average,
``Gamma_i = <exp(-z_i psi / phi_T)>``, which enters the axial Nernst-Planck
as a species-dependent potential ``w_i = -phi_T ln(Gamma_i) / z_i``. For a
uniform potential ``w_i = psi`` for every species, which is local Donnan;
that limit (``R << lambda_D``, i.e. ``eps -> infinity`` at fixed salt) is
what the calibration tests pin. The opposite limit (``eps -> 0``) is a
Gouy-Chapman layer at the wall around a core at the reservoir's potential. The reduction to 1-D assumes the slice is in radial equilibrium
(radial relaxation fast against axial transport) and that the axial
variation of the wall charge is slow against the radius; the second is the
weaker one, and is stated in ``docs/SCIENCE_PERM.md``.

Ions are points in this equation, as they are in the Donnan closure it
replaces, so the difference between the two is the closure alone.

The discretisation is finite volume on cells of equal width in ``r / R``:
the discrete Gauss law holds exactly (the cell weights sum to 1/2), so the
slice's electroneutrality is exact at every Newton iterate. Every slice's
Newton system is tridiagonal, and all slices are solved in one banded
LAPACK call.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded

from ..parameters import PARAMETERS as _P
from ._pnp_kernels import F_FARADAY, R_GAS

__all__ = ["RadialPartition", "radial_partition", "screening_ratio",
           "EPS0"]

#: Vacuum permittivity, F/m (CODATA 2018).
EPS0 = 8.8541878128e-12
#: As in the Donnan kernel: exponent arguments beyond 40 thermal voltages
#: are clipped so a runaway iterate stays finite.
_EXP_CLIP = 40.0


@dataclass
class RadialPartition:
    """The radial solution in every slice."""

    offset: np.ndarray        # V, (n_species, n): w_i = -phi_T ln Gamma_i / z_i
    psi: np.ndarray           # V, (n, cells): potential at the cell centres
    rho: np.ndarray           # cell centres, r / R
    converged: bool
    iterations: int

    @property
    def axis(self) -> np.ndarray:
        return self.psi[:, 0]

    @property
    def wall(self) -> np.ndarray:
        return self.psi[:, -1]


def _grid(cells: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    h = 1.0 / cells
    rho = (np.arange(cells) + 0.5) * h           # cell centres
    faces = np.arange(1, cells) * h              # interior faces
    weight = rho * h                             # disc weights, sum = 1/2
    return rho, faces, weight


def radial_partition(valences, reservoir, fixed, radius_m, thermal,
                     permittivity=None, cells=None, initial=None,
                     max_iterations: int = 80,
                     tol: float = 1e-9) -> RadialPartition:
    """Solve the radial PB equation in every slice.

    ``reservoir`` is (n_species, n) in mol/m^3, ``fixed`` (n,) in mol/m^3 as
    the Donnan closure takes it, ``radius_m`` (n,) the slice radius the
    charge sits on. ``initial`` (n, cells) is a starting potential in volts
    (the previous solution); without it Newton starts from zero. Newton
    stops when its step is below ``tol`` thermal voltages.
    """
    eps = (_P.value("permeation.permittivity_pore") if permittivity is None
           else permittivity) * EPS0
    cells = int(_P.value("permeation.radial_cells") if cells is None else cells)
    z = np.asarray(valences, dtype=float)
    cbar = np.asarray(reservoir, dtype=float)
    X = np.asarray(fixed, dtype=float)
    R = np.asarray(radius_m, dtype=float)
    n = len(X)
    rho, faces, weight = _grid(cells)
    h = 1.0 / cells

    # Dimensionless u = psi / phi_T. beta = F^2 R^2 / (eps R_gas T), m^3/mol.
    beta = F_FARADAY * R ** 2 / (eps * thermal)
    wall_flux = 0.5 * beta * X                   # rho du/drho at rho = 1

    u = (np.zeros((n, cells)) if initial is None
         else np.asarray(initial, dtype=float) / thermal)
    coupling = faces / h                         # face rho over spacing
    lower = np.zeros((n, cells))
    upper = np.zeros((n, cells))
    upper[:, :-1] = coupling
    lower[:, 1:] = coupling
    diag_lap = -(np.pad(coupling, (1, 0)) + np.pad(coupling, (0, 1)))

    converged, used = False, 0
    for used in range(1, max_iterations + 1):  # noqa: B007 - read after
        arg = np.clip(-z[:, None, None] * u[None], -_EXP_CLIP, _EXP_CLIP)
        boltz = cbar[:, :, None] * np.exp(arg)             # (s, n, cells)
        charge = (z[:, None, None] * boltz).sum(axis=0)     # sum z c
        dcharge = -(z[:, None, None] ** 2 * boltz).sum(axis=0)
        lap = diag_lap[None] * u
        lap[:, :-1] += coupling * u[:, 1:]
        lap[:, 1:] += coupling * u[:, :-1]
        residual = lap + beta[:, None] * weight[None] * charge
        residual[:, -1] += wall_flux
        jac_diag = diag_lap[None] + beta[:, None] * weight[None] * dcharge
        step = _tridiagonal(lower, jac_diag, upper, -residual)
        u = u + np.clip(step, -4.0, 4.0)
        if float(np.max(np.abs(step))) < tol:
            converged = True
            break

    arg = np.clip(-z[:, None, None] * u[None], -_EXP_CLIP, _EXP_CLIP)
    gamma = 2.0 * (np.exp(arg) * weight).sum(axis=2)        # (s, n)
    offset = -thermal * np.log(gamma) / z[:, None]
    return RadialPartition(offset, u * thermal, rho, converged, used)


def _tridiagonal(lower, diag, upper, rhs) -> np.ndarray:
    """Solve every row's tridiagonal system at once (rows uncoupled): the
    systems are laid end to end as one banded matrix with zero couplings
    between them, so a single LAPACK call does them all."""
    n, m = diag.shape
    ab = np.zeros((3, n * m))
    up = upper.copy()
    up[:, -1] = 0.0
    lo = lower.copy()
    lo[:, 0] = 0.0
    ab[0, 1:] = up.ravel()[:-1]
    ab[1] = diag.ravel()
    ab[2, :-1] = lo.ravel()[1:]
    return solve_banded((1, 1), ab, rhs.ravel()).reshape(n, m)


def screening_ratio(radius_m, species_conc, temperature,
                    permittivity=None) -> np.ndarray:
    """R / lambda_D of each slice at bath ionic strength: where it is well
    below 1, radial PB and local Donnan must agree; above 1 they need not."""
    eps = (_P.value("permeation.permittivity_pore") if permittivity is None
           else permittivity) * EPS0
    strength = 0.5 * sum(zz ** 2 * c for zz, c in species_conc)   # mol/m^3
    lam = np.sqrt(eps * R_GAS * temperature / (2.0 * F_FARADAY ** 2 * strength))
    return np.asarray(radius_m, dtype=float) / lam
