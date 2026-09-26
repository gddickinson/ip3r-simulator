"""The wall charge in the lumen's own shape: where it sits, and the ions'
equilibrium around it (Round 7.11).

The 1-D model spreads each lining group's charge along the axis (a Gaussian
of ``pore_charge.smoothing``) and over the slice's inscribed circle, then
neutralises every slice on its own (local Donnan). Round 7.6 found the real
lumen has corners the circle leaves out; the charges sit in some of them.
This module places the same charges on the voxelised lumen
(:mod:`ip3r.structure.pore_volume`) three ways, one closure each:

* ``slice`` — the 1-D model's charge per length (the same axial
  Gaussians) spread over each plane's real lumen region instead of the
  inscribed circle, and local Donnan per voxel: the 1-D closure, in the
  real cross-section.
* ``local`` — each group a 3-D Gaussian of the same width about its own
  charge centre, normalised over the lumen voxels it reaches (so the charge
  is conserved), and local Donnan per voxel. A carboxylate in a corner now
  concentrates its counter-ions in that corner.
* ``pb`` — the ``local`` density with the full nonlinear Poisson–Boltzmann
  equation in the lumen, ``ε∇²ψ = −F(Σ z c e^{−zψ/φT} + X)``, the potential
  zero on the bath voxels and no field into the protein (its permittivity
  taken as zero against the water's: the bound that confines every field
  line to the lumen). Donnan is its limit when the lumen is narrow against
  the Debye length; in a vestibule wider than that the charge is screened
  at the wall and the core stays near bulk.

The potential ``u`` (units of kT/e) is returned per voxel; a species of
valence z has equilibrium concentration ``c_bath e^{−z u}`` there, which is
all a linear-response conductance needs (:mod:`.charged3d`). Ions are points,
as in the 1-D model, and every species is counted in the neutrality of every
electrostatic voxel, as there (its own steric exclusion enters only the
conduction).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg
from scipy.spatial import cKDTree

from ..parameters import PARAMETERS as _P
from ..structure.pore_volume import PoreVolume
from ._pnp_kernels import F_FARADAY, R_GAS, _donnan_potential
from .ohmic3d import face_pairs
from .pore_charge import (AVOGADRO, ChargedGroup, PoreCharge, _centres,
                          charge_per_length)
from .radial_pb import EPS0

__all__ = ["CLOSURES_3D", "WallField", "voxel_centres", "slice_density",
           "local_density", "donnan_field", "poisson_boltzmann", "wall_field",
           "group_positions"]

#: The three placements, in the order they depart from the 1-D model.
CLOSURES_3D = ("slice", "local", "pb")

_EXP_CLIP = 40.0


@dataclass
class WallField:
    """A fixed-charge density on the lumen and the ions' potential in it."""

    closure: str
    fixed: np.ndarray            # grid shape, mol/m³ signed (ρ_fixed / F)
    potential: np.ndarray        # grid shape, kT/e; 0 outside ``mask``
    mask: np.ndarray             # the electrostatic volume
    placed: float                # e on the grid (Σ X h³)
    unreached: list[str] = field(default_factory=list)
    converged: bool = True
    iterations: int = 0
    self_energy: np.ndarray | None = None   # kT per z² (Round 7.15), or none

    def energy(self, valence: int) -> np.ndarray:
        """Boltzmann energy of a species, kT, on the grid (with its image
        cost z²W when the field carries one)."""
        e = valence * self.potential
        if self.self_energy is not None:
            e = e + valence ** 2 * np.where(self.mask, self.self_energy, 0.0)
        return np.clip(e, -_EXP_CLIP, _EXP_CLIP)

    def peak(self, valence: int, bath: float) -> float:
        """Highest equilibrium concentration of a species in the lumen, M."""
        e = self.energy(valence)[self.mask]
        return float(bath * np.exp(-e.min())) if e.size else 0.0


def voxel_centres(vol: PoreVolume, mask: np.ndarray | None = None) -> np.ndarray:
    """(n, 3) frame coordinates, Å, of ``mask``'s voxels in C order."""
    mask = vol.mask if mask is None else mask
    i, j, k = np.nonzero(mask)
    return np.column_stack([vol.xs[i], vol.xs[j], vol.zs[k]])


def _molar(e_per_A3: np.ndarray) -> np.ndarray:
    return e_per_A3 * 1e30 / AVOGADRO                  # e/Å³ → mol/m³


def slice_density(vol: PoreVolume, groups: list[ChargedGroup],
                  width: float | None = None) -> np.ndarray:
    """The 1-D model's charge per length (the same axial Gaussians), spread
    uniformly over each plane's lumen region (:meth:`PoreVolume.axis_regions`):
    the 1-D closure in the real cross-section, charge conserved."""
    width = _P.value("pore_charge.smoothing") if width is None else width
    per_length = charge_per_length(groups, vol.zs, width)       # e/Å
    x = np.zeros(vol.mask.shape)
    for k, region in enumerate(vol.axis_regions()):
        if region is None or per_length[k] == 0.0:
            continue
        area = np.count_nonzero(region) * vol.spacing ** 2
        x[:, :, k][region] = _molar(np.array(per_length[k] / area))
    return x


def group_positions(st, frame, groups: list[ChargedGroup]) -> np.ndarray:
    """Each lining group's charge centre in frame coordinates, Å."""
    centres = _centres(st)
    xyz = np.array([centres[(g.chain, g.res_seq)][1] for g in groups]
                   ).reshape(-1, 3)
    return frame.to_frame(xyz) if len(xyz) else xyz


def local_density(vol: PoreVolume, positions: np.ndarray,
                  charges: np.ndarray, labels: list[str] | None = None,
                  width: float | None = None
                  ) -> tuple[np.ndarray, list[str]]:
    """Each charge a 3-D Gaussian (``pore_charge.smoothing``) over the lumen
    voxels within ``charge3d.gaussian_reach`` widths of it, normalised over
    those voxels. Returns the density and the labels of charges that reach
    none (left out, and reported)."""
    width = _P.value("pore_charge.smoothing") if width is None else width
    reach = _P.value("charge3d.gaussian_reach") * width
    labels = labels or [str(i) for i in range(len(charges))]
    centres = voxel_centres(vol)
    tree = cKDTree(centres)
    flat = np.zeros(len(centres))
    unreached = []
    volume = vol.spacing ** 3
    for p, q, name in zip(np.asarray(positions, float), charges, labels):
        idx = np.asarray(tree.query_ball_point(p, reach), dtype=np.int64)
        if idx.size == 0:
            unreached.append(name)
            continue
        w = np.exp(-0.5 * np.sum((centres[idx] - p) ** 2, axis=1) / width ** 2)
        flat[idx] += q * w / (w.sum() * volume)
    x = np.zeros(vol.mask.shape)
    x[vol.mask] = _molar(flat)
    return x, unreached


def _bath(species) -> tuple[np.ndarray, np.ndarray]:
    valences = np.array([s.valence for s in species], dtype=float)
    conc = np.array([s.concentration * 1000.0 for s in species])  # mol/m³
    return valences, conc


def donnan_field(vol: PoreVolume, fixed: np.ndarray, species) -> np.ndarray:
    """Local Donnan potential per voxel, kT/e (the 1-D closure, voxel by voxel)."""
    valences, conc = _bath(species)
    x = fixed[vol.mask]
    u = np.zeros(vol.mask.shape)
    if x.size:
        u[vol.mask] = _donnan_potential(valences, conc[:, None] * np.ones(x.size),
                                        x, 1.0)
    return u


def poisson_boltzmann(vol: PoreVolume, fixed: np.ndarray, species,
                      permittivity: float | None = None,
                      initial: np.ndarray | None = None
                      ) -> tuple[np.ndarray, bool, int]:
    """Nonlinear Poisson–Boltzmann on the lumen, Newton with a capped step.

    Seven-point finite volumes on ``vol.mask`` (no flux into protein), u = 0
    on the bath voxels. The Jacobian ``L + h²A Σ z² c e^{−zu}`` is symmetric
    positive definite, solved by preconditioned conjugate gradients.
    Returns ``(u, converged, iterations)``."""
    eps = (_P.value("permeation.permittivity_pore") if permittivity is None
           else permittivity)
    temperature = _P.value("permeation.temperature")
    valences, conc = _bath(species)
    mask = vol.mask
    rows, cols = face_pairs(mask)
    n = int(mask.sum())
    free = ~(vol.top | vol.bottom)[mask]
    k = -np.ones(n, np.int64)
    k[free] = np.arange(int(free.sum()))
    ff = free[rows] & free[cols]
    degree = np.bincount(rows, minlength=n).astype(float)[free]
    offdiag = sparse.csr_matrix((-np.ones(int(ff.sum())), (k[rows[ff]], k[cols[ff]])),
                                shape=(int(free.sum()),) * 2)
    h = vol.spacing * 1e-10
    scale = h * h * F_FARADAY ** 2 / (EPS0 * eps * R_GAS * temperature)
    x = fixed[mask][free]
    u_full = np.zeros(n) if initial is None else initial[mask].astype(float)
    u_full[~free] = 0.0
    u = u_full[free]
    # L u for the free voxels: fixed (bath) neighbours are at u = 0.
    tol = _P.value("charge3d.newton_tolerance")
    cap = _P.value("charge3d.max_step")
    cg_tol = _P.value("pore3d.cg_tolerance")
    cg_max = int(_P.value("pore3d.cg_max_iterations"))
    converged, used = False, 0
    for used in range(1, int(_P.value("charge3d.newton_max_iterations")) + 1):  # noqa: B007
        arg = np.clip(-valences[:, None] * u[None, :], -_EXP_CLIP, _EXP_CLIP)
        boltz = conc[:, None] * np.exp(arg)
        g = degree * u + offdiag @ u - scale * ((valences[:, None] * boltz).sum(0) + x)
        diag = degree + scale * (valences[:, None] ** 2 * boltz).sum(0)
        jac = offdiag + sparse.diags(diag)
        step, info = cg(jac, -g, M=sparse.diags(1.0 / diag), rtol=cg_tol,
                        maxiter=cg_max)
        biggest = float(np.max(np.abs(step))) if step.size else 0.0
        if biggest > cap:
            step *= cap / biggest
        u = u + step
        if biggest < tol and info == 0:
            converged = True
            break
    out = np.zeros(mask.shape)
    u_full[free] = u
    out[mask] = u_full
    return out, converged, used


def wall_field(vol: PoreVolume, closure: str, species, charge: PoreCharge,
               positions: np.ndarray | None = None,
               width: float | None = None,
               permittivity: float | None = None) -> WallField:
    """The fixed charge of ``charge`` on ``vol`` and its equilibrium
    potential under ``closure`` (one of :data:`CLOSURES_3D`). ``positions``
    are the groups' charge centres (frame, Å), needed by ``local``/``pb``."""
    if closure not in CLOSURES_3D:
        raise ValueError(f"closure must be one of {CLOSURES_3D}, not {closure!r}")
    unreached: list[str] = []
    if closure == "slice":
        fixed = slice_density(vol, charge.groups, width=width)
    else:
        q = np.array([g.charge for g in charge.groups])
        fixed, unreached = local_density(
            vol, positions, q, [f"{g.label()}/{g.chain}" for g in charge.groups],
            width=width)
    placed = float(fixed.sum() * vol.spacing ** 3 * 1e-30 * AVOGADRO)
    u = donnan_field(vol, fixed, species)
    converged, used = True, 0
    if closure == "pb":
        u, converged, used = poisson_boltzmann(vol, fixed, species,
                                               permittivity=permittivity,
                                               initial=u)
    return WallField(closure, fixed, u, vol.mask, placed, unreached,
                     converged, used)
