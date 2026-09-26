"""The charge behind the wall: Poisson–Boltzmann with the protein in it
(Round 7.13).

Round 7.11's closures place only *lining* charges, and only in the lumen:
``local`` and ``pb`` normalise each group's Gaussian over the lumen voxels it
reaches, and ``pb`` lets no field into the protein. That leaves a salt
bridge with two readings, both wrong. *Full* keeps D2478 and ignores R2471′,
whose guanidinium sits a few Å behind it. *Paired* removes both, as if the
pair had no field. Both routes of Round 7.4 put D2478's pKa far below 7.3
(network 2.0–2.2 with the bridge counted, PROPKA 5.2), so the pair is an
ion pair, +1 and −1, at two places. Its field in the lumen is a dipole's,
seen through a low-permittivity protein.

This closure (``dielectric``) solves

    −∇·(ε ∇u) = (F² h² / ε0 R T) (Σ_s z_s c_s e^{−z_s u} [lumen] + X)

on the **whole box**, not only on the lumen: ε is the pore water's
(``permeation.permittivity_pore``) on every voxel open to the ion probe
(membrane sealed, as in :mod:`ip3r.structure.pore_volume`), the protein's
(``dielectric.eps_protein``) elsewhere, harmonic mean across a face. Mobile
ions live only in the conducting lumen, and the potential is zero on its
bath faces. Every other box face is insulating. X is every modelled
ionisable group in the box (formal charge, or a protonation reading), each a
Gaussian of ``dielectric.charge_width`` about its own charge centre over all
voxels (charge conserved): a charge behind the wall stays behind it.

``scope="lining"`` keeps only the lining groups (Round 7.11's set, now at
their own positions with the protein's permittivity); ``"all"`` adds every
other group, the bridge partners among them.

**Left out**: the Born (image) cost of an ion near the low-ε wall, which
every closure here omits alike. A calibration against :func:`.charge3d.poisson_boltzmann`
(ε_protein → 0 with the charge in the lumen) is in the tests.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg

from ..parameters import PARAMETERS as _P
from ..structure.pore import _heavy
from ..structure.pore_volume import PoreVolume, open_voxels
from ._pnp_kernels import F_FARADAY, R_GAS
from .charge3d import _EXP_CLIP, WallField, _bath, _molar
from .ohmic3d import face_pairs
from .pore_charge import AVOGADRO, CHARGE, PoreCharge, _centres
from .radial_pb import EPS0

__all__ = ["DIELECTRIC", "SCOPES", "solvent_mask", "permittivity_map", "point_density",
           "box_charges", "dielectric_pb", "dielectric_field"]

#: The closure's name beside :data:`.charge3d.CLOSURES_3D`.
DIELECTRIC = "dielectric"

#: Which groups carry charge: the lining set, or every group in the box.
SCOPES = ("lining", "all")


def solvent_mask(points: np.ndarray, vdw: np.ndarray, vol: PoreVolume) -> np.ndarray:
    """Voxels open to ``vol.probe`` (atoms in frame coordinates), the
    membrane sealed beyond ``pore3d.seal_radius`` as the lumen is: water."""
    solvent = open_voxels(points, vdw, vol.xs, vol.zs, vol.probe)
    x, y = np.meshgrid(vol.xs, vol.xs, indexing="ij")
    inside = (vol.zs >= vol.membrane[0]) & (vol.zs <= vol.membrane[1])
    solvent[:, :, inside] &= (np.hypot(x, y) <= _P.value("pore3d.seal_radius")
                              )[:, :, None]
    return solvent | vol.mask


def permittivity_map(solvent: np.ndarray, eps_water: float | None = None,
                     eps_protein: float | None = None) -> np.ndarray:
    """Relative permittivity per voxel: water where ``solvent``, protein
    (and sealed membrane) elsewhere."""
    w = _P.value("permeation.permittivity_pore") if eps_water is None else eps_water
    p = _P.value("dielectric.eps_protein") if eps_protein is None else eps_protein
    return np.where(solvent, float(w), float(p))


def point_density(vol: PoreVolume, positions: np.ndarray, charges: np.ndarray,
                  width: float | None = None) -> np.ndarray:
    """Each charge a Gaussian of ``width`` about its position, over every
    voxel within ``charge3d.gaussian_reach`` widths (lumen or protein),
    normalised over them: mol/m³ on the whole grid."""
    width = _P.value("dielectric.charge_width") if width is None else width
    reach = _P.value("charge3d.gaussian_reach") * width
    x = np.zeros(vol.mask.shape)
    h = vol.spacing
    for p, q in zip(np.asarray(positions, float).reshape(-1, 3), charges):
        lo = [np.searchsorted(a, c - reach) for a, c in ((vol.xs, p[0]), (vol.xs, p[1]), (vol.zs, p[2]))]
        hi = [np.searchsorted(a, c + reach, side="right") for a, c in ((vol.xs, p[0]), (vol.xs, p[1]), (vol.zs, p[2]))]
        if any(b <= a for a, b in zip(lo, hi)):
            continue                              # outside the box
        gx, gy, gz = np.meshgrid(vol.xs[lo[0]:hi[0]], vol.xs[lo[1]:hi[1]],
                                 vol.zs[lo[2]:hi[2]], indexing="ij")
        d2 = (gx - p[0]) ** 2 + (gy - p[1]) ** 2 + (gz - p[2]) ** 2
        w = np.where(d2 <= reach ** 2, np.exp(-0.5 * d2 / width ** 2), 0.0)
        if w.sum() == 0.0:
            continue
        x[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] += q * w / (w.sum() * h ** 3)
    return _molar(x)


def box_charges(st, frame, vol: PoreVolume, charge: PoreCharge, scope: str = "all",
                neutralise: frozenset[int] = frozenset(),
                charges: dict | None = None
                ) -> tuple[list[str], np.ndarray, np.ndarray]:
    """``(labels, positions, charges)`` of the groups carrying charge.
    ``lining`` = ``charge.groups`` as they stand (paired, neutralised,
    re-charged there); ``all`` adds every other modelled group whose centre
    lies in the box, except residue numbers in ``neutralise`` and partners
    of the bridges ``charge`` dropped."""
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES}, not {scope!r}")
    centres = _centres(st)
    keys = [(g.chain, g.res_seq) for g in charge.groups]
    q = [g.charge for g in charge.groups]
    if scope == "all":
        dropped = {k for b in charge.bridged for k in b.members()}
        lo = np.array([vol.xs[0], vol.xs[0], vol.zs[0]])
        hi = np.array([vol.xs[-1], vol.xs[-1], vol.zs[-1]])
        taken = set(keys)
        for k, (name, xyz) in sorted(centres.items()):
            if k in taken or k in dropped or k[1] in neutralise:
                continue
            f = frame.to_frame(xyz[None, :])[0]
            if np.all(f >= lo) and np.all(f <= hi):
                keys.append(k)
                q.append(CHARGE[name] if charges is None
                         else float(charges.get(k, CHARGE[name])))
    pos = frame.to_frame(np.array([centres[k][1] for k in keys]).reshape(-1, 3)) \
        if keys else np.zeros((0, 3))
    labels = [f"{centres[k][0]}{k[1]}/{k[0]}" for k in keys]
    return labels, pos, np.asarray(q, float)


def dielectric_pb(vol: PoreVolume, eps: np.ndarray, fixed: np.ndarray, species,
                  initial: np.ndarray | None = None
                  ) -> tuple[np.ndarray, bool, int]:
    """Nonlinear Poisson–Boltzmann over the whole grid with permittivity
    ``eps`` per voxel; ions only on ``vol.mask``; u = 0 on ``vol.top`` and
    ``vol.bottom``, every other face insulating. Newton; the step is capped
    by its largest change *in the lumen* (the protein's part is linear).
    Returns ``(u, converged, iterations)``, u in kT/e."""
    temperature = _P.value("permeation.temperature")
    valences, conc = _bath(species)
    shape = vol.mask.shape
    every = np.ones(shape, bool)
    rows, cols = face_pairs(every)
    flat_eps = eps.ravel()
    w = 2.0 * flat_eps[rows] * flat_eps[cols] / (flat_eps[rows] + flat_eps[cols])
    n = flat_eps.size
    free = ~(vol.top | vol.bottom).ravel()
    k = -np.ones(n, np.int64)
    k[free] = np.arange(int(free.sum()))
    degree = np.bincount(rows, weights=w, minlength=n)[free]
    ff = free[rows] & free[cols]
    offdiag = sparse.csr_matrix((-w[ff], (k[rows[ff]], k[cols[ff]])),
                                shape=(int(free.sum()),) * 2)
    h = vol.spacing * 1e-10
    scale = h * h * F_FARADAY ** 2 / (EPS0 * R_GAS * temperature)
    ions = vol.mask.ravel()[free]
    x = fixed.ravel()[free]
    u = np.zeros(int(free.sum())) if initial is None else initial.ravel()[free].astype(float)
    tol = _P.value("charge3d.newton_tolerance")
    cap = _P.value("charge3d.max_step")
    cg_tol = _P.value("pore3d.cg_tolerance")
    cg_max = int(_P.value("pore3d.cg_max_iterations"))
    converged, used = False, 0
    for used in range(1, int(_P.value("charge3d.newton_max_iterations")) + 1):  # noqa: B007
        arg = np.clip(-valences[:, None] * u[None, ions], -_EXP_CLIP, _EXP_CLIP)
        boltz = conc[:, None] * np.exp(arg)
        rho = x.copy()
        rho[ions] += (valences[:, None] * boltz).sum(0)
        g = degree * u + offdiag @ u - scale * rho
        diag = degree.copy()
        diag[ions] += scale * (valences[:, None] ** 2 * boltz).sum(0)
        jac = offdiag + sparse.diags(diag)
        step, info = cg(jac, -g, x0=None, M=sparse.diags(1.0 / diag),
                        rtol=cg_tol, maxiter=cg_max)
        biggest = float(np.max(np.abs(step[ions]))) if ions.any() else 0.0
        if biggest > cap and used > 1:
            step *= cap / biggest
        u = u + step
        if biggest < tol and info == 0:
            converged = True
            break
    out = np.zeros(n)
    out[free] = u
    return out.reshape(shape), converged, used


def dielectric_field(st, frame, vol: PoreVolume, species, charge: PoreCharge,
                     scope: str = "all", neutralise: frozenset[int] = frozenset(),
                     charges: dict | None = None,
                     eps_protein: float | None = None,
                     width: float | None = None) -> WallField:
    """The ``dielectric`` closure's :class:`.charge3d.WallField` on ``vol``:
    ``potential`` is the lumen's (zero elsewhere, as the other closures),
    ``fixed`` the whole grid's charge; ``unreached`` names nothing (no
    charge needs a lumen voxel). ``placed`` counts the charge in the box."""
    m = _heavy(st, include_hetero=False)
    solvent = solvent_mask(frame.to_frame(st.xyz[m]), st.vdw_radii()[m], vol)
    eps = permittivity_map(solvent, eps_protein=eps_protein)
    labels, pos, q = box_charges(st, frame, vol, charge, scope, neutralise, charges)
    fixed = point_density(vol, pos, q, width)
    u, ok, used = dielectric_pb(vol, eps, fixed, species)
    placed = float(fixed.sum() * vol.spacing ** 3 * 1e-30 * AVOGADRO)
    return WallField(DIELECTRIC, fixed, np.where(vol.mask, u, 0.0), vol.mask,
                     placed, [], ok, used)
