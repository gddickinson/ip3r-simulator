"""The image (Born) cost of an ion near the low-permittivity wall (Round 7.15).

Every closure so far treats an ion as a point in a mean field. A real ion
also polarises its surroundings, and near a wall of lower permittivity than
the water that polarisation pushes it away: the image force of Parsegian
1969's ion in a pore. Its energy is the ion's own reaction potential,

    W(r) = ½ z² [u_het(r) − u_bath(r)]        (kT, per unit charge squared)

where ``u_het`` is the potential at ``r`` of a unit charge at ``r`` in the
deposit's permittivity map, screened by the salt where ions reach, and
``u_bath`` the same charge in bulk water with the bath's screening. The
self-term of the charge on the grid is the same in both and cancels.

**Local solves.** The ions screen the reaction field, so each voxel's W is
solved on a cube of ``born.box_half_width`` about it, the potential zero on
its faces (linear Poisson–Boltzmann, seven-point finite volumes, harmonic
mean ε across a face). A voxel farther than ``born.reach`` from every
low-ε voxel takes W = 0 unsolved. Both are calibrated on the planar
Debye–Hückel wall (:func:`planar_image`) and scanned on 8TKF.

**The map.** Water is the region an ion's own sphere sweeps: every voxel
within the ion's radius of a voxel its centre can reach (the membrane
sealed as for the lumen). So no ion centre is closer than its radius to
the protein's dielectric. Round 7.13's closure puts the boundary at the
ion *centres* instead, which would put every wall-contact ion inside the
low-ε region (``surface="swept"`` in :mod:`.dielectric3d` is this map).
Screening charge lives where an ion centre can reach.

W depends on the geometry and the bath alone, so a deposit's field is
cached on disk (``data/cache/born``, keyed by the inputs' hash).
"""

from __future__ import annotations

import hashlib
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
from scipy import integrate, ndimage, sparse
from scipy.sparse.linalg import cg

from ..parameters import PARAMETERS as _P
from ..structure.pore import _heavy
from ..structure.pore_volume import PoreVolume
from ._pnp_kernels import F_FARADAY, R_GAS
from .dielectric3d import permittivity_map, solvent_mask, swept_mask
from .pore_charge import AVOGADRO
from .radial_pb import EPS0

__all__ = ["BornField", "ion_maps", "screening", "operator", "targets",
           "self_energy", "planar_image", "born_field", "bjerrum_vacuum"]

#: Bumped when the solve changes, so an old cache is not read.
_VERSION = "1"


def bjerrum_vacuum(temperature: float | None = None) -> float:
    """e² / (4π ε0 kT), Å."""
    t = _P.value("permeation.temperature") if temperature is None else temperature
    return F_FARADAY ** 2 / (4 * np.pi * EPS0 * R_GAS * t * AVOGADRO) * 1e10


def ion_maps(st, frame, vol: PoreVolume, eps_water: float | None = None,
             eps_protein: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """``(eps, ions)`` on ``vol``'s grid: the permittivity with water where
    an ion's sphere sweeps, and the voxels an ion centre reaches."""
    m = _heavy(st, include_hetero=False)
    ions = solvent_mask(frame.to_frame(st.xyz[m]), st.vdw_radii()[m], vol)
    water = swept_mask(ions, vol.probe, vol.spacing)
    return permittivity_map(water, eps_water, eps_protein), ions


def screening(vol: PoreVolume, ions: np.ndarray, species,
              potential: np.ndarray | None = None) -> np.ndarray:
    """The linearised screening term per voxel, h²F²Σz²c/(ε0RT): the bath's
    wherever an ion reaches, or on the lumen the mean field's
    ``c e^{−z u}`` when ``potential`` (kT/e) is given."""
    t = _P.value("permeation.temperature")
    h = vol.spacing * 1e-10
    scale = h * h * F_FARADAY ** 2 / (EPS0 * R_GAS * t)
    bulk = sum(s.valence ** 2 * s.concentration * 1000.0 for s in species)
    out = np.where(ions, scale * bulk, 0.0)
    if potential is not None:
        u = potential[vol.mask]
        local = sum(s.valence ** 2 * s.concentration * 1000.0
                    * np.exp(np.clip(-s.valence * u, -40, 40)) for s in species)
        out[vol.mask] = scale * local
    return out


def operator(eps: np.ndarray, k2: np.ndarray) -> sparse.csr_matrix:
    """−∇·(ε∇) + k² on a box, seven-point, harmonic-mean faces, potential
    zero one voxel beyond each face (that voxel's ε taken as its neighbour's)."""
    shape = eps.shape
    idx = np.arange(eps.size).reshape(shape)
    diag = k2.astype(float).ravel().copy()
    rows, cols, vals = [], [], []
    for ax in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[ax], b[ax] = slice(0, -1), slice(1, None)
        ea, eb = eps[tuple(a)], eps[tuple(b)]
        w = (2 * ea * eb / (ea + eb)).ravel()
        ia, ib = idx[tuple(a)].ravel(), idx[tuple(b)].ravel()
        rows += [ia, ib]
        cols += [ib, ia]
        vals += [-w, -w]
        diag += np.bincount(ia, w, eps.size) + np.bincount(ib, w, eps.size)
        for end in (0, -1):
            s = [slice(None)] * 3
            s[ax] = end
            diag += np.bincount(idx[tuple(s)].ravel(), eps[tuple(s)].ravel(),
                                eps.size)
    off = sparse.csr_matrix((np.concatenate(vals),
                             (np.concatenate(rows), np.concatenate(cols))),
                            shape=(eps.size,) * 2)
    return (off + sparse.diags(diag)).tocsr()


def targets(vol: PoreVolume, eps: np.ndarray, reach: float | None = None,
            eps_water: float | None = None) -> np.ndarray:
    """(n, 3) indices of the lumen voxels within ``born.reach`` of a voxel
    whose ε is not the water's (the rest take W = 0)."""
    reach = _P.value("born.reach") if reach is None else reach
    water = (_P.value("permeation.permittivity_pore") if eps_water is None
             else eps_water)
    low = eps != water
    if not low.any():
        return np.zeros((0, 3), int)
    d = ndimage.distance_transform_edt(~low, sampling=vol.spacing)
    return np.argwhere(vol.mask & (d <= reach))


# ---- the solves (module level, so worker processes can run them) ----------

_STATE: dict = {}


class _Template:
    """The sparsity of :func:`operator` on a cube of edge ``size``, built
    once: each voxel's operator is then its values poured into it."""

    def __init__(self, size: int):
        shape = (size,) * 3
        idx = np.arange(size ** 3).reshape(shape)
        ia, ib, edge = [], [], np.zeros(shape)
        for ax in range(3):
            a = [slice(None)] * 3
            b = [slice(None)] * 3
            a[ax], b[ax] = slice(0, -1), slice(1, None)
            ia.append(idx[tuple(a)].ravel())
            ib.append(idx[tuple(b)].ravel())
            for end in (0, -1):
                s = [slice(None)] * 3
                s[ax] = end
                edge[tuple(s)] += 1.0
        self.ia, self.ib = np.concatenate(ia), np.concatenate(ib)
        self.edge = edge.ravel()
        self.n = size ** 3
        diag = np.arange(self.n)
        rows = np.concatenate([self.ia, self.ib, diag])
        cols = np.concatenate([self.ib, self.ia, diag])
        order = sparse.coo_matrix((np.arange(1, len(rows) + 1, dtype=float),
                                   (rows, cols)), shape=(self.n,) * 2).tocsr()
        self.perm = order.data.astype(np.int64) - 1
        self.matrix = order

    def build(self, eps: np.ndarray, k2: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray]:
        """The operator of ``eps``/``k2`` (flat, this cube's) and its diagonal."""
        ea, eb = eps[self.ia], eps[self.ib]
        w = 2 * ea * eb / (ea + eb)
        diag = (k2 + eps * self.edge + np.bincount(self.ia, w, self.n)
                + np.bincount(self.ib, w, self.n))
        a = self.matrix.copy()
        a.data = np.concatenate([-w, -w, diag])[self.perm]
        return a, diag


def _setup(eps_pad, k2_pad, n, eps_ref, k2_ref, rhs, tol, maxiter):
    size = 2 * n + 1
    centre = np.ravel_multi_index((n, n, n), (size,) * 3)
    b = np.zeros(size ** 3)
    b[centre] = rhs
    t = _Template(size)
    _STATE.update(eps=eps_pad, k2=k2_pad, n=n, b=b, centre=centre, tol=tol,
                  maxiter=maxiter, template=t, x0=None)
    ref, diag = t.build(np.full(t.n, float(eps_ref)), np.full(t.n, float(k2_ref)))
    x = _solve(ref, diag)
    _STATE.update(ref=float(x[centre]), x0=x)


def _solve(a, diag) -> np.ndarray:
    s = _STATE
    x, _ = cg(a, s["b"], x0=s["x0"], M=sparse.diags(1.0 / diag), rtol=s["tol"],
              maxiter=s["maxiter"])
    return x


def _chunk(points: np.ndarray) -> np.ndarray:
    s = _STATE
    size = 2 * s["n"] + 1
    out = np.empty(len(points))
    for m, (i, j, k) in enumerate(points):
        cube = (slice(i, i + size), slice(j, j + size), slice(k, k + size))
        a, diag = s["template"].build(s["eps"][cube].ravel(), s["k2"][cube].ravel())
        out[m] = 0.5 * (_solve(a, diag)[s["centre"]] - s["ref"])
    return out


def _init(*args):
    _setup(*args)


def self_energy(eps: np.ndarray, k2: np.ndarray, points: np.ndarray,
                spacing: float, half: float | None = None,
                eps_ref: float | None = None, k2_ref: float | None = None,
                workers: int | None = None) -> np.ndarray:
    """W (kT per unit charge squared) at each voxel index in ``points``;
    the reference is uniform ``eps_ref`` (default the registered pore water) screened
    by ``k2_ref`` (default the map's largest). ``workers`` > 1 splits the
    points over processes."""
    half = _P.value("born.box_half_width") if half is None else half
    n = max(1, int(round(half / spacing)))
    eps_ref = (_P.value("permeation.permittivity_pore") if eps_ref is None
               else eps_ref)
    k2_ref = float(k2.max()) if k2_ref is None else k2_ref
    rhs = 4 * np.pi * bjerrum_vacuum() / spacing
    args = (np.pad(eps, n, mode="edge"), np.pad(k2, n, mode="edge"), n,
            eps_ref, k2_ref, rhs, _P.value("born.cg_tolerance"),
            int(_P.value("pore3d.cg_max_iterations")))
    points = np.asarray(points, int).reshape(-1, 3)
    workers = workers or 1
    if workers <= 1 or len(points) < 200:
        _setup(*args)
        return _chunk(points)
    chunks = np.array_split(points, workers * 8)
    with ProcessPoolExecutor(workers, initializer=_init, initargs=args) as ex:
        return np.concatenate(list(ex.map(_chunk, chunks)))


def planar_image(distance: float, eps_water: float, eps_wall: float,
                 debye: float = np.inf) -> float:
    """W (kT per z²) of a charge ``distance`` Å from a planar wall of
    ``eps_wall`` with Debye–Hückel salt (length ``debye``) on the water's
    side only: the reflected part of the half-space Green's function."""
    kappa = 0.0 if not np.isfinite(debye) else 1.0 / debye

    def f(k):
        p = np.hypot(k, kappa)
        return (k / p * np.exp(-2 * p * distance)
                * (eps_water * p - eps_wall * k) / (eps_water * p + eps_wall * k))
    return bjerrum_vacuum() / (2 * eps_water) * integrate.quad(f, 0, np.inf)[0]


@dataclass
class BornField:
    """A deposit's image self-energy on its electrostatic volume."""

    energy: np.ndarray          # grid, kT per unit z²; 0 off the lumen
    eps: np.ndarray             # the permittivity map it was solved in
    solved: int                 # voxels solved (the rest are beyond reach)
    seconds: float
    cached: bool = False
    screening: str = "bath"

    def on_axis(self, vol: PoreVolume) -> tuple[np.ndarray, np.ndarray]:
        """``(z, W)`` at the lumen voxel nearest the axis in each plane (NaN
        where the plane has none)."""
        r = vol.radius()
        out = np.full(len(vol.zs), np.nan)
        for k in range(len(vol.zs)):
            m = vol.mask[:, :, k]
            if m.any():
                i = np.argmin(np.where(m, r, np.inf))
                out[k] = self.energy[:, :, k].ravel()[i]
        return vol.zs, out


def _key(*parts) -> str:
    h = hashlib.sha256(_VERSION.encode())
    for p in parts:
        h.update(np.ascontiguousarray(p).tobytes() if isinstance(p, np.ndarray)
                 else repr(p).encode())
    return h.hexdigest()[:20]


def born_field(st, frame, vol: PoreVolume, species,
               potential: np.ndarray | None = None,
               workers: int | None = None, cache: bool = True,
               eps_protein: float | None = None) -> BornField:
    """The image self-energy on ``vol`` (bath screening, or the mean field's
    around ``potential``). Cached on disk by the inputs' hash."""
    from ..config import CACHE_DIR
    eps, ions = ion_maps(st, frame, vol, eps_protein=eps_protein)
    k2 = screening(vol, ions, species, potential)
    k2_bath = float(screening(vol, ions, species).max())
    pts = targets(vol, eps)
    half = _P.value("born.box_half_width")
    tol = _P.value("born.cg_tolerance")
    kind = "bath" if potential is None else "local"
    path = CACHE_DIR / "born" / f"{st.name}_{_key(eps, k2, pts, half, tol, vol.spacing, _P.value('permeation.temperature'))}.npz"
    if cache and path.exists():
        with np.load(path) as f:
            return BornField(f["energy"], eps, int(f["solved"]), 0.0, True, kind)
    t = time.time()
    w = np.zeros(vol.mask.shape)
    if len(pts):
        w[tuple(pts.T)] = self_energy(eps, k2, pts, vol.spacing, half,
                                      k2_ref=k2_bath,
                                      workers=workers or os.cpu_count())
    out = BornField(w, eps, len(pts), time.time() - t, False, kind)
    if cache:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, energy=w, solved=len(pts))
    return out
