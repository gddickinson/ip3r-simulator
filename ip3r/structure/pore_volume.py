"""The ion-accessible volume around the pore, as voxels (Round 7.6).

The 1-D profile (:mod:`.pore`) keeps one number per slice: the radius of the
largest circle about the axis that clears every atom within ±``pore.slab``.
That is the right quantity for "is it shut?", but a conductance integral
wants the *area* an ion can use, and a C4 pore is not a circle: the lumen
between four helices has corners the inscribed circle leaves out. This module
keeps the shape.

A voxel is open to an ion of radius ``probe`` when its centre lies at least
``vdW + probe`` from every protein heavy atom (hard spheres, the same
exclusion the 1-D model applies to the inscribed radius). The grid is a box
of half-width ``pore3d.box_half_width`` around the axis, reaching
``pore3d.bath_margin`` beyond each end of the pore-domain span.

**The membrane.** Inside the span the atom model has no lipid: outside the
pore-domain helices there is empty space that in a cell is bilayer. So there
only voxels within ``pore3d.seal_radius`` of the axis are kept. Beyond the
span, the box faces are bath: an ion that leaves the pore sideways, through a
window between cytosolic domains, reaches the bath as surely as one that
leaves along the axis. The conducting volume is every connected component
(face neighbours) that touches both baths; one touching a single bath, or
neither, is dropped. A dead end joined to the conducting volume stays in it
and carries no current (its potential is uniform).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .pore import _heavy
from .symmetry import Frame

__all__ = ["PoreVolume", "open_voxels", "pore_volume", "volume_from_atoms",
           "volume_from_mask"]

_NEIGHBOURS = 8          # nearest atoms tested per voxel (radii differ)


@dataclass
class PoreVolume:
    mask: np.ndarray          # (nx, nx, nz) bool: the conducting voxels
    xs: np.ndarray            # x and y grid coordinates, Å (frame)
    zs: np.ndarray            # z grid coordinates, Å (frame; cytosol +z)
    spacing: float            # Å
    top: np.ndarray           # bool, same shape: voxels held at the cytosolic bath
    bottom: np.ndarray        # bool: voxels held at the luminal bath
    membrane: tuple[float, float]
    probe: float              # Å, the ion radius the volume was cut for
    accessible: int = 0       # open voxels before the connectivity step

    @property
    def conducting(self) -> int:
        return int(self.mask.sum())

    def slice_area(self) -> np.ndarray:
        """The lumen's area in every z plane, Å²: the in-plane connected
        region of conducting voxels holding the open voxel nearest the axis
        (crevices reached only through another plane are left out). Zero
        where no conducting voxel lies within ``pore3d.seal_radius``."""
        return np.array([0.0 if r is None else np.count_nonzero(r) * self.spacing ** 2
                         for r in self.axis_regions()])

    def axis_regions(self) -> list[np.ndarray | None]:
        """Per z plane, the (nx, nx) mask of :meth:`slice_area`'s region, or
        None where the plane has no conducting voxel within the seal."""
        rho = self.radius()
        seal = _P.value("pore3d.seal_radius")
        out: list[np.ndarray | None] = []
        for k in range(len(self.zs)):
            plane = self.mask[:, :, k]
            if not plane.any():
                out.append(None)
                continue
            where = np.where(plane, rho, np.inf)
            i = np.unravel_index(np.argmin(where), where.shape)
            if where[i] > seal:
                out.append(None)
                continue
            labels, _ = ndimage.label(plane)
            out.append(labels == labels[i])
        return out

    def radius(self) -> np.ndarray:
        """Cylindrical radius of every (x, y) column, Å."""
        x, y = np.meshgrid(self.xs, self.xs, indexing="ij")
        return np.hypot(x, y)


def open_voxels(points: np.ndarray, vdw: np.ndarray, xs: np.ndarray,
                zs: np.ndarray, probe: float) -> np.ndarray:
    """Voxel centres (``xs`` × ``xs`` × ``zs``) clearing every atom by
    ``vdW + probe``. ``points`` are in frame coordinates, Å."""
    reach = float(vdw.max()) + probe if len(vdw) else 0.0
    lo = np.array([xs[0], xs[0], zs[0]]) - reach
    hi = np.array([xs[-1], xs[-1], zs[-1]]) + reach
    near = np.all((points >= lo) & (points <= hi), axis=1)
    points, vdw = points[near], vdw[near]
    x, y, z = np.meshgrid(xs, xs, zs, indexing="ij")
    if len(points) == 0:
        return np.ones(x.shape, bool)
    centres = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    k = min(_NEIGHBOURS, len(points))
    d, i = cKDTree(points).query(centres, k=k, distance_upper_bound=reach)
    d, i = d.reshape(len(centres), k), i.reshape(len(centres), k)
    radii = np.append(vdw, 0.0)[np.minimum(i, len(vdw))]
    clash = np.isfinite(d) & (d < radii + probe)
    return ~clash.any(axis=1).reshape(x.shape)


def volume_from_mask(accessible: np.ndarray, xs: np.ndarray, zs: np.ndarray,
                     membrane: tuple[float, float], probe: float,
                     seal: float | None = None) -> PoreVolume:
    """Seal the membrane, mark the baths, keep what connects them.

    Separate from the atom step so a test can hand it a drawn shape.
    """
    seal = _P.value("pore3d.seal_radius") if seal is None else seal
    spacing = float(xs[1] - xs[0])
    ok = accessible.copy()
    x, y = np.meshgrid(xs, xs, indexing="ij")
    in_membrane = (zs >= membrane[0]) & (zs <= membrane[1])
    ok[:, :, in_membrane] &= (np.hypot(x, y) <= seal)[:, :, None]
    face = np.zeros(ok.shape[:2], bool)
    face[[0, -1], :] = True
    face[:, [0, -1]] = True
    top = np.zeros(ok.shape, bool)
    top[:, :, -1] = True
    top[:, :, zs > membrane[1]] |= face[:, :, None]
    bottom = np.zeros(ok.shape, bool)
    bottom[:, :, 0] = True
    bottom[:, :, zs < membrane[0]] |= face[:, :, None]
    labels, _ = ndimage.label(ok)
    both = np.intersect1d(np.unique(labels[top & ok]), np.unique(labels[bottom & ok]))
    mask = np.isin(labels, both[both > 0])
    return PoreVolume(mask=mask, xs=xs, zs=zs, spacing=spacing, top=top & mask,
                      bottom=bottom & mask, membrane=membrane, probe=probe,
                      accessible=int(accessible.sum()))


def volume_from_atoms(points: np.ndarray, vdw: np.ndarray,
                      membrane: tuple[float, float], probe: float,
                      spacing: float | None = None, seal: float | None = None,
                      half_width: float | None = None,
                      margin: float | None = None) -> PoreVolume:
    """The conducting volume of atoms already in frame coordinates."""
    spacing = _P.value("pore3d.spacing") if spacing is None else spacing
    half = _P.value("pore3d.box_half_width") if half_width is None else half_width
    margin = _P.value("pore3d.bath_margin") if margin is None else margin
    n = int(round(half / spacing))
    xs = np.arange(-n, n + 1) * spacing
    zs = np.arange(membrane[0] - margin, membrane[1] + margin + 1e-9, spacing)
    acc = open_voxels(points, vdw, xs, zs, probe)
    return volume_from_mask(acc, xs, zs, membrane, probe, seal)


def pore_volume(st: Structure, frame: Frame, span: tuple[float, float],
                probe: float, spacing: float | None = None,
                seal: float | None = None) -> PoreVolume:
    """A deposit's conducting volume for an ion of radius ``probe``; the
    membrane is the pore-domain ``span``. Protein atoms only, as the 1-D
    permeation profile (a bound lipid is the preparation's, not a wall)."""
    m = _heavy(st, include_hetero=False)
    return volume_from_atoms(frame.to_frame(st.xyz[m]), st.vdw_radii()[m],
                             span, probe, spacing=spacing, seal=seal)
