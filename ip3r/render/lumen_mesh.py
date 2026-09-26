"""The lumen as a surface, coloured by the potential (Round 7.10).

Contours the conducting voxels of a :class:`~ip3r.physics.lumen_field.LumenField`
inside S0's window and within ``display.lumen_radius`` of the axis. The mask
is smoothed by ``display.lumen_smoothing`` for drawing only; every number the
panel prints comes from the unsmoothed mask. Each vertex takes φ from its
nearest conducting voxel, on the fixed ramp: 0 (luminal bath) blue, 1
(cytosolic bath) red. φ is a fraction of the applied voltage, so the scale
needs no auto-ranging.

Round 7.12 adds the charged readings (:mod:`ip3r.physics.lumen_charge`).
Each vertex keeps the voxel it reads from (``source``), so a colouring is
swapped by :meth:`LumenMesh.sample` without contouring again: the K+ drop
on the same 0–1 ramp, or the wall potential on a fixed diverging scale of
± ``display.lumen_potential_range`` kT/e, cation wells blue and repulsive
regions red (higher potential red, as for the drop).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from ..parameters import PARAMETERS as _P
from .colormaps import ramp

__all__ = ["LumenMesh", "lumen_mesh", "drawn_mask", "wall_colors",
           "COLOURINGS"]

#: What the surface can be coloured by: key -> label.
COLOURINGS = {"drop": "voltage drop (0 lumen, 1 cytosol)",
              "wall": "wall potential at equilibrium (kT/e)"}


class LumenMesh:
    """Triangles in the deposit's (lab) coordinates."""

    def __init__(self, positions, normals, colors, indices, phi, source=None):
        self.positions, self.normals = positions, normals
        self.colors, self.indices, self.phi = colors, indices, phi
        #: Per vertex, the grid index of the voxel its value is read from.
        self.source = source

    def sample(self, grid: np.ndarray) -> np.ndarray:
        """``grid``'s value at each vertex (grid shaped like the volume)."""
        return np.asarray(grid)[self.source]

    @property
    def n_triangles(self) -> int:
        return len(self.indices)


def drawn_mask(field, radius: float | None = None) -> np.ndarray:
    """The conducting voxels the surface encloses."""
    radius = _P.value("display.lumen_radius") if radius is None else radius
    vol = field.volume
    lo, hi = field.window
    inside = (vol.zs >= lo) & (vol.zs <= hi)
    return vol.mask & (vol.radius() <= radius)[:, :, None] & inside[None, None, :]


def lumen_mesh(field, frame, radius: float | None = None,
               smoothing: float | None = None) -> LumenMesh | None:
    """The drawn surface, or None when nothing conducts."""
    from skimage.measure import marching_cubes
    smoothing = (_P.value("display.lumen_smoothing") if smoothing is None
                 else smoothing)
    show = drawn_mask(field, radius)
    if not show.any():
        return None
    vol = field.volume
    h = vol.spacing
    grid = np.pad(show, 1).astype(np.float32)
    if smoothing > 0:
        grid = ndimage.gaussian_filter(grid, smoothing)
    verts, faces, normals, _ = marching_cubes(grid, 0.5)
    index = verts - 1.0                               # undo the padding
    origin = np.array([vol.xs[0], vol.xs[0], vol.zs[0]])
    local = origin + index * h
    # φ of the nearest drawn voxel (the surface lies half a voxel outside).
    near = ndimage.distance_transform_edt(~show, return_distances=False,
                                          return_indices=True)
    ijk = np.clip(np.rint(index).astype(int), 0, np.array(show.shape) - 1)
    src = tuple(near[a][ijk[:, 0], ijk[:, 1], ijk[:, 2]] for a in range(3))
    phi = field.phi[src]
    positions = frame.from_frame(local).astype(np.float32)
    lab_normals = normals @ frame.basis.T            # outward from the lumen
    return LumenMesh(positions, lab_normals.astype(np.float32), ramp(phi),
                     faces.astype(np.int32), phi, src)


def wall_colors(u: np.ndarray, span: float | None = None) -> np.ndarray:
    """The wall potential (kT/e) on the fixed diverging ramp: −span blue,
    0 pale, +span red; NaN grey."""
    span = _P.value("display.lumen_potential_range") if span is None else span
    return ramp(0.5 + 0.5 * np.asarray(u, float) / span)
