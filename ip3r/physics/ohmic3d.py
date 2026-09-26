"""Ohmic conductance of the pore's real shape, in three dimensions (Round 7.6).

The 1-D model (:mod:`.permeation`) integrates ``dz / (σ A(z))`` with
``A = π (r_free − r_ion)²`` from the inscribed radius. Here the same
electrolyte (same σ per species, same hard-sphere exclusion of each ion's
centre) fills the voxelised volume of :mod:`ip3r.structure.pore_volume`, and
Laplace's equation is solved in it: ``∇·(σ∇φ) = 0``, φ = 1 on the cytosolic
bath voxels and 0 on the luminal ones, no flux into protein or membrane.
Seven-point finite volumes: every pair of open face neighbours is joined by a
conductance ``σ h`` (area h², length h). Preconditioned conjugate gradients.

What it adds to the 1-D model, and nothing else: the non-circular lumen,
current spreading between wide and narrow parts, and any exit that does not
run along the axis. No charge, no concentration change: this is the neutral
reading, and a species' conductance is ``σ_s × g_s`` where ``g_s`` (metres)
depends on the shape alone.

Calibrated in ``tests/test_ohmic3d.py`` on a cylinder through a slab (the
series resistance plus Hall's two access terms), a dead-end pocket (carries
nothing, though a slice-area integral counts it), and a pore whose only exit
is sideways (shut to a profile along the axis, open here).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.channel import ChannelSummary, measure_channel
from ..structure.pore import PROFILE_MARGIN
from ..structure.pore_volume import PoreVolume, pore_volume
from .permeation import _species_conductivity, potassium_species

__all__ = ["Laplace", "geometric_conductance", "bernoulli_weight", "face_pairs", "Ohmic3D", "conductance_3d",
           "extrapolate", "cylinder_volume", "hall_cylinder"]


@dataclass
class Laplace:
    g: float                  # geometric conductance, m (conductance / σ)
    converged: bool
    iterations: int
    potential: np.ndarray = field(repr=False)   # per conducting voxel


def face_pairs(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Every ordered pair of face-neighbouring voxels of ``mask``, as indices
    into ``mask``'s voxels in C order (each pair appears both ways)."""
    index = -np.ones(mask.shape, np.int64)
    index[mask] = np.arange(int(mask.sum()))
    rows, cols = [], []
    for axis in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[axis], b[axis] = slice(0, -1), slice(1, None)
        pair = mask[tuple(a)] & mask[tuple(b)]
        ia, ib = index[tuple(a)][pair], index[tuple(b)][pair]
        rows += [ia, ib]
        cols += [ib, ia]
    return np.concatenate(rows), np.concatenate(cols)


def bernoulli_weight(ea: np.ndarray, eb: np.ndarray) -> np.ndarray:
    """Face weight between two voxels whose Boltzmann energies (kT) are
    ``ea`` and ``eb``: ``1 / mean(e^{+E})`` along the face with E linear,
    ``(eb − ea) / (e^{eb} − e^{ea})``. The Scharfetter–Gummel weight at zero
    current; ``e^{−ea}`` when the two are equal."""
    d = eb - ea
    small = np.abs(d) < 1e-9
    safe = np.where(small, 1.0, d)
    ratio = np.where(small, 1.0 - 0.5 * d, safe / np.expm1(np.where(small, 1.0, safe)))
    return np.exp(-ea) * ratio


def geometric_conductance(vol: PoreVolume, tol: float | None = None,
                          energy: np.ndarray | None = None) -> Laplace:
    """Solve Laplace's equation in ``vol`` and return ``g = I / (σ ΔV)``.

    ``energy`` (the grid's shape, kT) makes it ``∇·(e^{−E}∇μ) = 0``: the
    linear response of a species whose equilibrium concentration is the
    bath's times ``e^{−E}`` (:mod:`.charged3d`); ``g`` is then per bulk σ."""
    tol = _P.value("pore3d.cg_tolerance") if tol is None else tol
    mask = vol.mask
    n = int(mask.sum())
    if n == 0:
        return Laplace(0.0, True, 0, np.zeros(0))
    fixed = (vol.top | vol.bottom)[mask]
    value = vol.top[mask].astype(float)
    rows, cols = face_pairs(mask)
    if energy is None:
        w = np.ones(len(rows))
    else:
        e = np.asarray(energy, dtype=float)[mask]
        w = bernoulli_weight(e[rows], e[cols])
    degree = np.bincount(rows, weights=w, minlength=n)
    free = ~fixed
    k = -np.ones(n, np.int64)
    k[free] = np.arange(int(free.sum()))
    ff = free[rows] & free[cols]
    lap = sparse.csr_matrix((-w[ff], (k[rows[ff]], k[cols[ff]])),
                            shape=(int(free.sum()),) * 2)
    lap = lap + sparse.diags(degree[free])
    fx = free[rows] & fixed[cols]
    rhs = np.bincount(k[rows[fx]], weights=w[fx] * value[cols[fx]],
                      minlength=int(free.sum()))
    counter = {"n": 0}

    def _count(_):
        counter["n"] += 1

    phi = value.copy()
    info = 0
    if free.any():
        x, info = cg(lap, rhs, M=sparse.diags(1.0 / np.maximum(degree[free], 1e-300)),
                     rtol=tol, maxiter=int(_P.value("pore3d.cg_max_iterations")),
                     callback=_count)
        phi[free] = x
    # Current leaving the φ = 1 voxels into their free (or φ = 0) neighbours.
    out = vol.top[mask][rows]
    current = float(np.sum(w[out] * (phi[rows[out]] - phi[cols[out]])))
    return Laplace(g=current * vol.spacing * 1e-10, converged=info == 0,
                   iterations=counter["n"], potential=phi)


@dataclass
class Ohmic3D:
    name: str
    conductance: float                 # S, all species
    per_species: dict[str, float]      # geometric conductance g_s, m
    slice_reading: float               # S, ∫dz/(σ A_lumen) over the span ± S0's margin
    spacing: float                     # Å
    voxels: dict[str, int]             # conducting voxels per species
    converged: bool

    @property
    def conductance_pS(self) -> float:
        return self.conductance * 1e12

    @property
    def slice_pS(self) -> float:
        return self.slice_reading * 1e12


def conductance_3d(st: Structure, summary: ChannelSummary | None = None,
                   species=None, spacing: float | None = None,
                   seal: float | None = None) -> Ohmic3D:
    """A deposit's neutral conductance through its voxelised pore, in its
    family's recording bath unless ``species`` is given."""
    from .unitary import bath_for
    summary = summary or measure_channel(st)
    paralog = summary.numbering.paralog if summary.numbering else None
    species = species or potassium_species(bath=bath_for(paralog))
    temperature = _P.value("permeation.temperature")
    total, slice_total, per, voxels, ok = 0.0, 0.0, {}, {}, True
    h = _P.value("pore3d.spacing") if spacing is None else spacing
    for s in species:
        vol = pore_volume(st, summary.frame, summary.span, s.radius,
                          spacing=h, seal=seal)
        lap = geometric_conductance(vol)
        sigma = _species_conductivity(s, temperature)
        total += sigma * lap.g
        lo, hi = summary.span
        window = (vol.zs >= lo - PROFILE_MARGIN) & (vol.zs <= hi + PROFILE_MARGIN)
        area = vol.slice_area()[window] * 1e-20
        if np.all(area > 0):
            slice_total += sigma / np.sum(vol.spacing * 1e-10 / area)
        per[s.name], voxels[s.name] = lap.g, vol.conducting
        ok &= lap.converged
    return Ohmic3D(name=st.name, conductance=total, per_species=per,
                   slice_reading=slice_total, spacing=h, voxels=voxels,
                   converged=ok)


def extrapolate(coarse: Ohmic3D, fine: Ohmic3D) -> float:
    """Grid-extrapolated conductance, S: linear in the spacing (a voxel
    surface pinches narrow necks by up to half a voxel, an error of first
    order in h), through the two readings."""
    h1, h2 = coarse.spacing, fine.spacing
    return fine.conductance + (fine.conductance - coarse.conductance) * h2 / (h1 - h2)


# ---------------------------------------------------------------- calibration
def cylinder_volume(radius: float, length: float, spacing: float,
                    half_width: float, margin: float) -> PoreVolume:
    """A cylindrical hole of ``radius`` through a slab of ``length`` (Å),
    baths of ``margin`` either side: the shape with a known answer."""
    from ..structure.pore_volume import volume_from_mask
    n = int(round(half_width / spacing))
    xs = np.arange(-n, n + 1) * spacing
    zs = np.arange(-length / 2 - margin, length / 2 + margin + 1e-9, spacing)
    x, y, z = np.meshgrid(xs, xs, zs, indexing="ij")
    inside = np.abs(z) <= length / 2
    acc = ~inside | (np.hypot(x, y) <= radius)
    return volume_from_mask(acc, xs, zs, (-length / 2, length / 2), 0.0,
                            seal=radius)


def hall_cylinder(radius: float, length: float) -> float:
    """``g`` (m) of a cylinder in an infinite insulating slab between two
    half-space baths: ``1 / (L/(π a²) + 2 · 1/(4a))`` (Hall 1975)."""
    a, L = radius * 1e-10, length * 1e-10
    return 1.0 / (L / (np.pi * a * a) + 2.0 / (4.0 * a))
