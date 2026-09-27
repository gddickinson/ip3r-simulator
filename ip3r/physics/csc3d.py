"""Charge–space competition in the 3-D lumen (Round 7.19).

Round 7.17 gave each species a local excess chemical potential mu_i (hard
spheres + MSA, :mod:`.csc`) in the 1-D pore, under local neutrality, and
found RyR1's filter binding Ca2+ as Gillespie 2008's does while the pore's
P_Ca:P_K stayed below 1: the uncharged gate, in series, capped it, and a
1-D local-neutrality closure carries no field into an uncharged slice.
Here the same fluid lives on the voxelised lumen, so the wall field of
Rounds 7.11–7.12 can reach the gate:

* the wall's groups become fluid species on the voxels
  (:func:`wall_fluid`): two half-charged oxygens per acid and one sphere
  per base, placed exactly as the fixed charge is (``slice``: the 1-D
  model's charge per length over each plane's lumen; ``local``: each
  group's 3-D Gaussian about its own centre), so their charge *is* the
  fixed-charge map (tested);
* ``local`` closure (:func:`local_csc`): every voxel in equilibrium with
  the bath under local neutrality, :func:`.csc.partition` voxel by voxel —
  the 1-D closure, voxel by voxel;
* ``pb`` closure (:func:`pb_csc`): Poisson's equation on the lumen, the
  mean field reaching past the charges, with each species' excess over the
  bath taken from that locally neutral fluid and held in its Boltzmann
  factor. The reference is neutral on purpose: the MSA is a theory of a
  neutral mixture, and at a held potential with no neutrality its
  screening term (−l_B Γ z², Γ growing as √ρ) outgrows ln c, so a voxel's
  cations run away (tried first; the fixed-potential solve found no root
  below u ≈ −5 kT/e). Gillespie's density functional evaluates its
  screening on a reference fluid for the same reason. Outside the wall
  groups' reach the reference is the bath, so the gate feels the mean
  field alone, which is the question.

Each species' equilibrium enters the conduction as an energy per voxel,
``E_i = −ln(c_i / c_i,bath)`` = z_i u + mu_i − mu_i,bath (kT), exactly
where Round 7.15's image cost entered: :class:`SpeciesField`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from ..parameters import PARAMETERS as _P
from ..structure.pore_volume import PoreVolume
from . import csc
from .charge3d import _EXP_CLIP, local_density, poisson_boltzmann, slice_density
from .selectivity import thermal_voltage

__all__ = ["PLACEMENTS", "wall_fluid", "fluid_for", "SpeciesField",
           "local_csc", "pb_csc"]

#: How the wall's groups are placed, as :mod:`.charge3d` places the charge.
PLACEMENTS = ("slice", "local")

_ACIDS, _BASES = ("ASP", "GLU"), ("LYS", "ARG")


def _placed(vol: PoreVolume, groups, positions, placement, charges):
    """mol/m³ of ``charges`` (one per group) on the lumen voxels, flat."""
    if placement == "slice":
        grid = slice_density(vol, [replace(g, charge=float(q))
                                   for g, q in zip(groups, charges)])
    else:
        grid, _ = local_density(vol, positions, np.asarray(charges, float))
    return grid[vol.mask]


def wall_fluid(vol: PoreVolume, groups, positions=None,
               placement: str = "local") -> csc.Structural:
    """The lining groups as fluid species on ``vol``'s lumen voxels (flat,
    C order): acid oxygens (two per acid) and bases, each placed as
    ``placement`` places the charge, so ``charge()`` is the fixed map."""
    if placement not in PLACEMENTS:
        raise ValueError(f"placement must be one of {PLACEMENTS}, not {placement!r}")
    m = int(vol.mask.sum())
    rows, names = [], []
    pos = np.zeros((0, 3)) if positions is None else np.asarray(positions, float)
    for name, kinds, per_group, key in (
            ("O", _ACIDS, 2.0, "csc.oxygen_diameter"),
            ("N", _BASES, 1.0, "csc.base_diameter")):
        pick = [i for i, g in enumerate(groups) if g.res_name in kinds]
        if not pick:
            continue
        mine = [groups[i] for i in pick]
        where = pos[pick] if len(pos) else None
        count = per_group * _placed(vol, mine, where, placement,
                                    np.ones(len(mine)))
        charge = _placed(vol, mine, where, placement,
                         [g.charge for g in mine])
        valence = np.divide(charge, count, out=np.zeros_like(count),
                            where=count > 0)
        rows.append((count, valence, _P.value(key)))
        names.append(name)
    if not rows:
        return csc.Structural.empty(m)
    return csc.Structural(np.array([r[0] for r in rows]),
                          np.array([r[1] for r in rows]),
                          np.array([r[2] for r in rows]), tuple(names))


def fluid_for(species, wall: csc.Structural, **kw):
    """(Fluid, activities (n,) mol/m³, water activity, bath (n,) mol/m³)."""
    fl = csc.Fluid.of(species, wall, **kw)
    bath = np.array([s.concentration * 1000.0 for s in species])
    a, aw = csc.bath_activity(bath, fl)
    return fl, a, aw, bath


@dataclass
class SpeciesField:
    """Each species' equilibrium on the grid, as an energy (kT) relative
    to its bath: c_i = c_i,bath e^{−E_i}. What a linear-response
    conductance needs, per species rather than per valence."""

    closure: str
    energies: dict[str, np.ndarray]     # name -> grid (0 off the mask)
    potential: np.ndarray               # kT/e on the grid
    mask: np.ndarray
    excess: dict[str, np.ndarray] = field(default_factory=dict)  # mu - mu_bath, kT
    converged: bool = True
    iterations: int = 0

    def energy(self, name: str) -> np.ndarray:
        return np.clip(self.energies[name], -_EXP_CLIP, _EXP_CLIP)


def _energies(vol, species, c, bath):
    out = {}
    for i, s in enumerate(species):
        g = np.zeros(vol.mask.shape)
        g[vol.mask] = -np.log(np.maximum(c[i], 1e-300) / bath[i])
        out[s.name] = g
    return out


def _to_grid(vol, flat):
    g = np.zeros(vol.mask.shape)
    g[vol.mask] = flat
    return g


def _walled(wall: csc.Structural, fixed) -> np.ndarray:
    """Voxels a wall group or a fixed charge reaches; elsewhere a neutral
    voxel with no wall is the bath, and is not solved."""
    return (wall.density > 0).any(axis=0) | (np.asarray(fixed) != 0)


def local_csc(vol: PoreVolume, fixed: np.ndarray, species,
              wall: csc.Structural, **fluid_kw) -> SpeciesField:
    """Every lumen voxel in equilibrium with the bath under local
    neutrality with the fluid's excess (:func:`.csc.partition`)."""
    fl, a, aw, bath = fluid_for(species, wall, **fluid_kw)
    x = fixed[vol.mask]
    on = np.flatnonzero(_walled(wall, x))
    alone = csc.packing(np.zeros((len(species), on.size)), np.zeros(on.size), fl,
                        wall.take(on))
    if np.any(alone >= csc._PACKING_CAP):
        raise ValueError(
            f"the wall's groups alone pack {alone.max():.2f} in "
            f"{int(np.sum(alone >= csc._PACKING_CAP))} voxels (no fluid fits): "
            "pore_charge.smoothing is too narrow for csc.structural_volume "
            f"{fl.wall_volume:g}")
    # Identical voxels (the slice placement fills a plane with one state)
    # are solved once.
    state = np.vstack([x[on], wall.density[:, on], wall.valence[:, on]]).T
    uniq, first, back = np.unique(state, axis=0, return_index=True,
                                  return_inverse=True)
    back = back.ravel()
    m = len(uniq)
    th = thermal_voltage()
    c = np.repeat(bath[:, None], x.size, axis=1)     # neutral, wall-free: bath
    psi = np.zeros(x.size)
    p = csc.partition(np.repeat(a[:, None], m, axis=1), x[on[first]], fl, th, aw,
                      wall=wall.take(on[first]),
                      initial=(np.repeat(bath[:, None], m, axis=1),
                               np.full(m, fl.water)))
    c[:, on], psi[on] = p.c[:, back], p.psi[back] / th
    mu_bath = np.log(a / bath)
    excess = {}
    for i, s in enumerate(species):
        e = np.zeros(x.size)
        e[on] = p.mu[i, back] - mu_bath[i]
        excess[s.name] = _to_grid(vol, e)
    return SpeciesField("local", _energies(vol, species, c, bath),
                        _to_grid(vol, psi), vol.mask, excess, p.converged,
                        p.iterations)


def pb_csc(vol: PoreVolume, fixed: np.ndarray, species, wall: csc.Structural,
           permittivity: float | None = None, reference: SpeciesField | None = None,
           **fluid_kw) -> SpeciesField:
    """Poisson's equation on the lumen (no field into protein, u = 0 on the
    baths, as the ``pb`` closure) with each species' excess over the bath,
    mu_i − mu_i,bath, taken from the locally neutral fluid (``reference``,
    :func:`local_csc` unless given) and held fixed in its Boltzmann factor.
    In the Donnan limit it is that closure again (tested)."""
    reference = reference or local_csc(vol, fixed, species, wall, **fluid_kw)
    offset = np.stack([reference.excess[s.name] for s in species])
    u, ok, used = poisson_boltzmann(vol, fixed, species, permittivity=permittivity,
                                    initial=reference.potential, offset=offset)
    energies = {s.name: np.where(vol.mask, s.valence * u + reference.excess[s.name], 0.0)
                for s in species}
    return SpeciesField("pb", energies, u, vol.mask, reference.excess,
                        bool(ok and reference.converged), used)
