"""The Gummel loop shared by the closures that give each species its own
potential offset: ``radial`` (Round 4) and ``csc`` (Round 7.17).

Both close a charged slice with a local equilibrium that is not one
Donnan potential. Each species i carries an offset ``w_i`` (V), and the
Nernst-Planck solve runs in ``applied + w_i``. One pass:

1. Nernst-Planck for every species in its own potential;
2. the reservoir each slice is in equilibrium with,
   ``a_i = c_i exp(z_i w_i / phi_T)``;
3. the closure's local equilibrium on that reservoir, giving new offsets;
4. ``w`` relaxed towards them; the ends stay pinned by their own baths.

The radial closure's reservoir is a concentration and its offset the
cross-section partition (:mod:`.radial_pb`). The csc closure's reservoir
is an activity and its offset psi + mu_i kT/(z_i e), with mu_i the excess
chemical potential (:mod:`.csc`). The loop was the radial one; it moved
here unchanged when the second closure needed it.
"""

from __future__ import annotations

import numpy as np

from ..parameters import PARAMETERS as _P
from ._pnp_kernels import _nernst_planck

__all__ = ["solve_offsets", "RadialClosure", "CscClosure"]


class RadialClosure:
    """Cylindrical Poisson-Boltzmann across each slice."""

    name = "radial"

    def __init__(self, species, radius, fixed, thermal):
        from .radial_pb import radial_partition
        self._solve = radial_partition
        self.valences = np.array([s.valence for s in species], dtype=float)
        self.wall = np.zeros_like(radius) if fixed is None else fixed
        # The slice the charge sits on: the same floor map_charge divides by.
        self.slice_r = np.maximum(
            radius, _P.value("permeation.radius_potassium") * 1e-10)
        self.thermal = thermal
        self.species = species

    def start(self):
        """(offsets (n_sp, n), left reservoir, right reservoir)."""
        n = len(self.wall)
        reference = np.array([[0.5 * (s.concentration + s.right) * 1000.0] * n
                              for s in self.species])
        part = self._solve(self.valences, reference, self.wall, self.slice_r,
                           self.thermal)
        offset = part.offset.copy()
        self.warm, self.ok, self.part = part.psi, part.converged, part
        for end, bath_of in ((0, lambda s: s.concentration),
                             (-1, lambda s: s.right)):
            own = np.array([[bath_of(s) * 1000.0] for s in self.species])
            offset[:, end] = self._solve(self.valences, own, self.wall[[end]],
                                         self.slice_r[[end]],
                                         self.thermal).offset[:, 0]
        left = [s.concentration * 1000.0 for s in self.species]
        right = [s.right * 1000.0 for s in self.species]
        return offset, left, right

    def update(self, reservoir):
        part = self._solve(self.valences, reservoir, self.wall, self.slice_r,
                           self.thermal, initial=self.warm)
        self.warm, self.ok, self.part = part.psi, self.ok and part.converged, part
        return part.offset

    def extra(self, offset) -> dict:
        return {"species_offset_mV": {s.name: (offset[i] * 1e3).tolist()
                                      for i, s in enumerate(self.species)},
                "radial_axis_mV": [float(self.part.axis.min() * 1e3),
                                   float(self.part.axis.max() * 1e3)],
                "radial_wall_mV": [float(self.part.wall.min() * 1e3),
                                   float(self.part.wall.max() * 1e3)],
                "radial_converged": self.ok}


class CscClosure:
    """Local electroneutrality with the hard-sphere + MSA excess."""

    name = "csc"

    def __init__(self, species, fixed, thermal, wall=None, **fluid_kw):
        from . import csc
        self.csc = csc
        n = len(fixed)
        self.fixed = np.asarray(fixed, float)
        self.wall = (csc.Structural.from_fixed(self.fixed) if wall is None
                     else wall)
        self.fluid = csc.Fluid.of(species, self.wall, **fluid_kw)
        self.species, self.thermal, self.n = species, thermal, n
        self.valences = self.fluid.valence

    def _equilibrium(self, activity, fixed, wall, warm=None):
        return self.csc.partition(activity, fixed, self.fluid, self.thermal,
                                  self.water, wall=wall, initial=warm)

    def _guess(self, bath, m):
        return (np.repeat(np.asarray(bath, float)[:, None], m, axis=1),
                np.full(m, self.fluid.water))

    def start(self):
        mean = np.array([0.5 * (s.concentration + s.right) * 1000.0
                         for s in self.species])
        act, self.water = self.csc.bath_activity(mean, self.fluid)
        part = self._equilibrium(np.repeat(act[:, None], self.n, axis=1),
                                 self.fixed, self.wall,
                                 self._guess(mean, self.n))
        offset = part.offset(self.valences, self.thermal)
        self.warm, self.ok = part, part.converged
        ends = []
        for end, bath_of in ((0, lambda s: s.concentration),
                             (-1, lambda s: s.right)):
            own = np.array([bath_of(s) * 1000.0 for s in self.species])
            a, _ = self.csc.bath_activity(own, self.fluid)
            p = self._equilibrium(a[:, None], self.fixed[[end]],
                                  self.wall.take([end]), self._guess(own, 1))
            offset[:, end] = p.offset(self.valences, self.thermal)[:, 0]
            self.ok = self.ok and p.converged
            ends.append(a.tolist())
        return offset, ends[0], ends[1]

    def update(self, reservoir):
        part = self._equilibrium(reservoir, self.fixed, self.wall, self.warm)
        self.warm, self.ok = part, self.ok and part.converged
        return part.offset(self.valences, self.thermal)

    def extra(self, offset) -> dict:
        p = self.warm
        return {"species_offset_mV": {s.name: (offset[i] * 1e3).tolist()
                                      for i, s in enumerate(self.species)},
                "excess_kT": {s.name: p.mu[i].tolist()
                              for i, s in enumerate(self.species)},
                "water_M": (p.water / 1000.0).tolist(),
                "csc_converged": self.ok,
                "csc_terms": {"hs": self.fluid.hs, "msa": self.fluid.msa,
                              "wall_volume": self.fluid.wall_volume,
                              "bjerrum_A": self.fluid.lb}}


def solve_offsets(closure, z, radius, voltage, species, areas, fixed, thermal,
                  temperature, max_iterations, tol, relaxation, charged,
                  symmetric):
    """The Gummel loop with a per-species offset closure (see the module)."""
    from .permeation import _assemble
    offset, left_res, right_res = closure.start()
    left = {s.name: left_res[i]
            * float(np.exp(-s.valence * offset[i, 0] / thermal))
            for i, s in enumerate(species)}
    right = {s.name: right_res[i]
             * float(np.exp(-s.valence * offset[i, -1] / thermal))
             for i, s in enumerate(species)}
    applied = np.linspace(0.0, voltage, len(z))
    excluded = [s.name for s in species if float(areas[s.name].min()) <= 0.0]
    threshold = tol * max(abs(voltage), thermal)
    concentrations, fluxes = {}, {}
    converged, used = False, 0
    for used in range(1, max_iterations + 1):  # noqa: B007 - read after the loop
        for i, s in enumerate(species):
            if s.name in excluded:
                concentrations[s.name] = np.full_like(z, left[s.name])
                fluxes[s.name] = 0.0
                continue
            concentrations[s.name], fluxes[s.name] = _nernst_planck(
                z, areas[s.name], applied + offset[i], s.valence,
                s.diffusivity, thermal, left[s.name], right[s.name])
        reservoir = np.array([
            np.maximum(concentrations[s.name], 0.0)
            * np.exp(np.clip(s.valence * offset[i] / thermal, -40.0, 40.0))
            for i, s in enumerate(species)])
        step = closure.update(reservoir) - offset
        step[:, [0, -1]] = 0.0
        offset = offset + relaxation * step
        change = float(np.max(np.abs(step)))
        if change < threshold:
            converged = True
            break
    extra = closure.extra(offset)
    ok = extra.get("radial_converged", extra.get("csc_converged", True))
    return _assemble(z, radius, voltage, species, temperature,
                     applied + offset[0], concentrations, fluxes,
                     converged and ok, used, charged, symmetric,
                     excluded, fixed, offset[0], extra, closure=closure.name)
