"""Steady Poisson–Nernst–Planck on the voxelised lumen (Round 7.23).

Round 7.19 read permeability from linear response: identical baths, each
species' equilibrium held fixed, one Boltzmann-weighted Laplace solve per
species. The measured ratios (Xu 2006, Vais 2010) are bi-ionic reversal
potentials instead. There the baths differ and the net current is zero,
but each species still carries a current, so neither the concentrations
nor the potential are at equilibrium. This module solves that state.

**Nernst–Planck in Slotboom form.** With a species' energy
``E_i = z_i ψ + μ_i`` (kT; ψ the potential in kT/e, μ_i a fixed excess
per voxel), the flux is ``J_i = −D_i e^{−E_i} ∇n_i`` with
``n_i = c_i e^{E_i}``. Steady state is ``∇·(e^{−E_i}∇n_i) = 0``, which is
exactly Round 7.11's weighted Laplace problem
(:func:`.ohmic3d.geometric_conductance`, Scharfetter–Gummel faces). With
n_i held at each bath's value (luminal ``c_lum``, cytosolic
``c_cyt e^{z_i v}``), ``n_i = n_lum + (n_cyt − n_lum) φ_i`` and the
particle current is ``D_i g_i (n_cyt − n_lum)``, where φ_i and g_i are
that solve's potential and conductance. Each species moves in its own
volume (its own radius); elsewhere in the electrostatic volume it takes
n from its nearest voxel.

**Poisson.** ψ = v φ_0 + u, where φ_0 is the uncharged lumen's Laplace
potential (1 at the cytosolic bath, 0 at the luminal) and u = 0 on both
baths. Laplace's operator annihilates v φ_0, so u obeys Round 7.11's
Poisson–Boltzmann with each species' charge ``n_i e^{−z_i(vφ_0 + u) − μ_i}``,
that is :func:`.charge3d.poisson_boltzmann` with the applied field,
the excess and ln n_i all in its ``offset``. Impermeant species (Vais's
NMDG⁺) are Boltzmann-distributed in their own region and carry no flux.

**Gummel.** Alternate the two with the other held until u stops moving.
At v = 0 between identical baths every n_i is uniform and the solve is
Round 7.11's ``pb`` equilibrium exactly (tested). A small v then gives
Round 7.19's linear-response conductances (tested).

**Round 7.27** adds two kinds of energy. ``hidden`` energies act on
transport but not on Poisson: a Ca²⁺ site whose −2e is fixed per bound
Ca²⁺ (compensated), so what it gathers carries no net charge. A
``coupling`` returns energies that depend on the state itself (a saturable
site's, K⁺ blocked by the Ca²⁺ it holds). They are recomputed from each
Gummel iterate's n and ψ, damped, and must settle with u.

Frame and signs follow :mod:`.selectivity`: the lumen is the ``bottom``
bath, v = V_cyt − V_lumen, and current is positive from lumen to cytosol
for a cation.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy import ndimage

from ..parameters import PARAMETERS as _P
from ..structure.pore_volume import PoreVolume
from ._pnp_kernels import F_FARADAY
from .charge3d import _EXP_CLIP, poisson_boltzmann
from .ohmic3d import geometric_conductance
from .selectivity import thermal_voltage

__all__ = ["Domain", "domain", "Impermeant", "Steady", "steady_state"]

#: Floor on n (mol/m³) before its log enters the Poisson offset; a species
#: absent from both ends of a region is absent there, not −∞.
_N_FLOOR = 1e-12


@dataclass
class Domain:
    """The volumes one deposit's solves share: the electrostatic one (the
    smallest ion's), each species' own, and the applied-field profile."""

    elec: PoreVolume
    vols: dict[str, PoreVolume]
    phi0: np.ndarray                    # grid, 1 cytosol -> 0 lumen on elec
    nearest: dict[str, tuple] = field(repr=False, default_factory=dict)


def domain(elec: PoreVolume, vols: dict[str, PoreVolume]) -> Domain:
    lap = geometric_conductance(elec)
    phi0 = np.zeros(elec.mask.shape)
    phi0[elec.mask] = lap.potential
    nearest = {}
    for name, vol in vols.items():
        if not np.any(elec.mask & ~vol.mask):
            continue
        # every voxel -> the nearest voxel of this species' own volume
        _, idx = ndimage.distance_transform_edt(~vol.mask, return_indices=True)
        nearest[name] = tuple(idx)
    return Domain(elec, vols, phi0, nearest)


@dataclass
class Impermeant:
    """A species confined to ``region`` (a grid mask): Boltzmann there with
    ``concentration`` (M) at ψ = 0, absent elsewhere, no flux."""

    name: str
    valence: int
    concentration: float
    region: np.ndarray


@dataclass
class Steady:
    v: float                            # V, cytosol - lumen
    u: np.ndarray                       # kT/e on the grid (psi = v phi0 + u)
    n: dict[str, np.ndarray]            # Slotboom n per species, mol/m³ (grid)
    flux: dict[str, float]              # mol/s, cytosol -> lumen
    phi: dict[str, np.ndarray]          # each species' solve, per own voxel
    converged: bool
    iterations: int
    # Round 7.27: the coupling's energies at the end, (visible, hidden)
    coupled: tuple = field(default=({}, {}), repr=False)

    def energy(self, name: str, hidden: bool = True):
        """The coupling's energy on species ``name`` (kT grid, or 0):
        what Poisson sees, plus the hidden part when ``hidden``."""
        vis, hid = self.coupled
        e = vis.get(name, 0.0)
        return e + hid.get(name, 0.0) if hidden else e

    @property
    def currents(self) -> dict[str, float]:
        """A per species, positive from lumen to cytosol for a cation."""
        return {k: -F_FARADAY * self._z[k] * j for k, j in self.flux.items()}

    @property
    def current(self) -> float:
        return float(sum(self.currents.values()))

    _z: dict = field(default_factory=dict, repr=False)


def _ends(s, vhat: float) -> tuple[float, float]:
    """(n at the lumen, n at the cytosol), mol/m³."""
    return s.concentration * 1000.0, s.right * 1000.0 * np.exp(s.valence * vhat)


def _transport(dom: Domain, s, psi, excess, vhat, initial):
    """One species' steady Nernst–Planck in its own volume given ψ."""
    vol = dom.vols[s.name]
    e = s.valence * psi
    if excess is not None:
        e = e + excess
    e = np.clip(e, -_EXP_CLIP, _EXP_CLIP)
    lap = geometric_conductance(vol, energy=e, initial=initial)
    lo, hi = _ends(s, vhat)
    n = np.zeros(vol.mask.shape)
    n[vol.mask] = lo + (hi - lo) * lap.potential
    if s.name in dom.nearest:                       # elec voxels it cannot reach
        n = n[dom.nearest[s.name]]
    flux = s.diffusivity * lap.g * (hi - lo)
    return n, flux, lap.potential, lap.converged


def _merge(*dicts) -> dict:
    """Energies per species summed over several dicts (None skipped)."""
    out: dict = {}
    for d in dicts:
        for k, v in (d or {}).items():
            out[k] = out[k] + v if k in out else v
    return out


def _offsets(dom, species, n, excess, vhat, impermeant):
    """Poisson's offset per species: ψ's applied part, the excess and
    −ln(n / c_ref), so that c = c_ref e^{−z u − offset} = n e^{−E}."""
    mask = dom.elec.mask
    rows, refs = [], []
    for s in species:
        ref = max(s.concentration, s.right)
        o = s.valence * vhat * dom.phi0 - np.log(
            np.maximum(n[s.name], _N_FLOOR) / (ref * 1000.0))
        if excess and s.name in excess:
            o = o + excess[s.name]
        rows.append(np.where(mask, o, 0.0))
        refs.append(replace(s, concentration=ref, concentration_right=None))
    for m in impermeant:
        o = np.where(m.region, m.valence * vhat * dom.phi0, 2 * _EXP_CLIP)
        rows.append(np.where(mask, o, 0.0))
        refs.append(m)
    return np.stack(rows), refs


def steady_state(dom: Domain, species, v: float, fixed: np.ndarray,
                 excess: dict[str, np.ndarray] | None = None,
                 impermeant=(), initial: Steady | None = None,
                 permittivity: float | None = None,
                 hidden: dict[str, np.ndarray] | None = None,
                 coupling=None) -> Steady:
    """The steady state at applied voltage ``v`` (V) between the baths the
    species carry (``concentration`` = lumen, ``right`` = cytosol).
    ``fixed`` is the wall's charge (mol/m³ on the grid), ``excess`` each
    species' held excess chemical potential (kT on the grid, 0 = the
    reference bath's), ``initial`` a previous solution to start from.
    ``hidden``: energies transport sees and Poisson does not. ``coupling``:
    ``f(psi, n) -> (visible, hidden)``, energies recomputed each iteration,
    starting from ``coupling.start()``."""
    vhat = v / thermal_voltage()
    u = (np.zeros(dom.elec.mask.shape) if initial is None
         else initial.u.copy())
    phis = {} if initial is None else dict(initial.phi)
    vis, hid = ({}, {}) if coupling is None else (
        initial.coupled if initial is not None else coupling.start())
    tol = _P.value("pnp3d.gummel_tolerance")
    damping = _P.value("pnp3d.gummel_damping")
    ok, used = False, 0
    for used in range(1, int(_P.value("pnp3d.gummel_max_iterations")) + 1):  # noqa: B007
        psi = vhat * dom.phi0 + u
        seen = _merge(excess, vis)
        moved = _merge(seen, hidden, hid)
        n, flux, inner = {}, {}, True
        for s in species:
            n[s.name], flux[s.name], phis[s.name], c_ok = _transport(
                dom, s, psi, moved.get(s.name), vhat, phis.get(s.name))
            inner &= c_ok
        offset, refs = _offsets(dom, species, n, seen, vhat, impermeant)
        new, p_ok, _ = poisson_boltzmann(dom.elec, fixed, refs,
                                         permittivity=permittivity,
                                         initial=u, offset=offset)
        change = float(np.max(np.abs(new - u)))
        u = u + damping * (new - u)
        if coupling is not None:
            vis, hid, moved_c = _relax(coupling(vhat * dom.phi0 + u, n),
                                       (vis, hid))
            change = max(change, moved_c)
        if change < tol and inner and p_ok:
            ok = True
            break
    psi = vhat * dom.phi0 + u
    moved = _merge(excess, vis, hidden, hid)
    for s in species:                     # the currents in the final field
        n[s.name], flux[s.name], phis[s.name], _ = _transport(
            dom, s, psi, moved.get(s.name), vhat, phis.get(s.name))
    return Steady(v, u, n, flux, phis, ok, used, coupled=(vis, hid),
                  _z={s.name: s.valence for s in species})


def _relax(new: tuple, old: tuple) -> tuple:
    """The coupling's energies moved ``pnp3d.coupling_damping`` of the way
    to their new values; and the largest move (kT) on the lumen."""
    w = _P.value("pnp3d.coupling_damping")
    out, most = [], 0.0
    for fresh, prev in zip(new, old):
        d = {}
        for k, e in fresh.items():
            p = prev.get(k, 0.0)
            d[k] = p + w * (e - p)
            most = max(most, float(np.max(np.abs(d[k] - p))))
        out.append(d)
    return out[0], out[1], most
