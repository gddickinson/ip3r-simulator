"""Ion permeation: 1-D drift-diffusion over the measured pore profile.

Ported from the PIEZO1 simulator (``piezo1/physics/permeation.py``), whose
docstring records how the method was arrived at; the physics is unchanged and
the equations are in ``docs/SCIENCE.md``. In short:

* steady Nernst-Planck for each species down a channel of varying
  cross-section (the accessible area is the free radius less the ion's own),
  Scharfetter-Gummel discretised so drift-dominated slices stay stable;
* the potential closed in the **electroneutral limit**, because full Poisson
  coupling does not converge when the screening length exceeds the pore
  radius (PIEZO1 measured the divergence). A neutral pore between identical
  baths is closed by ohmic current continuity; a charged one by local
  electroneutrality against the wall charge (a local Donnan partition);
* Hall's access resistance in series at both mouths;
* :func:`series_conductance`, the closed-form resistor sum that a neutral
  pore must reduce to, derived without the solver, as its independent check.

What IP3R adds is only the input: the profile is this project's ``r_free``
(protein heavy atoms less their van der Waals radii), and the fixed charge
comes from the deposit's own side-chain atoms
(:mod:`ip3r.physics.pore_charge`). Every transport constant is registered
under ``permeation.*``; the in-pore diffusivity and the ion radius are the
two nobody has measured, and the answer is swept over them, not tuned.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from ._pnp_kernels import (F_FARADAY, R_GAS, _charge_diagnostics,
                           _donnan_potential, _face_conductance,
                           _nernst_planck, _neutrality_step, _ohmic_potential)

__all__ = ["IonSpecies", "PermeationResult", "solve_pnp", "series_conductance",
           "access_resistance", "potassium_species", "blocking_mechanisms",
           "debye_length", "bulk_conductivity", "F_FARADAY", "R_GAS"]


@dataclass(frozen=True)
class IonSpecies:
    """One permeant species. ``concentration`` is the bath at ``z[0]`` (M);
    ``concentration_right`` the bath at ``z[-1]`` (None = symmetric)."""

    name: str
    valence: int
    diffusivity: float                 # m^2/s, already scaled for the pore
    radius: float                      # A, crystal
    concentration: float               # M
    concentration_right: float | None = None

    @property
    def right(self) -> float:
        return (self.concentration if self.concentration_right is None
                else self.concentration_right)

    @property
    def symmetric(self) -> bool:
        return self.right == self.concentration


def potassium_species(bath: float | None = None,
                      diffusion_scale: float | None = None,
                      ion_radius: float | None = None) -> list[IonSpecies]:
    """Symmetric KCl: the condition every IP3R unitary conductance is quoted in.

    ``diffusion_scale`` and ``ion_radius`` override the registered in-pore
    diffusivity fraction and K+ radius, for the sensitivity sweep.
    """
    bath = _P.value("permeation.bath_concentration") if bath is None else bath
    scale = (_P.value("permeation.diffusion_scale") if diffusion_scale is None
             else diffusion_scale)
    k_radius = (_P.value("permeation.radius_potassium") if ion_radius is None
                else ion_radius)
    return [
        IonSpecies("K+", +1, _P.value("permeation.diffusion_potassium") * scale,
                   k_radius, bath),
        IonSpecies("Cl-", -1, _P.value("permeation.diffusion_chloride") * scale,
                   _P.value("permeation.radius_chloride"), bath),
    ]


@dataclass
class PermeationResult:
    """A current, and everything needed to disbelieve it."""

    current: float                     # A, at `voltage`, series-corrected
    conductance: float                 # S
    voltage: float                     # V
    z: np.ndarray                      # m
    radius: np.ndarray                 # m
    potential: np.ndarray              # V
    concentrations: dict = field(default_factory=dict)   # mol/m^3
    fluxes: dict = field(default_factory=dict)           # mol/s
    blocked_by: str | None = None
    access_ohm: float = 0.0
    pore_ohm: float = 0.0
    pore_current: float = 0.0          # A, signed, before the access correction
    converged: bool = True
    iterations: int = 0
    meta: dict = field(default_factory=dict)

    @property
    def conductance_pS(self) -> float:
        return self.conductance * 1e12

    @property
    def is_conducting(self) -> bool:
        return self.blocked_by is None and self.conductance > 0.0

    def summary(self) -> str:
        if self.blocked_by:
            return f"no current: {self.blocked_by}"
        return (f"{self.conductance_pS:.0f} pS ({self.current * 1e12:.2f} pA at "
                f"{self.voltage * 1e3:.0f} mV); pore {self.pore_ohm / 1e9:.2f} "
                f"GOhm + access {self.access_ohm / 1e9:.2f} GOhm")


def access_resistance(mouth_radius: float, conductivity: float) -> float:
    """Hall's spreading resistance of one circular mouth, ``1/(4 sigma a)``, ohm."""
    if mouth_radius <= 0 or conductivity <= 0:
        return np.inf
    return 1.0 / (4.0 * conductivity * mouth_radius)


def _species_conductivity(s: IonSpecies, temperature: float,
                          concentration: float | None = None) -> float:
    """Nernst-Einstein conductivity of one species, S/m (c in M)."""
    c = s.concentration if concentration is None else concentration
    return (s.valence ** 2 * F_FARADAY ** 2 * s.diffusivity * c * 1000.0
            / (R_GAS * temperature))


def bulk_conductivity(species: list[IonSpecies], temperature: float) -> float:
    """Ohmic conductivity of the bath, S/m."""
    return sum(_species_conductivity(s, temperature) for s in species)


def _accessible_area(radius_m: np.ndarray, ion_radius_A: float) -> np.ndarray:
    """Cross-section an ion's centre can occupy (hard-sphere exclusion), m^2."""
    usable = np.maximum(radius_m - ion_radius_A * 1e-10, 0.0)
    return np.pi * usable ** 2


def series_conductance(z_A: np.ndarray, radius_A: np.ndarray,
                       species: list[IonSpecies] | None = None,
                       temperature: float | None = None) -> dict:
    """``R = integral dz / (sigma A(z)) + 2 R_access``, without the solver.

    What a neutral pore between identical baths must reduce to. Derived here
    with no Scharfetter-Gummel, no Gummel loop and no Donnan step, so a
    disagreement with :func:`solve_pnp` means one of them is wrong.
    """
    species = species or potassium_species()
    temperature = temperature or _P.value("permeation.temperature")
    z = np.asarray(z_A, dtype=float) * 1e-10
    radius = np.asarray(radius_A, dtype=float) * 1e-10
    per_length = np.zeros_like(z)
    for s in species:
        per_length += (_species_conductivity(s, temperature)
                       * _accessible_area(radius, s.radius))
    if np.any(per_length <= 0):
        return {"conductance": 0.0, "pore_ohm": np.inf, "access_ohm": 0.0,
                "blocked": "a slice is too narrow for any ion to enter"}
    pore_ohm = float(np.trapezoid(1.0 / per_length, z))
    sigma = bulk_conductivity(species, temperature)
    access_ohm = 2.0 * access_resistance(float(max(radius[0], radius[-1])), sigma)
    return {"conductance": 1.0 / (pore_ohm + access_ohm), "pore_ohm": pore_ohm,
            "access_ohm": access_ohm, "conductivity": sigma}


def debye_length(species: list[IonSpecies], temperature: float,
                 permittivity_relative: float) -> float:
    """Screening length of the bath, m — reported because it decides whether
    the electroneutral limit (and a continuum model at all) applies."""
    ionic_strength = 0.5 * sum(s.valence ** 2 * s.concentration * 1000.0
                               for s in species)
    if ionic_strength <= 0:
        return np.inf
    return float(np.sqrt(permittivity_relative * 8.8541878128e-12 * R_GAS
                         * temperature / (2.0 * F_FARADAY ** 2 * ionic_strength)))


def blocking_mechanisms(radius_m: np.ndarray,
                        species: list[IonSpecies]) -> list[str]:
    """Every reason no current flows. IP3R has no wetting model, so this is
    the steric test alone: the narrowest slice against the smallest cation."""
    smallest = min(s.radius for s in species if s.valence > 0)
    if float(radius_m.min()) <= smallest * 1e-10:
        return [f"sterically occluded: the narrowest slice is "
                f"{radius_m.min() * 1e10:.2f} A, below the {smallest:.2f} A "
                f"radius of the smallest permeant cation"]
    return []


def solve_pnp(z_A: np.ndarray, radius_A: np.ndarray,
              voltage: float | None = None,
              species: list[IonSpecies] | None = None,
              fixed_charge: np.ndarray | None = None,
              max_iterations: int = 400, tol: float = 1e-10,
              relaxation: float = 0.4) -> PermeationResult:
    """Solve 1-D drift-diffusion over a profile (``z_A``, ``radius_A`` in A).

    ``fixed_charge`` is the wall charge per slice as a signed molar-equivalent
    density (mol/m^3, rho_fixed / F), from
    :func:`ip3r.physics.pore_charge.map_charge`. With none, and identical
    baths, the arithmetic is the ohmic closure and must agree with
    :func:`series_conductance`; with charge the closure is local
    electroneutrality (the ohmic operator has no term the charge could enter
    through).
    """
    voltage = _P.value("permeation.test_voltage") if voltage is None else voltage
    species = species or potassium_species()
    temperature = _P.value("permeation.temperature")
    thermal = R_GAS * temperature / F_FARADAY

    z = np.asarray(z_A, dtype=float) * 1e-10
    radius = np.asarray(radius_A, dtype=float) * 1e-10
    order = np.argsort(z)
    z, radius = z[order], radius[order]

    reasons = blocking_mechanisms(radius, species)
    if reasons:
        return PermeationResult(
            current=0.0, conductance=0.0, voltage=voltage, z=z, radius=radius,
            potential=np.linspace(0.0, voltage, len(z)),
            blocked_by=" AND ".join(reasons),
            meta={"n_slices": len(z), "min_radius_A": float(radius.min()) * 1e10})

    areas = {s.name: _accessible_area(radius, s.radius) for s in species}
    fixed = (None if fixed_charge is None
             else np.asarray(fixed_charge, dtype=float)[order])
    charged = fixed is not None and bool(np.any(fixed != 0.0))
    symmetric = all(s.symmetric for s in species)

    if charged or not symmetric:
        reference = np.array([[0.5 * (s.concentration + s.right) * 1000.0] * len(z)
                              for s in species])
        psi = _donnan_potential([s.valence for s in species], reference,
                                np.zeros_like(z) if fixed is None else fixed,
                                thermal)
    else:
        psi = None

    def _bath(value, valence, end):
        if psi is None:
            return value * 1000.0
        return value * 1000.0 * float(np.exp(-valence * psi[end] / thermal))

    left = {s.name: _bath(s.concentration, s.valence, 0) for s in species}
    right = {s.name: _bath(s.right, s.valence, -1) for s in species}

    applied = np.linspace(0.0, voltage, len(z))
    potential = applied if psi is None else applied + psi
    valences = [s.valence for s in species]
    zero_charge = np.zeros_like(z)
    threshold = (tol * max(abs(voltage), 1e-12) if psi is None
                 else tol * max(abs(voltage), thermal))
    concentrations, fluxes = {}, {}
    converged, used = False, 0
    # A species excluded from any slice carries exactly zero flux (slices are
    # in series), and its zero-area rows would make the matrix singular.
    excluded = [s.name for s in species if float(areas[s.name].min()) <= 0.0]

    for used in range(1, max_iterations + 1):  # noqa: B007 - read after the loop
        for s in species:
            if s.name in excluded:
                concentrations[s.name] = np.full_like(z, left[s.name])
                fluxes[s.name] = 0.0
                continue
            concentrations[s.name], fluxes[s.name] = _nernst_planck(
                z, areas[s.name], potential, s.valence, s.diffusivity, thermal,
                left[s.name], right[s.name])
        if psi is None:
            face = _face_conductance(z, areas, concentrations, species,
                                     temperature)
            updated = _ohmic_potential(z, face, 0.0, voltage)
            change = float(np.max(np.abs(updated - potential)))
            potential = (1.0 - relaxation) * potential + relaxation * updated
        else:
            step = _neutrality_step(
                valences, [concentrations[s.name] for s in species],
                zero_charge if fixed is None else fixed, thermal)
            change = float(np.max(np.abs(step)))
            potential = potential + relaxation * step
        if change < threshold:
            converged = True
            break

    current = sum(s.valence * F_FARADAY * fluxes[s.name] for s in species)
    sigma = bulk_conductivity(species, temperature)
    lam = debye_length(species, temperature,
                       _P.value("permeation.permittivity_pore"))
    access_ohm = 2.0 * access_resistance(float(max(radius[0], radius[-1])), sigma)
    # A chord conductance at zero drive is 0/0 (PIEZO1 once reported 1/R_access
    # from a 1e-20 A solver residue at 0 V), so it is gated explicitly.
    pore_ohm = (abs(voltage / current)
                if current != 0 and voltage != 0 else np.inf)
    total_ohm = pore_ohm + access_ohm
    meta = {"n_slices": len(z), "conductivity_S_per_m": sigma,
            "species": [s.name for s in species],
            "debye_length_A": lam * 1e10,
            "min_radius_A": float(radius.min()) * 1e10,
            "double_layers_overlap": bool(lam > float(radius.min())),
            "fixed_charge": charged, "symmetric_baths": symmetric,
            "excluded": excluded}
    if charged or not symmetric:
        meta.update(_charge_diagnostics(concentrations, species, fixed, psi))
    return PermeationResult(
        current=voltage / total_ohm if np.isfinite(total_ohm) else 0.0,
        conductance=1.0 / total_ohm if np.isfinite(total_ohm) else 0.0,
        voltage=voltage, z=z, radius=radius, potential=potential,
        concentrations=concentrations, fluxes=fluxes, access_ohm=access_ohm,
        pore_ohm=pore_ohm, pore_current=float(current), converged=converged,
        iterations=used, meta=meta)
