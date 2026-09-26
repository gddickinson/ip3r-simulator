"""The charged pore's conductance in three dimensions (Round 7.11).

**Why no 3-D drift-diffusion solve is needed.** Between identical baths at
small voltage, a steady drift-diffusion pore is a perturbation of its
equilibrium: each species' flux is ``J = −(D c_eq / kT) ∇μ``, divergence
free, with the electrochemical potential μ fixed at the two baths. So each
species is one Laplace problem whose conductivity is its bulk σ times its
equilibrium partition ``c_eq / c_bath = e^{−z u}``, and the species add in
parallel. Neutrality, Poisson or Donnan, enters only through the equilibrium
``u``. :func:`linear_response_1d` is that statement in 1-D, and the tests
hold the existing 1-D solver (:func:`.permeation.solve_pnp`, Gummel loop
at a finite voltage) to it. In 3-D the equilibrium is
:mod:`.charge3d`'s and the Laplace solve :func:`.ohmic3d.geometric_conductance`
with an energy per voxel.

This is also why rings of alternating sign lower the 1-D conductance: each
species' resistance is ``∫ dz / (σ A e^{−z u})``, dominated by the zone where
it is the excluded co-ion, and in a series pore every carrier crosses one.
In 3-D a co-ion excluded from a corner may still pass along the axis.

Readings, in order of departure from the 1-D model: ``slice`` (its charge
per length over the real cross-section), ``local`` (each group at its own
centre), ``pb`` (screened by Poisson–Boltzmann). Each is set beside the same
deposit's neutral 3-D reading; the ratio is the charge's effect.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.channel import ChannelSummary, measure_channel
from ..structure.pore_volume import pore_volume
from ._pnp_kernels import _donnan_potential
from .charge3d import CLOSURES_3D, WallField, group_positions, wall_field
from .dielectric3d import DIELECTRIC, dielectric_field
from .ohmic3d import geometric_conductance
from .permeation import (_accessible_area, _species_conductivity,
                         access_resistance, bulk_conductivity,
                         potassium_species)
from .pore_charge import PoreCharge, pore_charge

__all__ = ["linear_response_1d", "Wall3D", "wall_3d", "mutant_panel_3d", "wall_scan", "SCAN_KEYS",
           "MutantRow3D"]


def linear_response_1d(z_A: np.ndarray, radius_A: np.ndarray, species,
                       fixed: np.ndarray | None = None) -> float:
    """Zero-voltage conductance of a 1-D charged pore, S: local Donnan
    equilibrium, each species a series resistor ``∫ dz / (σ A e^{−z u})``
    (the Scharfetter–Gummel face weight between slices), the species in
    parallel, Hall's access in series. No solver loop."""
    from .ohmic3d import bernoulli_weight
    temperature = _P.value("permeation.temperature")
    z = np.asarray(z_A, float) * 1e-10
    radius = np.asarray(radius_A, float) * 1e-10
    valences = [s.valence for s in species]
    bath = np.array([[s.concentration * 1000.0] * len(z) for s in species])
    wall = np.zeros_like(z) if fixed is None else np.asarray(fixed, float)
    u = _donnan_potential(valences, bath, wall, 1.0)
    g = 0.0
    for s in species:
        area = _accessible_area(radius, s.radius)
        if np.any(area <= 0):
            continue
        e = s.valence * u
        # Face between slices i, i+1: area and Boltzmann factor both vary.
        w = bernoulli_weight(e[:-1], e[1:])
        a_face = 0.5 * (area[:-1] + area[1:])
        r = np.sum(np.diff(z) / (a_face * w)) / _species_conductivity(s, temperature)
        g += 1.0 / r
    sigma = bulk_conductivity(species, temperature)
    access = 2.0 * access_resistance(float(max(radius[0], radius[-1])), sigma)
    return 1.0 / (1.0 / g + access) if g > 0 else 0.0


@dataclass
class Wall3D:
    """One deposit's 3-D conductance, neutral and under each closure."""

    name: str
    neutral: float                          # S
    charged: dict[str, float]               # closure -> S
    fields: dict[str, WallField] = field(default_factory=dict, repr=False)
    per_species: dict[str, dict[str, float]] = field(default_factory=dict)
    charge: PoreCharge | None = field(default=None, repr=False)
    spacing: float = 0.0
    converged: bool = True
    bath: float = 0.0                       # M, the cation's

    def ratio(self, closure: str) -> float:
        return self.charged[closure] / self.neutral if self.neutral else float("nan")

    def row(self) -> str:
        parts = "  ".join(f"{c} {self.charged[c] * 1e12:6.1f} (x{self.ratio(c):4.2f})"
                          for c in self.charged)
        return (f"{self.name:5s} neutral {self.neutral * 1e12:6.1f}  {parts} pS"
                f"{'' if self.converged else '  [n.c.]'}")

    def peaks(self) -> dict[str, float]:
        """Highest counter-ion (cation) concentration per closure, M."""
        return {c: f.peak(+1, self.bath) for c, f in self.fields.items()}


def _profile(st: Structure, summary: ChannelSummary):
    from .unitary import permeation_profile
    return permeation_profile(st, summary)


def wall_3d(st: Structure, summary: ChannelSummary | None = None,
            closures=CLOSURES_3D, neutralise: frozenset[int] = frozenset(),
            pair_bridges: bool = False, spacing: float | None = None,
            species=None, width: float | None = None,
            permittivity: float | None = None,
            keep_fields: bool = False, scope: str = "all",
            eps_protein: float | None = None) -> Wall3D:
    """A deposit's K+ conductance in 3-D, neutral and under each closure,
    in its family's bath. The electrostatics live on the smallest ion's
    volume; each species conducts on its own. A closure may also be
    ``dielectric`` (:mod:`.dielectric3d`, with ``scope`` and
    ``eps_protein``)."""
    from .unitary import bath_for
    summary = summary or measure_channel(st)
    paralog = summary.numbering.paralog if summary.numbering else None
    species = species or potassium_species(bath=bath_for(paralog))
    temperature = _P.value("permeation.temperature")
    h = _P.value("pore3d.spacing") if spacing is None else spacing
    charge = pore_charge(st, summary.frame, _profile(st, summary),
                         pair_bridges=pair_bridges, neutralise=neutralise)
    vols = {s.name: pore_volume(st, summary.frame, summary.span, s.radius,
                                spacing=h) for s in species}
    smallest = min(species, key=lambda s: s.radius)
    elec = vols[smallest.name]
    positions = group_positions(st, summary.frame, charge.groups)
    sigma = {s.name: _species_conductivity(s, temperature) for s in species}
    per: dict[str, dict[str, float]] = {"neutral": {}}
    ok = True
    neutral = 0.0
    for s in species:
        lap = geometric_conductance(vols[s.name])
        per["neutral"][s.name] = lap.g
        neutral += sigma[s.name] * lap.g
        ok &= lap.converged
    charged, fields = {}, {}
    for closure in closures:
        if closure == DIELECTRIC:
            wf = dielectric_field(st, summary.frame, elec, species, charge,
                                  scope=scope, neutralise=neutralise,
                                  eps_protein=eps_protein)
        else:
            wf = wall_field(elec, closure, species, charge, positions=positions,
                            width=width, permittivity=permittivity)
        ok &= wf.converged
        total, per[closure] = 0.0, {}
        for s in species:
            lap = geometric_conductance(vols[s.name], energy=wf.energy(s.valence))
            per[closure][s.name] = lap.g
            total += sigma[s.name] * lap.g
            ok &= lap.converged
        charged[closure] = total
        if keep_fields:
            fields[closure] = wf
        else:
            fields[closure] = WallField(closure, np.zeros(0), wf.potential,
                                        wf.mask, wf.placed, wf.unreached,
                                        wf.converged, wf.iterations)
    return Wall3D(st.name, neutral, charged, fields, per, charge, h, ok,
                  bath=max(s.concentration for s in species if s.valence > 0))


@dataclass(frozen=True)
class MutantRow3D:
    name: str
    measured_ratio: float
    one_d: float                  # the 1-D charged ratio (Round 6)
    ratios: dict[str, float]      # closure -> mutant / wild type, 3-D


def mutant_panel_3d(st: Structure | None = None, closures=CLOSURES_3D,
                    spacing: float | None = None, progress=None
                    ) -> tuple[Wall3D, list[MutantRow3D]]:
    """Xu 2006's RyR1 charge mutants in 3-D: each residue neutralised on all
    four subunits, the ratio to the wild type under every closure."""
    from ..io import loader
    from .ryr_mutants import mutant_panel, mutants, open_deposit
    st = st or loader.load(open_deposit())
    summary = measure_channel(st)
    wt = wall_3d(st, summary, closures=closures, spacing=spacing)
    wt_1d, rows_1d = mutant_panel(st)
    one_d = {r.name: r.charged_ratio for r in rows_1d}
    wt_pS = _P.value("permeation.published_ryr1")
    out = []
    for i, (name, (res, pS)) in enumerate(mutants().items()):
        if progress:
            progress(i, len(mutants()), name)
        m = wall_3d(st, summary, closures=closures, spacing=spacing,
                    neutralise=frozenset({res}))
        out.append(MutantRow3D(name, pS / wt_pS, one_d.get(name, float("nan")),
                               {c: m.charged[c] / wt.charged[c] for c in closures}))
    return wt, out



#: What :func:`wall_scan` can move: the placement's width, the lining rule,
#: the screening permittivity.
SCAN_KEYS = ("pore_charge.smoothing", "pore_charge.lining_margin",
             "permeation.permittivity_pore")


def wall_scan(st: Structure, key: str, values, spacing: float | None = None,
              pair_bridges: bool = False, closures=CLOSURES_3D
              ) -> list[tuple[float, Wall3D]]:
    """:func:`wall_3d` with one registered constant moved through ``values``
    (restored afterwards, whatever happens)."""
    if key not in SCAN_KEYS:
        raise ValueError(f"key must be one of {SCAN_KEYS}, not {key!r}")
    summary = measure_channel(st)
    before = _P.overrides().get(key)
    out = []
    try:
        for v in values:
            _P.set_value(key, v)
            out.append((float(v), wall_3d(st, summary, closures=closures,
                                          spacing=spacing,
                                          pair_bridges=pair_bridges)))
    finally:
        if before is None:
            _P.reset(key)
        else:
            _P.set_value(key, before)
    return out
