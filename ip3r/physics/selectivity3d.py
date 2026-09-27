"""Selectivity in 3-D (Round 7.19): P_Ca:P_K from the charged lumen's
linear response, with the charge–space excess per voxel.

**The ruler.** Between identical baths at small voltage, each species'
conductance is its bulk σ_i times a Boltzmann-weighted Laplace solve on its
own volume (:mod:`.charged3d`), and σ_i = z_i² e² D_i c_i / kT. GHK's
permeability at 0 mV in the same bath is g_i = z_i² e² P_i c_i / kT, so
P_i = D_i × g_i/σ_i: the ratio P_Ca:P_K is D_Ca g_Ca / (D_K g_K), both
weighted Laplace conductances in metres. It is the permeability of the
independence regime, the same integral Round 7.17 took in 1-D
(:func:`.csc_readings.shares`, ``1/P = ∫ e^{E}/(D A) dz``), and it is not
a reversal potential: in 1-D the two differ (9HEO csc 0.87 in linear
response, 0.64 at Xu's reversal), and both 1-D readings are printed beside
the 3-D ones so the difference is carried, not hidden.

**The bath.** Symmetric: the family's KCl (RyR1 250 mM, IP3R the Vais bath)
with ``selectivity.ryr1_cacl2_lumen`` CaCl2 on both sides, Cl⁻ balancing.

**Readings** (:data:`READINGS`), each on the same volumes:
``neutral`` (shape and diffusivity alone); ``slice`` / ``local`` / ``pb``
(Round 7.11's closures, point ions); ``slice + csc`` / ``local + csc``
(:func:`.csc3d.local_csc`, every voxel a locally neutral charge–space
fluid); ``pb + csc`` (:func:`.csc3d.pb_csc`, Poisson's mean field with that
fluid's excess held per voxel).

Where each ion's resistance lies: the drop of its electrochemical
potential (the Laplace solution, 1 cytosol → 0 lumen, averaged over each
plane's lumen region) across ± ``lumen.constriction_half_width`` about the
filter and the gate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.channel import ChannelSummary, measure_channel
from ..structure.pore_volume import pore_volume
from .charge3d import group_positions, wall_field
from .csc3d import local_csc, pb_csc, wall_fluid
from .lumen_field import plane_means
from .ohmic3d import geometric_conductance
from .pore_charge import pore_charge

__all__ = ["READINGS", "CSC_READINGS", "Reading3D", "Selectivity3D", "Pore",
           "prepare", "selectivity_3d", "mutant_panel", "bath_species",
           "one_d_readings", "measured_ratio"]

#: Point-ion readings, then the charge–space ones (label -> placement).
READINGS = ("neutral", "slice", "local", "pb", "slice + csc", "local + csc",
            "pb + csc")
CSC_READINGS = {"slice + csc": "slice", "local + csc": "local",
                "pb + csc": "local"}
_CATION, _DIVALENT, _ANION = "K+", "Ca2+", "Cl-"


def bath_species(ryr: bool):
    """The symmetric mixed bath: the family's KCl + CaCl2 on both sides."""
    from .csc_readings import _mixed
    kcl = _P.value("permeation.ryr1_bath_concentration" if ryr
                   else "permeation.bath_concentration")
    return _mixed(kcl, _P.value("selectivity.ryr1_cacl2_lumen"))


def measured_ratio(ryr: bool) -> float:
    return _P.value("selectivity.published_ryr1_pca_pk" if ryr
                    else "selectivity.published_pca_pk")


@dataclass
class Reading3D:
    label: str
    g: dict[str, float]                 # weighted geometric conductance, m
    p: dict[str, float]                 # D x g, m³/s
    shares: dict[str, dict[str, float]]  # "filter"/"gate" -> species -> share
    peak: dict[str, float]              # M, highest in the lumen
    converged: bool = True

    @property
    def ratio(self) -> float:
        return self.p[_DIVALENT] / self.p[_CATION] if self.p[_CATION] else float("nan")

    @property
    def pcl_pk(self) -> float:
        return self.p[_ANION] / self.p[_CATION] if self.p[_CATION] else float("nan")


@dataclass
class Selectivity3D:
    name: str
    ryr: bool
    spacing: float
    readings: dict[str, Reading3D]
    measured: float
    one_d: dict[str, float] = field(default_factory=dict)
    unreached: list[str] = field(default_factory=list)
    net: float = 0.0

    def ratio(self, label: str) -> float:
        return self.readings[label].ratio

    @property
    def converged(self) -> bool:
        return all(r.converged for r in self.readings.values())


@dataclass
class Pore:
    """What every reading of one deposit shares: its channel, bath and the
    species' volumes (the electrostatics on the smallest ion's)."""

    st: Structure
    summary: ChannelSummary
    species: list
    vols: dict
    elec: object
    ryr: bool
    spacing: float


def prepare(st: Structure, summary: ChannelSummary | None = None,
            spacing: float | None = None) -> Pore:
    from ..core.annotations import is_ryr
    summary = summary or measure_channel(st)
    ryr = bool(summary.numbering and is_ryr(summary.numbering.paralog))
    species = bath_species(ryr)
    h = _P.value("pore3d.spacing") if spacing is None else spacing
    vols = {s.name: pore_volume(st, summary.frame, summary.span, s.radius,
                                spacing=h) for s in species}
    smallest = min(species, key=lambda s: s.radius)
    return Pore(st, summary, species, vols, vols[smallest.name], ryr, h)


def _shares(pore: Pore, name: str, potential: np.ndarray) -> dict[str, float]:
    vol = pore.vols[name]
    grid = np.full(vol.mask.shape, np.nan)
    grid[vol.mask] = potential
    means = plane_means(vol, grid, np.ones(len(vol.zs), bool))
    ok = np.isfinite(means)
    half = _P.value("lumen.constriction_half_width")
    out = {}
    for where in ("filter", "gate"):
        zc = pore.summary.constrictions[where].z
        ends = np.interp([zc - half, zc + half], vol.zs[ok], means[ok])
        out[where] = float(abs(ends[1] - ends[0]))
    return out


def _read(pore: Pore, label: str, energies) -> Reading3D:
    """``energies``: species name -> grid (kT), or None for the neutral pore."""
    g, p, peak, ok = {}, {}, {}, True
    shares: dict[str, dict[str, float]] = {"filter": {}, "gate": {}}
    for s in pore.species:
        e = None if energies is None else energies[s.name]
        lap = geometric_conductance(pore.vols[s.name], energy=e)
        g[s.name], p[s.name] = lap.g, s.diffusivity * lap.g
        ok &= lap.converged
        for where, v in _shares(pore, s.name, lap.potential).items():
            shares[where][s.name] = v
        m = pore.elec.mask
        peak[s.name] = (s.concentration if e is None
                        else float(s.concentration * np.exp(-e[m].min())))
    return Reading3D(label, g, p, shares, peak, bool(ok))


def one_d_readings(st: Structure, ryr: bool,
                   neutralise: frozenset[int] = frozenset()) -> dict[str, float]:
    """The 1-D pore in the same bath: linear response under donnan and csc
    (Round 7.17's integral), and csc at the family's reversal protocol."""
    from .csc_readings import csc_wall, shares
    from .selectivity import ryr1_calcium_ratio, selectivity
    w = csc_wall(st, neutralise)
    kcl = next(s.concentration for s in bath_species(ryr) if s.name == _CATION)
    out = {f"1-D {c}": shares(w, c, kcl=kcl).ratio for c in ("donnan", "csc")}
    if ryr:
        out["1-D csc, reversal"] = float(ryr1_calcium_ratio(
            w.z, w.radius, w.fixed, closure="csc", structural=w.wall)[1])
    else:
        out["1-D csc, reversal"] = float(selectivity(
            w.z, w.radius, w.fixed, closure="csc", structural=w.wall).pca_pk)
    return out


def selectivity_3d(st_or_pore, readings=READINGS,
                   neutralise: frozenset[int] = frozenset(),
                   pair_bridges: bool = False, spacing: float | None = None,
                   one_d: bool = True, **fluid_kw) -> Selectivity3D:
    """Every reading of one deposit (``fluid_kw`` to the csc fluid, e.g.
    ``wall_volume``)."""
    pore = (st_or_pore if isinstance(st_or_pore, Pore)
            else prepare(st_or_pore, spacing=spacing))
    st, summary, sp = pore.st, pore.summary, pore.species
    from .charged3d import _profile
    charge = pore_charge(st, summary.frame, _profile(st, summary),
                         pair_bridges=pair_bridges, neutralise=neutralise)
    positions = group_positions(st, summary.frame, charge.groups)
    out, fields, unreached = {}, {}, set()

    def field_of(closure):
        if closure not in fields:
            fields[closure] = wall_field(pore.elec, closure, sp, charge,
                                         positions=positions)
            unreached.update(fields[closure].unreached)
        return fields[closure]

    for label in readings:
        if label == "neutral":
            out[label] = _read(pore, label, None)
        elif label not in CSC_READINGS:
            wf = field_of(label)
            out[label] = _read(pore, label, {s.name: wf.energy(s.valence) for s in sp})
            out[label].converged &= wf.converged
        else:
            placement = CSC_READINGS[label]
            key = f"{placement} + csc"
            if key not in fields:
                fixed = field_of(placement).fixed
                wall = wall_fluid(pore.elec, charge.groups, positions, placement)
                fields[key] = (local_csc(pore.elec, fixed, sp, wall, **fluid_kw), fixed, wall)
            f, fixed, wall = fields[key]
            if label.startswith("pb"):
                f = pb_csc(pore.elec, fixed, sp, wall, reference=f, **fluid_kw)
            out[label] = _read(pore, label, {s.name: f.energy(s.name) for s in sp})
            out[label].converged &= f.converged
    return Selectivity3D(st.name, pore.ryr, pore.spacing, out,
                         measured_ratio(pore.ryr),
                         one_d_readings(st, pore.ryr, neutralise) if one_d else {},
                         sorted(unreached), charge.net_charge)


def mutant_panel(st: Structure | None = None, readings=READINGS,
                 spacing: float | None = None, progress=None, **kw
                 ) -> list[tuple[str, float, Selectivity3D]]:
    """Xu 2006's wild type and five mutants on the open RyR1 deposit:
    (name, measured P_Ca:P_K, reading). The volumes are shared."""
    from ..io import loader
    from .ryr_mutants import mutants, open_deposit
    st = st or loader.load(open_deposit())
    pore = prepare(st, spacing=spacing)
    key = "selectivity.published_ryr1_pca_pk"
    todo = [("wild type", frozenset(), _P.value(key))]
    for name, (res, _) in sorted(mutants().items(), key=lambda kv: kv[1][0]):
        todo.append((name, frozenset({res}), _P.value(f"{key}_{name.lower()}")))
    out = []
    for i, (name, off, pca) in enumerate(todo):
        if progress:
            progress(i, len(todo), name)
        out.append((name, pca, selectivity_3d(pore, readings, neutralise=off, **kw)))
    return out
