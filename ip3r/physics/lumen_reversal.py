"""The lumen at reversal: each ion's concentration and electrochemical drop
in Round 7.23's steady state (Round 7.24).

Rounds 7.10–7.18 coloured the lumen at equilibrium or in linear response
between identical baths, and only for K+ (Round 7.10's open item: Cl⁻ was
never drawn). The measured selectivities are read where the baths differ and
no net current flows. There every ion still carries a current, and this
module reads that state on the voxels, one experiment of the family's
protocol (:func:`.reversal3d.experiments`) under one reading
(:data:`.reversal3d.READINGS`):

- **concentration** ``c_i = n_i e^{−(z_i ψ + μ_i)}`` (M): what Poisson
  counts, so on the whole electrostatic volume (the smallest ion's). Where a
  larger ion's centre cannot reach, it takes its Slotboom n from its nearest
  own voxel, as :mod:`.pnp3d` does; the value drawn is the model's.
- **electrochemical drop** ``(n_i − n_lumen) / (n_cytosol − n_lumen)``:
  0 at the luminal bath, 1 at the cytosolic. The species' flux is
  ``D_i e^{−E_i} ∇n_i``, so where this rises steeply is where that ion's
  resistance lies. Each ion has its own; at reversal their currents cancel,
  not vanish.

The surface is cut from the reversal's own grid (``reversal3d.spacing``,
the electrostatic volume), so it is not the equilibrium box's 0.5 Å K+
lumen; the panel says so. Plane means over S0's window give the plot.

Round 7.26: ``reading`` may instead name one of Round 7.25's candidate walls
(:data:`.wall_candidates.CANDIDATES`), solved as the search solved it: point
ions under Poisson with that wall's charge and Ca²⁺-only energy.

Round 7.28: the ``Ca2+ site`` candidate is Round 7.27's compensated,
K⁺-blocking site, solved with its coupling; the reading then also carries
the site's occupancy θ and K⁺'s block energy −ln(1 − fθ) on the grid
(zero outside the band, NaN off the lumen), and the Ca²⁺ ions it holds.

Round 7.30: ``ca`` replaces the Ca²⁺ experiment's luminal CaCl₂ (IP3R
only) as Round 7.29's sweep does (:func:`.mole_fraction.experiment`), so
the site can be drawn filling; at Vais's own 10 mM it is his experiment.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.channel import ChannelSummary, measure_channel
from ..structure.pore import PROFILE_MARGIN
from .lumen_field import LumenField, field_from_volume, plane_means
from .pnp3d import domain
from .reversal3d import (READINGS, _impermeant, concentration, experiments,
                         reversal, steady_at_reversal, wall)
from .selectivity import thermal_voltage
from .selectivity3d import prepare
from .wall_candidates import CANDIDATES, candidate_wall

__all__ = ["ReversalLumen", "reversal_lumen", "ion_grids", "SPECIES",
           "READINGS", "CANDIDATES"]

#: Every ion a protocol can carry, in the order the panel lists them.
SPECIES = ("K+", "Cl-", "Ca2+")


@dataclass
class ReversalLumen:
    """One experiment under one reading at its reversal, on the lumen."""

    lumen: LumenField = field(repr=False)   # the reversal grid's own lumen
    reading: str
    experiment: str
    v: float                                # V_rev, V (cytosol − lumen)
    conc: dict[str, np.ndarray] = field(repr=False)   # M; NaN off the lumen
    drop: dict[str, np.ndarray] = field(repr=False)   # 0..1; NaN off the lumen
    currents: dict[str, float] = field(default_factory=dict)  # A at V_rev
    baths: dict[str, tuple[float, float]] = field(default_factory=dict)  # M
    converged: bool = True
    solves: int = 0
    wall: str = ""                          # a candidate wall's description
    #: Round 7.28, a site's reading: θ and K+'s block energy (kT) on the
    #: grid, and the Ca2+ ions held; None / nan without a site or Ca2+.
    occupancy: np.ndarray | None = field(default=None, repr=False)
    k_block: np.ndarray | None = field(default=None, repr=False)
    held: float = float("nan")

    @property
    def ca(self) -> float:
        """Luminal CaCl₂ (M) of the experiment; 0 without Ca²⁺."""
        return self.baths.get("Ca2+", (0.0, 0.0))[0]

    @property
    def candidate(self) -> bool:
        """A Round 7.25 candidate wall, not the deposit's own."""
        return self.reading in CANDIDATES

    @property
    def z(self) -> np.ndarray:
        return self.lumen.z

    @property
    def species(self) -> list[str]:
        return [s for s in SPECIES if s in self.conc]

    def _means(self, grid: np.ndarray) -> np.ndarray:
        vol = self.lumen.volume
        lo, hi = self.lumen.window
        return plane_means(vol, grid, (vol.zs >= lo) & (vol.zs <= hi))

    def conc_3d(self, name: str) -> np.ndarray:
        """Mean concentration (M) over each window plane's lumen region."""
        return self._means(self.conc[name])

    def peak(self, name: str) -> tuple[float, float]:
        """The highest concentration (M) in the window's lumen, and its z."""
        vol = self.lumen.volume
        lo, hi = self.lumen.window
        c = np.where((vol.zs >= lo) & (vol.zs <= hi), self.conc[name], np.nan)
        k = np.unravel_index(int(np.nanargmax(c)), c.shape)
        return float(c[k]), float(vol.zs[k[2]])

    def drop_3d(self, name: str) -> np.ndarray:
        """The ion's drop over the window, 0 at its luminal end and 1 at its
        cytosolic (NaN if it does not rise across the window)."""
        f = self._means(self.drop[name])
        lo, hi = f[0], f[-1]
        if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
            return np.full(len(f), np.nan)
        return (f - lo) / (hi - lo)

    def drop_across(self, name: str, z0: float, half_width: float) -> float:
        """Share of the ion's window drop within ``z0 ± half_width``."""
        f = self.drop_3d(name)
        if not np.all(np.isfinite(f)):
            return float("nan")
        return float(np.interp(z0 + half_width, self.z, f)
                     - np.interp(z0 - half_width, self.z, f))

    def steepest_z(self, name: str) -> float:
        """z where the ion's drop rises fastest (Å): its resistance."""
        f = self.drop_3d(name)
        if not np.all(np.isfinite(f)) or len(f) < 2:
            return float("nan")
        return float(self.z[int(np.argmax(np.gradient(f, self.z)))])

    def summary(self) -> str:
        cur = ", ".join(f"{k} {self.currents[k] * 1e12:+.2f} pA"
                        for k in self.species)
        lum = (f" (luminal CaCl2 {self.ca * 1e3:.3g} mM)" if self.ca > 0
               else "")
        text = (f"{self.lumen.name}, {self.experiment} experiment{lum}, "
                f"{self.reading}{f' ({self.wall})' if self.wall else ''}, "
                f"at reversal V = {self.v * 1e3:+.2f} mV"
                f"{'' if self.converged else ' (n.c.)'}; each ion's current "
                f"there: {cur}")
        for name in self.species:
            c, z = self.peak(name)
            text += (f"; {name} peaks at {c:.3g} M (z {z:+.1f} Å), its drop "
                     f"steepest at z {self.steepest_z(name):+.1f} Å")
        if self.occupancy is not None:
            text += (f"; the site holds {self.held:.2f} Ca2+ (θ up to "
                     f"{np.nanmax(self.occupancy):.2f}), K+'s block up to "
                     f"{np.nanmax(self.k_block):.2f} kT")
        return text

    def occupancy_3d(self) -> np.ndarray | None:
        """Plane-mean θ over the window's lumen (None without a site)."""
        return None if self.occupancy is None else self._means(self.occupancy)


def ion_grids(dom, steady, species, excess=None):
    """Each species' concentration (M) and electrochemical drop (0 at the
    luminal bath, 1 at the cytosolic) on the grid of ``steady``, NaN off the
    electrostatic volume; and its baths (lumen, cytosol) in M. The drop is
    NaN throughout for a species whose n is the same at both baths."""
    mask = dom.elec.mask
    conc, drop, baths = {}, {}, {}
    for s in species:
        c = concentration(dom, steady, s.name, s.valence, excess) / 1000.0
        conc[s.name] = np.where(mask, c, np.nan)
        n_lum = s.concentration * 1000.0
        n_cyt = s.right * 1000.0 * np.exp(s.valence * steady.v / thermal_voltage())
        d = ((steady.n[s.name] - n_lum) / (n_cyt - n_lum) if n_cyt != n_lum
             else np.full(mask.shape, np.nan))
        drop[s.name] = np.where(mask, d, np.nan)
        baths[s.name] = (s.concentration, s.right)
    return conc, drop, baths


def reversal_lumen(st: Structure, reading: str = "pb + csc",
                   experiment: str = "Ca2+",
                   summary: ChannelSummary | None = None,
                   spacing: float | None = None, ca: float | None = None,
                   **fluid_kw) -> ReversalLumen:
    """``experiment`` of ``st``'s family protocol under ``reading`` (or one
    of the candidate walls), solved to its reversal on
    ``reversal3d.spacing``'s grid (or ``spacing``) and read on the lumen.
    ``ca`` (M) replaces the luminal CaCl₂ of IP3R's Ca²⁺ experiment."""
    if reading not in READINGS + CANDIDATES:
        raise ValueError(f"reading must be one of {READINGS + CANDIDATES}, "
                         f"not {reading!r}")
    summary = summary or measure_channel(st)
    h = _P.value("reversal3d.spacing") if spacing is None else spacing
    pore = prepare(st, summary, spacing=h)
    exps = {e.name: e for e in experiments(pore.ryr)}
    if experiment not in exps:
        raise ValueError(f"{'RyR1' if pore.ryr else 'IP3R'}'s protocol has no "
                         f"{experiment} experiment (it has {sorted(exps)})")
    exp = exps[experiment]
    if ca is not None:
        exp = _with_calcium(pore, exp, ca)
    lo, hi = summary.span
    p = summary.profile
    smallest = min(pore.species, key=lambda s: s.radius)
    f = field_from_volume(pore.elec, (lo - PROFILE_MARGIN, hi + PROFILE_MARGIN),
                          r_free=(p.z, p.r_free), name=st.name,
                          species=smallest.name)
    shut = [s.name for s in exp.species() if not pore.vols[s.name].conducting]
    if not f.conducts or shut:
        raise ValueError(f"{st.name} is shut: no {', '.join(shut) or smallest.name}"
                         " path joins the two baths, so there is no reversal "
                         "to solve")
    dom = domain(pore.elec, pore.vols)
    w = wall(pore)
    note, site = "", {}
    if reading in CANDIDATES:
        cw = candidate_wall(pore, reading, w.fixed)
        names = {s.name for s in exp.species()}
        excess = {k: e for k, e in (cw.excess or {}).items() if k in names} or None
        coupling = cw.coupling(dom)
        v, steady, n = reversal(dom, exp.species(), cw.fixed, excess,
                                _impermeant(pore, exp), coupling=coupling)
        ok, note = steady.converged, cw.description
        if coupling is not None and "Ca2+" in names:
            site = _site_grids(dom, steady, coupling, pore.elec.mask)
    else:
        v, steady, excess, ok, n = steady_at_reversal(pore, dom, w, reading,
                                                      exp, **fluid_kw)
    conc, drop, baths = ion_grids(dom, steady, exp.species(), excess)
    return ReversalLumen(f, reading, experiment, float(v), conc, drop,
                         dict(steady.currents), baths,
                         bool(ok and f.converged), n, note, **site)


def _with_calcium(pore, exp, ca: float):
    """Vais's Ca²⁺ experiment with ``ca`` M luminal CaCl₂ (Round 7.29's)."""
    from .mole_fraction import experiment
    if pore.ryr or exp.name != "Ca2+":
        raise ValueError("luminal Ca2+ is swept in Vais 2010's Ca2+ experiment "
                         "only: give an IP3R deposit and that experiment")
    if not ca > 0:
        raise ValueError(f"luminal CaCl2 must be positive, not {ca!r}")
    return experiment(float(ca))


def _site_grids(dom, steady, coupling, mask) -> dict:
    """θ, K+'s block energy (NaN off the lumen) and the ions held, from the
    solved state as :func:`.ca_site.read` takes them."""
    from .ca_site import N_AVOGADRO
    psi = steady.v / thermal_voltage() * dom.phi0 + steady.u
    c_free = coupling.free(psi, steady.n["Ca2+"])
    theta = coupling.occupancy(c_free)
    visible, _ = coupling.energies(c_free)
    block = visible.get("K+", np.zeros(mask.shape))
    held = float(theta[coupling.mask].sum() * coupling.s
                 * (dom.elec.spacing * 1e-10) ** 3 * N_AVOGADRO)
    return {"occupancy": np.where(mask, theta, np.nan),
            "k_block": np.where(mask, block, np.nan), "held": held}

