"""Round 7.25's candidate walls, built on one deposit's reversal grid so the
lumen box can draw them at reversal (Round 7.26).

Each is a charge map (mol/m³) and a Ca²⁺-only energy (kT) on the pore's
electrostatic voxels, read as Round 7.25's search read it: point ions under
Poisson, no charge–space fluid (the csc fluid is the deposit's own groups,
which a candidate replaces). The candidates:

* ``span well``: no charge, and a Ca²⁺-only well of
  ``wallsearch.gui_well_depth`` kT over the membrane span, the search's best
  uncharged wall on both 8TKF and 7T3T;
* ``span well + charge``: that well with the deposit's own lining charge
  (Round 7.23's ``local`` placement) added, as a binding site carries charge;
* ``ring pair``: no deposit charge, and the opposite-charge C4 ring pair of
  :func:`.wall_search.ring_search` with the highest B, the parallel-path
  route round the series bound (its search is two minutes of linear
  response, memoised per deposit and grid until a parameter changes);
* ``Ca2+ site`` (Round 7.28): no charge, and Round 7.27's compensated,
  K⁺-blocking saturable site over the span at ``casite.gui_depth`` kT (the
  depth at which it gives Vais's 15.2 on 8TKF). Its energies depend on the
  Ca²⁺ it holds, so it is a :mod:`.pnp3d` coupling, not a fixed excess:
  :meth:`CandidateWall.coupling` builds it on the reversal's domain.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from .charge3d import local_density
from .selectivity_bound import band_mask

__all__ = ["CANDIDATES", "CANDIDATE_LABELS", "CandidateWall",
           "candidate_wall", "best_ring"]

CANDIDATES = ("span well", "span well + charge", "ring pair", "Ca2+ site")
CANDIDATE_LABELS = {
    "span well": "Ca2+-only well over the span, uncharged (7.25's best)",
    "span well + charge": "that well + the deposit's charge",
    "ring pair": "opposite C4 ring pair at the highest B (7.25)",
    "Ca2+ site": "compensated Ca2+ site blocking K+, at 15.2 (7.27)",
}


@dataclass
class CandidateWall:
    name: str
    fixed: np.ndarray = field(repr=False)          # mol/m³ on pore.elec
    well: np.ndarray | None = field(repr=False)    # Ca2+ energy, kT; None
    description: str = ""
    site: object | None = None                     # ca_site.Site; None
    mask: np.ndarray | None = field(default=None, repr=False)  # its band

    def coupling(self, dom):
        """Round 7.27's :class:`.ca_site.SiteCoupling` on ``dom`` (None for
        a wall without a site)."""
        if self.site is None:
            return None
        from .ca_site import SiteCoupling, site_density
        density = site_density(self.mask, dom.elec.spacing,
                               _P.value("casite.sites"))
        return SiteCoupling(dom, self.site, self.mask, density)

    @property
    def excess(self) -> dict[str, np.ndarray] | None:
        """The species-specific energy :func:`.pnp3d.steady_state` holds."""
        return None if self.well is None else {"Ca2+": self.well}


_RINGS: dict[tuple[str, float], object] = {}
_P.subscribe(_RINGS.clear)


def best_ring(pore):
    """The converged ring pair of :func:`.wall_search.ring_search` with the
    highest B on ``pore``'s deposit and grid (memoised)."""
    from .wall_search import ring_search
    key = (pore.st.name, float(pore.spacing))
    if key not in _RINGS:
        rings = [r for r in ring_search(pore.st, spacing=pore.spacing)
                 if r.converged]
        if not rings:
            raise ValueError(f"no ring pair converged on {pore.st.name}")
        _RINGS[key] = max(rings, key=lambda r: r.bound)
    return _RINGS[key]


def candidate_wall(pore, name: str, deposit_fixed: np.ndarray) -> CandidateWall:
    """The candidate ``name`` on ``pore``; ``deposit_fixed`` is the deposit's
    own charge map (:func:`.reversal3d.wall`'s ``fixed``)."""
    from .wall_search import regions, ring_positions
    if name not in CANDIDATES:
        raise ValueError(f"candidate must be one of {CANDIDATES}, not {name!r}")
    if pore.ryr:
        raise ValueError("the candidate walls were searched against Vais "
                         "2010's pair: give an IP3R deposit")
    zero = np.zeros_like(deposit_fixed)
    if name == "ring pair":
        r = best_ring(pore)
        charges = np.r_[np.full(4, r.charge), np.full(4, -r.charge)]
        fixed, _ = local_density(pore.elec, ring_positions(pore.summary, r.z),
                                 charges)
        return CandidateWall(name, fixed, None,
                             f"opposite C4 rings at z {r.z:+.1f} Å, ±{r.charge:g} "
                             f"e per site (linear-response B {r.bound:.2f})")
    lo, hi = regions(pore.summary)["span"]
    if name == "Ca2+ site":
        from .ca_site import Site
        site = Site(_P.value("casite.gui_depth"), True,
                    _P.value("casite.block"), "span")
        return CandidateWall(name, zero, None,
                             f"{site.label} (z {lo:+.1f} to {hi:+.1f} Å, "
                             f"{_P.value('casite.sites'):g} sites), uncharged",
                             site, band_mask(pore, lo, hi))
    d = _P.value("wallsearch.gui_well_depth")
    well = -d * band_mask(pore, lo, hi)
    charged = name.endswith("charge")
    return CandidateWall(name, deposit_fixed if charged else zero, well,
                         f"Ca2+-only well of {d:g} kT over the span (z {lo:+.1f}"
                         f" to {hi:+.1f} Å)"
                         + (", with the deposit's charge" if charged else
                            ", uncharged"))
