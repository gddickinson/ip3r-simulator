"""The conduction pathway: radius profile along the four-fold axis.

Two radii are reported, and they answer different questions:

* ``r_min`` — the minimum distance from the axis to any protein heavy-atom
  **centre** within a slab. This is the quantity ``ip3r_genes`` S0 committed
  (``structure_pore.tsv``), a lower bound that models no atom size, the honest
  number for a 3.3 Å map. Computing it here with the same definition is what
  lets the published profile be re-derived.
* ``r_free`` — the same minimum with each atom's van der Waals radius
  subtracted: the largest sphere centred on the axis that fits. This is the
  quantity an ion "sees", and the one the viewer draws as the pore surface.

The narrowest luminal point is the selectivity filter and the narrowest
cytosolic point the gate — S0's rule, applied to the half of the
pore-domain span either side of its midpoint.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .symmetry import Frame, subunit_ca

__all__ = ["PoreProfile", "Constriction", "pore_profile", "tm_span",
           "find_constrictions", "lining_residues"]


@dataclass
class PoreProfile:
    z: np.ndarray            # axial position, Å (frame coordinates)
    r_min: np.ndarray        # min heavy-atom centre distance to the axis, Å
    r_free: np.ndarray       # the same less the atom's vdW radius, Å
    meta: dict = field(default_factory=dict)

    def at(self, z: float) -> float:
        return float(np.interp(z, self.z, self.r_min))


@dataclass
class Constriction:
    name: str
    z: float
    radius: float
    residues: list[str]


def _heavy(st: Structure, include_hetero: bool = False) -> np.ndarray:
    """Heavy atoms that can line the pore.

    ``include_hetero=True`` keeps every non-water heavy atom, which is what S0
    measured (its parser kept HETATM records); the default keeps protein only,
    which is what an ion sees once bound lipids and detergent are discounted.
    """
    m = st.element != "H"
    return m if include_hetero else m & st.mask_protein() & ~st.hetero


def tm_span(st: Structure, frame: Frame, residues: tuple[int, int]) -> tuple[float, float]:
    """Axial extent (2nd-98th percentile) of a residue range in subunit 1."""
    ca = subunit_ca(st)[frame.chains[0]]
    pts = np.array([xyz for r, xyz in ca.items()
                    if residues[0] <= r <= residues[1]])
    if len(pts) == 0:
        raise ValueError(f"no resolved residues in {residues}")
    z = frame.to_frame(pts)[:, 2]
    return float(np.percentile(z, 2)), float(np.percentile(z, 98))


def pore_profile(st: Structure, frame: Frame, z0: float, z1: float,
                 step: float | None = None, slab: float | None = None,
                 include_hetero: bool = False) -> PoreProfile:
    step = _P.value("pore.step") if step is None else step
    slab = _P.value("pore.slab") if slab is None else slab
    m = _heavy(st, include_hetero)
    f = frame.to_frame(st.xyz[m])
    vdw = st.vdw_radii()[m]
    r = np.hypot(f[:, 0], f[:, 1])
    z = f[:, 2]
    order = np.argsort(z)
    z_sorted = z[order]
    zs, rmin, rfree = [], [], []
    # Grid anchored the way S0 anchored it (np.arange from z0), so the two
    # profiles land on comparable sample points.
    for zi in np.arange(z0, z1 + step, step):
        lo, hi = np.searchsorted(z_sorted, [zi - slab, zi + slab + 1e-9])
        idx = order[lo:hi]
        if len(idx) < 4:
            continue
        zs.append(float(zi))
        rmin.append(float(r[idx].min()))
        rfree.append(float((r[idx] - vdw[idx]).min()))
    return PoreProfile(np.array(zs), np.array(rmin), np.array(rfree),
                       meta={"step": step, "slab": slab, "z0": z0, "z1": z1,
                             "include_hetero": include_hetero})


def lining_residues(st: Structure, frame: Frame, z: float, radius: float,
                    slab: float | None = None, tol: float | None = None,
                    include_hetero: bool = False) -> list[str]:
    """Residues with a heavy atom at the constriction, as ``RESNUM`` labels."""
    slab = _P.value("pore.lining_slab") if slab is None else slab
    tol = _P.value("pore.lining_tol") if tol is None else tol
    m = _heavy(st, include_hetero)
    f = frame.to_frame(st.xyz[m])
    r = np.hypot(f[:, 0], f[:, 1])
    sel = (np.abs(f[:, 2] - z) <= slab) & (r <= radius + tol)
    labels = {f"{st.res_name[m][i]}{st.res_seq[m][i]}" for i in np.flatnonzero(sel)}
    return sorted(labels, key=lambda s: (int("".join(c for c in s if c.isdigit())), s))


def find_constrictions(st: Structure, frame: Frame, span: tuple[float, float],
                       profile: PoreProfile | None = None,
                       include_hetero: bool = False) -> dict[str, Constriction]:
    """Filter (narrowest luminal point) and gate (narrowest cytosolic point).

    ``span`` is the pore-domain axial extent; the profile is taken 12 Å beyond
    each end and the gate searched up to 6 Å past the cytosolic end — S0's
    windows, so the same constriction is found.
    """
    lo, hi = span
    prof = profile or pore_profile(st, frame, lo - 12.0, hi + 12.0,
                                   include_hetero=include_hetero)
    mid = 0.5 * (lo + hi)
    out = {}
    for name, sel in (("filter", prof.z < mid),
                      ("gate", (prof.z >= mid) & (prof.z <= hi + 6.0))):
        if not sel.any():
            continue
        i = np.flatnonzero(sel)[np.argmin(prof.r_min[sel])]
        out[name] = Constriction(name, float(prof.z[i]), float(prof.r_min[i]),
                                 lining_residues(st, frame, float(prof.z[i]),
                                                 float(prof.r_min[i]),
                                                 include_hetero=include_hetero))
    return out
