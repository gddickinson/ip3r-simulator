"""The lumen and the potential in it: where the applied voltage falls
(Round 7.10).

Round 7.6's 3-D ohmic solve (:mod:`.ohmic3d`) returns a conductance, but the
potential it solves for says more: *where* along the pore the voltage drops.
A 1-D model puts the drop wherever the inscribed circle is smallest,
because it integrates ``dz / A(z)`` with ``A = π (r_free − r_ion)²``. The
real lumen has corners the circle misses, and current spreads between
narrow and wide parts. The two can disagree about which constriction holds
the field, and that is what a voltage-dependent blocker or a charged ring
would feel.

Both routes are read over S0's window (the pore-domain span ±
``PROFILE_MARGIN``), for one ion (the family bath's cation) and with the
same hard-sphere exclusion:

- **3-D**: φ from Laplace's equation in the voxelised lumen, averaged over
  each plane's lumen region (:meth:`PoreVolume.axis_regions`). φ = 1 at the
  cytosolic bath and 0 at the luminal one.
- **1-D**: φ from the cumulative ``∫ dz / A_1d`` across the window.

Each is normalised to its own drop across the window, so the curves say
*where inside the window* the voltage falls. ``window_share`` is the 3-D
drop's share of the whole applied voltage; the rest falls in the vestibules
and the baths.

No charge and no concentration change: the neutral reading, as in
:mod:`.ohmic3d`. The drawn surface is :mod:`ip3r.render.lumen_mesh`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..structure.channel import ChannelSummary, measure_channel
from ..structure.pore import PROFILE_MARGIN
from ..structure.pore_volume import PoreVolume, pore_volume
from .ohmic3d import geometric_conductance
from .permeation import potassium_species

__all__ = ["LumenField", "field_from_volume", "lumen_field", "plane_means"]


@dataclass
class LumenField:
    name: str
    species: str
    probe: float                  # Å, the ion radius the volume was cut for
    volume: PoreVolume = field(repr=False)
    phi: np.ndarray = field(repr=False)   # grid; NaN outside the conducting voxels
    g: float                      # m, geometric conductance (this ion's volume)
    window: tuple[float, float]   # Å, S0's window along z
    z: np.ndarray                 # plane centres inside the window, Å
    area_3d: np.ndarray           # Å², the lumen region of each plane
    area_1d: np.ndarray           # Å², π (r_free − probe)²
    phi_3d: np.ndarray            # raw φ, the mean over each plane's region
    window_share: float           # 3-D drop across the window / applied voltage
    converged: bool = True

    @property
    def conducts(self) -> bool:
        return self.volume.conducting > 0 and self.g > 0

    @property
    def shut_1d(self) -> bool:
        """The inscribed circle admits no ion somewhere in the window."""
        return not np.all(self.area_1d > 0)

    @property
    def drop_3d(self) -> np.ndarray:
        """3-D φ in the window, 0 at its luminal end and 1 at its cytosolic end."""
        lo, hi = self.phi_3d[0], self.phi_3d[-1]
        if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
            return np.full(len(self.z), np.nan)
        return (self.phi_3d - lo) / (hi - lo)

    @property
    def drop_1d(self) -> np.ndarray:
        """The 1-D model's φ, normalised the same way (NaN if it is shut)."""
        if self.shut_1d or len(self.z) < 2:
            return np.full(len(self.z), np.nan)
        step = np.gradient(self.z) / self.area_1d
        r = np.cumsum(step) - step / 2 - step[0] / 2
        return r / r[-1]

    def half_z(self, route: str = "3d") -> float:
        """z at which half the window's drop has fallen, Å (NaN if none)."""
        f = self.drop_3d if route == "3d" else self.drop_1d
        if not np.all(np.isfinite(f)):
            return float("nan")
        return float(np.interp(0.5, f, self.z))

    def drop_across(self, z0: float, half_width: float) -> tuple[float, float]:
        """Share of the window's drop that falls within ``z0 ± half_width``,
        (3-D, 1-D)."""
        def share(f):
            if not np.all(np.isfinite(f)):
                return float("nan")
            return float(np.interp(z0 + half_width, self.z, f)
                         - np.interp(z0 - half_width, self.z, f))
        return share(self.drop_3d), share(self.drop_1d)

    def summary(self) -> str:
        if not self.conducts:
            return f"{self.name}: no {self.species} path joins the two baths"
        text = (f"{self.name}, {self.species} (r {self.probe:.2f} Å): "
                f"{self.window_share:.0%} of the applied voltage falls in "
                f"the window; half of that is gone at z = {self.half_z('3d'):+.1f} Å")
        if self.shut_1d:
            return text + " (the inscribed circle is shut: no 1-D reading)"
        return text + f" (1-D: {self.half_z('1d'):+.1f} Å)"


def plane_means(vol: PoreVolume, grid: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """Mean of ``grid`` over each kept plane's lumen region (NaN where a
    plane has none)."""
    regions = [r for r, k in zip(vol.axis_regions(), keep) if k]
    plane = grid[:, :, keep]
    return np.array([np.nan if r is None else float(plane[:, :, k][r].mean())
                     for k, r in enumerate(regions)])


def field_from_volume(vol: PoreVolume, window: tuple[float, float],
                      r_free=None, name: str = "", species: str = "") -> LumenField:
    """Solve for φ in ``vol`` and read it over ``window``.

    ``r_free`` is the 1-D profile as ``(z, r_free)`` arrays (Å); without it
    the 1-D area is left out (NaN).
    """
    lap = geometric_conductance(vol)
    phi = np.full(vol.mask.shape, np.nan)
    phi[vol.mask] = lap.potential
    keep = (vol.zs >= window[0]) & (vol.zs <= window[1])
    regions = [r for r, k in zip(vol.axis_regions(), keep) if k]
    z = vol.zs[keep]
    area3 = np.array([0.0 if r is None else r.sum() * vol.spacing ** 2
                      for r in regions])
    mean = plane_means(vol, phi, keep)
    if r_free is None:
        area1 = np.full(len(z), np.nan)
    else:
        rf = np.interp(z, *r_free)
        area1 = np.pi * np.maximum(rf - vol.probe, 0.0) ** 2
    share = (float(mean[-1] - mean[0]) if len(mean) and
             np.isfinite(mean[[0, -1]]).all() else float("nan"))
    return LumenField(name=name, species=species, probe=vol.probe, volume=vol,
                      phi=phi, g=lap.g, window=window, z=z, area_3d=area3,
                      area_1d=area1, phi_3d=mean, window_share=share,
                      converged=lap.converged)


def lumen_field(st: Structure, summary: ChannelSummary | None = None,
                species=None, spacing: float | None = None) -> LumenField:
    """A deposit's lumen and potential for its family bath's cation (or
    ``species``), on the registered ``pore3d`` grid."""
    from .unitary import bath_for
    summary = summary or measure_channel(st)
    if species is None:
        paralog = summary.numbering.paralog if summary.numbering else None
        species = next(s for s in potassium_species(bath=bath_for(paralog))
                       if s.valence > 0)
    vol = pore_volume(st, summary.frame, summary.span, species.radius,
                      spacing=spacing)
    lo, hi = summary.span
    p = summary.profile
    return field_from_volume(vol, (lo - PROFILE_MARGIN, hi + PROFILE_MARGIN),
                             r_free=(p.z, p.r_free), name=st.name,
                             species=species.name)
