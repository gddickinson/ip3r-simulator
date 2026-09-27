"""The gate's own geometry (Round 7.21): widen a deposit's gate and read
the conductance and P_Ca:P_K again.

Round 7.19 found 9HEO's uncharged gate holding over half of Ca²⁺'s
resistance and read that as the cap on P_Ca:P_K. The test is to take the
gate away. :func:`widen_gate` moves every atom within
``gate.widen_half_width`` of the gate constriction radially outward by Δ
times a cos² taper (full at the gate, zero at the window's edge). The move
is the same at every subunit, so the four-fold axis is unchanged, and it
is the same for every atom at one height, so the gate's radius grows by Δ
while its charges ride with their residues. The filter, ``10 Å`` luminal
of 9HEO's gate, lies outside the window.

:func:`gate_scan` then reads the neutral K⁺ conductance (Round 7.6's 3-D
solve) and the charged, charge–space P_Ca:P_K (Round 7.19's linear
response) at each Δ, with the filter and gate shares read at the
deposited constrictions, so "the gate" stays the same place as it opens.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.channel import ChannelSummary, measure_channel
from .ohmic3d import conductance_3d
from .selectivity3d import prepare, selectivity_3d

__all__ = ["widen_gate", "widening", "GateRow", "GateScan", "gate_scan", "READINGS",
           "deltas"]

#: The Round 7.19 readings the scan repeats.
READINGS = ("neutral", "pb", "local + csc", "pb + csc")


def widening(z: np.ndarray, z_gate: float, delta: float,
             half_width: float | None = None) -> np.ndarray:
    """Radial displacement (Å) at each height: Δ cos²(π(z − z_gate)/2w)
    inside ±w, zero outside."""
    w = _P.value("gate.widen_half_width") if half_width is None else half_width
    dz = np.asarray(z, float) - z_gate
    return np.where(np.abs(dz) < w, delta * np.cos(np.pi * dz / (2 * w)) ** 2, 0.0)


def widen_gate(st: Structure, summary: ChannelSummary, delta: float,
               half_width: float | None = None) -> Structure:
    """A copy of ``st`` with its gate widened by ``delta`` Å (see module)."""
    x = summary.frame.to_frame(st.xyz)
    r = np.hypot(x[:, 0], x[:, 1])
    shift = widening(x[:, 2], summary.constrictions["gate"].z, delta, half_width)
    scale = np.ones_like(r)
    off = r > 1e-9                       # an atom on the axis has no outward
    scale[off] = (r[off] + shift[off]) / r[off]
    x[:, 0] *= scale
    x[:, 1] *= scale
    xyz = summary.frame.from_frame(x).astype(st.xyz.dtype)
    return st.copy_with_coords(xyz, name=f"{st.name} gate +{delta:g} A")


def deltas() -> list[float]:
    top, step = _P.value("gate.widen_max"), _P.value("gate.widen_step")
    return [float(d) for d in np.arange(0.0, top + step / 2, step)]


@dataclass
class GateRow:
    """One widening: the gate's radius as measured, the neutral K⁺
    conductance (pS) and each reading's ratio and shares."""

    delta: float
    gate_radius: float
    narrowest: str                       # which constriction is now narrowest
    conductance: float
    ratios: dict = field(default_factory=dict)
    gate_share: dict = field(default_factory=dict)      # reading -> Ca2+ share
    filter_share: dict = field(default_factory=dict)
    converged: bool = True


@dataclass
class GateScan:
    name: str
    measured: float                      # the family's measured P_Ca:P_K
    rows: list


def gate_scan(st: Structure, widths=None, readings=READINGS,
              spacing: float | None = None, progress=None) -> GateScan:
    """:func:`widen_gate` at each Δ, measured and read (see module)."""
    summary = measure_channel(st)
    rows, measured = [], float("nan")
    widths = deltas() if widths is None else widths
    for i, d in enumerate(widths):
        if progress:
            progress(i, len(widths), d)
        wide = widen_gate(st, summary, d)
        seen = measure_channel(wide)
        radius = seen.profile.at(summary.constrictions["gate"].z)
        narrow = min(seen.constrictions, key=lambda k: seen.constrictions[k].radius)
        g = conductance_3d(wide, summary, spacing=spacing)
        # The deposited frame and constrictions: the gate stays where it was.
        sel = selectivity_3d(prepare(wide, summary, spacing=spacing),
                             readings=readings, one_d=False)
        measured = sel.measured
        rows.append(GateRow(
            d, radius, narrow, g.conductance_pS,
            {k: x.ratio for k, x in sel.readings.items()},
            {k: x.shares["gate"]["Ca2+"] for k, x in sel.readings.items()},
            {k: x.shares["filter"]["Ca2+"] for k, x in sel.readings.items()},
            g.converged and sel.converged))
    return GateScan(st.name, measured, rows)
