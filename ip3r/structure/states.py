"""The ITPR3 state panel: how the pore changes between gating states.

``ip3r_genes`` S11 chose one ITPR3 deposition per conformational state (apo,
resting, labile resting, preactivated, activated, higher-order inhibited) by
stated rules. Measuring each with :func:`measure_channel` — same axis
method, same profile, same constriction rule — turns that panel into the
gating transition seen at the pore: in the activated state 8TKF the gate
widens from ~2.5 Å to ~5.9 Å while the filter barely moves.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..io import loader
from ..io.registry import load_registry
from .channel import ChannelSummary, measure_channel

__all__ = ["StateRow", "state_panel"]


@dataclass
class StateRow:
    pdb_id: str
    state: str
    resolution: float
    ip3_bound: bool
    summary: ChannelSummary

    def radius(self, name: str) -> float:
        c = self.summary.constrictions.get(name)
        return float("nan") if c is None else c.radius


def state_panel(paralog: str = "ITPR3", progress=None) -> list[StateRow]:
    """Measure every human deposit of ``paralog`` in the registry."""
    entries = [e for e in load_registry() if e.paralog == paralog and e.human]
    rows = []
    for i, e in enumerate(entries):
        if progress:
            progress(i, len(entries), e.pdb_id)
        s = measure_channel(loader.load(e.pdb_id))
        rows.append(StateRow(e.pdb_id, e.state, e.resolution, e.ip3_bound, s))
    return sorted(rows, key=lambda r: r.radius("gate"))
