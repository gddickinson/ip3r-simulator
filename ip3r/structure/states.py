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

__all__ = ["StateRow", "panel_entries", "state_panel"]


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


def panel_entries(paralog: str = "ITPR3", extended: bool = False) -> list:
    """The registry entries :func:`state_panel` measures: ``paralog``'s
    deposits in its reference numbering (human for an IP3R, the curated
    rabbit panel for RyR1), controls never, ``extended`` ones on request."""
    return [e for e in load_registry() if e.paralog == paralog
            and (e.human or e.family == "RyR") and not e.is_control
            and (extended or not e.is_extended)]


def state_panel(paralog: str = "ITPR3", progress=None,
                extended: bool = False) -> list[StateRow]:
    """Measure every :func:`panel_entries` deposit: the publication's panel,
    or with ``extended`` every registered full-length deposit too."""
    entries = panel_entries(paralog, extended)
    rows = []
    for i, e in enumerate(entries):
        if progress:
            progress(i, len(entries), e.pdb_id)
        s = measure_channel(loader.load(e.pdb_id))
        rows.append(StateRow(e.pdb_id, e.state, e.resolution, e.ip3_bound, s))
    return sorted(rows, key=lambda r: r.radius("gate"))
