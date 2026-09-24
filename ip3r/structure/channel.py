"""One call that measures a tetrameric IP3R deposit the way S0 measured 6DQN.

:func:`measure_channel` chains the structural modules — four-fold axis (two
ways), C4 residual, numbering, pore-domain span, pore profile, filter and
gate, IP3 sites and their contacts — into a :class:`ChannelSummary`. The
findings checks read it, and so does the viewer's info panel, so the numbers
a user sees and the numbers that are checked are the same numbers.

The pore-domain span comes from the paralog's PF00520 element, which is only
meaningful if the deposit is in that paralog's human numbering; a deposit
that fails the numbering check (rat 7LHF) gets its span from the residues
nearest the axis in the membrane half of the protein instead, and says so.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.annotations import elements
from ..core.structure import Structure
from .ligand import contacts, ligand_sites
from .numbering import NumberingCheck, best_numbering
from .pore import (PROFILE_MARGIN, Constriction, PoreProfile, find_constrictions,
                   pore_profile, tm_span)
from .symmetry import (Frame, axis_by_centroids, axis_by_superposition,
                       c4_residual, subunit_ca, tetramer_frame)

__all__ = ["ChannelSummary", "measure_channel"]


@dataclass
class ChannelSummary:
    name: str
    frame: Frame
    superposition_angle: float
    axis_centroid: np.ndarray
    axis_disagreement_deg: float
    c4_residual: float
    numbering: NumberingCheck | None
    span: tuple[float, float]
    span_source: str
    profile: PoreProfile
    constrictions: dict[str, Constriction]
    ip3_contacts: dict[str, list] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    @property
    def n_ip3(self) -> int:
        return len(self.ip3_contacts)


def _geometric_span(st: Structure, frame: Frame) -> tuple[float, float]:
    """Axial extent of CA atoms within 30 Å of the axis on the -z half."""
    ca = frame.to_frame(st.xyz[st.mask_ca()])
    r = np.hypot(ca[:, 0], ca[:, 1])
    sel = (r < 30.0) & (ca[:, 2] < np.median(ca[:, 2]))
    z = ca[sel, 2]
    return float(np.percentile(z, 2)), float(np.percentile(z, 60))


def measure_channel(st: Structure, include_hetero: bool = True) -> ChannelSummary:
    """Measure a tetramer. ``include_hetero`` as S0 (keep HETATM heavy atoms)."""
    ca = subunit_ca(st)
    _, _, angle = axis_by_superposition(ca)
    ax_c, _, _ = axis_by_centroids(ca)
    frame = tetramer_frame(st)
    disagreement = float(np.degrees(np.arccos(min(1.0, abs(frame.axis @ ax_c)))))
    notes = []
    num = best_numbering(st)
    if num is not None:
        pore = [e for e in elements(num.paralog) if e.name == "channel"]
        span = tm_span(st, frame, (pore[0].start, pore[0].end))
        source = f"{num.paralog} PF00520 {pore[0].start}-{pore[0].end}"
    else:
        span = _geometric_span(st, frame)
        source = "geometric (no human numbering fits)"
        notes.append("numbering matches no human paralog; residue-keyed "
                     "annotation is not painted on this deposit")
    profile = pore_profile(st, frame, span[0] - PROFILE_MARGIN, span[1] + PROFILE_MARGIN,
                           include_hetero=include_hetero)
    cons = find_constrictions(st, frame, span, profile, include_hetero)
    ip3 = {}
    for site in ligand_sites(st):
        c = contacts(st, site)
        ip3[site.subunit] = c
    return ChannelSummary(st.name, frame, angle, ax_c, disagreement,
                          c4_residual(st, frame), num, span, source, profile,
                          cons, ip3, notes)
