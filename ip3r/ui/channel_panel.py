"""The channel tab: what the loaded deposit measures as, and its pore.

Shows the :class:`~ip3r.structure.channel.ChannelSummary` of the current
structure — axis found two ways, C4 residual, numbering verdict, filter and
gate — and plots its pore profile. When the deposit is 6DQN the profile
``ip3r_genes`` S0 committed is drawn on the same axes, so the published and
the recomputed curves can be compared by eye as well as by the check.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QCheckBox, QLabel, QPushButton, QVBoxLayout, QWidget

from ..core import genes_data as G
from .plot_canvas import PALETTE, PlotCanvas

__all__ = ["ChannelPanel"]


class ChannelPanel(QWidget):
    pore_toggled = pyqtSignal(bool)
    states_requested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        self.text = QLabel("Load a structure to measure it.")
        self.text.setWordWrap(True)
        lay.addWidget(self.text)
        self.show_pore = QCheckBox("Draw the pore in the viewport "
                                   "(probe radius = r_free)")
        self.show_pore.toggled.connect(self.pore_toggled.emit)
        lay.addWidget(self.show_pore)
        self.states_btn = QPushButton("Compare the ITPR3 gating states "
                                      "(measures 7 deposits)")
        self.states_btn.clicked.connect(self.states_requested.emit)
        lay.addWidget(self.states_btn)
        self.canvas = PlotCanvas(self, height=3.6)
        lay.addWidget(self.canvas, 1)

    def show_summary(self, s) -> None:
        num = s.numbering
        num_txt = (f"{num.paralog} numbering ({num.identity:.1%} of {num.n_compared} "
                   f"side-chain residues identical; {num.n_stubbed} stubbed "
                   f"residues not scored"
                   + (f"; mismatch segments {', '.join(f'{a}-{b}' for a, b in num.mismatch_segments[:4])}"
                      if num.mismatch_segments else "") + ")"
                   if num else "matches no human paralog numbering")
        cons = "<br>".join(
            f"<b>{c.name}</b>: r = {c.radius:.2f} Å at z = {c.z:+.1f} Å, lined by "
            f"{', '.join(c.residues)}" for c in s.constrictions.values())
        ip3 = (f"{s.n_ip3} IP3 site(s); contacts on subunit "
               + "; ".join(f"{k}: {len(v['same_subunit'])}" for k, v in s.ip3_contacts.items())
               if s.n_ip3 else "no IP3 bound")
        self.text.setText(
            f"<b>{s.name}</b><br>"
            f"Subunit→subunit rotation {s.superposition_angle:.2f}° (C4 expects 90°); "
            f"superposition and centroid axes differ by {s.axis_disagreement_deg:.3f}°; "
            f"C4 residual {s.c4_residual:.3f} Å.<br>"
            f"Numbering: {num_txt}.<br>"
            f"Pore-domain span z = {s.span[0]:.1f} to {s.span[1]:.1f} Å ({s.span_source}).<br>"
            f"{cons}<br>{ip3}."
            + ("<br><i>" + "; ".join(s.notes) + "</i>" if s.notes else ""))
        self._plot(s)

    def _plot(self, s) -> None:
        ax = self.canvas.reset()
        p = s.profile
        ax.plot(p.z, p.r_min, color=PALETTE[0], lw=1.4,
                label="recomputed: min heavy-atom distance to axis")
        ax.plot(p.z, np.maximum(p.r_free, 0), color=PALETTE[2], lw=1.0,
                label="recomputed: free radius (less vdW)")
        if s.name.upper() == "6DQN" and G.available(
                "s0_baseline/review_figures/structure_pore.tsv"):
            rows = G.read_tsv("s0_baseline/review_figures/structure_pore.tsv")
            ax.plot([float(r["z_along_axis"]) for r in rows],
                    [float(r["min_heavy_atom_radius"]) for r in rows],
                    color=PALETTE[1], lw=1.0, ls="--", label="published (ip3r_genes S0)")
        for c in s.constrictions.values():
            ax.axvline(c.z, color="#8a8f99", lw=0.6, ls=":")
            ax.annotate(c.name, (c.z, c.radius), xytext=(4, 8),
                        textcoords="offset points", color="#d7dbe3", fontsize=7)
        ax.axvspan(*s.span, color="#5b9ef2", alpha=0.07, lw=0)
        ax.set_xlabel("z along the four-fold axis (Å; luminal ← → cytosolic)")
        ax.set_ylabel("radius (Å)")
        ax.set_ylim(0, 20)
        self.canvas.legend(ax, loc="upper left")
        self.canvas.draw_now()

    def show_states(self, rows) -> None:
        """Overlay the pore profiles of the ITPR3 state panel."""
        self.states_btn.setEnabled(True)
        ax = self.canvas.reset()
        lines = []
        for k, r in enumerate(rows):
            p = r.summary.profile
            ax.plot(p.z, p.r_min, lw=1.6 if r.state == "activated" else 0.9,
                    color=PALETTE[k % len(PALETTE)],
                    label=f"{r.pdb_id} {r.state} (gate {r.radius('gate'):.2f} Å)")
            lines.append(f"{r.pdb_id} {r.state}: gate {r.radius('gate'):.2f} Å, "
                         f"filter {r.radius('filter'):.2f} Å")
        ax.set_xlabel("z along the four-fold axis (Å)")
        ax.set_ylabel("min heavy-atom radius (Å)")
        ax.set_ylim(0, 16)
        ax.set_title("ITPR3 state panel (ip3r_genes S11): the gate opens only "
                     "in the activated state")
        self.canvas.legend(ax, loc="upper left")
        self.canvas.draw_now()
        self.text.setText("<br>".join(lines) + "<br><i>Each deposit measured "
                          "with the same axis, profile and constriction rules; "
                          "z origins differ slightly between deposits.</i>")
