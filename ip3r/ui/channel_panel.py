"""The channel tab: what the loaded deposit measures as, and its pore.

Shows the :class:`~ip3r.structure.channel.ChannelSummary` of the current
structure — axis found two ways, C4 residual, numbering verdict, filter and
gate — and plots its pore profile. The state comparison and the unitary
conductance follow the loaded deposit's family (:meth:`set_paralog`): the
ITPR3 panel, or the curated RyR1 panel with its charge mutants. When the
deposit is 6DQN the profile
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
    lumen_toggled = pyqtSignal(bool)
    states_requested = pyqtSignal()
    unitary_requested = pyqtSignal()
    mutants_requested = pyqtSignal()

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
        self.show_lumen = QCheckBox("Draw the lumen (3-D, coloured by the "
                                    "potential) and plot where the voltage falls")
        self.show_lumen.setToolTip(
            "The ion-accessible volume of Round 7.6's 3-D solve, K+ (Cl- not "
            "drawn), coloured by the potential with the cytosol at 1 (red) and "
            "the lumen at 0 (blue). Solved on the deposit, so it is hidden on "
            "a morph or mode frame.")
        self.show_lumen.toggled.connect(self.lumen_toggled.emit)
        lay.addWidget(self.show_lumen)
        self.lumen_info = QLabel("")
        self.lumen_info.setWordWrap(True)
        lay.addWidget(self.lumen_info)
        self.states_btn = QPushButton("")
        self.states_btn.clicked.connect(self.states_requested.emit)
        lay.addWidget(self.states_btn)
        self.unitary_btn = QPushButton("Unitary K+ conductance of each state "
                                       "(drift-diffusion over the pore)")
        self.unitary_btn.clicked.connect(self.unitary_requested.emit)
        lay.addWidget(self.unitary_btn)
        self.mutants_btn = QPushButton("RyR1 charge mutants: model vs Xu 2006 "
                                       "(on the open deposit)")
        self.mutants_btn.clicked.connect(self.mutants_requested.emit)
        lay.addWidget(self.mutants_btn)
        self.unitary_rows = None
        self.mutant_rows = None
        self.panel_paralog = "ITPR3"
        self.set_paralog(None)
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

    def set_lumen_info(self, html: str) -> None:
        self.lumen_info.setText(html)

    def show_lumen_field(self, f, s) -> None:
        """The lumen's area and the potential along the window, 3-D against
        the 1-D model's inscribed circle; constrictions marked."""
        from ..parameters import PARAMETERS as _P
        axes = self.canvas.reset(2, 1)
        top, bottom = axes[0, 0], axes[1, 0]
        top.plot(f.z, f.area_3d, color=PALETTE[0], lw=1.4, label="3-D: lumen region")
        top.plot(f.z, f.area_1d, color=PALETTE[2], lw=1.0,
                 label="1-D: π (r_free − r_ion)²")
        top.set_ylabel(f"{f.species} area (Å²)")
        top.set_title(f"{f.name}: where the voltage falls, S0's window", fontsize=8)
        bottom.plot(f.z, f.drop_3d, color=PALETTE[0], lw=1.4, label="3-D (Laplace)")
        bottom.plot(f.z, f.drop_1d, color=PALETTE[2], lw=1.0, label="1-D (∫dz/A)")
        bottom.set_ylabel("φ, share of the\nwindow's drop")
        bottom.set_xlabel("z along the four-fold axis (Å; luminal ← → cytosolic)")
        bottom.set_ylim(-0.02, 1.02)
        w = _P.value("lumen.constriction_half_width")
        lines = [f.summary() + "."]
        for c in s.constrictions.values():
            for ax in (top, bottom):
                ax.axvline(c.z, color="#8a8f99", lw=0.6, ls=":")
            d3, d1 = f.drop_across(c.z, w)
            bottom.annotate(c.name, (c.z, 0.05), xytext=(3, 0),
                            textcoords="offset points", color="#d7dbe3", fontsize=7)
            if f.conducts:
                lines.append(f"{c.name} (z {c.z:+.1f} ± {w:.0f} Å): {d3:.0%} of "
                             f"the drop in 3-D, " + (f"{d1:.0%} in 1-D" if
                                                     np.isfinite(d1) else "1-D shut"))
        self.canvas.legend(top, loc="upper left")
        self.canvas.legend(bottom, loc="upper left")
        self.canvas.draw_now()
        self.set_lumen_info("<br>".join(lines) + "<br><i>Neutral pore, the "
                            "family's bath; each curve normalised to its own "
                            "drop across the window.</i>")

    def set_paralog(self, paralog: str | None) -> None:
        """Point the state and conductance buttons at the loaded deposit's
        panel: RyR1's for a RyR1 deposit, otherwise ITPR3's."""
        from ..core.annotations import is_ryr
        from ..io.registry import load_registry
        self.panel_paralog = "RYR1" if is_ryr(paralog) else "ITPR3"
        n = sum(1 for e in load_registry() if e.paralog == self.panel_paralog
                and (e.human or e.family == "RyR"))
        self.states_btn.setText(f"Compare the {self.panel_paralog} gating states "
                                f"(measures {n} deposits)")
        self.mutants_btn.setVisible(self.panel_paralog == "RYR1")

    def show_states(self, rows) -> None:
        """Overlay the pore profiles of the state panel."""
        self.states_btn.setEnabled(True)
        ax = self.canvas.reset()
        lines = []
        for k, r in enumerate(rows):
            p = r.summary.profile
            ax.plot(p.z, p.r_min, lw=1.6 if r.state in ("activated", "open") else 0.9,
                    color=PALETTE[k % len(PALETTE)],
                    label=f"{r.pdb_id} {r.state} (gate {r.radius('gate'):.2f} Å)")
            lines.append(f"{r.pdb_id} {r.state}: gate {r.radius('gate'):.2f} Å, "
                         f"filter {r.radius('filter'):.2f} Å")
        ax.set_xlabel("z along the four-fold axis (Å)")
        ax.set_ylabel("min heavy-atom radius (Å)")
        ax.set_ylim(0, 16)
        ax.set_title("ITPR3 state panel (ip3r_genes S11)"
                     if self.panel_paralog == "ITPR3" else
                     "RyR1 state panel (curated: scripts/curate_ryr.py)")
        self.canvas.legend(ax, loc="upper left")
        self.canvas.draw_now()
        self.text.setText("<br>".join(lines) + "<br><i>Each deposit measured "
                          "with the same axis, profile and constriction rules; "
                          "z origins differ slightly between deposits.</i>")

    def show_unitary(self, rows) -> None:
        """Bars of conductance per state against the measured values."""
        from ..physics.unitary import published
        self.unitary_btn.setEnabled(True)
        bath = next((u.bath for u in rows if u.bath), None)
        self.unitary_rows = rows
        ax = self.canvas.reset()
        x = np.arange(len(rows))
        series = (("neutral", "no wall charge"),
                  ("charged", "lining side chains charged"),
                  ("paired", "charged, salt bridges cancelled"))
        for k, (attr, label) in enumerate(series):
            vals = [getattr(u, attr).conductance_pS if getattr(u, attr).converged
                    else np.nan for u in rows]
            err = [[v - u.sweep[attr][0], u.sweep[attr][1] - v]
                   if u.sweep and v > 0 and np.all(np.isfinite(u.sweep[attr]))
                   else [0.0, 0.0] for u, v in zip(rows, vals)]
            ax.bar(x + (k - 1) * 0.28, vals, 0.28, color=PALETTE[k], label=label,
                   yerr=np.array(err).T, ecolor="#8a8f99", capsize=2)
        for k, (name, g) in enumerate(published(self.panel_paralog).items()):
            ax.axhline(g, color=PALETTE[len(series) + k], lw=1.0, ls="--",
                       label=f"measured: {name} {g:.0f} pS")
        for i, u in enumerate(rows):
            if not u.neutral.is_conducting:
                ax.annotate("shut", (i, 8), ha="center", color="#8a8f99", fontsize=7)
        ax.set_xticks(x, [f"{u.name}\n{u.state}" for u in rows], fontsize=6)
        from ..parameters import PARAMETERS as _P
        mM = 1000 * (bath or _P.value("permeation.bath_concentration"))
        ax.set_ylabel(f"K+ conductance, symmetric {mM:.0f} mM KCl (pS)")
        ax.set_title("Continuum model of each pore; whiskers: diffusivity "
                     "0.25-1x bulk, ion radius 1-2 Å", fontsize=8)
        self.canvas.legend(ax, loc="center left")
        self.canvas.draw_now()
        self.text.setText("<br>".join(u.row() for u in rows) + "<br><i>" + "; ".join(
            f"{u.name}: {u.charge.summary()}; paired: {u.paired_charge.summary()}"
            for u in rows if u.neutral.is_conducting)
            + ". A continuum of point ions in a pore a few ions wide: the "
            "comparison is of magnitude, not a fit.</i>")


    def show_mutants(self, result) -> None:
        """Measured vs modelled conductance ratio of each RyR1 charge mutant."""
        wt, rows = result
        self.mutants_btn.setEnabled(True)
        self.mutant_rows = rows
        ax = self.canvas.reset()
        x = np.arange(len(rows))
        series = (("measured_ratio", "measured (Xu 2006)"),
                  ("charged_ratio", "model: lining charges"),
                  ("paired_ratio", "model: salt bridges cancelled"))
        for k, (attr, label) in enumerate(series):
            ax.bar(x + (k - 1) * 0.28, [getattr(r, attr) for r in rows], 0.28,
                   color=PALETTE[(k + 3) % len(PALETTE)] if k == 0 else PALETTE[k],
                   label=label)
        ax.axhline(1.0, color="#8a8f99", lw=0.8, ls=":")
        ax.set_xticks(x, [r.name for r in rows])
        ax.set_ylabel("conductance, mutant / wild type")
        ax.set_ylim(0, 1.45)                    # room for the legend above 1
        ax.set_title(f"RyR1 charge neutralisation, modelled on {wt.name}; "
                     "a ratio cancels the unmeasured transport constants",
                     fontsize=8)
        self.canvas.legend(ax, loc="upper left")
        self.canvas.draw_now()
        self.text.setText(wt.row() + "<br>" + "<br>".join(r.row() for r in rows))
