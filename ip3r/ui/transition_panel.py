"""The transition tab: one paralog's two states, morphed and scored.

Choose an end state for the loaded deposit (or press the family's preset:
ITPR3 resting → activated 8TKG → 8TKF, or RyR1 primed → open from the
curated panel's ``morph_start``/``morph_end`` roles), build, then scrub or play the morph, colour the
receptor by how far each residue moved, and read how well the elastic
network of the start state predicts the move.

Everything shown says what it is: the morph is an interpolation (side
chains too: each atom both deposits resolve is interpolated to its end
position, so the gate plotted along the path is the gate drawn, and the
rigid-side-chain gate beside it shows what the shortcut would have said),
and the overlap is printed beside the value a random direction of the same
symmetry would get. The headline plot is the collective A modes taken
together (``TransitionOverlap.subspace``): which single A mode carries the
move is a cutoff choice, the A subspace is not.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QHBoxLayout,
                             QLabel, QPushButton, QSlider, QVBoxLayout, QWidget)

from ..core.annotations import ELEMENT_LABELS
from ..io.registry import load_registry
from ..structure.morph import NOTE
from .plot_canvas import PALETTE, PlotCanvas

__all__ = ["TransitionPanel", "PRESET", "preset_for", "LOWER_PLOTS"]

#: The resting -> activated pair of ITPR3 (S11's state panel; Round 2).
PRESET = ("8TKG", "8TKF")


def preset_for(family: str) -> tuple[str, str, str] | None:
    """``(start, end, label)``: ITPR3's pair, or RyR1's same-preparation
    primed -> open pair as ``scripts/curate_ryr.py`` marked it."""
    if family != "RyR":
        return (*PRESET, "Resting → activated")
    reg = load_registry()
    start = next((e.pdb_id for e in reg if "morph_start" in e.roles), None)
    end = next((e.pdb_id for e in reg if "morph_end" in e.roles), None)
    return (start, end, "Primed → open") if start and end else None

#: What the lower plot can show (the A subspace is always on top).
LOWER_PLOTS = {"gate": "the pore along the morph",
               "modes": "every mode, by symmetry",
               "elements": "mean displacement by element"}

_IRREP_COLOUR = {"A": PALETTE[1], "B": PALETTE[0], "E": PALETTE[2], "mixed": "#8a8f99"}


class TransitionPanel(QWidget):
    build_requested = pyqtSignal(str, str, str)     # end id, fit, method
    preset_requested = pyqtSignal(str, str)         # start id, end id
    frame_requested = pyqtSignal(int)
    play_requested = pyqtSignal()
    stop_requested = pyqtSignal()
    paint_toggled = pyqtSignal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        note = QLabel(NOTE + " C-alphas follow the path, and every atom both "
                      "deposits resolve is interpolated to its end position; "
                      "ligands and unmatched atoms ride their C-alpha rigidly.")
        note.setWordWrap(True)
        lay.addWidget(note)

        form = QFormLayout()
        self.start = QLabel("—")
        self.end = QComboBox()
        self.fit = QComboBox()
        self.fit.addItem("pore domain (membrane held still)", "pore")
        self.fit.addItem("whole tetramer (least RMSD)", "global")
        self.method = QComboBox()
        self.method.addItem("restrained (Cα–Cα distances kept)", "restrained")
        self.method.addItem("linear (chords; shows the artefact)", "linear")
        form.addRow("Start (loaded)", self.start)
        form.addRow("End", self.end)
        form.addRow("Superpose on", self.fit)
        form.addRow("Interpolation", self.method)
        lay.addLayout(form)

        row = QHBoxLayout()
        self.build = QPushButton("Build")
        self.build.clicked.connect(lambda: self.build_requested.emit(
            self.end.currentData() or "", self.fit.currentData(), self.method.currentData()))
        self._preset = (*PRESET, "Resting → activated")
        self.preset = QPushButton("")
        self.preset.clicked.connect(lambda: self.preset_requested.emit(*self._preset[:2]))
        self._set_preset("IP3R")
        row.addWidget(self.build)
        row.addWidget(self.preset)
        lay.addLayout(row)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)

        row = QHBoxLayout()
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self._slid)
        self.frame_label = QLabel("")
        self.frame_label.setMinimumWidth(
            self.frame_label.fontMetrics().horizontalAdvance("29/29 (start) gate 5.85 Å") + 8)
        self.play = QPushButton("Play")
        self.play.clicked.connect(self.play_requested.emit)
        self.stop = QPushButton("Stop")
        self.stop.clicked.connect(self.stop_requested.emit)
        for w in (self.slider, self.frame_label, self.play, self.stop):
            row.addWidget(w, 1 if w is self.slider else 0)
        lay.addLayout(row)
        self.paint = QCheckBox("Colour by displacement (fixed scale)")
        self.paint.toggled.connect(self.paint_toggled.emit)
        lay.addWidget(self.paint)

        row = QHBoxLayout()
        row.addWidget(QLabel("Lower plot"))
        self.lower = QComboBox()
        for key, label in LOWER_PLOTS.items():
            self.lower.addItem(label, key)
        self.lower.currentIndexChanged.connect(self._redraw)
        row.addWidget(self.lower, 1)
        lay.addLayout(row)
        self.canvas = PlotCanvas(self, height=5.6, rows=2)
        lay.addWidget(self.canvas, 1)
        self.report = QLabel("")
        self.report.setWordWrap(True)
        self.report.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.report)
        self.result = None
        self._marker = None
        self._following = False
        self.set_enabled_controls(False)

    # --------------------------------------------------------------- state

    def _set_preset(self, family: str) -> None:
        p = preset_for(family)
        self.preset.setVisible(p is not None)
        if p:
            self._preset = p
            self.preset.setText(f"{p[2]} ({p[0]} → {p[1]})")

    def set_start(self, pdb_id: str | None, paralog: str | None) -> None:
        """The loaded deposit; the end list is its paralog's other deposits."""
        self.clear()
        self.start.setText(pdb_id or "—")
        self.end.clear()
        if pdb_id and paralog:
            entries = [e for e in load_registry() if e.paralog == paralog]
            self._set_preset(entries[0].family if entries else "IP3R")
            for e in entries:
                if (e.human or e.family == "RyR") and e.pdb_id != pdb_id:
                    self.end.addItem(f"{e.pdb_id} — {e.state}", e.pdb_id)
            i = self.end.findData(self._preset[1])
            self.end.setCurrentIndex(max(i, 0))
        self.build.setEnabled(self.end.count() > 0)
        if pdb_id and not paralog:
            self.status.setText(f"{pdb_id} is in no human paralog numbering, so its "
                                "residues cannot be matched to another deposit.")

    def set_enabled_controls(self, on: bool) -> None:
        for w in (self.slider, self.play, self.stop, self.paint):
            w.setEnabled(on)

    def set_busy(self, text: str) -> None:
        self.status.setText(text)
        self.build.setEnabled(False)

    def failed(self, text: str) -> None:
        self.status.setText(f"not built: {text.splitlines()[0]}")
        self.build.setEnabled(True)

    def show_result(self, result) -> None:
        self.result = result
        tr, traj, ov = result.transition, result.trajectory, result.overlap
        self.build.setEnabled(True)
        corr = ("" if tr.meta["correspondence_determined"] else
                " Every cyclic subunit correspondence fits alike (a C4-symmetric "
                "pair); the deposited chain labels were kept.")
        self.status.setText(
            f"{tr.start_id} → {tr.end_id}: {len(tr.residues):,} residues × "
            f"{tr.n_subunits} subunits in the basis; {tr.fit} fit RMSD "
            f"{tr.fit_rmsd:.2f} Å, overall {tr.rmsd:.2f} Å. {traj.summary()}." + corr)
        self._following = True
        self.slider.setRange(0, len(traj) - 1)
        self.slider.setValue(0)
        self._following = False
        self._frame_text(0)
        self.set_enabled_controls(True)
        self.report.setText("<br>".join(ov.report()))
        self._plot(tr, ov, result.gate)

    def follow(self, i: int) -> None:
        """Move the slider with playback without requesting a frame."""
        self._following = True
        self.slider.setValue(i)
        self._following = False
        self._frame_text(i)

    def clear(self) -> None:
        self.result = None
        self.status.setText("")
        self.report.setText("")
        self.frame_label.setText("")
        self.paint.setChecked(False)
        self.set_enabled_controls(False)
        self._marker = None
        self.canvas.reset(2, 1)
        self.canvas.draw_idle()

    # ------------------------------------------------------------ internals

    def _frame_text(self, i: int) -> None:
        n = self.slider.maximum()
        gate = ""
        if self.result is not None and n:
            g = self.result.gate
            gate = f" gate {g.gate[i]:.2f} Å"
            if self._marker is not None:
                self._marker.set_xdata([g.fraction[i]])
                self.canvas.draw_idle()
        self.frame_label.setText(f"{i}/{n}" + (" (start)" if i == 0 else
                                               " (end)" if i == n else "") + gate)

    def _slid(self, i: int) -> None:
        self._frame_text(i)
        if not self._following:
            self.frame_requested.emit(i)

    def _redraw(self) -> None:
        r = self.result
        if r is not None:
            self._plot(r.transition, r.overlap, r.gate)

    def _plot(self, tr, ov, gate) -> None:
        """The A subspace on top, always; the lower plot as chosen."""
        top, ax = self.canvas.reset(2, 1)[:, 0]
        self._plot_subspace(top, tr, ov)
        self._marker = None
        which = self.lower.currentData()
        if which == "modes":
            self._plot_modes(ax, tr, ov)
        elif which == "elements":
            means = sorted(tr.element_means().items(), key=lambda kv: kv[1])
            ax.barh([ELEMENT_LABELS.get(k, k) for k, _ in means],
                    [v for _, v in means], color=PALETTE[0])
            ax.set_title(f"mean Cα displacement by element (Å, {tr.fit} fit)")
            ax.tick_params(axis="y", labelsize=7)
        else:
            self._plot_gate(ax, tr, gate)
        self.canvas.draw_now()

    def _plot_subspace(self, ax, tr, ov) -> None:
        """The headline: the collective A modes taken together, lowest
        first. The split among them moves with the ANM cutoff; the total
        does not (Round 2, ``network_checks.cutoff_scan``)."""
        sub = ov.subspace("A")
        if not len(sub.index):
            ax.text(0.5, 0.5, "no collective A modes", transform=ax.transAxes,
                    ha="center", color="grey")
            return
        k = np.arange(1, len(sub.index) + 1)
        ax.bar(k, ov.overlap[sub.index], color=_IRREP_COLOUR["A"], alpha=0.55,
               label="each A mode")
        ax.plot(k, sub.cumulative, "o-", color="#d7dbe3", lw=1.4, ms=3,
                label="A modes together")
        ax.plot(k, sub.null, color="#d7dbe3", lw=1, ls=":",
                label="random, same symmetry")
        ax.axhline(sub.ceiling, color=_IRREP_COLOUR["A"], lw=0.8, ls="--",
                   label=f"ceiling √(A share) {sub.ceiling:.2f}")
        ax.set_xticks(k)
        ax.set_xticklabels([f"#{i + 1}" for i in sub.index], fontsize=7)
        ax.set_ylim(0, 1.75)                       # legend above the ceiling
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel(f"collective A modes of {ov.reference}, lowest first")
        ax.set_ylabel("overlap")
        ax.set_title(f"A subspace: {sub.total:.2f} of the move "
                     f"(random {sub.null[-1]:.3f})")
        self.canvas.legend(ax, loc="upper left", ncol=2)

    def _plot_modes(self, ax, tr, ov) -> None:
        k = np.arange(1, len(ov.overlap) + 1)
        local = ~ov.modes.is_collective()
        colours = [_IRREP_COLOUR.get(s, "#8a8f99") for s in ov.modes.symmetry]
        bars = ax.bar(k, ov.overlap, color=colours)
        for b, lo in zip(bars, local):
            if lo:
                b.set_alpha(0.3)
                b.set_hatch("//")
        ax.plot(k, ov.cumulative, color="#d7dbe3", lw=1.2, label="cumulative")
        ax.plot(k, ov.null_cumulative, color="#d7dbe3", lw=1, ls=":",
                label="random, same symmetry")
        for name in ("A", "B", "E"):
            ax.bar([0], [0], color=_IRREP_COLOUR[name], label=f"{name} mode")
        ax.set_xlim(0.4, len(k) + 0.6)
        ax.set_xticks(k[::2] if len(k) > 12 else k)
        ax.set_ylim(0, 1.45)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel(f"ANM mode of {ov.reference} (hatched: local artefact)")
        ax.set_ylabel("overlap |cos|")
        ax.set_title("every mode, by symmetry")
        self.canvas.legend(ax, loc="upper left", ncol=3)

    def _plot_gate(self, ax, tr, g) -> None:
        ax.plot(g.fraction, g.gate, color=PALETTE[0], lw=1.6,
                label="gate, interpolated")
        ax.plot(g.fraction, g.rigid_gate, color=PALETTE[0], lw=1, ls="--",
                label="gate, rigid side chains")
        ax.plot(g.fraction, g.filter, color=PALETTE[2], lw=1.2, label="filter")
        self._marker = ax.axvline(g.fraction[self.slider.value()], color="#d7dbe3", lw=0.8)
        ax.set_xlim(0, 1)
        ax.set_xlabel(f"path fraction (gate half-way at {g.half_open():.2f})")
        ax.set_ylabel("radius r_min (Å)")
        ax.set_title("the pore along the morph")
        self.canvas.legend(ax, loc="lower right", ncol=1)
