"""The modes tab: elastic-network normal modes of the loaded tetramer.

Computes the ANM on a worker thread, lists the modes with their C4 irrep,
and animates the selected one in the viewport. The animation is the mode's
displacement field added to the coordinates with a sinusoidal amplitude —
**an illustration of a direction of motion**, not a trajectory: an ANM
gives shapes and relative stiffnesses, not amplitudes or time scales, and
the panel says so.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (QHBoxLayout, QLabel, QPushButton, QSlider,
                             QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)
from PyQt6.QtCore import Qt

__all__ = ["ModesPanel", "IRREP_TEXT"]

IRREP_TEXT = {
    "A": "A — all four subunits move alike; can couple to IP3 binding at all "
         "four sites and to a symmetric pore opening",
    "B": "B — neighbouring subunits move in opposition",
    "E": "E — degenerate pair; tilts or shears the tetramer, cannot open a "
         "four-fold pore at first order",
    "mixed": "not cleanly one irrep (deposit asymmetry)",
}


class ModesPanel(QWidget):
    compute_requested = pyqtSignal()
    animate_requested = pyqtSignal(int, float)     # mode index, amplitude Å
    stop_requested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        note = QLabel("Anisotropic network model over C-alpha atoms (every "
                      "residue all four subunits resolve, strided). Modes are "
                      "directions of collective motion with relative "
                      "stiffness; amplitudes and time scales are illustrative. "
                      "κ is the fraction of sites a mode moves (Brüschweiler "
                      "1995); low-κ modes are network artefacts.")
        note.setWordWrap(True)
        lay.addWidget(note)
        row = QHBoxLayout()
        self.compute = QPushButton("Compute modes")
        self.compute.clicked.connect(self.compute_requested.emit)
        row.addWidget(self.compute)
        self.status = QLabel("")
        row.addWidget(self.status, 1)
        lay.addLayout(row)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["#", "eigenvalue", "irrep", "χ(C4)", "κ"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.itemSelectionChanged.connect(self._selected)
        lay.addWidget(self.table, 1)
        self.explain = QLabel("")
        self.explain.setWordWrap(True)
        lay.addWidget(self.explain)
        row = QHBoxLayout()
        row.addWidget(QLabel("Amplitude"))
        self.amp = QSlider(Qt.Orientation.Horizontal)
        self.amp.setRange(1, 40)
        self.amp.setValue(12)
        self.amp.valueChanged.connect(self._selected)
        row.addWidget(self.amp, 1)
        self.amp_label = QLabel("12 Å")
        row.addWidget(self.amp_label)
        self.stop = QPushButton("Stop")
        self.stop.clicked.connect(self.stop_requested.emit)
        row.addWidget(self.stop)
        lay.addLayout(row)
        self.modes = None
        self.local = {}

    def set_busy(self, text: str) -> None:
        self.status.setText(text)
        self.compute.setEnabled(False)

    def show_modes(self, modes, meta: str, local: dict | None = None) -> None:
        """``local`` maps a local mode's index to where it sits ("on residue 86, ...")."""
        self.modes = modes
        self.local = local or {}
        self.compute.setEnabled(True)
        a = modes.first("A")
        where = sorted(set(self.local.values()))
        self.status.setText(meta + (f"; lowest collective A mode is #{a + 1}"
                                    if a is not None else "")
                            + (f"; local artefacts {'; '.join(where)}" if where else ""))
        kappa = modes.collectivity()
        self.table.setRowCount(modes.n_modes)
        for i in range(modes.n_modes):
            for j, text in enumerate((str(i + 1), f"{modes.eigenvalues[i]:.3e}",
                                      modes.symmetry[i], f"{modes.character[i]:+.3f}",
                                      f"{kappa[i]:.2f}")):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()

    def _selected(self) -> None:
        self.amp_label.setText(f"{self.amp.value()} Å")
        rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        if not rows or self.modes is None:
            return
        i = rows[0].row()
        local = "" if self.modes.is_collective()[i] else (
            " Collectivity κ is below threshold: a weakly attached fragment "
            f"of the network ({self.local.get(i, 'unlocated')}), not a collective "
            "motion.")
        self.explain.setText(f"Mode {i + 1}: {IRREP_TEXT.get(self.modes.symmetry[i], '')}."
                             + local)
        self.animate_requested.emit(i, float(self.amp.value()))

    def select(self, index: int, amplitude: float) -> bool:
        """Animate mode ``index`` at ``amplitude`` Å, as a click would (a
        session restore); False if the modes are not computed or too few."""
        if self.modes is None or not 0 <= index < self.modes.n_modes:
            return False
        for w in (self.amp, self.table):
            w.blockSignals(True)
        self.amp.setValue(int(round(amplitude)))
        self.table.selectRow(index)
        for w in (self.amp, self.table):
            w.blockSignals(False)
        self._selected()                          # one animation, however reached
        return True

    def clear(self) -> None:
        self.modes = None
        self.table.setRowCount(0)
        self.status.setText("")
        self.explain.setText("")
        self.compute.setEnabled(True)


def mode_frame(base: np.ndarray, disp: np.ndarray, phase: float) -> np.ndarray:
    """Coordinates at ``phase`` (radians) of a sinusoidal mode animation."""
    return base + np.sin(phase) * disp
