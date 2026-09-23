"""Help → Parameters: every registered number with its provenance.

Read-only by design in this version: the findings checks refuse to run
against a modified registry, so editing a value belongs with a visible
"parameters modified" banner, which is a roadmap item.
"""

from __future__ import annotations

from PyQt6.QtWidgets import QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout

from ..parameters import PARAMETERS

__all__ = ["ParametersDialog"]


class ParametersDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Registered parameters")
        self.resize(1100, 600)
        lay = QVBoxLayout(self)
        rows = PARAMETERS.provenance_rows()
        cols = ["key", "value", "unit", "kind", "citation", "source_note"]
        t = QTableWidget(len(rows), len(cols))
        t.setHorizontalHeaderLabels(cols)
        for i, r in enumerate(rows):
            for j, c in enumerate(cols):
                t.setItem(i, j, QTableWidgetItem(f"{r[c]:g}" if c == "value" else str(r[c])))
        t.resizeColumnsToContents()
        lay.addWidget(t)
