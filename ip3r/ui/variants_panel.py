"""The variants tab: the S17 variant harvest, placed on the structure.

A variant is keyed on its paralog's canonical human residue number. It can
be highlighted on the loaded deposit only if that deposit is in the same
paralog's numbering (the S24 rule) — otherwise the panel says why nothing is
highlighted rather than lighting up whatever residue carries that number.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from ..config import PARALOGS
from ..core.annotations import constraint_at, element_of, variants

__all__ = ["VariantsPanel"]

_BUCKETS = ("P/LP", "B/LB", "VUS", "conflicting", "other")


class VariantsPanel(QWidget):
    highlight_residue = pyqtSignal(str, int)      # paralog, residue

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.paralog = QComboBox()
        self.paralog.addItems(PARALOGS)
        self.bucket = QComboBox()
        self.bucket.addItems(["all"] + list(_BUCKETS))
        self.bucket.setCurrentText("P/LP")
        for w in (self.paralog, self.bucket):
            w.currentIndexChanged.connect(self.refresh)
        row.addWidget(QLabel("Paralog"))
        row.addWidget(self.paralog)
        row.addWidget(QLabel("Class"))
        row.addWidget(self.bucket)
        lay.addLayout(row)
        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(
            ["residue", "change", "class", "element", "deep JSD", "family JSD", "condition"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.itemSelectionChanged.connect(self._selected)
        lay.addWidget(self.table, 1)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)
        self.refresh()

    def refresh(self) -> None:
        gene = self.paralog.currentText()
        want = self.bucket.currentText()
        rows = [v for v in variants(gene) if want == "all" or v["class_bucket"] == want]
        resi = np.array([v["resi"] for v in rows], int)
        deep = constraint_at(gene, resi, "deep") if len(rows) else []
        fam = constraint_at(gene, resi, "family") if len(rows) else []
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for i, v in enumerate(rows):
            vals = [v["resi"], f"{v['ref_aa']}{v['resi']}{v['alt_aa']}", v["class_bucket"],
                    element_of(gene, v["resi"]) or "", deep[i], fam[i], v["condition"]]
            for j, x in enumerate(vals):
                item = QTableWidgetItem()
                if isinstance(x, (int, float, np.floating)):
                    item.setData(Qt.ItemDataRole.DisplayRole,
                                 float(x) if np.isfinite(x) else float("nan"))
                else:
                    item.setText(str(x))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(i, j, item)
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        self.status.setText(f"{len(rows)} variant(s) in {gene} "
                            f"({'all classes' if want == 'all' else want}), "
                            "from the ip3r_genes S17 harvest (ClinVar + curated UniProt).")

    def _selected(self) -> None:
        items = self.table.selectedItems()
        if items:
            resi = int(self.table.item(items[0].row(), 0).data(Qt.ItemDataRole.DisplayRole))
            self.highlight_residue.emit(self.paralog.currentText(), resi)

    def report(self, text: str) -> None:
        self.status.setText(text)
