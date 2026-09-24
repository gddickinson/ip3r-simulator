"""The variants tab: the S17 variant harvest, placed on the structure.

A variant is keyed on its paralog's canonical human residue number. It can
be drawn on the loaded deposit only if that deposit is in the same
paralog's numbering (the S24 rule) — otherwise the panel says why nothing is
drawn rather than lighting up whatever residue carries that number.

"Draw on structure" puts a sphere on every variant residue of the chosen
class on all four subunits. With a layer chosen under "VUS by layer", each
VUS is placed against its gene's labelled medians on that layer (Paper 5
§8): the table gains a stratum column, the plot shows the three classes and
both medians, and the VUS spheres take the stratum's colour.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel,
                             QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from ..analysis.vus_figure import draw_strip
from ..analysis.vus_strata import (CLASS_COLORS, STRATUM_COLORS, STRATUM_LABELS,
                                   stratum_of)
from ..config import PARALOGS
from ..core.annotations import LAYERS, constraint_at, element_of, variants
from ..render.variant_spheres import CLASS_ORDER, resource_stratification
from .plot_canvas import PlotCanvas

__all__ = ["VariantsPanel"]

_BUCKETS = CLASS_ORDER
_COLUMNS = ["residue", "change", "class", "element", "deep JSD", "family JSD",
            "stratum", "condition"]


def _swatch(rgb) -> str:
    r, g, b = (int(255 * c) for c in rgb)
    return f'<span style="color:#{r:02x}{g:02x}{b:02x}">●</span>'


class VariantsPanel(QWidget):
    highlight_residue = pyqtSignal(str, int)            # paralog, residue
    draw_requested = pyqtSignal(str, tuple, object)     # paralog, classes, layer | None

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.paralog = QComboBox()
        self.paralog.addItems(PARALOGS)
        self.bucket = QComboBox()
        self.bucket.addItems(["all"] + list(_BUCKETS))
        self.bucket.setCurrentText("P/LP")
        self.layer = QComboBox()
        self.layer.addItem("none", None)
        for la in LAYERS:
            self.layer.addItem(la, la)
        self.draw = QCheckBox("Draw on structure")
        for w in (self.paralog, self.bucket, self.layer):
            w.currentIndexChanged.connect(self.refresh)
        self.draw.toggled.connect(lambda *_: self.emit_draw())
        for text, w in (("Paralog", self.paralog), ("Class", self.bucket),
                        ("VUS by layer", self.layer)):
            row.addWidget(QLabel(text))
            row.addWidget(w)
        row.addWidget(self.draw)
        row.addStretch(1)
        lay.addLayout(row)
        self.legend = QLabel("")
        self.legend.setWordWrap(True)
        lay.addWidget(self.legend)
        self.table = QTableWidget(0, len(_COLUMNS))
        self.table.setHorizontalHeaderLabels(_COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.itemSelectionChanged.connect(self._selected)
        lay.addWidget(self.table, 1)
        self.canvas = PlotCanvas(self, height=2.8)
        lay.addWidget(self.canvas)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)
        self.strat = None
        self.refresh()

    def classes(self) -> tuple[str, ...]:
        want = self.bucket.currentText()
        return tuple(_BUCKETS) if want == "all" else (want,)

    def refresh(self) -> None:
        gene = self.paralog.currentText()
        want = self.bucket.currentText()
        layer = self.layer.currentData()
        self.strat = resource_stratification(gene, layer) if layer else None
        rows = [v for v in variants(gene) if want == "all" or v["class_bucket"] == want]
        resi = np.array([v["resi"] for v in rows], int)
        deep = constraint_at(gene, resi, "deep") if len(rows) else []
        fam = constraint_at(gene, resi, "family") if len(rows) else []
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for i, v in enumerate(rows):
            vals = [v["resi"], f"{v['ref_aa']}{v['resi']}{v['alt_aa']}", v["class_bucket"],
                    element_of(gene, v["resi"]) or "", deep[i], fam[i],
                    self._stratum(v), v["condition"]]
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
        self._legend()
        self._plot()
        self.status.setText(f"{len(rows)} variant(s) in {gene} "
                            f"({'all classes' if want == 'all' else want}), "
                            "from the ip3r_genes S17 harvest (ClinVar + curated UniProt).")
        self.emit_draw()

    def _stratum(self, v: dict) -> str:
        if self.strat is None or v["class_bucket"] != "VUS":
            return ""
        s = self.strat
        return stratum_of(s.vus.get(int(v["resi"]), np.nan), s.median_pathogenic,
                          s.median_benign)

    def _legend(self) -> None:
        parts = [f"{_swatch(CLASS_COLORS[c])} {c}" for c in self.classes()
                 if not (c == "VUS" and self.strat is not None)]
        if self.strat is not None and "VUS" in self.classes():
            parts += [f"{_swatch(STRATUM_COLORS[k])} {STRATUM_LABELS[k]}"
                      for k in ("pathogenic-like", "between", "benign-like", "not scored")]
        self.legend.setText("&nbsp;&nbsp;".join(parts))

    def _plot(self) -> None:
        ax = self.canvas.reset()
        if self.strat is None:
            text = ("choose a layer under “VUS by layer” to place each VUS\n"
                    "against its gene's labelled medians (Paper 5 §8)"
                    if self.layer.currentData() is None else
                    "no stratification: a labelled class has no scored position")
            ax.text(0.5, 0.5, text, ha="center", va="center", color="#8a8f99",
                    fontsize=8, transform=ax.transAxes)
            ax.set_axis_off()
        else:
            draw_strip(ax, self.strat)
        self.canvas.draw_now()

    def emit_draw(self) -> None:
        self.draw_requested.emit(self.paralog.currentText(),
                                 self.classes() if self.draw.isChecked() else (),
                                 self.layer.currentData())

    def follow(self, paralog: str | None) -> None:
        """A deposit loaded: show its paralog's variants, and redraw them."""
        if paralog in PARALOGS and paralog != self.paralog.currentText():
            self.paralog.setCurrentText(paralog)          # refresh() redraws
        else:
            self.emit_draw()

    def _selected(self) -> None:
        items = self.table.selectedItems()
        if items:
            resi = int(self.table.item(items[0].row(), 0).data(Qt.ItemDataRole.DisplayRole))
            self.highlight_residue.emit(self.paralog.currentText(), resi)

    def report(self, text: str) -> None:
        self.status.setText(text)
