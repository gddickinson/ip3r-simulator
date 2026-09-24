"""Left dock: choose a deposition, and how it is drawn.

The list is the curated registry (``resources/structures.json``), marked by
whether each file is already downloaded. Loading happens on a worker thread;
the panel emits :attr:`load_requested` and the main window does the rest.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QListWidget, QListWidgetItem,
                             QPushButton, QVBoxLayout, QWidget)

from ..core.annotations import ELEMENT_COLORS, ELEMENT_LABELS, ELEMENT_ORDER, LAYER_LABELS
from ..core.modules import MODULES_KEY
from ..io.loader import is_local
from ..parameters import PARAMETERS as _P
from ..io.registry import load_registry
from ..render.representations import COLOR_LABELS, STYLE_LABELS, ColorBy, Style

__all__ = ["StructurePanel"]


def _swatch(rgb) -> str:
    r, g, b = (int(255 * c) for c in rgb)
    return f"<span style='color: rgb({r},{g},{b})'>■</span>"


class StructurePanel(QWidget):
    load_requested = pyqtSignal(str)            # pdb id
    style_changed = pyqtSignal()
    sites_toggled = pyqtSignal(str, bool)       # site class, on

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)

        box = QGroupBox("Depositions (from ip3r_genes S0/S11/S22)")
        bl = QVBoxLayout(box)
        self.list = QListWidget()
        self.list.itemDoubleClicked.connect(lambda it: self._load(it))
        bl.addWidget(self.list)
        row = QHBoxLayout()
        self.load_btn = QPushButton("Load")
        self.load_btn.clicked.connect(lambda: self._load(self.list.currentItem()))
        row.addWidget(self.load_btn)
        row.addStretch(1)
        bl.addLayout(row)
        self.info = QLabel("")
        self.info.setWordWrap(True)
        self.info.setTextFormat(Qt.TextFormat.RichText)
        bl.addWidget(self.info)
        lay.addWidget(box)

        box = QGroupBox("Representation")
        form = QFormLayout(box)
        self.style = QComboBox()
        for s in Style:
            self.style.addItem(STYLE_LABELS[s], s)
        self.color = QComboBox()
        for c in ColorBy:
            if c is not ColorBy.VALUE:
                self.color.addItem(COLOR_LABELS[c], c)
        self.layer = QComboBox()
        for key, label in LAYER_LABELS.items():
            self.layer.addItem(label, key)
        self.layer.setEnabled(False)
        self.ligands = QCheckBox("Show ligands (IP3 enlarged)")
        self.ligands.setChecked(True)
        form.addRow("Style", self.style)
        form.addRow("Colour", self.color)
        form.addRow("Layer", self.layer)
        form.addRow(self.ligands)
        for w in (self.style, self.color, self.layer):
            w.currentIndexChanged.connect(self._restyle)
        self.ligands.toggled.connect(lambda _: self.style_changed.emit())
        lay.addWidget(box)

        box = QGroupBox("Subunits")
        self.chain_row = QHBoxLayout(box)
        self.chain_boxes: dict[str, QCheckBox] = {}
        lay.addWidget(box)

        box = QGroupBox("Measured sites (human numbering only)")
        sl = QVBoxLayout(box)
        self.site_boxes = {}
        for key, label in (("ip3_contact", "Ten IP3 contacts (S0)"),
                           ("filter_lining", "Filter lining (S0)"),
                           ("gate_lining", "Gate lining (S0)"),
                           (MODULES_KEY, "Paper 6 modules: ligand core (green) "
                                         "and pore less the loop (magenta), Cα")):
            cb = QCheckBox(label)
            cb.toggled.connect(lambda on, k=key: self.sites_toggled.emit(k, on))
            sl.addWidget(cb)
            self.site_boxes[key] = cb
        lay.addWidget(box)

        self.legend = QLabel()
        self.legend.setTextFormat(Qt.TextFormat.RichText)
        self.legend.setWordWrap(True)
        lay.addWidget(self.legend)
        lay.addStretch(1)
        self.refresh_list()
        self._update_legend()

    # ------------------------------------------------------------ registry

    def refresh_list(self) -> None:
        current = self.current_id()
        self.list.clear()
        for e in load_registry():
            mark = "●" if is_local(e.pdb_id) else "○"
            item = QListWidgetItem(f"{mark} {e.label}")
            item.setData(Qt.ItemDataRole.UserRole, e.pdb_id)
            item.setToolTip(f"{e.title}\nroles: {', '.join(e.roles)}"
                            f"\n{'IP3 bound' if e.ip3_bound else 'no IP3'}"
                            f"\n{'downloaded' if is_local(e.pdb_id) else 'will be fetched from RCSB'}")
            self.list.addItem(item)
            if e.pdb_id == current:
                self.list.setCurrentItem(item)

    def current_id(self) -> str | None:
        it = self.list.currentItem()
        return None if it is None else it.data(Qt.ItemDataRole.UserRole)

    def select(self, pdb_id: str) -> None:
        for i in range(self.list.count()):
            it = self.list.item(i)
            if it.data(Qt.ItemDataRole.UserRole) == pdb_id.upper():
                self.list.setCurrentItem(it)
                self._load(it)
                return

    def _load(self, item) -> None:
        if item is not None:
            self.load_requested.emit(item.data(Qt.ItemDataRole.UserRole))

    # -------------------------------------------------------------- styling

    def current_style(self) -> Style:
        return self.style.currentData()

    def current_color(self) -> ColorBy:
        return self.color.currentData()

    def current_layer(self) -> str:
        return self.layer.currentData()

    def visible_chains(self) -> frozenset | None:
        if not self.chain_boxes:
            return None
        return frozenset(c for c, b in self.chain_boxes.items() if b.isChecked())

    def set_chains(self, chains: list[str]) -> None:
        for b in self.chain_boxes.values():
            b.deleteLater()
        self.chain_boxes = {}
        for c in chains:
            b = QCheckBox(c)
            b.setChecked(True)
            b.toggled.connect(lambda _: self.style_changed.emit())
            self.chain_row.addWidget(b)
            self.chain_boxes[c] = b

    def set_info(self, html: str) -> None:
        self.info.setText(html)

    def _restyle(self) -> None:
        self.layer.setEnabled(self.current_color() is ColorBy.CONSERVATION)
        self._update_legend()
        self.style_changed.emit()

    def _update_legend(self) -> None:
        if self.current_color() is ColorBy.ELEMENT_DOMAIN:
            rows = [f"{_swatch(ELEMENT_COLORS[k])} {ELEMENT_LABELS[k]}"
                    for k in ELEMENT_ORDER]
            self.legend.setText("<br>".join(rows))
        elif self.current_color() is ColorBy.CONSERVATION:
            self.legend.setText(
                "JSD on a fixed 0.50-0.95 scale: <span style='color:#3b4cc0'>■</span> "
                "variable → <span style='color:#b40426'>■</span> conserved; "
                "<span style='color:#6b6d75'>■</span> not scored (never the low end).")
        elif self.current_color() is ColorBy.DISPLACEMENT:
            top = _P.value("display.displacement_max")
            self.legend.setText(
                f"C-alpha displacement between two states on a fixed 0-{top:g} Å "
                "scale: <span style='color:#3b4cc0'>■</span> still → "
                "<span style='color:#b40426'>■</span> moves most; "
                "<span style='color:#6b6d75'>■</span> not measured (outside the "
                "common basis, or no transition built — use the Transition tab).")
        else:
            self.legend.setText("")
