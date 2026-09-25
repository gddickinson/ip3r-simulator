"""Left dock: choose a deposition, and how it is drawn.

The list is the curated registry (``resources/structures.json`` and
``ryr1.json``), grouped by family with RyR1 collapsed until one of its
deposits is chosen, and marked by whether each file is already downloaded. Loading happens on a worker thread;
the panel emits :attr:`load_requested` and the main window does the rest.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QPushButton, QTreeWidget,
                             QTreeWidgetItem, QVBoxLayout, QWidget)

from ..core.annotations import ELEMENT_COLORS, ELEMENT_LABELS, ELEMENT_ORDER, LAYER_LABELS
from ..core.modules import MODULES_KEY
from ..io.loader import is_local
from ..parameters import PARAMETERS as _P
from ..io.registry import load_registry
from ..render.colormaps import PLDDT_COLORS, SEAM_COLORS, SHELL_COLORS
from ..render.representations import COLOR_LABELS, STYLE_LABELS, ColorBy, Style
from ..structure.graft import FILL_MODES
from ..structure.shells import SHELLS
from .view_state import ParameterFollower

__all__ = ["StructurePanel", "FAMILY_LABELS"]

#: Headings of the deposition list, by registry family.
FAMILY_LABELS = {"IP3R": "IP3 receptors", "RyR": "RyR1 (rabbit)"}
_FAMILY_ROLE = Qt.ItemDataRole.UserRole + 1


def _swatch(rgb) -> str:
    r, g, b = (int(255 * c) for c in rgb)
    return f"<span style='color: rgb({r},{g},{b})'>■</span>"


class StructurePanel(QWidget):
    load_requested = pyqtSignal(str)            # pdb id
    style_changed = pyqtSignal()
    sites_toggled = pyqtSignal(str, bool)       # site class, on
    completeness_changed = pyqtSignal(str)      # a FILL_MODES key

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)

        box = QGroupBox("Depositions")
        box.setToolTip("IP3R: from ip3r_genes S0/S11/S22; RyR1: curated here")
        bl = QVBoxLayout(box)
        self.list = QTreeWidget()
        self.list.setHeaderHidden(True)
        self.list.setRootIsDecorated(True)
        self.list.itemDoubleClicked.connect(lambda it, _col: self._load(it))
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
            if c not in (ColorBy.VALUE, ColorBy.PLDDT):
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
        self.completeness = QComboBox()
        for key, label, tip in FILL_MODES:
            self.completeness.addItem(label, key)
            self.completeness.setItemData(self.completeness.count() - 1, tip,
                                          Qt.ItemDataRole.ToolTipRole)
        self.completeness.currentIndexChanged.connect(self._completeness)
        form.addRow("Completeness", self.completeness)
        self.fill_info = QLabel("")
        self.fill_info.setWordWrap(True)
        self.fill_info.setTextFormat(Qt.TextFormat.RichText)
        form.addRow(self.fill_info)
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
                           (MODULES_KEY, "Paper 6 modules (Cα)")):
            cb = QCheckBox(label)
            if key == MODULES_KEY:
                cb.setToolTip("Ligand core (green) and pore less the luminal loop (magenta)")
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
        self.follower = ParameterFollower(self)
        self.follower.changed.connect(self._parameters_changed)

    # ------------------------------------------------------------ registry

    def refresh_list(self) -> None:
        current = self.current_id()
        expanded = {self.list.topLevelItem(i).data(0, _FAMILY_ROLE)
                    for i in range(self.list.topLevelItemCount())
                    if self.list.topLevelItem(i).isExpanded()}
        first = not self.list.topLevelItemCount()
        self.list.clear()
        groups: dict[str, QTreeWidgetItem] = {}
        for e in load_registry():
            fam = e.family
            if fam not in groups:
                groups[fam] = QTreeWidgetItem(self.list, [fam])
                groups[fam].setFlags(Qt.ItemFlag.ItemIsEnabled)   # a heading, not a deposit
                groups[fam].setData(0, _FAMILY_ROLE, fam)
            mark = "●" if is_local(e.pdb_id) else "○"
            item = QTreeWidgetItem(groups[fam], [f"{mark} {e.label}"])
            item.setData(0, Qt.ItemDataRole.UserRole, e.pdb_id)
            item.setToolTip(0, f"{e.title}\nroles: {', '.join(e.roles)}"
                               f"\n{'IP3 bound' if e.ip3_bound else 'no IP3'}"
                               f"\n{'downloaded' if is_local(e.pdb_id) else 'will be fetched from RCSB'}")
            if e.pdb_id == current:
                self.list.setCurrentItem(item)
        for fam, g in groups.items():
            g.setText(0, f"{FAMILY_LABELS.get(fam, fam)} ({g.childCount()})")
            # IP3R open on first build, RyR1 collapsed; afterwards as the user left them.
            g.setExpanded(fam == "IP3R" if first else fam in expanded)
        if current is not None:
            it = self._item(current)
            if it is not None:
                it.parent().setExpanded(True)

    def _item(self, pdb_id: str):
        """The list item of a deposit, or None."""
        for i in range(self.list.topLevelItemCount()):
            g = self.list.topLevelItem(i)
            for j in range(g.childCount()):
                if g.child(j).data(0, Qt.ItemDataRole.UserRole) == pdb_id.upper():
                    return g.child(j)
        return None

    def ids(self) -> list[str]:
        """Every deposit in the list, in order."""
        return [self.list.topLevelItem(i).child(j).data(0, Qt.ItemDataRole.UserRole)
                for i in range(self.list.topLevelItemCount())
                for j in range(self.list.topLevelItem(i).childCount())]

    def current_id(self) -> str | None:
        it = self.list.currentItem()
        return None if it is None else it.data(0, Qt.ItemDataRole.UserRole)

    def select(self, pdb_id: str) -> None:
        it = self._item(pdb_id)
        if it is not None:
            it.parent().setExpanded(True)
            self.list.setCurrentItem(it)
            self._load(it)

    def _load(self, item) -> None:
        pdb_id = None if item is None else item.data(0, Qt.ItemDataRole.UserRole)
        if pdb_id:                                  # a family heading has none
            self.load_requested.emit(pdb_id)

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

    # --------------------------------------------------------- completeness

    def current_completeness(self) -> str:
        return self.completeness.currentData()

    def set_completeness(self, key: str) -> None:
        i = self.completeness.findData(key)
        if i >= 0:
            self.completeness.setCurrentIndex(i)

    def _completeness(self) -> None:
        self._update_legend()
        self.completeness_changed.emit(self.current_completeness())

    def set_fill_info(self, html: str) -> None:
        self.fill_info.setText(html)

    def _fill_legend(self) -> str:
        if self.current_completeness() == "none":
            return ""
        edges = [_P.value(k) for k in ("display.plddt_low", "graft.plddt_confident",
                                       "display.plddt_very_high")]
        names = [f"&lt; {edges[0]:g}", f"{edges[0]:g}-{edges[1]:g}",
                 f"{edges[1]:g}-{edges[2]:g}", f"&ge; {edges[2]:g}"]
        rows = " ".join(f"{_swatch(c)} {n}" for c, n in zip(PLDDT_COLORS, names))
        tol = _P.value("graft.join_tolerance")
        return (f"<br><b>AlphaFold fill</b>, pLDDT: {rows}<br>"
                f"{_swatch(SEAM_COLORS[True])} seam closes (≤ {tol:g} Å) "
                f"{_swatch(SEAM_COLORS[False])} seam broken")

    def _parameters_changed(self) -> None:
        """The colour scales are registered values: when one the drawing
        uses is edited, legend and colours move together (``_restyle``).
        Otherwise only the legend is re-read."""
        scaled = self.current_color() in (ColorBy.DISPLACEMENT, ColorBy.LIGAND_SHELL)
        if scaled or self.current_completeness() != "none":
            self._restyle()
        else:
            self._update_legend()

    def _restyle(self) -> None:
        self.layer.setEnabled(self.current_color() is ColorBy.CONSERVATION)
        self._update_legend()
        self.style_changed.emit()

    def _update_legend(self) -> None:
        self._colour_legend()
        self.legend.setText(self.legend.text() + self._fill_legend())

    def _colour_legend(self) -> None:
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
        elif self.current_color() is ColorBy.LIGAND_SHELL:
            edges = [_P.value(k) for k in ("ligand.contact_cutoff",
                                           "ligand.shell_second_edge",
                                           "ligand.shell_third_edge",
                                           "ligand.shell_radius")]
            lo, rows = 0.0, []
            for rgb, name, hi in zip(SHELL_COLORS, SHELLS, edges):
                rows.append(f"{_swatch(rgb)} {name} ({lo:g}-{hi:g} Å)")
                lo = hi
            self.legend.setText(
                "All-atom distance of each residue to the IP3 on its own "
                "subunit:<br>" + "<br>".join(rows)
                + f"<br><span style='color:#6b6d75'>■</span> beyond {lo:g} Å, "
                "or no IP3 bound on that subunit.")
        else:
            self.legend.setText("")
