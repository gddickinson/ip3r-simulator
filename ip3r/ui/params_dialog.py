"""Help → Parameters: every registered number, editable, with its provenance.

One row per parameter, grouped by category: current value, documented
default, unit, bounds and citation; the full reference is on the tooltip.
Ported from the PIEZO1 editor, with import/export added.

**An edit is never silent.** Overridden rows are amber and bold, the main
window shows a banner for as long as anything differs from its default
(:mod:`ip3r.ui.params_banner`), and the findings checks refuse to confirm
against a modified registry. Values outside the declared bounds are clamped
and the clamp is reported.

**Nothing typed here persists.** The GUI does not remember overrides between
sessions: a value typed once must not silently change what a later run
computes. "Export…" writes the set in the ``IP3R_PARAMETERS`` format, which
is the one route to reproducing it (headless or in the GUI, via "Import…").
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (QAbstractItemView, QCheckBox, QDialog, QFileDialog,
                             QHBoxLayout, QLabel, QLineEdit, QPushButton,
                             QTreeWidget, QTreeWidgetItem, QVBoxLayout)

from ..parameters import PARAMETERS, references

__all__ = ["ParametersDialog", "COLUMNS", "OVERRIDDEN"]

COLUMNS = ["Parameter", "Value", "Default", "Unit", "Bounds", "Kind", "Source"]
VALUE = COLUMNS.index("Value")
OVERRIDDEN = QColor("#f2a65a")
CITED = QColor("#7ed67e")
UNCITED = QColor("#8a919e")
KEY = Qt.ItemDataRole.UserRole


def _bounds(p) -> str:
    lo = "" if p.minimum is None else f"{p.minimum:g}"
    hi = "" if p.maximum is None else f"{p.maximum:g}"
    return f"{lo} – {hi}" if lo or hi else ""


class ParametersDialog(QDialog):
    """Browse and edit the registry."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Parameters — values, defaults and sources")
        self.resize(1180, 680)
        self.references = references()
        self._items: dict[str, QTreeWidgetItem] = {}

        lay = QVBoxLayout(self)
        intro = QLabel(
            "Every number the calculations use. Double-click a <b>Value</b> to "
            "change it; values outside the bounds are clamped. "
            "<span style='color:#f2a65a'>Amber</span> marks a value that differs "
            "from its documented default: the findings checks will not confirm "
            "while any does, and nothing typed here is remembered after you quit.")
        intro.setWordWrap(True)
        lay.addWidget(intro)

        row = QHBoxLayout()
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Filter by key, name, unit, citation or category…")
        self.filter_edit.textChanged.connect(self.apply_filter)
        row.addWidget(self.filter_edit, 1)
        self.modified_only = QCheckBox("Show only modified")
        self.modified_only.toggled.connect(self.apply_filter)
        row.addWidget(self.modified_only)
        lay.addLayout(row)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(len(COLUMNS))
        self.tree.setHeaderLabels(COLUMNS)
        self.tree.setAlternatingRowColors(True)
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked
                                  | QAbstractItemView.EditTrigger.SelectedClicked
                                  | QAbstractItemView.EditTrigger.EditKeyPressed)
        self.tree.itemChanged.connect(self._on_edited)
        lay.addWidget(self.tree, 1)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)

        buttons = QHBoxLayout()
        for text, slot in (("Reset selected", self.reset_selected),
                           ("Reset all to documented defaults", self.reset_all),
                           ("Import…", self._import), ("Export…", self._export)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            buttons.addWidget(b)
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.setDefault(True)
        close.clicked.connect(self.accept)
        buttons.addWidget(close)
        lay.addLayout(buttons)
        self._populate()

    # ------------------------------------------------------------- building

    def _populate(self) -> None:
        self.tree.blockSignals(True)
        self.tree.clear()
        self._items.clear()
        for category in PARAMETERS.categories():
            parent = QTreeWidgetItem([category])
            parent.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self.tree.addTopLevelItem(parent)
            parent.setFirstColumnSpanned(True)
            for p in PARAMETERS.in_category(category):
                item = self._row(p)
                parent.addChild(item)
                self._items[p.key] = item
            parent.setExpanded(True)
        for column in range(len(COLUMNS)):
            self.tree.resizeColumnToContents(column)
        self.tree.blockSignals(False)
        self.apply_filter()
        self._refresh_status()

    def _row(self, p) -> QTreeWidgetItem:
        item = QTreeWidgetItem([p.name, f"{PARAMETERS.value(p.key):g}", f"{p.default:g}",
                                p.unit, _bounds(p), p.kind, p.citation])
        item.setData(0, KEY, p.key)
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
        source = self.references.get(p.citation, p.citation)
        tip = (f"{p.key}\n{p.description}\nkind: {p.kind}"
               + (f"\nallowed: {_bounds(p)}" if _bounds(p) else "")
               + f"\n\nSource: {source}"
               + (f"\n{p.source_note}" if p.source_note else ""))
        for column in range(len(COLUMNS)):
            item.setToolTip(column, tip)
        self._paint(item, p)
        return item

    def _paint(self, item: QTreeWidgetItem, p) -> None:
        overridden = not PARAMETERS.is_default(p.key)
        colour = OVERRIDDEN if overridden else (CITED if p.cited else UNCITED)
        font = item.font(0)
        font.setBold(overridden)
        for column in range(len(COLUMNS)):
            item.setForeground(column, QBrush(colour))
            item.setFont(column, font)

    def _rewrite(self, key: str) -> None:
        item = self._items[key]
        self.tree.blockSignals(True)
        item.setText(VALUE, f"{PARAMETERS.value(key):g}")
        self._paint(item, PARAMETERS.get(key))
        self.tree.blockSignals(False)

    # -------------------------------------------------------------- editing

    def edit(self, key: str, text: str) -> bool:
        """Type ``text`` into ``key``'s Value cell (the smoke test's route);
        returns whether it was accepted as a number."""
        self._items[key].setText(VALUE, text)            # routes through _on_edited
        return self._last_ok

    _last_ok = True

    def _on_edited(self, item: QTreeWidgetItem, column: int) -> None:
        key = item.data(0, KEY)
        if column != VALUE or not key:
            return
        p = PARAMETERS.get(key)
        text = item.text(VALUE).strip()
        try:
            requested = float(text)
        except ValueError:
            requested = float("nan")
        if requested != requested:                       # not a number (or "nan")
            self._last_ok = False
            self._rewrite(key)
            return self._say(f"<span style='color:#f26d6d'>{text!r} is not a number; "
                             f"{p.name} left at {PARAMETERS.value(key):g}</span>")
        self._last_ok = True
        applied = PARAMETERS.set_value(key, requested)
        self._rewrite(key)
        if applied != requested:
            self._say(f"<span style='color:#f2a65a'>{requested:g} is outside "
                      f"{_bounds(p)} for {p.name}; clamped to {applied:g} {p.unit}</span>")
        else:
            self._refresh_status()

    def reset_selected(self) -> None:
        for item in self.tree.selectedItems():
            key = item.data(0, KEY)
            if key:
                PARAMETERS.reset(key)
                self._rewrite(key)
        self._refresh_status()

    def reset_all(self) -> None:
        PARAMETERS.reset()
        for key in self._items:
            self._rewrite(key)
        self.apply_filter()
        self._refresh_status()

    def _import(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Import parameter overrides", "",
                                              "JSON (*.json)")
        if not path:
            return
        try:
            unknown = PARAMETERS.read_overrides(path)
        except (OSError, ValueError) as exc:
            return self._say(f"<span style='color:#f26d6d'>not imported: {exc}</span>")
        self._populate()
        if unknown:
            self._say(f"<span style='color:#f2a65a'>imported; ignored unknown keys: "
                      f"{', '.join(unknown)}</span><br>{PARAMETERS.override_summary()}")

    def _export(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Export parameter overrides",
                                              "ip3r_parameters.json", "JSON (*.json)")
        if path:
            PARAMETERS.write_overrides(path)
            self._say(f"{len(PARAMETERS.overrides())} override(s) written to {path}; "
                      f"reproduce with <tt>IP3R_PARAMETERS={path}</tt>")

    # ------------------------------------------------------------- filtering

    def apply_filter(self) -> None:
        needle, only = self.filter_edit.text(), self.modified_only.isChecked()
        for i in range(self.tree.topLevelItemCount()):
            parent = self.tree.topLevelItem(i)
            shown = 0
            for j in range(parent.childCount()):
                child = parent.child(j)
                ok = PARAMETERS.matches(child.data(0, KEY), needle, only)
                child.setHidden(not ok)
                shown += ok
            parent.setHidden(shown == 0)

    def _say(self, html: str) -> None:
        self.status.setText(html)

    def _refresh_status(self) -> None:
        if PARAMETERS.modified:
            self._say(f"<b style='color:#f2a65a'>{PARAMETERS.override_summary()}</b><br>"
                      "The findings checks will report <i>not run</i> until every "
                      "value is back at its default.")
        else:
            self._say(f"{len(PARAMETERS)} parameters, all at their documented defaults.")
