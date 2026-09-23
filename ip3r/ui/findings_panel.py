"""The findings tab: every ip3r_genes result this application re-derives.

Checks are grouped by the paper they belong to. Running them happens on a
worker thread (the structural ones parse six tetramers). Selecting a check
shows the claim as published, how this application re-derives it, both
numbers, the verdict, and — where there is one — an exhibit drawn from the
check's own numbers. Checks about a structure offer to show it: the main
window loads the deposit and highlights what the check measured.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (QHBoxLayout, QLabel, QPushButton, QSplitter,
                             QTextBrowser, QTreeWidget, QTreeWidgetItem,
                             QVBoxLayout, QWidget)

from ..analysis import exhibits
from ..analysis.checks import PAPERS, all_checks, run_checks
from ..config import genes_results
from .plot_canvas import PlotCanvas
from .theme import STATUS_COLORS
from .workers import run_async

__all__ = ["FindingsPanel"]

_KIND_TEXT = {"recomputed": "recomputed from primary data (structures)",
              "rederived": "re-derived from the committed input tables",
              "read": "read from the published table"}


class FindingsPanel(QWidget):
    show_structure = pyqtSignal(str, str)           # pdb id, check id

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        head = QLabel(f"Re-deriving the results of <b>ip3r_genes</b> "
                      f"(<code>{genes_results()}</code>).")
        head.setWordWrap(True)
        lay.addWidget(head)
        row = QHBoxLayout()
        self.run_btn = QPushButton("Run all checks")
        self.run_btn.clicked.connect(self.run_all)
        row.addWidget(self.run_btn)
        self.summary = QLabel("")
        row.addWidget(self.summary, 1)
        lay.addLayout(row)
        split = QSplitter(Qt.Orientation.Vertical)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["check", "verdict", "kind"])
        self.tree.itemSelectionChanged.connect(self._selected)
        split.addWidget(self.tree)
        low = QWidget()
        ll = QVBoxLayout(low)
        ll.setContentsMargins(0, 0, 0, 0)
        self.detail = QTextBrowser()
        self.detail.setOpenExternalLinks(False)
        ll.addWidget(self.detail, 2)
        self.show_btn = QPushButton("Show on structure")
        self.show_btn.clicked.connect(self._show)
        self.show_btn.setEnabled(False)
        ll.addWidget(self.show_btn)
        self.canvas = PlotCanvas(low, height=2.8)
        self.canvas.reset().set_axis_off()
        ll.addWidget(self.canvas, 2)
        split.addWidget(low)
        lay.addWidget(split, 1)
        self.results: dict = {}
        self._items: dict = {}
        self._populate()

    def _populate(self) -> None:
        self.tree.clear()
        groups = {}
        for paper, title in PAPERS.items():
            g = QTreeWidgetItem([title, "", ""])
            g.setFirstColumnSpanned(True)
            groups[paper] = g
        for c in all_checks():
            it = QTreeWidgetItem([c.id, "—", c.kind])
            it.setData(0, Qt.ItemDataRole.UserRole, c.id)
            it.setToolTip(0, c.claim)
            groups[c.paper].addChild(it)
            self._items[c.id] = it
        for g in groups.values():
            if g.childCount():
                self.tree.addTopLevelItem(g)
                g.setExpanded(True)
        self.tree.resizeColumnToContents(0)

    def run_all(self) -> None:
        self.run_btn.setEnabled(False)
        self.summary.setText("running…")
        run_async(run_checks, on_done=self._done,
                  on_error=lambda e: (self.summary.setText(e),
                                      self.run_btn.setEnabled(True)))

    def _done(self, results) -> None:
        self.run_btn.setEnabled(True)
        counts: dict[str, int] = {}
        for r in results:
            self.results[r.check.id] = r
            counts[r.status] = counts.get(r.status, 0) + 1
            it = self._items[r.check.id]
            it.setText(1, r.status)
            it.setForeground(1, QColor(STATUS_COLORS.get(r.status, "#d7dbe3")))
        self.summary.setText(", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
        self._selected()

    def _current(self):
        items = self.tree.selectedItems()
        if not items:
            return None
        return items[0].data(0, Qt.ItemDataRole.UserRole)

    def _selected(self) -> None:
        cid = self._current()
        if cid is None:
            return
        c = next(x for x in all_checks() if x.id == cid)
        r = self.results.get(cid)
        colour = STATUS_COLORS.get(r.status, "#d7dbe3") if r else "#8a8f99"
        html = [f"<h3>{c.id}</h3><p><b>Claim.</b> {c.claim}</p>",
                f"<p><b>How it is re-derived</b> ({_KIND_TEXT[c.kind]}). {c.method}</p>",
                f"<p><b>Sources:</b> {', '.join(c.sources) or '—'}"
                + (f"; structures {', '.join(c.structures)}" if c.structures else "") + "</p>"]
        if r:
            o = r.outcome
            html.append(f"<p><b>Verdict:</b> <span style='color:{colour}'>{o.status}"
                        f"</span> ({r.seconds:.1f} s)</p>")
            html.append(f"<p><b>Published:</b> {o.published}<br><b>Found:</b> {o.found}</p>")
            if o.detail:
                html.append(f"<p>{o.detail}</p>")
        else:
            html.append("<p><i>Not run yet.</i></p>")
        self.detail.setHtml("".join(html))
        self.show_btn.setEnabled(bool(c.structures))
        ax = self.canvas.reset()
        if r is None or not exhibits.draw(ax, cid, r.outcome):
            ax.set_axis_off()
        self.canvas.draw_now()

    def _show(self) -> None:
        cid = self._current()
        c = next((x for x in all_checks() if x.id == cid), None)
        if c and c.structures:
            self.show_structure.emit(c.structures[0], cid)
