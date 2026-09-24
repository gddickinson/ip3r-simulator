"""The Tree tab: Paper 2's RyR-rooted maximum-likelihood tree.

Drawn by :func:`ip3r.analysis.tree_figure.draw_tree` from the committed
``rooted.nwk`` with this project's reader: the paralogue clades and the
outgroup boxed, hagfish and lamprey tips marked, supported nodes dotted.
The matplotlib toolbar zooms and pans; clicking near a tip names it.
"""

from __future__ import annotations

from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PyQt6.QtWidgets import (QCheckBox, QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
                             QWidget)

from ..analysis.tree import group_of, load_tree
from ..analysis.tree_figure import draw_tree, layout
from .plot_canvas import PlotCanvas
from .workers import run_async

__all__ = ["TreePanel"]

_VERT = ("ITPR1", "ITPR2", "ITPR3", "vertebrate_basal")

_KEY = ("Boxes: each paralogue's whole clade (the MRCA of every record naming it; "
        "size and SH-aLRT/UFBoot beside it) and the RyR outgroup. Filled dot: the "
        "record names the paralogue. Open circle: a vertebrate record naming none. "
        "Yellow diamond: hagfish or lamprey, with the support of its clade and of "
        "where it joins. Small white dot: a node clearing both registered support bars.")


class TreePanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.root = None
        self._tips: list = []
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.labels = QCheckBox("Tip labels")
        self.supports = QCheckBox("Support dots")
        self.supports.setChecked(True)
        for b in (self.labels, self.supports):
            b.toggled.connect(self.redraw)
        self.zoom_btn = QPushButton("Vertebrates")
        self.zoom_btn.setToolTip("Zoom to the vertebrate tips (toolbar Home restores)")
        self.zoom_btn.clicked.connect(self.zoom_vertebrates)
        self.reload_btn = QPushButton("Reload")
        self.reload_btn.clicked.connect(self.load)
        row.addWidget(self.labels)
        row.addWidget(self.supports)
        row.addStretch(1)
        row.addWidget(self.zoom_btn)
        row.addWidget(self.reload_btn)
        lay.addLayout(row)
        self.canvas = PlotCanvas(self, width=5.5, height=8.0)
        lay.addWidget(NavigationToolbar2QT(self.canvas, self))
        lay.addWidget(self.canvas, 1)
        key = QLabel(_KEY)
        key.setWordWrap(True)
        lay.addWidget(key)
        self.status = QLabel("Not loaded.")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)
        self.canvas.mpl_connect("button_press_event", self._clicked)

    def load(self) -> None:
        """Parse the committed tree on a worker, then draw it."""
        self.status.setText("Reading phylogeny/rooted.nwk …")
        run_async(load_tree, on_done=self._loaded, on_error=self._failed)

    def ensure_loaded(self) -> None:
        if self.root is None:
            self.load()

    def _loaded(self, root) -> None:
        self.root = root
        self.redraw()

    def _failed(self, err: str) -> None:
        msg = err.splitlines()[0] if err else "failed"
        if msg.startswith("GenesDataMissing"):
            msg = "ip3r_genes results are not available (" + msg + ")"
        self.status.setText(msg)

    def redraw(self, *_, keep_view=None) -> None:
        if self.root is None:
            return
        ax = self.canvas.reset()
        info = draw_tree(ax, self.root, labels=self.labels.isChecked(),
                         supports=self.supports.isChecked())
        if keep_view is not None:
            ax.set_xlim(*keep_view[0])
            ax.set_ylim(*keep_view[1])
        self._tips = layout(self.root).tips
        self.canvas.draw_now()
        clades = ", ".join(f"{g} {n}" for g, n in info["clades"].items())
        cyc = " + ".join(map(str, info["cyclostome_clades"]))
        text = (f"{info['tips']} tips; clades {clades}; cyclostome-only clades of {cyc}")
        if self.supports.isChecked():
            text += f"; {info['supported_nodes']} bipartitions clear both support bars"
        self.status.setText(text + ".")
        self.info = info

    def zoom_vertebrates(self) -> None:
        rows = [i for i, t in enumerate(self._tips) if group_of(t.label) in _VERT]
        if not rows:
            return
        lay = layout(self.root)
        xs = [lay.x[id(self._tips[i])] for i in rows]
        box_edge = max(lay.x.values()) * 1.02          # where draw_tree ends its boxes
        view = ((min(xs) - 0.4, box_edge + 0.35 * (box_edge - min(xs))), (max(rows) + 1.5, min(rows) - 1.5))
        self.labels.blockSignals(True)
        self.labels.setChecked(True)                 # readable at this scale
        self.labels.blockSignals(False)
        self.redraw(keep_view=view)

    def _clicked(self, event) -> None:
        if event.inaxes is None or event.ydata is None or not self._tips:
            return
        i = int(round(event.ydata))
        if 0 <= i < len(self._tips):
            self.status.setText(f"Row {i}: {self._tips[i].label}")
