"""The Genomes tab: Papers 3 and 4 as a genome × paralog grid.

One row per assembly in the retention sweep (309), one column per cell
(ITPR1-3 and the RyR control), coloured by a chosen layer: what the sweep
found, which known genes it missed, the S15a evidence state, how a
protein-database search could reach the gene, or whether the cell carries
more lesions than its own genome's identity-matched siblings (S15b §8). A contig-N50 strip on a fixed
log scale runs beside the grid, and sorting by N50 draws the contiguity bar.
Clicking a row names the genome and lists its four cells.
"""

from __future__ import annotations

from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel, QPushButton,
                             QVBoxLayout, QWidget)

from ..analysis import genome_grid as GG
from ..analysis.grid_figure import LAYER_TITLES, draw_grid
from .plot_canvas import PlotCanvas
from .workers import run_async

__all__ = ["GenomesPanel"]

_ALL = "all classes"
#: The class a layer's published statement is about: a check's "Show"
#: opens the layer on it (the lesion excess is a bird result).
LAYER_FOCUS = {"lesion": "Aves"}


class GenomesPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.grid: GG.GenomeGrid | None = None
        self.shown: GG.GenomeGrid | None = None
        self.info: dict = {}
        self._pending_layer: str | None = None
        self._loading = False
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.layer = QComboBox()
        for key, title in LAYER_TITLES.items():
            self.layer.addItem(title, key)
        self.order = QComboBox()
        self.order.addItems(list(GG.ORDERS))
        self.vclass = QComboBox()
        self.vclass.addItem(_ALL)
        self.above = QCheckBox("Above the bar only")
        self.above.setToolTip("Only assemblies whose contig N50 reaches the "
                              "registered contiguity bar (genomes.contiguity_bar_bp)")
        self.reload_btn = QPushButton("Reload")
        self.reload_btn.clicked.connect(self.load)
        for w in (self.layer, self.order, self.vclass):
            w.currentIndexChanged.connect(self.redraw)
        self.above.toggled.connect(self.redraw)
        for w in (QLabel("Colour"), self.layer, QLabel("Sort"), self.order):
            row.addWidget(w)
        row.addStretch(1)
        lay.addLayout(row)
        row = QHBoxLayout()                # two rows: the dock stays narrow
        for w in (self.vclass, self.above):
            row.addWidget(w)
        row.addStretch(1)
        row.addWidget(self.reload_btn)
        lay.addLayout(row)
        self.canvas = PlotCanvas(self, width=5.5, height=8.0)
        lay.addWidget(NavigationToolbar2QT(self.canvas, self))
        lay.addWidget(self.canvas, 1)
        self.status = QLabel("Not loaded.")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)
        self.detail = QLabel("")
        self.detail.setWordWrap(True)
        lay.addWidget(self.detail)
        self.canvas.mpl_connect("button_press_event", self._clicked)

    # ------------------------------------------------------------ loading
    def load(self) -> None:
        self._loading = True
        self.status.setText("Reading the sweep's per-cell tables …")
        run_async(GG.load_grid, on_done=self._loaded, on_error=self._failed)

    def ensure_loaded(self) -> None:
        if self.grid is None and not self._loading:
            self.load()

    def show_layer(self, layer: str) -> None:
        """Select a layer (from a check's "Show"), loading first if needed."""
        self._pending_layer = layer
        if self.grid is None:
            self.ensure_loaded()
        else:
            self._apply_pending()

    def _apply_pending(self) -> None:
        if self._pending_layer is not None:
            self.order.setCurrentText("contig N50")
            focus = LAYER_FOCUS.get(self._pending_layer)
            self.vclass.setCurrentText(focus if focus and self.vclass.findText(focus) >= 0
                                       else _ALL)
            self.layer.setCurrentIndex(self.layer.findData(self._pending_layer))
            self._pending_layer = None
            self.redraw()

    def _loaded(self, grid) -> None:
        self._loading = False
        self.grid = grid
        self.vclass.blockSignals(True)
        self.vclass.clear()
        self.vclass.addItem(_ALL)
        self.vclass.addItems(sorted({g.vclass for g in grid.genomes}))
        self.vclass.blockSignals(False)
        self.redraw()
        self._apply_pending()

    def _failed(self, err: str) -> None:
        self._loading = False
        msg = err.splitlines()[0] if err else "failed"
        if msg.startswith("GenesDataMissing"):
            msg = "ip3r_genes results are not available (" + msg + ")"
        self.status.setText(msg)

    # ------------------------------------------------------------ drawing
    def redraw(self, *_) -> None:
        if self.grid is None:
            return
        vclass = self.vclass.currentText()
        rows = [i for i, g in enumerate(self.grid.genomes)
                if (vclass == _ALL or g.vclass == vclass)
                and (not self.above.isChecked() or g.above_bar)]
        grid = GG.order(self.grid.subset(rows), self.order.currentText())
        self.shown = grid
        layer = self.layer.currentData()
        ax = self.canvas.reset()
        if grid.n == 0:
            self.status.setText("No genome matches the filter.")
            self.canvas.draw_now()
            return
        self.info = draw_grid(ax, grid, layer)
        ax.set_title(LAYER_TITLES[layer])
        self.canvas.draw_now()
        itpr = grid.counts(layer, ["ITPR1", "ITPR2", "ITPR3"])
        parts = ", ".join(f"{k} {v}" for k, v in sorted(itpr.items(), key=lambda x: -x[1]))
        self.status.setText(f"{grid.n} genomes ({self.info['above_bar']} above the "
                            f"contiguity bar). ITPR cells: {parts}.")

    def _clicked(self, event) -> None:
        if event.inaxes is None or event.ydata is None or self.shown is None:
            return
        i = int(round(event.ydata))
        if not 0 <= i < self.shown.n:
            return
        g = self.shown.genomes[i]
        cells = self.shown.row(i)
        text = "; ".join(f"{c}: " + ", ".join(v for v in cells[c].values() if v)
                         for c in self.shown.cells)
        self.detail.setText(f"{g.organism} ({g.accession}, {g.vclass}); contig N50 "
                            f"{g.contig_n50:,.0f} bp, {g.level}, {g.source} — {text}")
