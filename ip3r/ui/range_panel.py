"""The Range tab: Paper 1's presence and absence across the tree of life.

One bar per clade of the S20 proteome sweep (the prokaryotes collapsed to one
row per domain unless expanded): the fraction of swept proteomes carrying an
IP3 receptor call on a fixed 0–1 scale, coloured by supergroup. A cross marks
a clade whose absence was confirmed in controlled genome assemblies (S23),
rebuilt here from the per-genome ledgers. Clicking a row names the clade and
lists its genome-level absences.
"""

from __future__ import annotations

from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PyQt6.QtWidgets import (QCheckBox, QHBoxLayout, QLabel, QPushButton, QSpinBox,
                             QVBoxLayout, QWidget)

from ..analysis import range_table as RT
from ..analysis.range_figure import draw_range
from .plot_canvas import PlotCanvas
from .workers import run_async

__all__ = ["RangePanel"]


class RangePanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.data: RT.Range | None = None
        self.rows: list[tuple] = []
        self.info: dict = {}
        self._loading = False
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.min_n = QSpinBox()
        self.min_n.setRange(1, 100)
        self.min_n.setValue(1)
        self.min_n.setToolTip("Leave out clades with fewer swept proteomes than this")
        self.collapse = QCheckBox("Collapse prokaryotes")
        self.collapse.setChecked(True)
        self.marks = QCheckBox("Genome absences")
        self.marks.setChecked(True)
        self.marks.setToolTip("Mark clades whose absence held in controlled "
                              "genome assemblies (S23)")
        self.reload_btn = QPushButton("Reload")
        self.reload_btn.clicked.connect(self.load)
        self.min_n.valueChanged.connect(self.redraw)
        self.collapse.toggled.connect(self.redraw)
        self.marks.toggled.connect(self.redraw)
        for w in (QLabel("Min proteomes"), self.min_n, self.collapse, self.marks):
            row.addWidget(w)
        row.addStretch(1)
        row.addWidget(self.reload_btn)
        lay.addLayout(row)
        self.canvas = PlotCanvas(self, width=5.5, height=9.0)
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
        self.status.setText("Reading the proteome and genome sweeps …")
        run_async(RT.load_range, on_done=self._loaded, on_error=self._failed)

    def ensure_loaded(self) -> None:
        if self.data is None and not self._loading:
            self.load()

    def _loaded(self, data) -> None:
        self._loading = False
        self.data = data
        self.redraw()

    def _failed(self, err: str) -> None:
        self._loading = False
        msg = err.splitlines()[0] if err else "failed"
        if msg.startswith("GenesDataMissing"):
            msg = "ip3r_genes results are not available (" + msg + ")"
        self.status.setText(msg)

    # ------------------------------------------------------------ drawing
    def _absence_marks(self) -> dict:
        out = {}
        for c in self.data.clades:
            held = [a for a in self.data.absences_in(c.clade) if a.holds]
            out[c.clade] = (sum(a.clade == c.clade for a in held),
                            sum(a.clade != c.clade for a in held))
        return out

    def redraw(self, *_) -> None:
        if self.data is None:
            return
        ax = self.canvas.reset()
        clades = [(c.clade, c.supergroup, c.n, c.present) for c in self.data.clades]
        info = draw_range(ax, clades, self._absence_marks() if self.marks.isChecked() else None,
                          min_n=self.min_n.value(), collapse=self.collapse.isChecked())
        self.info, self.rows = info, info["rows"]
        ax.set_title("Where the IP3 receptor is: S20 proteomes, S23 genomes")
        self.canvas.draw_now()
        held = sum(a.holds for a in self.data.absences)
        self.status.setText(
            f"{info['clades']} rows shown; {info['with_call']} with a call; "
            f"{info['present']:,} of {info['proteomes']:,} proteomes carry one. "
            f"{held} of {len(self.data.absences)} absence targets hold in "
            f"controlled genomes.")

    def _clicked(self, event) -> None:
        if event.inaxes is None or event.ydata is None or not self.rows:
            return
        i = int(round(event.ydata))
        if not 0 <= i < len(self.rows):
            return
        label, group, n, k, names = self.rows[i]
        text = f"{label} ({group}): {k:,} of {n:,} proteomes carry an IP3R call."
        if len(names) == 1 and self.data is not None:
            parts = [f"{a.clade} ({a.rank}): {a.controlled}/{a.genomes} genomes "
                     f"controlled, {a.with_full} with a gene"
                     + (f", {a.trace_only} trace only" if a.trace_only else "")
                     for a in self.data.absences_in(names[0])]
            if parts:
                text += " Genome absences — " + "; ".join(parts) + "."
        self.detail.setText(text)
