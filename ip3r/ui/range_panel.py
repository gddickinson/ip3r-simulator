"""The Range tab: Paper 1's presence and absence across the tree of life.

One bar per clade of the S20 proteome sweep (the prokaryotes collapsed to one
row per domain unless expanded): the fraction of swept proteomes carrying an
IP3 receptor call on a fixed 0–1 scale, coloured by supergroup. A cross marks
a clade whose absence was confirmed in controlled genome assemblies (S23),
rebuilt here from the per-genome ledgers. Clicking a row names the clade and
lists its genome-level absences; double-clicking it opens the clade's S23
genomes one by one (the "Genomes (S23)" view: control verdict, copy-ledger
status, contiguity against the genome's own bar, complete gene models).
"""

from __future__ import annotations

from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel, QPushButton,
                             QSpinBox, QVBoxLayout, QWidget)

from ..analysis import range_genomes as RG
from ..analysis import range_table as RT
from ..analysis.range_figure import draw_genomes, draw_range
from .plot_canvas import PlotCanvas
from .workers import run_async

__all__ = ["RangePanel", "CLADES", "GENOMES"]

CLADES, GENOMES = "Clades (S20)", "Genomes (S23)"


def _load():
    return RT.load_range(), RG.genome_rows()


class RangePanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.data: RT.Range | None = None
        self.genomes: list[RG.GenomeRow] = []
        self.shown_genomes: list[RG.GenomeRow] = []
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
        lay.addLayout(row)
        row = QHBoxLayout()                # two rows: the dock stays narrow
        self.view = QComboBox()
        self.view.addItems([CLADES, GENOMES])
        self.clade = QComboBox()
        self.clade.setToolTip("The S20 clade whose S23 genomes are drawn")
        self.view.currentIndexChanged.connect(self.redraw)
        self.clade.currentIndexChanged.connect(self.redraw)
        for w in (QLabel("View"), self.view, self.clade):
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
        run_async(_load, on_done=self._loaded, on_error=self._failed)

    def ensure_loaded(self) -> None:
        if self.data is None and not self._loading:
            self.load()

    def _loaded(self, loaded) -> None:
        self._loading = False
        self.data, self.genomes = loaded
        keep = self.clade.currentData()
        self.clade.blockSignals(True)
        self.clade.clear()
        for c, n in RG.clades_with_genomes(self.genomes):
            self.clade.addItem(f"{c} ({n})", c)
        if keep is not None and self.clade.findData(keep) >= 0:
            self.clade.setCurrentIndex(self.clade.findData(keep))
        self.clade.blockSignals(False)
        self.redraw()

    def show_genomes(self, clade: str) -> bool:
        """Open ``clade``'s S23 genomes; False if it has none."""
        i = self.clade.findData(clade)
        if i < 0:
            return False
        self.clade.blockSignals(True)
        self.clade.setCurrentIndex(i)
        self.clade.blockSignals(False)
        if self.view.currentText() == GENOMES:
            self.redraw()
        else:
            self.view.setCurrentText(GENOMES)
        return True

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
        genomes = self.view.currentText() == GENOMES
        for w in (self.min_n, self.collapse, self.marks):
            w.setEnabled(not genomes)
        self.clade.setEnabled(genomes)
        if genomes:
            self._draw_genomes()
            return
        self.shown_genomes = []
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

    def _draw_genomes(self) -> None:
        clade = self.clade.currentData()
        rows = RG.in_clade(self.genomes, clade) if clade else []
        self.shown_genomes, self.rows = rows, []
        ax = self.canvas.reset()
        if not rows:
            self.status.setText("No S23 genome in this clade.")
            self.canvas.draw_now()
            return
        info = draw_genomes(ax, rows, f"{clade}: S23 genomes")
        self.info = info
        self.canvas.draw_now()
        self.status.setText(
            f"{clade}: {info['genomes']} genomes, {info['controlled']} controlled, "
            f"{info['with_gene']} with a complete gene model ({info['copies']} models). "
            "Grey: not in the ledger.")

    def _genome_clicked(self, event) -> None:
        i = int(round(event.ydata))
        if not 0 <= i < len(self.shown_genomes):
            return
        r = self.shown_genomes[i]
        placed = f"placed by {r.placed_by}" if r.placed_by else "phylum not swept in S20"
        self.detail.setText(
            f"{r.organism} ({r.accession}; {r.phylum or '—'} / {r.klass or '—'}, "
            f"{placed}): control {r.verdict.replace('_', ' ') or '—'}; "
            f"{r.status.replace('_', ' ') or '—'}; {r.copies} complete gene models; "
            f"contig N50 {r.contig_n50:,.0f} bp against its bar {r.bar_bp:,.0f} bp.")

    def _clicked(self, event) -> None:
        if event.inaxes is None or event.ydata is None:
            return
        if self.shown_genomes:
            self._genome_clicked(event)
            return
        if not self.rows:
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
            if self.clade.findData(names[0]) >= 0:
                if event.dblclick:
                    self.show_genomes(names[0])
                    return
                text += " Double-click for its S23 genomes."
        self.detail.setText(text)
