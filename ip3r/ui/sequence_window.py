"""The sequence window (View → Sequence…, Ctrl+Shift+Q; Round 7.18, after
PIEZO1's).

One chain's construct at a time, in the deposit's author numbers, with a
track painted behind the letters (element, deep JSD, or which residues are
resolved) and the measured sites underlined. A drag selects onto the 3-D
model, on this chain or, ticked, on every subunit (the C4 copies); a
selection made on the model shows here. It reads and writes the one
:class:`~ip3r.ui.selection.SelectionController`, so the two cannot drift.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel, QScrollArea,
                             QVBoxLayout, QWidget)

from .sequence_model import SITE_COLOURS, TRACKS, chain_sequence, track
from .sequence_view import SequenceView

__all__ = ["SequenceWindow"]


class SequenceWindow(QWidget):
    def __init__(self, selection, parent=None) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Sequence")
        self.resize(820, 560)
        self.selection = selection
        self.structure = None
        self.paralog: str | None = None
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel("Chain"))
        self.chain = QComboBox()
        self.chain.currentIndexChanged.connect(lambda _: self._show_chain())
        row.addWidget(self.chain)
        row.addWidget(QLabel("Track"))
        self.track = QComboBox()
        for k, label in TRACKS.items():
            self.track.addItem(label, k)
        self.track.currentIndexChanged.connect(lambda _: self._paint())
        row.addWidget(self.track, 1)
        self.everywhere = QCheckBox("select on every subunit")
        self.everywhere.setToolTip("A drag selects the same residue numbers on "
                                   "each chain that resolves them.")
        row.addWidget(self.everywhere)
        lay.addLayout(row)
        keys = "  ".join(f"<span style='color:{c}'>▁ {k.replace('_', ' ')}</span>"
                         for k, c in SITE_COLOURS.items())
        legend = QLabel(f"underlined: {keys}; dim letters are unresolved")
        lay.addWidget(legend)
        self.info = QLabel("")
        self.info.setWordWrap(True)
        lay.addWidget(self.info)
        self.view = SequenceView()
        self.view.selection_changed.connect(self._dragged)
        self.scroll = QScrollArea()
        self.scroll.setWidget(self.view)
        self.scroll.setWidgetResizable(True)
        lay.addWidget(self.scroll, 1)
        selection.changed.connect(self._follow)

    # --------------------------------------------------------------- state

    def set_structure(self, st, paralog: str | None) -> None:
        self.structure, self.paralog = st, paralog
        self.chain.blockSignals(True)
        self.chain.clear()
        for c in st.chains:
            if ((st.chain == c) & ~st.hetero).any():
                self.chain.addItem(c, c)
        self.chain.blockSignals(False)
        self._show_chain()

    def show_residue(self, chain: str | None = None, residue: int | None = None) -> None:
        """Bring ``chain`` up and scroll to ``residue``."""
        if chain is not None and self.chain.findData(chain) >= 0:
            self.chain.setCurrentIndex(self.chain.findData(chain))
        if residue is not None:
            self.scroll.verticalScrollBar().setValue(
                max(0, self.view.row_top(residue) - 40))
        self.show()
        self.raise_()

    @property
    def current_chain(self) -> str | None:
        return self.chain.currentData()

    # -------------------------------------------------------------- drawing

    def _show_chain(self) -> None:
        st, ch = self.structure, self.current_chain
        if st is None or ch is None:
            self.view.set_sequence(None)
            return
        self.seq = chain_sequence(st, ch)
        self.view.set_sequence(self.seq)
        self._paint()
        self._follow()

    def _paint(self) -> None:
        if self.structure is None or self.current_chain is None:
            return
        key = self.track.currentData()
        self.view.set_decorations(track(self.seq, key, self.paralog))
        s = self.seq
        numbering = (f"in human {self.paralog} numbering" if self.paralog
                     else "in no human numbering: element, JSD and sites are grey")
        what = ("the construct (poly-seq scheme)" if s.source == "construct"
                else "the built residues (no poly-seq scheme in the file)")
        self.info.setText(f"{self.structure.name} chain {s.chain}: {len(s):,} residues of "
                          f"{what}, {s.n_resolved:,} resolved; {numbering}.")

    # ------------------------------------------------------------ selection

    def _dragged(self, lo: int, hi: int) -> None:
        st, ch = self.structure, self.current_chain
        chains = ([c for c in (self.chain.itemData(i) for i in range(self.chain.count()))]
                  if self.everywhere.isChecked() else [ch])
        picked = set()
        for c in chains:
            built = set(st.res_seq[(st.chain == c) & ~st.hetero].tolist())
            picked |= {(c, r) for r in range(lo, hi + 1) if r in built}
        self.selection.select(picked)

    def _follow(self) -> None:
        ch = self.current_chain
        self.view.set_selected(r for c, r in self.selection.residues if c == ch)
