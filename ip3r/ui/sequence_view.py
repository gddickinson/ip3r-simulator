"""A painted, wrapped sequence grid with drag selection (Round 7.18,
PIEZO1's ``sequence_view`` adapted).

Painted with QPainter in one pass: ITPR subunits are ~2,700 residues, and a
widget per residue would be thousands to lay out on every resize.

Selection is reported in **residue numbers**, never string offsets: a
construct starts wherever its author numbering starts and may skip, and an
offset handed to something that expects a residue number is the numbering
bug this project keeps catching elsewhere.
"""

from __future__ import annotations

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPen
from PyQt6.QtWidgets import QSizePolicy, QToolTip, QWidget

__all__ = ["SequenceView", "CHEMISTRY"]

#: Letter colours by broad chemistry.
CHEMISTRY = {**{a: "#6fb1ff" for a in "AVLIMFWPC"}, **{a: "#7ed67e" for a in "STNQGY"},
             **{a: "#f26d6d" for a in "DE"}, **{a: "#f2a65a" for a in "KRH"}}


class SequenceView(QWidget):
    #: The selected residue numbers changed (first, last), from a drag.
    selection_changed = pyqtSignal(int, int)

    GUTTER, PAD, RULER = 64, 8, 16

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.sequence = None
        self.decorations: dict = {}
        #: Residue numbers shown as selected (set by the window, or a drag).
        self.selected: set[int] = set()
        self._anchor: int | None = None
        self._columns = 50
        font = QFont("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self.setFont(font)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

    # ------------------------------------------------------------------ data

    def set_sequence(self, sequence, decorations: dict | None = None) -> None:
        self.sequence, self.decorations = sequence, decorations or {}
        self.selected, self._anchor = set(), None
        self._relayout()
        self.update()

    def set_decorations(self, decorations: dict) -> None:
        self.decorations = decorations or {}
        self.update()

    def set_selected(self, residues) -> None:
        self.selected = {int(r) for r in residues}
        self.update()

    # ---------------------------------------------------------------- layout

    def _cell(self) -> float:
        return max(QFontMetricsF(self.font()).horizontalAdvance("W"), 7.0) + 2.0

    def _row_height(self) -> float:
        return QFontMetricsF(self.font()).height() + self.RULER

    def _relayout(self) -> None:
        usable = max(self.width() - self.GUTTER - 2 * self.PAD, self._cell() * 10)
        self._columns = max(10, int(usable / self._cell()) // 10 * 10)   # tens
        if self.sequence is not None:
            rows = -(-len(self.sequence) // self._columns)
            self.setMinimumHeight(int(rows * self._row_height() + 2 * self.PAD))

    def resizeEvent(self, event) -> None:            # noqa: N802
        super().resizeEvent(event)
        self._relayout()

    def index_at(self, x: float, y: float) -> int | None:
        if self.sequence is None:
            return None
        row = int((y - self.PAD - self.RULER) // self._row_height())
        col = int((x - self.GUTTER - self.PAD) // self._cell())
        if row < 0 or not 0 <= col < self._columns:
            return None
        i = row * self._columns + col
        return i if 0 <= i < len(self.sequence) else None

    def row_top(self, residue: int) -> int:
        """The y of the row holding ``residue`` (for scrolling to it)."""
        try:
            i = self.sequence.positions.index(int(residue))
        except (AttributeError, ValueError):
            return 0
        return int(self.PAD + (i // self._columns) * self._row_height())

    # ----------------------------------------------------------------- paint

    def paintEvent(self, event) -> None:             # noqa: N802
        p = QPainter(self)
        p.fillRect(self.rect(), QColor("#14171c"))
        if self.sequence is None:
            p.setPen(QColor("#6f7684"))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "no structure loaded")
            p.end()
            return
        p.setFont(self.font())
        cell, line, row_h = self._cell(), QFontMetricsF(self.font()).height(), self._row_height()
        small = QFont(self.font())
        small.setPointSizeF(self.font().pointSizeF() * 0.7)
        visible = event.rect()
        seq = self.sequence
        for i, (letter, residue) in enumerate(zip(seq.letters, seq.positions)):
            row, col = divmod(i, self._columns)
            y = self.PAD + self.RULER + row * row_h
            if y + line < visible.top() or y - self.RULER > visible.bottom():
                continue
            x = self.GUTTER + self.PAD + col * cell
            box = QRectF(x, y, cell, line)
            if col == 0:
                p.setPen(QColor("#6f7684"))
                p.drawText(QRectF(0, y, self.GUTTER, line),
                           Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                           str(residue))
            dec = self.decorations.get(residue)
            if residue in self.selected:
                p.fillRect(box, QColor(255, 210, 64, 150))
            elif dec is not None and dec.background:
                c = QColor(dec.background)
                c.setAlpha(130)
                p.fillRect(box, c)
            colour = QColor(CHEMISTRY.get(letter, "#c8ccd4"))
            if not seq.resolved[i]:
                colour.setAlpha(80)                   # unresolved: dimmed, kept
            p.setPen(colour)
            p.drawText(box, Qt.AlignmentFlag.AlignCenter, letter)
            if dec is not None and dec.underline:
                p.setPen(QPen(QColor(dec.underline), 2.0))
                p.drawLine(QPointF(x + 1, y + line), QPointF(x + cell - 1, y + line))
            if residue % 10 == 0:
                p.setFont(small)
                p.setPen(QColor("#4a515c"))
                p.drawText(QRectF(x - cell, y - self.RULER, cell * 3, self.RULER - 2),
                           Qt.AlignmentFlag.AlignCenter, str(residue))
                p.setFont(self.font())
        p.end()

    # ----------------------------------------------------------------- input

    def _residue(self, event) -> int | None:
        i = self.index_at(event.position().x(), event.position().y())
        return None if i is None else self.sequence.positions[i]

    def mousePressEvent(self, event) -> None:        # noqa: N802
        r = self._residue(event)
        if r is not None:
            self._anchor = r
            self._drag_to(r)

    def mouseMoveEvent(self, event) -> None:         # noqa: N802
        r = self._residue(event)
        if r is None:
            return
        dec = self.decorations.get(r)
        if dec is not None and dec.tooltip:
            QToolTip.showText(event.globalPosition().toPoint(), dec.tooltip, self)
        if self._anchor is not None and event.buttons():
            self._drag_to(r)

    def mouseReleaseEvent(self, event) -> None:      # noqa: N802
        self._anchor = None

    def _drag_to(self, r: int) -> None:
        lo, hi = sorted((self._anchor, r))
        self.selected = {q for q in self.sequence.positions if lo <= q <= hi}
        self.update()
        self.selection_changed.emit(lo, hi)
