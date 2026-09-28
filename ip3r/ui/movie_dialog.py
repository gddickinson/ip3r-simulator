"""File → Record movie…: what to record and how, then where to save it.

The defaults come from the ``movie.*`` parameters. What cannot be recorded
now (no mode animating, no transition built, no MP4 encoder) is listed but
disabled, with the reason as its tooltip, so the dialog never offers a movie
that would fail.
"""

from __future__ import annotations

from PyQt6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                             QFileDialog, QFormLayout, QLabel, QSpinBox)

from ..parameters import PARAMETERS as _P
from .movie_model import FORMATS, SOURCE_LABELS, SOURCES, MovieSpec, formats_available

__all__ = ["MovieDialog", "record_movie"]


def _disable(combo: QComboBox, reasons: dict) -> None:
    model = combo.model()
    for i in range(combo.count()):
        why = reasons.get(combo.itemData(i))
        if why:
            model.item(i).setEnabled(False)
            model.item(i).setToolTip(why)


class MovieDialog(QDialog):
    def __init__(self, recorder, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Record movie")
        self.recorder = recorder
        form = QFormLayout(self)
        self.source = QComboBox()
        for s in SOURCES:
            self.source.addItem(SOURCE_LABELS[s], s)
        reasons = recorder.unavailable()
        _disable(self.source, reasons)
        usable = [s for s in SOURCES if not reasons[s]]
        if usable:                           # the most specific thing on screen
            self.source.setCurrentIndex(SOURCES.index(usable[-1]))
        self.frames = QSpinBox()
        self.frames.setRange(2, 720)
        self.fps = QSpinBox()
        self.fps.setRange(1, 60)
        self.fps.setValue(int(_P.value("movie.fps")))
        self.fps.setSuffix(" frames/s")
        self.turn = QCheckBox("Turn the camera once while it plays")
        self.hud = QCheckBox("Include the HUD (title, scale bar, axis)")
        self.hud.setChecked(True)
        self.width = QSpinBox()
        self.width.setRange(64, 3840)
        self.width.setValue(int(_P.value("movie.max_width")))
        self.width.setSuffix(" px wide at most")
        self.fmt = QComboBox()
        for f in FORMATS:
            self.fmt.addItem(f.upper(), f)
        _disable(self.fmt, formats_available())
        self.note = QLabel()
        self.note.setWordWrap(True)
        form.addRow("Show", self.source)
        form.addRow("Frames", self.frames)
        form.addRow("Rate", self.fps)
        form.addRow("", self.turn)
        form.addRow("", self.hud)
        form.addRow("Size", self.width)
        form.addRow("Format", self.fmt)
        form.addRow(self.note)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        self.ok = buttons.addButton("Record…", QDialogButtonBox.ButtonRole.AcceptRole)
        self.ok.setEnabled(bool(usable))
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)
        self.source.currentIndexChanged.connect(self._follow)
        self._follow()

    def _follow(self) -> None:
        """Frames and the camera turn mean different things per source."""
        s = self.source.currentData()
        morph = self.recorder.win.morph.result
        self.frames.setEnabled(s != "transition")
        self.frames.setValue(int(_P.value("movie.turntable_frames" if s == "turntable"
                                          else "movie.mode_frames")))
        self.turn.setEnabled(s != "turntable")
        self.turn.setChecked(s == "turntable")
        if s == "transition" and morph is not None:
            n = len(morph.trajectory)
            hold = int(_P.value("movie.hold_frames"))
            self.note.setText(f"The morph's {n} solved frames there and back, each end "
                              f"held {hold} frames: {2 * n - 2 + 2 * hold} frames.")
        elif s == "mode":
            self.note.setText("One cycle of the animating mode; its shape is the "
                              "model's, its amplitude and speed are illustrative.")
        else:
            self.note.setText("One full turn about the screen's vertical axis.")
        why = self.recorder.unavailable().get(s)
        self.ok.setEnabled(not why)
        if why:
            self.note.setText(why)

    def spec(self) -> MovieSpec:
        return MovieSpec(source=self.source.currentData(), frames=self.frames.value(),
                         fps=float(self.fps.value()), turn=self.turn.isChecked(),
                         hold=int(_P.value("movie.hold_frames")),
                         hud=self.hud.isChecked(), max_width=self.width.value(),
                         fmt=self.fmt.currentData())


def record_movie(win) -> None:
    """The File menu's entry: the dialog, a file name, then the recording."""
    rec = win.movies
    if rec.busy:
        win.statusBar().showMessage("a movie is already being recorded")
        return
    dlg = MovieDialog(rec, win)
    if dlg.exec() != QDialog.DialogCode.Accepted:
        return
    spec = dlg.spec()
    name = f"{win.scene.structure.name}_{spec.source}{FORMATS[spec.fmt]}"
    path, _ = QFileDialog.getSaveFileName(win, "Save movie", name,
                                          f"{spec.fmt.upper()} (*{FORMATS[spec.fmt]})")
    if path:
        rec.start(spec, path)
