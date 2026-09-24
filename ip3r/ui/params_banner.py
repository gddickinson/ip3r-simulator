"""The "parameters modified" banner above the viewport.

Shown for as long as any registered value differs from its documented
default, whatever set it (the Parameters dialog, or an ``IP3R_PARAMETERS``
file at start-up), and hidden otherwise. It says which values, and that the
findings checks will not confirm, so no number on screen can be mistaken for
one computed at the published defaults.
"""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel, QPushButton

from ..parameters import PARAMETERS

__all__ = ["ParametersBanner"]

STYLE = ("QFrame#paramsBanner { background: #5a3a12; border-bottom: 2px solid #f2a65a; }"
         "QFrame#paramsBanner QLabel { color: #ffe2c0; }")


class ParametersBanner(QFrame):
    """Amber strip; ``edit_requested`` opens the dialog. The host shows it
    while ``visibility_changed`` says so (a toolbar, in the main window)."""

    edit_requested = pyqtSignal()
    visibility_changed = pyqtSignal(bool)   # the host strip follows
    _changed = pyqtSignal()              # registry callback -> GUI thread

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setObjectName("paramsBanner")
        self.setStyleSheet(STYLE)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(10, 4, 10, 4)
        self.text = QLabel()
        lay.addWidget(self.text, 1)
        edit = QPushButton("Edit…")
        edit.clicked.connect(self.edit_requested)
        lay.addWidget(edit)
        reset = QPushButton("Reset all")
        reset.clicked.connect(lambda: PARAMETERS.reset())
        lay.addWidget(reset)
        self._changed.connect(self.refresh)
        notify = self._changed.emit              # one object, so unsubscribe matches
        PARAMETERS.subscribe(notify)
        self.destroyed.connect(lambda *_: PARAMETERS.unsubscribe(notify))
        self.refresh()

    def refresh(self) -> None:
        overrides = PARAMETERS.overrides()
        self.visibility_changed.emit(bool(overrides))
        if overrides:
            listed = ", ".join(f"{k} = {v:g} (default {PARAMETERS.default(k):g})"
                               for k, v in sorted(overrides.items())[:3])
            more = f", and {len(overrides) - 3} more" if len(overrides) > 3 else ""
            self.text.setText(f"<b>Parameters modified</b>: {listed}{more}. "
                              "Findings checks will not confirm.")
            self.setToolTip(PARAMETERS.override_summary() + "\n\nResults computed now "
                            "are not the published defaults'. Results already on "
                            "screen keep the values they were computed with.")
