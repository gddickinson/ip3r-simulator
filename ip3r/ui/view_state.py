"""Panel view state for sessions, and controls that follow the parameters.

Two small tools the panels share:

* **Restoring a saved control.** ``set_combo``, ``set_spin`` and
  ``set_check`` put a value from a session file onto a widget. A value that
  cannot go there (an unknown key, a wrong type, a number outside the
  widget's range) is not forced: it is reported in ``notes``, and the
  widget keeps what it had or takes the clamped value, which is said.
* **Seeded controls.** A spin box whose starting value came from a
  registered parameter (``Seeded``) follows an edit of that parameter while
  it still shows the value it was seeded with. A value the user typed is
  theirs and is kept. Without this, editing ``puff.n_channels`` changed
  nothing in the GUI: the panel simulated the number its spin box had been
  given when it was built.

``ParameterFollower`` carries the registry's change callback onto the GUI
thread as a Qt signal, and unsubscribes when its owner goes.
"""

from __future__ import annotations

import math

from PyQt6.QtCore import QObject, pyqtSignal

from ..parameters import PARAMETERS

__all__ = ["ParameterFollower", "Seeded", "set_combo", "set_spin", "set_check",
           "is_number"]


def is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def set_combo(combo, value, what: str, notes: list[str]) -> bool:
    """Select the item whose data is ``value``; note it if there is none."""
    i = combo.findData(value)
    if i < 0:
        notes.append(f"unknown {what} {value!r} kept as it was")
        return False
    combo.setCurrentIndex(i)
    return True


def set_spin(spin, value, what: str, notes: list[str]) -> bool:
    """Set a spin box or slider, reporting a clamp or a refusal."""
    if not is_number(value):
        notes.append(f"{what} {value!r} is not a number; kept as it was")
        return False
    lo, hi = spin.minimum(), spin.maximum()
    if not lo <= value <= hi:
        notes.append(f"{what} {value:g} outside {lo:g}-{hi:g}; clamped")
    spin.setValue(type(spin.value())(min(max(value, lo), hi)))
    return True


def set_check(box, value, what: str, notes: list[str]) -> bool:
    if not isinstance(value, bool):
        notes.append(f"{what} {value!r} is not true/false; kept as it was")
        return False
    box.setChecked(value)
    return True


class ParameterFollower(QObject):
    """``changed`` fires on the GUI thread after every effective change of a
    registered value, whoever made it (dialog, banner, session restore)."""

    changed = pyqtSignal()

    def __init__(self, owner: QObject) -> None:
        super().__init__(owner)
        notify = self.changed.emit              # one object, so unsubscribe matches
        PARAMETERS.subscribe(notify)
        self.destroyed.connect(lambda *_: PARAMETERS.unsubscribe(notify))


class Seeded:
    """Spin boxes seeded from parameters, following an edit until the user
    types their own value.

    ``getter`` is called at seeding and at every change, so it may depend on
    other controls (the puff cluster's size depends on the receptor chosen).
    """

    def __init__(self) -> None:
        self._rows: dict[int, list] = {}           # id(spin) -> [spin, getter, seeded]

    def seed(self, spin, getter) -> None:
        """Give ``spin`` its parameter value now (overwriting what it shows)."""
        spin.setValue(type(spin.value())(getter()))
        self._rows[id(spin)] = [spin, getter, spin.value()]

    def is_seeded(self, spin) -> bool:
        """True while ``spin`` still shows the value it was last seeded with."""
        row = self._rows.get(id(spin))
        return row is not None and spin.value() == row[2]

    def follow(self) -> list[str]:
        """Re-seed every spin still at its seeded value; returns their names
        (``objectName``) for a status line."""
        moved = []
        for spin, getter, _ in list(self._rows.values()):
            if self.is_seeded(spin):
                before = spin.value()
                self.seed(spin, getter)
                if spin.value() != before:
                    moved.append(spin.objectName() or "a control")
        return moved
