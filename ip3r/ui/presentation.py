"""Full-screen viewing (View → Full screen, F11 or the platform's own key;
after PIEZO1's ``presentation.py``).

Hides the panels, the menu and status bars so the 3-D view fills the
screen, for looking, demonstrating or recording. Leaving puts everything
back *as it was*: a panel closed before entering stays closed, a floating
one floats again, and the window returns to its previous state (normal or
maximised). Three things deliberately stay:

* the HUD, whose scale bar and prediction banner are what make a
  full-screen image readable (View → HUD switches its parts);
* the amber "parameters modified" strip, when parameters are modified: it
  is a warning about what is on screen, not furniture;
* every shortcut: the menus' actions are also the window's own, so they
  fire with the menu bar hidden (Esc leaves full screen).

The camera refits on the resize by itself while it has not been moved.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, QTimer

__all__ = ["PresentationController", "HINT_MS"]

#: How long the "how to leave" hint shows on entering (ms).
HINT_MS = 3000


class PresentationController:
    def __init__(self, win) -> None:
        self.win = win
        self._restore: dict | None = None

    @property
    def active(self) -> bool:
        return self._restore is not None

    def toggle(self, on: bool | None = None) -> None:
        on = (not self.active) if on is None else bool(on)
        self.enter() if on else self.leave()

    def enter(self) -> None:
        if self.active:
            return
        win = self.win
        self._restore = {
            "state": win.windowState() & ~Qt.WindowState.WindowFullScreen,
            "docks": [(d, d.isVisible(), d.isFloating()) for d in win.docks.docks],
            "menu": win.menuBar().isVisible(),
            "status": win.statusBar().isVisible(),
        }
        for d in win.docks.docks:
            d.hide()
        win.menuBar().setVisible(False)
        win.statusBar().setVisible(False)
        win.showFullScreen()
        win.viewport.setFocus()
        self._sync_action(True)
        win.hud.set_readout("presentation", "full screen: Esc or F11 to leave")
        QTimer.singleShot(HINT_MS, lambda: win.hud.set_readout("presentation", ""))

    def leave(self) -> None:
        if not self.active:
            return
        win, r = self.win, self._restore
        self._restore = None
        win.setWindowState(r["state"])
        if not r["state"] & (Qt.WindowState.WindowMaximized | Qt.WindowState.WindowMinimized):
            win.showNormal()
        win.menuBar().setVisible(r["menu"])
        win.statusBar().setVisible(r["status"])
        for dock, visible, floating in r["docks"]:
            dock.setFloating(floating)
            dock.setVisible(visible)
        win.hud.set_readout("presentation", "")
        self._sync_action(False)

    def _sync_action(self, on: bool) -> None:
        a = getattr(self.win, "fullscreen_action", None)
        if a is not None and a.isChecked() != on:
            a.blockSignals(True)
            a.setChecked(on)
            a.blockSignals(False)
