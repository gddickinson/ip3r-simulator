"""The panels as docks, and the layout they are in (Round 7.18, after
PIEZO1's ``docks.py``).

Every panel is a movable, floatable, closable dock with an object name (Qt
saves a layout by name). The shipped layout is captured once the window
has been shown and sized, so **View → Reset layout** always has somewhere
to return to. With ``remember`` on (the application; never the smoke test,
which must start from the shipped layout) the layout is saved through
:class:`QSettings` on close and restored on the next start. A saved layout
from a larger screen is clamped: the window is moved back onto the screen
it opens on, so it can never restore off every monitor.
"""

from __future__ import annotations

from PyQt6.QtCore import QByteArray, QSettings, Qt
from PyQt6.QtWidgets import QDockWidget, QScrollArea

__all__ = ["DockManager", "SETTINGS_ORG", "SETTINGS_APP", "LAYOUT_VERSION"]

SETTINGS_ORG, SETTINGS_APP = "ip3r_simulation", "ip3r"
#: Bumped when the set of docks changes: an older saved layout is ignored.
LAYOUT_VERSION = 1


class DockManager:
    def __init__(self, win, remember: bool = False) -> None:
        self.win = win
        self.remember = remember
        self.docks: list[QDockWidget] = []
        self._default: QByteArray | None = None

    def add(self, title: str, widget, area, scroll: bool = True,
            min_width: int = 300) -> QDockWidget:
        dock = QDockWidget(title, self.win)
        dock.setObjectName(f"dock_{title.lower()}")
        if scroll:
            holder = QScrollArea()
            holder.setWidget(widget)
            holder.setWidgetResizable(True)
            # Fit the width, scroll only vertically: a sideways scroll hid subunit D.
            holder.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            widget = holder
        dock.setWidget(widget)
        dock.setMinimumWidth(min_width)
        self.win.addDockWidget(area, dock)
        self.docks.append(dock)
        return dock

    def panel_actions(self) -> list:
        """One checkable action per dock (View → Panels)."""
        return [d.toggleViewAction() for d in self.docks]

    # --------------------------------------------------------------- layout

    def capture_default(self) -> None:
        """Remember the shipped layout (called once, after the first sizing)."""
        if self._default is None:
            self._default = self.win.saveState(LAYOUT_VERSION)

    def reset(self) -> None:
        for d in self.docks:
            d.setFloating(False)
            d.show()
        if self._default is not None:
            self.win.restoreState(self._default, LAYOUT_VERSION)

    def _settings(self) -> QSettings:
        return QSettings(SETTINGS_ORG, SETTINGS_APP)

    def save(self) -> None:
        if not self.remember:
            return
        s = self._settings()
        s.setValue("layout/geometry", self.win.saveGeometry())
        s.setValue("layout/state", self.win.saveState(LAYOUT_VERSION))

    def restore(self) -> bool:
        """Restore a saved layout; False when there is none (or it is stale)."""
        if not self.remember:
            return False
        s = self._settings()
        geom, state = s.value("layout/geometry"), s.value("layout/state")
        if geom is None or state is None:
            return False
        self.win.restoreGeometry(geom)
        ok = self.win.restoreState(state, LAYOUT_VERSION)
        self._clamp()
        return bool(ok)

    def _clamp(self) -> None:
        screen = self.win.screen()
        if screen is None:
            return
        avail = screen.availableGeometry()
        g = self.win.frameGeometry()
        w, h = min(g.width(), avail.width()), min(g.height(), avail.height())
        x = min(max(g.x(), avail.x()), avail.right() - w + 1)
        y = min(max(g.y(), avail.y()), avail.bottom() - h + 1)
        if (w, h) != (g.width(), g.height()):
            self.win.resize(w, h)
        self.win.move(x, y)
