"""The main window's menus (split from ``main_window.py`` in Round 7.18,
when View gained the sequence window, measuring, the HUD and the panels,
and the Analyses and Help menus arrived).

Every shortcut here is listed in :data:`~ip3r.ui.help_content.SHORTCUTS`
(a test holds the two together).
"""

from __future__ import annotations

from PyQt6.QtGui import QAction, QKeySequence
from PyQt6.QtWidgets import QFileDialog, QInputDialog, QMessageBox

from .. import __version__
from ..config import genes_results
from .analyses import ANALYSES, GROUPS
from .help_content import DOCS
from .help_dialog import open_document

__all__ = ["build_menus", "add_action", "HUD_LABELS"]

HUD_LABELS = {"title": "Deposit title", "scale_bar": "Scale bar",
              "gnomon": "Orientation (four-fold axis)", "readouts": "Selection and distance"}


def add_action(win, menu, text, slot, shortcut=None, checkable=False, tip="") -> QAction:
    a = QAction(text, win)
    if shortcut:
        a.setShortcut(QKeySequence(shortcut))
    if checkable:
        a.setCheckable(True)
    if tip:
        a.setToolTip(tip)
    a.triggered.connect(slot)
    menu.addAction(a)
    return a


def build_menus(win) -> None:
    mb = win.menuBar()
    _file(win, mb.addMenu("&File"))
    _view(win, mb.addMenu("&View"))
    _analyses(win, mb.addMenu("&Analyses"))
    _help(win, mb.addMenu("&Help"))


def _file(win, f) -> None:
    add_action(win, f, "Open session…", win.sessions.open, "Ctrl+O")
    add_action(win, f, "Save session…", win.sessions.save, "Ctrl+Shift+S")
    f.addSeparator()
    add_action(win, f, "Fetch all registry structures", win.fetch_all)
    add_action(win, f, "Save screenshot…", lambda: save_screenshot(win, hud=True), "Ctrl+S",
               tip="The viewport as shown, HUD included.")
    add_action(win, f, "Save screenshot without the HUD…",
               lambda: save_screenshot(win, hud=False))
    f.addSeparator()
    add_action(win, f, "Quit", win.close, QKeySequence.StandardKey.Quit)


def _view(win, v) -> None:
    sc = win.scene
    add_action(win, v, "Side view (cytosol up)", sc.side_view, "Ctrl+1")
    add_action(win, v, "Top view (down the pore)", sc.top_view, "Ctrl+2")
    add_action(win, v, "IP3 site", win._site_view, "Ctrl+3")
    add_action(win, v, "Fit to view", sc.fit_view, "Ctrl+0")
    add_action(win, v, "Toggle spin", lambda: win.viewport.set_spin(
        0.0 if win.viewport._spin_speed else 20.0), "Space")
    v.addSeparator()
    add_action(win, v, "Sequence…", lambda: win.show_sequence(), "Ctrl+Shift+Q")
    win.measure_action = add_action(
        win, v, "Measure distances", win.selection.arm, "Ctrl+M", checkable=True,
        tip="Every two clicked atoms close a distance.")
    add_action(win, v, "Clear selection and distances", win.clear_selection, "Esc")
    v.addSeparator()
    hud = v.addMenu("HUD")
    win.hud_actions = {}
    for key, label in HUD_LABELS.items():
        a = add_action(win, hud, label, lambda on, k=key: win.set_hud(k, on), checkable=True)
        a.setChecked(getattr(win.hud.settings, key))
        win.hud_actions[key] = a
    panels = v.addMenu("Panels")
    for a in win.docks.panel_actions():
        panels.addAction(a)
    add_action(win, v, "Reset layout", win.docks.reset)


def _analyses(win, menu) -> None:
    menu.setToolTipsVisible(True)
    win.analysis_actions = {}
    for group in GROUPS:
        sub = menu.addMenu(group)
        sub.setToolTipsVisible(True)
        for a in ANALYSES:
            if a.group == group:
                act = add_action(win, sub, f"{a.label}…", lambda _=False, a=a: win.open_analysis(a),
                                 tip=f"{a.about} Takes {a.duration}.")
                win.analysis_actions[a.key] = act
    menu.addSeparator()
    add_action(win, menu, "Command…", lambda: _custom(win),
               tip="Any python -m ip3r command, run the same way.")


def _custom(win) -> None:
    text, ok = QInputDialog.getText(win, "Run a command", "python -m ip3r",
                                    text="info " + (win.scene.structure.name
                                                    if win.scene.structure else "6DQN"))
    if ok and text.strip():
        win.open_command(text.split(), title=text.strip())


def _help(win, h) -> None:
    add_action(win, h, "Guide", win.show_help, "F1")
    docs = h.addMenu("Documents")
    for label, path in DOCS:
        add_action(win, docs, label, lambda _=False, p=path: open_document(p) or
                   win.statusBar().showMessage(f"{p} is not in this checkout"))
    h.addSeparator()
    add_action(win, h, "Parameters…", win.edit_parameters, "Ctrl+Shift+P")
    add_action(win, h, "About", lambda: _about(win))


def save_screenshot(win, hud: bool = True, path: str | None = None) -> str | None:
    """The viewport to a PNG; with ``hud``, the overlay's scale bar and
    labels are in it (a figure without a scale states nothing)."""
    if path is None:
        path, _ = QFileDialog.getSaveFileName(win, "Save screenshot", "ip3r.png", "PNG (*.png)")
    if not path:
        return None
    image = win.viewport.grab() if hud else win.viewport.grabFramebuffer()
    image.save(path)
    win.statusBar().showMessage(f"saved {path}")
    return path


def _about(win) -> None:
    QMessageBox.about(win, "About", (
        f"<b>IP3R Structural Simulator {__version__}</b><p>Physics-driven "
        "3-D model of the IP3 receptor, and a re-derivation of the "
        f"ip3r_genes results found at<br><code>{genes_results()}</code>.</p>"
        "<p>Ported from the PIEZO1 simulator. Help → Guide (F1) walks "
        "through the interface; README.md and INTERFACE.md map the code.</p>"))
