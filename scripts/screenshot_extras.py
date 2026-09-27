"""The smoke test's Round 7.18 steps: the PIEZO1 features ported to this GUI.

HUD (title, scale bar, the four-fold axis), a click selection drawn and
shown in the sequence window, a drag in the sequence window selecting on
every subunit, the right-click menu, a measured distance equal to the
coordinates', a screenshot that carries the HUD, a closed dock restored by
Reset layout, the guide's shortcut table, an Analyses-menu command run
in its own process and stamped, and full screen: panels hidden, F11 still
firing with the menu bar gone, and each panel restored as it was.
"""

from __future__ import annotations

import numpy as np


def _view(win, app, out) -> bool:
    from ip3r.ui.context_menu import build_context_menu
    from ip3r.ui.menus import save_screenshot
    st, hud, sel = win.scene.structure, win.hud, win.selection
    if st is None:
        return True
    if not hud.title.startswith(st.name) or hud.axis is None:
        raise RuntimeError(f"HUD title / axis: {hud.title!r}")
    if hud.scale_bar() is None:
        raise RuntimeError("no scale bar fits the view")
    ca = np.flatnonzero((st.atom_name == "CA") & ~st.hetero & (st.chain == st.chains[0]))
    i = int(ca[len(ca) // 2])
    ch, r = str(st.chain[i]), int(st.res_seq[i])
    win._picked(i)
    if sel.residues != {(ch, r)} or win.scene.scene.get("selection") is None:
        raise RuntimeError(f"a click selected {sel.residues}")
    if "selected" not in hud.readouts.get("selection", ""):
        raise RuntimeError("the HUD does not name the selection")
    win.show_sequence(chain=ch, residue=r)
    sw = win.sequence_window
    if not sw.isVisible() or r not in sw.view.selected:
        raise RuntimeError("the sequence window does not show the model's selection")
    sw.everywhere.setChecked(True)
    sw._dragged(r, r + 4)
    chains = {c for c, _ in sel.residues}
    if len(chains) < 2 or not all(r <= q <= r + 4 for _, q in sel.residues):
        raise RuntimeError(f"a drag on every subunit selected {sorted(sel.residues)[:6]}")
    menu = build_context_menu(win, i)
    texts = [a.text() for a in menu.actions()]
    if not any(t.startswith(f"Select {r} on every subunit") for t in texts):
        raise RuntimeError(f"context menu: {texts[:6]}")
    win.measure_action.trigger()
    j = int(ca[len(ca) // 2 + 10])
    win._picked(i)
    win._picked(j)
    d = float(np.linalg.norm(st.xyz[i] - st.xyz[j]))
    if not sel.distances or f"{d:.2f} Å" not in hud.readouts.get("distance", ""):
        raise RuntimeError(f"distance readout: {hud.readouts}")
    win.measure_action.trigger()
    app.processEvents()
    with_hud = save_screenshot(win, hud=True, path=str(out / "gui_selection.png"))
    bare = win.viewport.grabFramebuffer()
    from PyQt6.QtGui import QImage
    if QImage(with_hud).convertToFormat(bare.format()) == bare:
        raise RuntimeError("the screenshot does not carry the HUD")
    win.grab().save(str(out / "gui_selection_window.png"))
    sw.close()
    win.clear_selection()
    if sel.residues or sel.distances or win.scene.scene.get("selection") is not None:
        raise RuntimeError("Esc left a selection drawn")
    return False


def _layout(win, app, out) -> bool:
    dock = win.docks.docks[-1]
    dock.toggleViewAction().trigger()
    app.processEvents()
    if dock.isVisible():
        raise RuntimeError("closing the Analysis dock left it shown")
    win.docks.reset()
    app.processEvents()
    if not dock.isVisible():
        raise RuntimeError("Reset layout did not bring the Analysis dock back")
    win.show_help("Keyboard shortcuts")
    body = win.help_dialog.body.toPlainText()
    if "Ctrl+M" not in body or "F1" not in body:
        raise RuntimeError(f"guide shortcuts: {body[:120]}")
    win.help_dialog.show_topic("Channel and the lumen")
    if "u + W" not in win.help_dialog.body.toPlainText():
        raise RuntimeError("the guide does not describe the K+ energy colouring")
    win.help_dialog.close()
    return False


def _analysis_start(win, app, out) -> bool:
    from ip3r.ui.analyses import ANALYSES
    info = next(a for a in ANALYSES if a.key == "info")
    win._smoke_result = win.open_analysis(info)
    if win._smoke_result is None:
        raise RuntimeError("the Analyses menu refused 'info' with a deposit on screen")
    return False


def _analysis(win, app, out) -> bool:
    w, name = win._smoke_result, win.scene.structure.name
    if w.running:
        return True
    text = w.text()
    if f"python -m ip3r info {name}" not in text or "registered defaults" not in text:
        raise RuntimeError(f"no stamp: {text[:200]}")
    if "exit 0" not in w.state.text() or "numbering" not in text:
        raise RuntimeError(f"info failed: {w.state.text()} / {text[-300:]}")
    w.grab().save(str(out / "gui_analysis.png"))
    w.close()
    return False


def _fullscreen_start(win, app, out) -> bool:
    """Close the Structure dock first: leaving must not reopen it."""
    win.docks.docks[0].hide()
    win._smoke_fit = win.scene.fit_target
    win.fullscreen_action.trigger()
    return False


def _check_f11_binding(win):
    """The design, independent of focus: full screen hides the menu bar,
    so its shortcut must be the window's own action, in window context."""
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QKeySequence
    act = win.fullscreen_action
    if act not in win.actions() or not act.isEnabled():
        raise RuntimeError("the full-screen action is not the window's own, or is disabled")
    if QKeySequence("F11") not in act.shortcuts() or \
            act.shortcutContext() != Qt.ShortcutContext.WindowShortcut:
        raise RuntimeError(f"full screen's shortcuts {act.shortcuts()} / context")


def _fullscreen(win, app, out) -> bool:
    import sys
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QKeySequence, QShortcutEvent
    from PyQt6.QtTest import QTest
    if not win.isFullScreen():
        return True                              # the platform animates it
    tries = getattr(win, "_smoke_f11", 0)
    if tries == 0:
        app.processEvents()
        if any(d.isVisible() for d in win.docks.docks) or win.statusBar().isVisible():
            raise RuntimeError("full screen left a panel or the status bar shown")
        if win.viewport.width() < 0.95 * win.width():
            raise RuntimeError(f"the viewport fills {win.viewport.width()} of {win.width()} px")
        if "Esc" not in win.hud.readouts.get("presentation", ""):
            raise RuntimeError("no hint on how to leave full screen")
        win.grab().save(str(out / "gui_fullscreen.png"))
        _check_f11_binding(win)
    elif not win.presentation.active:
        print("  F11 left full screen as a key press", file=sys.stderr)
        return False                             # F11 left, menu bar hidden
    if tries >= 4:
        # A window shortcut fires only in the focused window, and a
        # background process cannot take focus from the app the user is in.
        # Deliver F11 as the shortcut map would, to the bound action.
        app.sendEvent(win.fullscreen_action, QShortcutEvent(QKeySequence("F11"), False))
        app.processEvents()
        if win.presentation.active:
            raise RuntimeError("F11 did not leave full screen with the menu bar hidden")
        print("  window not focused: F11 delivered to its action as a shortcut event",
              file=sys.stderr)
        return False
    win._smoke_f11 = tries + 1
    win.activateWindow()
    QTest.keyClick(win.viewport, Qt.Key.Key_F11)   # a shortcut, menu bar hidden
    app.processEvents()
    return True


def _fullscreen_left(win, app, out) -> bool:
    if win.isFullScreen():
        return True
    app.processEvents()
    structure, analysis = win.docks.docks[0], win.docks.docks[-1]
    if structure.isVisible() or not analysis.isVisible():
        raise RuntimeError("leaving full screen did not restore each panel as it was")
    if win.fullscreen_action.isChecked() or not win.statusBar().isVisible():
        raise RuntimeError("the menu check or the status bar was not restored")
    structure.show()
    return False


_STEPS = (_view, _layout, _analysis_start, _analysis, _fullscreen_start, _fullscreen,
          _fullscreen_left)
EXTRAS_STEPS = len(_STEPS)


def extras_step(k: int, win, app, out) -> bool:
    """Run step ``k``; True means "not ready yet, call again"."""
    return bool(_STEPS[k](win, app, out))
