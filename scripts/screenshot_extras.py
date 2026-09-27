"""The smoke test's Round 7.18 steps: the PIEZO1 features ported to this GUI.

HUD (title, scale bar, the four-fold axis), a click selection drawn and
shown in the sequence window, a drag in the sequence window selecting on
every subunit, the right-click menu, a measured distance equal to the
coordinates', a screenshot that carries the HUD, a closed dock restored by
Reset layout, the guide's shortcut table, and an Analyses-menu command run
in its own process and stamped.
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


_STEPS = (_view, _layout, _analysis_start, _analysis)
EXTRAS_STEPS = len(_STEPS)


def extras_step(k: int, win, app, out) -> bool:
    """Run step ``k``; True means "not ready yet, call again"."""
    return bool(_STEPS[k](win, app, out))
