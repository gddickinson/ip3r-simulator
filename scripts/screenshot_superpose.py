"""The smoke test's Round 7.20 steps: a second deposit drawn superposed.

8TKF is drawn on 8TKG by the Structure panel's Superpose selector and must
equal the headless :func:`ip3r.structure.superpose.superpose` atom for atom.
The legend names it, a session holds the choice, and the overlay follows a
hidden subunit. Then 8TKF itself is loaded: the choice (8TKF on 8TKF) is no
longer a candidate and is dropped, while 8TKG, now the candidate, is fitted
the other way round. Clearing the selector removes every overlay batch.
"""

from __future__ import annotations

import numpy as np

_BATCH = "overlay:ribbon"


def _start(win, app, out) -> None:
    sp = win.structure_panel
    sp.color.setCurrentIndex(sp.color.findText("Uniform"))
    sp.select("8TKG")


def _drawn(win, app, out) -> bool:
    from ip3r.io.loader import load
    from ip3r.structure.superpose import superpose
    sp, sc, ctl = win.structure_panel, win.scene, win.superposed
    if sc.structure is None or sc.structure.name != "8TKG":
        return True
    if sp.current_superpose()[0] != "8TKF":
        if not sp.set_superpose("8TKF", "pore"):
            raise RuntimeError("8TKF is not offered on 8TKG")
        return True
    if ctl.result is None:
        if ctl.message.startswith("not superposed"):
            raise RuntimeError(ctl.message)
        return True
    ref = superpose(sc.structure, load("8TKF"), "pore")
    if not np.allclose(ctl.view.structure.xyz, ref.structure.xyz, atol=1e-3):
        raise RuntimeError("the drawn overlay is not the headless superposition")
    if sc.scene.get(_BATCH) is None or "8TKF, superposed" not in sp.legend.text():
        raise RuntimeError("overlay not drawn, or not in the legend")
    if "8TKF on 8TKG (pore fit" not in sp.superpose_info.text():
        raise RuntimeError(f"no fit summary: {sp.superpose_info.text()}")
    held = win.sessions.capture().superpose
    if held != {"pdb": "8TKF", "fit": "pore"}:
        raise RuntimeError(f"session holds superpose {held!r}")
    app.processEvents()
    win.grab().save(str(out / "gui_superpose.png"))
    box = sp.chain_boxes[sorted(sp.chain_boxes)[0]]
    box.setChecked(False)
    if ctl.view.visible_chains != sp.visible_chains():
        raise RuntimeError("the overlay did not follow a hidden subunit")
    box.setChecked(True)
    print(f"  {ctl.message}")
    sp.select("8TKF")
    return False


def _reloaded(win, app, out) -> bool:
    sp, sc = win.structure_panel, win.scene
    if sc.structure is None or sc.structure.name != "8TKF":
        return True
    if sp.current_superpose()[0]:
        raise RuntimeError("8TKF kept as its own overlay")
    if sc.scene.get(_BATCH) is not None:
        raise RuntimeError("the old overlay survived the new deposit")
    sp.set_superpose("8TKG", "global")
    return False


def _reversed(win, app, out) -> bool:
    sp, sc, ctl = win.structure_panel, win.scene, win.superposed
    if ctl.result is None:
        if ctl.message.startswith("not superposed"):
            raise RuntimeError(ctl.message)
        return True
    if (ctl.result.shown_id, ctl.result.other_id, ctl.result.fit) != ("8TKF", "8TKG", "global"):
        raise RuntimeError(f"wrong pair: {ctl.message}")
    sp.set_superpose("", "global")
    if sc.scene.get(_BATCH) is not None or ctl.view is not None:
        raise RuntimeError("clearing the selector left the overlay drawn")
    sp.color.setCurrentIndex(sp.color.findText("Functional element"))
    return False


STEPS = (_start, _drawn, _reloaded, _reversed)
