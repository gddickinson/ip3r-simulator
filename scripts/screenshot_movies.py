"""The smoke test's movie steps (group ``movies``): File → Record movie…'s
recorder, driven through the real window, writes the README's movies into
``docs/img`` as animated WebP, and each is checked.

- ``movie_turntable.webp``: 6DQN turning once, by functional element;
- ``movie_transition.webp``: the 8TKG → 8TKF morph there and back,
  coloured by displacement;
- ``movie_transition_pore.webp``: the same morph with two opposite subunits
  and the ligands hidden and the camera on the pore domain, so the gate is
  seen opening;
- ``movie_mode.webp``: 8TKG's lowest collective A mode through one cycle.

The window is made compact for them (docks hidden, the viewport
:data:`VIEW` points), so the Retina capture is only halved to :data:`WIDTH`
and the HUD's text stays readable; full screen would shrink it five-fold.
Every movie must play exactly its plan's frames (a GIF stores identical
consecutive frames once, so the total duration is what is checked), must
move (the first frame differs from the one an eighth of the way in: a C4
tetramer turned by half a turn looks the same), and must leave the view as
it found it: the camera, the morph's frame, the mode animating.
"""

from __future__ import annotations

import numpy as np

import screenshot_core as core
from ip3r.parameters import PARAMETERS as _P
from ip3r.ui.movie_model import MovieSpec, delay_ms, plan

__all__ = ["STEPS"]

#: The viewport while recording, in points, and the movie's width in pixels.
VIEW, WIDTH = (720, 450), 720
#: Room around the molecule, as the camera's fit margin.
MARGIN = 1.18
_RUN: dict = {}


def _spec(source: str, **kw) -> MovieSpec:
    frames = int(_P.value("movie.turntable_frames" if source == "turntable"
                          else "movie.mode_frames"))
    return MovieSpec(source, frames=frames, fps=float(_P.value("movie.fps")),
                     hold=int(_P.value("movie.hold_frames")), max_width=WIDTH,
                     fmt="webp", **kw)


def _frames(path) -> tuple[list[np.ndarray], int]:
    """The movie's stored frames and its total duration in ms."""
    from PIL import Image
    out, total = [], 0
    with Image.open(path) as im:
        for k in range(im.n_frames):
            im.seek(k)
            out.append(np.asarray(im.convert("RGB"), dtype=np.int16))
            total += im.info["duration"]          # a WebP frame's delay is read on load
    return out, total


def _compact(win) -> bool:
    """Docks and status bar hidden, the viewport :data:`VIEW`; True once
    the viewport has that size (a resize lands a moment later)."""
    if "_docks" not in _RUN:
        _RUN["_docks"] = [(d, d.isVisible()) for d in win.docks.docks]
        for d, _ in _RUN["_docks"]:
            d.hide()
        win.statusBar().hide()
    extra = win.height() - win.viewport.height()
    if (win.viewport.width(), win.viewport.height()) != VIEW:
        win.resize(VIEW[0], VIEW[1] + extra)
        return False
    return True


def _restore_window(win) -> None:
    for d, visible in _RUN.pop("_docks", []):
        d.setVisible(visible)
    win.statusBar().show()
    win.resize(1560, 960)


def _visible_mask(sc) -> np.ndarray:
    mask = np.zeros(sc.structure.n_atoms, bool)
    mask[sc._visible_atoms()] = True
    return mask


def _frame_view(win, atoms: np.ndarray | None = None) -> None:
    """Side view, framed on ``atoms`` (all visible ones by default) with
    room around them, then left alone (no refit behind the movie's back)."""
    sc = win.scene
    sc.side_view()
    xyz = sc.view.structure.xyz
    everything = xyz[sc._visible_atoms()]
    sc.scene.camera.frame(everything if atoms is None else xyz[atoms],
                          margin=MARGIN, scene=everything)
    sc.navigated()
    win.viewport.update()


def _record(win, out, name: str, spec: MovieSpec, setup=None) -> bool:
    """Compact the window, run ``setup`` once, record once, then wait for the
    file and check it; True while any of that is still under way."""
    rec, cam = win.movies, win.scene.scene.camera
    if name not in _RUN:
        if rec.busy or not _compact(win):
            return True
        if setup is not None:
            setup()
        n_solved = len(win.morph.result.trajectory) if win.morph.result is not None else None
        _RUN[name] = {"rotation": cam.rotation.copy(), "frame": win.morph.frame,
                      "mode": win.scene.animated_mode, "n": len(plan(spec, n_solved)),
                      "delay": delay_ms(spec.fps, spec.fmt)}
        if not rec.start(spec, str(out / name)):
            raise RuntimeError(f"{name}: the recorder refused: "
                               f"{win.statusBar().currentMessage()}")
        return True
    if rec.busy:
        return True
    run, path = _RUN[name], out / name
    if not path.exists():
        raise RuntimeError(f"{name} was not written: {win.statusBar().currentMessage()}")
    frames, total = _frames(path)
    if total != run["n"] * run["delay"]:
        raise RuntimeError(f"{name}: plays {total} ms, its plan is {run['n']} frames "
                           f"of {run['delay']} ms")
    if np.abs(frames[0] - frames[len(frames) // 8]).mean() < 1.0:
        raise RuntimeError(f"{name}: nothing moves in the first eighth of the movie")
    if frames[0].shape[1] > WIDTH:
        raise RuntimeError(f"{name}: {frames[0].shape[1]} px wide, above {WIDTH}")
    if not np.allclose(cam.rotation, run["rotation"]):
        raise RuntimeError(f"{name}: the camera was not put back")
    if win.morph.frame != run["frame"] or win.scene.animated_mode != run["mode"]:
        raise RuntimeError(f"{name}: the morph frame or the mode was not put back")
    print(f"  {name}: {run['n']} frames ({len(frames)} stored), {frames[0].shape[1]} x "
          f"{frames[0].shape[0]} px, {path.stat().st_size / 1e6:.1f} MB")
    return False


def _colour(win, text: str) -> None:
    sp = win.structure_panel
    sp.color.setCurrentIndex(sp.color.findText(text))


def _turntable(win, app, out) -> bool:
    def setup():
        sp = win.structure_panel
        if win.scene.structure is None or win.scene.structure.name != "6DQN":
            raise RuntimeError("the movies group starts on 6DQN")
        sp.style.setCurrentIndex(sp.style.findText("Cartoon"))
        _colour(win, "Functional element")
        sp.site_boxes["ip3_contact"].setChecked(True)
        _frame_view(win)
    return _record(win, out, "movie_turntable.webp", _spec("turntable"), setup)


def _transition_start(win, app, out) -> bool:
    win.structure_panel.site_boxes["ip3_contact"].setChecked(False)
    core.enter_transition(win, app, out)
    return False


def _transition(win, app, out) -> bool:
    if win.morph.result is None:
        return True

    def setup():
        _colour(win, "Displacement (Transition tab)")
        _frame_view(win)
    return _record(win, out, "movie_transition.webp", _spec("transition"), setup)


def _pore(win, app, out) -> bool:
    """Subunits B and D hidden, the camera on the pore domain."""
    sp, sc, name = win.structure_panel, win.scene, "movie_transition_pore.webp"

    def setup():
        sp.ligands.setChecked(False)              # the bound lipids hide the gate
        for c in sorted(sp.chain_boxes)[1::2]:
            sp.chain_boxes[c].setChecked(False)
        lo, hi = sc.summary.span
        z = sc.summary.frame.to_frame(sc.view.structure.xyz)[:, 2]
        pore = _visible_mask(sc) & (z > lo - 5) & (z < hi + 5)
        _frame_view(win, np.flatnonzero(pore))
    if _record(win, out, name, _spec("transition"), setup):
        return True
    for b in sp.chain_boxes.values():
        b.setChecked(True)
    sp.ligands.setChecked(True)
    return False


def _check_dialog(win) -> None:
    """With a mode animating and no morph built, File → Record movie… offers
    the mode first, refuses the transition with its reason, and every
    format this machine can write."""
    from ip3r.ui.menus import leaf_actions
    from ip3r.ui.movie_dialog import MovieDialog
    if not any(a.text() == "Record movie…" for a in leaf_actions(win.menuBar())):
        raise RuntimeError("File → Record movie… is missing")
    dlg = MovieDialog(win.movies, win)
    try:
        items = dlg.source.model()
        state = {dlg.source.itemData(i): (items.item(i).isEnabled(), items.item(i).toolTip())
                 for i in range(dlg.source.count())}
        if dlg.source.currentData() != "mode" or not state["turntable"][0]:
            raise RuntimeError(f"the dialog offers {dlg.source.currentData()}: {state}")
        if state["transition"][0] or "Transition tab" not in state["transition"][1]:
            raise RuntimeError(f"the transition is offered without a morph: {state}")
        spec = dlg.spec()
        if spec.frames != _P.value("movie.mode_frames") or spec.turn:
            raise RuntimeError(f"the mode's defaults are not the parameters: {spec}")
    finally:
        dlg.deleteLater()


def _mode_start(win, app, out) -> bool:
    win.morph.reset()
    _colour(win, "Functional element")
    win.compute_modes()
    return False


def _mode(win, app, out) -> bool:
    modes = win.modes.modes
    if modes is None:
        return True
    if win.scene.animated_mode is None:
        win.modes.table.selectRow(modes.first("A"))
        return True
    if _record(win, out, "movie_mode.webp", _spec("mode"), lambda: _frame_view(win)):
        return True
    _check_dialog(win)
    win.scene.stop_animation()
    _restore_window(win)
    return False


STEPS = (_turntable, _transition_start, _transition, _pore, _mode_start, _mode)
