"""File-menu session save/open, and restoring one onto the window.

A session is the view, never the results (see :mod:`ip3r.io.session`).
Restoring is asynchronous in three places: the structure loads on a worker,
a transition is rebuilt on a worker, and a mode animation needs the modes
computed on a worker. So the controller holds the session as *pending* and
finishes it from :meth:`loaded`, :meth:`transition_built` and
:meth:`modes_computed`, which the main window calls, in that order: a mode
animation starts only after the transition's frame is shown, because
showing a frame stops it. A pending restore is
dropped if a different deposit arrives first. Style and camera applied to
whatever happened to be loaded would be a valid-looking wrong view.

Parameters: a session saved under edited values says so, and restoring asks
whether to apply them. Applying them always re-measures the deposit, even if
it is the one on screen, because the channel summary it shows was made under
the old values.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PyQt6.QtWidgets import QFileDialog, QMessageBox

from ..config import DATA_DIR
from ..io.registry import get_entry
from ..io.session import Session, load_session, parameter_differences, save_session
from ..parameters import PARAMETERS
from ..render.representations import ColorBy, Style

__all__ = ["SessionController"]


class SessionController:
    def __init__(self, window) -> None:
        self.win = window
        self.path: Path | None = None
        self.pending: Session | None = None
        self._frame: tuple[str, int] | None = None     # (end id, frame) awaited
        self._mode: dict | None = None                 # {index, amplitude} awaited
        self._modes_for: str | None = None             # deposit the ANM runs on
        self.notes: list[str] = []

    # ------------------------------------------------------------ capture

    def capture(self) -> Session:
        win, sp = self.win, self.win.structure_panel
        st = win.scene.structure
        cam = win.viewport.scene.camera if win.viewport.scene is not None else None
        chains = sp.visible_chains()
        s = Session(
            structure=st.name if st is not None else "",
            n_atoms=int(st.n_atoms) if st is not None else 0,
            style=sp.current_style().value, color_by=sp.current_color().value,
            layer=sp.current_layer() or "", show_ligands=sp.ligands.isChecked(),
            visible_chains=([] if chains is None or len(chains) == len(sp.chain_boxes)
                            else sorted(chains)),
            sites=[k for k, b in sp.site_boxes.items() if b.isChecked()],
            show_pore=win.channel.show_pore.isChecked(),
            completeness=sp.current_completeness(),
            tab=win.tabs.tabText(win.tabs.currentIndex()),
            dynamics=win.dynamics.view_state(), variants=win.variants.view_state(),
            parameters=PARAMETERS.overrides())
        if win.scene.animated_mode is not None:
            index, amplitude = win.scene.animated_mode
            s.modes = {"index": index, "amplitude": amplitude}
        if cam is not None:
            s.camera_rotation = [float(v) for v in cam.rotation]
            s.camera_pivot = [float(v) for v in cam.pivot]
            s.camera_distance = float(cam.distance)
            s.camera_pan = [float(v) for v in cam.pan]
            s.camera_slab = -1.0 if cam.slab_front is None else float(cam.slab_front)
            s.orthographic = bool(cam.orthographic)
        r = win.morph.result
        if r is not None:
            s.transition = {"end": r.transition.end_id, "fit": r.transition.fit,
                            "method": r.trajectory.method, "frame": int(win.morph.frame),
                            "paint": win.transition.paint.isChecked()}
        return s

    # -------------------------------------------------------------- files

    def save(self) -> None:
        if self.win.scene.structure is None:
            return self._status("nothing to save — load a structure first")
        start = str(self.path or (DATA_DIR / "sessions" / "session.json"))
        name, _ = QFileDialog.getSaveFileName(self.win, "Save session", start,
                                              "Session (*.json)")
        if name:
            self.save_to(Path(name))

    def save_to(self, path: Path) -> Session:
        s = self.capture()
        save_session(s, path)
        self.path = path
        self._status(f"session saved to {path} — {s.describe()}")
        return s

    def open(self) -> None:
        start = str(self.path or (DATA_DIR / "sessions"))
        name, _ = QFileDialog.getOpenFileName(self.win, "Open session", start,
                                              "Session (*.json)")
        if not name:
            return
        try:
            session = load_session(name)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self.win, "Could not read that session", str(exc))
            return
        self.path = Path(name)
        self.apply(session)

    # ------------------------------------------------------------ restore

    def apply(self, session: Session, parameters: str | None = None) -> bool:
        """Start restoring ``session``. ``parameters`` is ``"apply"`` or
        ``"keep"`` when the saved set differs from the current one; ``None``
        asks. Returns False if nothing was started (cancelled or unknown)."""
        win = self.win
        if get_entry(session.structure) is None:
            QMessageBox.warning(win, "Session not restored",
                                f"{session.structure or 'No structure'} is not "
                                "in the structure registry.")
            return False
        self.notes, self._frame, self._mode = [], None, None
        diffs = parameter_differences(session.parameters, PARAMETERS.overrides())
        changed = False
        if diffs:
            choice = parameters or self._ask(diffs)
            if choice is None:
                return False
            if choice == "apply":
                unknown = PARAMETERS.replace(session.parameters)
                changed = True
                if unknown:
                    self.notes.append(f"unregistered parameters ignored: "
                                      f"{', '.join(unknown)}")
            else:
                self.notes.append(f"kept current parameters ({len(diffs)} "
                                  "differ from the session's) — numbers may differ")
        self.pending = session
        st = win.scene.structure
        if st is not None and st.name == session.structure and not changed:
            self._restore_view(st)
        else:
            win.structure_panel.select(session.structure)   # list follows the load
        return True

    def loaded(self, st) -> None:
        """Called by the window after every load has been drawn."""
        s, self.pending = self.pending, None
        self._frame = self._mode = None
        if s is None:
            return
        if s.structure != st.name:
            return self._status(f"session restore abandoned: {st.name} was "
                                f"loaded instead of {s.structure}")
        self._restore_view(st, s)

    def transition_built(self, result) -> None:
        if self._frame is None or result.transition.end_id != self._frame[0]:
            return
        frame, self._frame = self._frame[1], None
        self.win.morph.show_frame(frame)
        self.win.transition.follow(self.win.morph.frame)
        self._restore_mode()

    def modes_computed(self) -> None:
        """Called by the window when an ANM finishes (for any reason)."""
        st = self.win.scene.structure
        if self._mode is None or st is None or st.name != self._modes_for:
            return
        m, self._mode = self._mode, None
        if not self.win.modes.select(m["index"], m["amplitude"]):
            self.notes.append(f"mode #{m['index'] + 1} not in this network; "
                              "not animated")
        self._finish()

    def _restore_view(self, st, s: Session | None = None) -> None:
        s = s or self.pending
        self.pending = None
        win, sp = self.win, self.win.structure_panel
        if s.n_atoms and s.n_atoms != st.n_atoms:
            self.notes.append(f"{st.name} now has {st.n_atoms:,} atoms, the session "
                              f"{s.n_atoms:,}: the file changed; check the camera")
        combos = [(sp.style, _enum(Style, s.style), "style"),
                  (sp.color, _enum(ColorBy, s.color_by), "colour"),
                  (sp.layer, s.layer or None, "layer"),
                  (sp.completeness, s.completeness, "completeness")]
        boxes = [sp.ligands, *sp.chain_boxes.values(), *sp.site_boxes.values()]
        for w in [c for c, _, _ in combos] + boxes:
            w.blockSignals(True)
        try:
            for combo, value, what in combos:
                i = combo.findData(value) if value is not None else -1
                if i >= 0:
                    combo.setCurrentIndex(i)
                elif value is not None or what != "layer":
                    self.notes.append(f"unknown {what} kept as it was")
            sp.ligands.setChecked(s.show_ligands)
            for c, b in sp.chain_boxes.items():
                b.setChecked(not s.visible_chains or c in s.visible_chains)
            for k, b in sp.site_boxes.items():
                b.setChecked(k in s.sites)
        finally:
            for w in [c for c, _, _ in combos] + boxes:
                w.blockSignals(False)
        missing = sorted(set(s.sites) - set(sp.site_boxes))
        if missing:
            self.notes.append(f"unknown site sets ignored: {', '.join(missing)}")
        sp._restyle()                    # legend, layer enablement, one redraw
        win._apply_sites()
        win.channel.show_pore.blockSignals(True)
        win.channel.show_pore.setChecked(s.show_pore)
        win.channel.show_pore.blockSignals(False)
        win.scene.show_pore(s.show_pore)
        sp._update_legend()
        win.fills.request(sp.current_completeness())
        self._set_camera(s)
        self.notes += win.dynamics.restore(s.dynamics) + win.variants.restore(s.variants)
        win.scene.stop_animation()          # a running mode belongs to the old view
        for i in range(win.tabs.count()):
            if win.tabs.tabText(i) == s.tab:
                win.tabs.setCurrentIndex(i)
        self._mode = dict(s.modes) or None
        self._restore_transition(s)

    def _set_camera(self, s: Session) -> None:
        if self.win.viewport.scene is None:
            return
        cam = self.win.viewport.scene.camera
        cam.rotation = np.array(s.camera_rotation, dtype=float)
        cam.pivot = np.array(s.camera_pivot, dtype=float)
        cam.distance = float(s.camera_distance)
        cam.pan = np.array(s.camera_pan, dtype=float)
        cam.slab_front = None if s.camera_slab < 0 else float(s.camera_slab)
        cam.orthographic = bool(s.orthographic)
        self.win.scene.fit_target = None      # a saved camera, not a fit to redo on resize
        self.win.viewport.update()

    def _restore_transition(self, s: Session) -> None:
        t = s.transition
        if not t:
            return self._restore_mode()
        tp = self.win.transition
        if tp.end.findData(t["end"]) < 0:
            self.notes.append(f"transition to {t['end']} not available from "
                              f"{s.structure}; not rebuilt")
            return self._restore_mode()
        for combo, key in ((tp.end, "end"), (tp.fit, "fit"), (tp.method, "method")):
            i = combo.findData(t.get(key))
            if i >= 0:
                combo.setCurrentIndex(i)
        tp.paint.blockSignals(True)
        tp.paint.setChecked(bool(t.get("paint")))
        tp.paint.blockSignals(False)
        self._frame = (t["end"], int(t.get("frame", 0)))
        self._status(f"session: rebuilding {s.structure} → {t['end']}…")
        self.win.build_transition(t["end"], tp.fit.currentData(), tp.method.currentData())

    def _restore_mode(self) -> None:
        """Last step: re-select the animating mode, computing modes first
        unless this deposit's are already in the panel."""
        if self._mode is None:
            return self._finish()
        st = self.win.scene.structure
        self._modes_for = st.name
        if self.win.modes.modes is not None:
            return self.modes_computed()
        self._status(f"session: computing {st.name}'s modes to animate "
                     f"#{self._mode['index'] + 1}…")
        self.win.compute_modes()

    # ------------------------------------------------------------ helpers

    def _finish(self) -> None:
        self._status("session restored" + ("; " + "; ".join(self.notes)
                                           if self.notes else ""))

    def _status(self, text: str) -> None:
        self.win.statusBar().showMessage(text)

    def _ask(self, diffs: dict) -> str | None:
        rows = "".join(f"<br>{k}: session {_fmt(a)}, now {_fmt(b)}"
                       for k, (a, b) in list(diffs.items())[:12])
        more = f"<br>… and {len(diffs) - 12} more" if len(diffs) > 12 else ""
        box = QMessageBox(QMessageBox.Icon.Question, "Session parameters",
                          "This session was saved under a different parameter "
                          f"set:{rows}{more}<br><br>Apply the session's set "
                          "(the structure is re-measured), or keep the current "
                          "one?", parent=self.win)
        apply_btn = box.addButton("Apply session's", QMessageBox.ButtonRole.AcceptRole)
        keep_btn = box.addButton("Keep current", QMessageBox.ButtonRole.RejectRole)
        box.addButton(QMessageBox.StandardButton.Cancel)
        box.exec()
        return ("apply" if box.clickedButton() is apply_btn else
                "keep" if box.clickedButton() is keep_btn else None)


def _enum(kind, value: str):
    try:
        return kind(value)
    except ValueError:
        return None


def _fmt(v) -> str:
    return "default" if v is None else f"{v:g}"
