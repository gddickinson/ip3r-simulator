"""Records a movie of the 3-D view (File → Record movie…).

The recorder walks a :func:`~ip3r.ui.movie_model.plan` one frame per turn of
the event loop: it sets the camera and the scene state, renders, grabs the
viewport (with the HUD, or the bare framebuffer) and scales it to the
movie's width. The window stays responsive, and Esc cancels. The live
animations are paused while it records. Afterwards the view is put back as
it was (the camera, the spin, the mode animating or the morph's frame), and
the frames are encoded on a worker, since encoding a GIF blocks for seconds.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QObject, Qt, QTimer, pyqtSignal
from PyQt6.QtGui import QImage

from ..render.camera import quat_from_axis_angle, quat_multiply
from ..render.representations import Style
from .movie_model import MovieSpec, fit_size, plan, write_movie
from .workers import run_async

__all__ = ["MovieRecorder", "image_array"]


def image_array(image: QImage) -> np.ndarray:
    """A QImage as an H x W x 3 uint8 array (row padding removed)."""
    image = image.convertToFormat(QImage.Format.Format_RGB888)
    w, h = image.width(), image.height()
    buf = np.frombuffer(image.constBits().asstring(image.sizeInBytes()), np.uint8)
    return buf.reshape(h, image.bytesPerLine())[:, :3 * w].reshape(h, w, 3).copy()


class MovieRecorder(QObject):
    #: frames captured, frames in the movie
    progress = pyqtSignal(int, int)
    #: the file written
    finished = pyqtSignal(str)
    #: why it stopped
    failed = pyqtSignal(str)

    def __init__(self, win) -> None:
        super().__init__(win)
        self.win = win
        self.busy = False
        self.frames: list[np.ndarray] = []
        self._cancel = False
        self.progress.connect(lambda k, n: win.statusBar().showMessage(
            f"recording frame {k} of {n} (Esc cancels)"))
        self.finished.connect(lambda p: win.statusBar().showMessage(f"saved {p}"))
        self.failed.connect(lambda why: win.statusBar().showMessage(f"no movie: {why}"))

    # ------------------------------------------------------------ what can move

    def unavailable(self) -> dict[str, str | None]:
        """Source → None when it can be recorded now, else why not."""
        sc = self.win.scene
        none = None if sc.structure is not None else "no structure is loaded"
        return {"turntable": none,
                "mode": none or (None if sc.animated_mode is not None else
                                 "animate a mode in the Modes tab first"),
                "transition": none or (None if self.win.morph.result is not None else
                                       "build a transition in the Transition tab first")}

    # ------------------------------------------------------------------ record

    def start(self, spec: MovieSpec, path: str) -> bool:
        why = self.unavailable()[spec.source]
        if self.busy or why:
            self.failed.emit("a movie is already being recorded" if self.busy else why)
            return False
        sc, vp = self.win.scene, self.win.viewport
        morph = self.win.morph
        n_solved = len(morph.result.trajectory) if morph.result is not None else None
        self._plan = plan(spec, n_solved)
        self._spec, self._path = spec, path
        self._camera = sc.scene.camera
        self._saved = {"rotation": self._camera.rotation.copy(), "spin": vp._spin_speed,
                       "mode": sc.animated_mode, "frame": morph.frame}
        if spec.source == "mode":
            self._base, self._disp = sc._base_xyz.copy(), sc._disp
        vp.clear_animations()
        vp.set_spin(0.0)
        if spec.source != "turntable" and sc.view.style is Style.CARTOON:
            sc.view.style = Style.TUBE          # as the live animations draw it: a
            sc.fill.restyle(Style.TUBE, sc.view.visible_chains)   # cartoon per frame
            sc.view.rebuild()                   # is too slow, and frame 0 must match
        vp.setEnabled(False)                 # a drag would move the camera mid-movie
        self.busy, self._cancel, self.frames = True, False, []
        QTimer.singleShot(0, self._step)
        return True

    def cancel(self) -> None:
        self._cancel = True

    def _apply(self, frame) -> None:
        sc, cam = self.win.scene, self._camera
        q = quat_from_axis_angle(np.array([0.0, 1.0, 0.0]), np.radians(frame.angle))
        rot = quat_multiply(q, self._saved["rotation"])
        cam.rotation = rot / np.linalg.norm(rot)
        if self._spec.source == "mode":
            xyz = self._base + np.sin(frame.state) * self._disp
            sc.view.update_coords(xyz)
            sc.move_overlays(xyz)
        elif self._spec.source == "transition":
            self.win.morph.show_frame(int(frame.state), stop=False)

    def _label(self, frame) -> str:
        """What the frame shows, for the HUD: a movie names what moves."""
        if self._spec.source == "transition":
            r = self.win.morph.result
            i, n = int(frame.state), len(r.trajectory)
            return (f"{r.transition.start_id} → {r.transition.end_id}, interpolated: "
                    f"frame {i}/{n - 1}, gate {r.gate.gate[i]:.2f} Å")
        if self._spec.source == "mode":
            idx, amp = self._saved["mode"]
            modes = self.win.scene.modes[0]
            return (f"normal mode #{idx + 1} ({modes.symmetry[idx]}): shape from the "
                    f"network, {amp:.0f} Å amplitude illustrative")
        return ""

    def _grab(self) -> np.ndarray:
        vp = self.win.viewport
        vp.repaint()                          # render this frame now, not later
        image = vp.grab().toImage() if self._spec.hud else vp.grabFramebuffer()
        if self.frames:                       # the window resized mid-movie: keep
            h, w = self.frames[0].shape[:2]   # the first frame's size
        else:
            w, h = fit_size(image.width(), image.height(), self._spec.max_width,
                            even=self._spec.fmt == "mp4")
        image = image.scaled(w, h, Qt.AspectRatioMode.IgnoreAspectRatio,
                             Qt.TransformationMode.SmoothTransformation)
        return image_array(image)

    def _step(self) -> None:
        if self._cancel:
            self._restore()
            self.failed.emit("cancelled")
            return
        k = len(self.frames)
        try:
            self._apply(self._plan[k])
            if self._spec.hud:
                self.win.hud.set_readout("movie", self._label(self._plan[k]))
            self.frames.append(self._grab())
        except Exception as exc:              # never leave the view half-recorded
            self._restore()
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        self.progress.emit(k + 1, len(self._plan))
        if k + 1 < len(self._plan):
            QTimer.singleShot(0, self._step)
            return
        self._restore()
        self.win.statusBar().showMessage(f"encoding {len(self.frames)} frames…")
        self.busy = True
        run_async(write_movie, self._path, self.frames, self._spec,
                  on_done=self._written, on_error=self._encode_failed)

    def _restore(self) -> None:
        """The view as it was before recording."""
        sc, vp, s = self.win.scene, self.win.viewport, self._saved
        self.win.hud.set_readout("movie", "")
        self._camera.rotation = s["rotation"]
        vp.set_spin(s["spin"])
        vp.setEnabled(True)
        if s["mode"] is not None:                 # it was animating: animate again
            sc.animate_mode(*s["mode"])
        elif self._spec.source == "transition" and self.win.morph.result is not None:
            self.win.morph.show_frame(s["frame"], stop=False)
        vp.update()
        self.busy = False

    def _written(self, path) -> None:
        self.busy = False
        self.finished.emit(str(path))

    def _encode_failed(self, why: str) -> None:
        self.busy = False
        self.failed.emit(why.splitlines()[0])
