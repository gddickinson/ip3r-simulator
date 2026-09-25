"""The lumen overlay: solve for the potential on a worker, draw the surface,
plot where the voltage falls (Round 7.10).

Like the fill (:mod:`.fill_controller`), the lumen belongs to one deposit:
a load rebuilds it if it is switched on, and a result arriving after a
newer request (another deposit, the box unticked, a session restore) is
dropped. The solve takes about ten seconds on an open deposit at the
registered 0.5 Å grid, so it never runs on the main thread.
"""

from __future__ import annotations

from ..physics.lumen_field import lumen_field
from ..render.lumen_mesh import lumen_mesh
from .workers import run_async

__all__ = ["LumenController"]


class LumenController:
    def __init__(self, scene, panel, status=None) -> None:
        self.scene, self.panel = scene, panel
        self.status = status or (lambda msg: None)
        self.field = None
        self.mesh = None
        self._token = 0
        #: The last outcome as text (for the smoke test and the status bar).
        self.message = ""
        panel.lumen_toggled.connect(self.request)

    @property
    def busy(self) -> bool:
        return self.panel.show_lumen.isChecked() and self.field is None \
            and not self.message

    def loaded(self, st) -> None:
        self.field = self.mesh = None
        self.request(self.panel.show_lumen.isChecked())

    def request(self, on: bool) -> None:
        self._token += 1
        token, st, summary = self._token, self.scene.structure, self.scene.summary
        self.field = self.mesh = None
        self.message = ""
        if not on or st is None:
            self.scene.show_lumen(None)
            self.panel.set_lumen_info("")
            return
        self.panel.set_lumen_info(f"solving for the potential in {st.name}'s "
                                  "lumen (3-D, about ten seconds)…")

        def work():
            field = lumen_field(st, summary)
            return token, field, lumen_mesh(field, summary.frame)
        run_async(work, on_done=self._built,
                  on_error=lambda e: token == self._token and self._failed(e))

    def _built(self, result) -> None:
        token, field, mesh = result
        if token != self._token:
            return
        self.field, self.mesh = field, mesh
        self.message = field.summary()
        self.scene.show_lumen(mesh)
        self.panel.show_lumen_field(field, self.scene.summary)
        self.status(self.message)

    def _failed(self, error) -> None:
        self.message = f"lumen not built: {error}"
        self.scene.show_lumen(None)
        self.panel.set_lumen_info(f"<span style='color:#e06c6c'>{self.message}</span>")
        self.status(self.message)
