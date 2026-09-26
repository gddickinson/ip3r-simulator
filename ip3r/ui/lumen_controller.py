"""The lumen overlay: solve for the potential on a worker, draw the surface,
plot where the voltage falls (Rounds 7.10, 7.12).

Like the fill (:mod:`.fill_controller`), the lumen belongs to one deposit:
a load rebuilds it if it is switched on, and a result arriving after a
newer request (another deposit, the box unticked, a session restore) is
dropped. The neutral solve takes about ten seconds on an open deposit at the
registered 0.5 Å grid, and a wall charge 10–30 s more, so neither runs on
the main thread. A new charge choice re-solves only the charge (the neutral
field is kept); a new colouring only recolours the surface.
"""

from __future__ import annotations

import numpy as np

from ..physics.lumen_charge import charged_lumen
from ..physics.lumen_field import lumen_field
from ..render.lumen_mesh import lumen_mesh, wall_colors
from ..render.colormaps import ramp
from .workers import run_async

__all__ = ["LumenController"]


class LumenController:
    def __init__(self, scene, panel, status=None) -> None:
        self.scene, self.panel = scene, panel
        self.box = panel.lumen_box
        self.status = status or (lambda msg: None)
        self.field = None
        self.mesh = None
        #: The wall-charge reading drawn, or None (the neutral pore).
        self.charged = None
        self._token = 0
        #: The last outcome as text (for the smoke test and the status bar).
        self.message = ""
        panel.lumen_toggled.connect(self.request)
        self.box.charge_changed.connect(self.recharge)
        self.box.colour_changed.connect(self.recolour)

    @property
    def busy(self) -> bool:
        if not self.box.show.isChecked() or self.message.startswith("lumen not"):
            return False
        return self.field is None or (self.box.closure is not None
                                      and self.charged is None)

    def loaded(self, st) -> None:
        self.field = self.mesh = self.charged = None
        self.request(self.box.show.isChecked())

    def request(self, on: bool) -> None:
        """Solve everything again (a load, a toggle, a session restore)."""
        self._token += 1
        token, st, summary = self._token, self.scene.structure, self.scene.summary
        self.field = self.mesh = self.charged = None
        self.message = ""
        if not on or st is None:
            self.scene.show_lumen(None)
            self.panel.set_lumen_info("")
            return
        closure, pairs = self.box.closure, self.box.pairs.isChecked()
        self.panel.set_lumen_info(
            f"solving for the potential in {st.name}'s lumen (3-D, about ten "
            "seconds" + (", and the wall charge half a minute more" if closure
                         else "") + ")…")

        def work():
            field = lumen_field(st, summary)
            mesh = lumen_mesh(field, summary.frame)
            charged = (charged_lumen(st, field, closure, summary, pairs)
                       if closure and field.conducts else None)
            return token, field, mesh, charged
        run_async(work, on_done=self._built,
                  on_error=lambda e: token == self._token and self._failed(e))

    def recharge(self) -> None:
        """The wall charge changed: keep the neutral field, re-solve the charge."""
        if self.field is None:
            return self.request(self.box.show.isChecked())
        self._token += 1
        token, st, summary = self._token, self.scene.structure, self.scene.summary
        field, mesh = self.field, self.mesh
        closure, pairs = self.box.closure, self.box.pairs.isChecked()
        self.charged = None
        self.message = ""
        if closure is None or not field.conducts:
            return self._built((token, field, mesh, None))
        self.panel.set_lumen_info(f"placing {st.name}'s wall charge ({closure}) "
                                  "on the lumen (10–30 seconds)…")
        run_async(lambda: (token, field, mesh,
                           charged_lumen(st, field, closure, summary, pairs)),
                  on_done=self._built,
                  on_error=lambda e: token == self._token and self._failed(e))

    def _built(self, result) -> None:
        token, field, mesh, charged = result
        if token != self._token:
            return
        self.field, self.mesh, self.charged = field, mesh, charged
        self.message = (charged or field).summary()
        self.recolour()
        self.panel.show_lumen_field(field, self.scene.summary, charged)
        self.status(self.message)

    def values(self) -> np.ndarray | None:
        """What the surface shows at each vertex, per the colour choice."""
        if self.mesh is None:
            return None
        if self.box.colouring == "wall":
            if self.charged is None:          # the neutral pore: no wall potential
                return np.zeros(len(self.mesh.positions))
            return self.mesh.sample(self.charged.u)
        if self.charged is None:
            return self.mesh.phi
        return self.mesh.sample(self.charged.mu)

    def recolour(self) -> None:
        if self.mesh is None:
            self.scene.show_lumen(None)
            return
        v = self.values()
        self.mesh.colors = wall_colors(v) if self.box.colouring == "wall" else ramp(v)
        self.scene.show_lumen(self.mesh)

    def _failed(self, error) -> None:
        self.message = f"lumen not built: {error}"
        self.scene.show_lumen(None)
        self.panel.set_lumen_info(f"<span style='color:#e06c6c'>{self.message}</span>")
        self.status(self.message)
