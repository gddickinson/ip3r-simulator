"""The lumen overlay: solve for the potential on a worker, draw the surface,
plot where the voltage falls (Rounds 7.10, 7.12, 7.16).

Like the fill (:mod:`.fill_controller`), the lumen belongs to one deposit:
a load rebuilds it if it is switched on, and a result arriving after a
newer request (another deposit, the box unticked, a session restore) is
dropped. The neutral solve takes about ten seconds on an open deposit at the
registered 0.5 Å grid, and a wall charge 10–30 s more (the dielectric
closure, solved over the whole box, longer), so neither runs on
the main thread. A new charge choice re-solves only the charge (the neutral
field is kept, unless the image cost changes the grid); a new colouring
only recolours the surface.

Round 7.16's image cost is Round 7.15's W on a 1 Å grid
(``born.lumen_spacing``): ticking it re-cuts the neutral field there too.
Its first solve on a deposit is a quarter hour on every core, on the
worker, and it is cached on disk after that.

Round 7.24's steady state at reversal replaces the whole reading: its own
lumen (the reversal grid's electrostatic volume) and surface, each ion's
concentration and drop; any change to it solves again from scratch.
Round 7.26: one of Round 7.25's candidate walls is read the same way, and
the deposit's own reading of the same experiment last drawn (``own``) is
plotted beside it.
"""

from __future__ import annotations

import numpy as np

from ..physics.dielectric3d import DIELECTRIC
from ..physics.lumen_charge import charged_lumen
from ..physics.lumen_field import lumen_field
from ..physics.lumen_reversal import reversal_lumen
from ..render.lumen_mesh import (conc_colors, image_colors, lumen_mesh,
                                 reversal_key, wall_colors)
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
        #: Round 7.24's reading at reversal drawn, or None.
        self.reversal = None
        #: The deposit's own wall at reversal last drawn on this deposit
        #: (Round 7.26: set beside a candidate wall of the same experiment).
        self.own = None
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
        if self.box.reversal is not None:
            return self.reversal is None
        return self.field is None or (self.box.closure is not None
                                      and self.charged is None)

    def loaded(self, st) -> None:
        self.field = self.mesh = self.charged = self.reversal = self.own = None
        self.request(self.box.show.isChecked())

    def request(self, on: bool) -> None:
        """Solve everything again (a load, a toggle, a session restore)."""
        self._token += 1
        token, st, summary = self._token, self.scene.structure, self.scene.summary
        self.field = self.mesh = self.charged = self.reversal = None
        self.message = ""
        if not on or st is None:
            self.scene.show_lumen(None)
            self.panel.set_lumen_info("")
            return
        if self.box.reversal is not None:
            return self._request_reversal(token, st, summary)
        closure, pairs = self.box.closure, self.box.pairs.isChecked()
        image, spacing = self.box.image, self.box.spacing
        self.panel.set_lumen_info(
            f"solving for the potential in {st.name}'s lumen (3-D, about ten "
            "seconds" + ("" if not closure else self._wait(closure, image))
            + ")…")

        def work():
            field = lumen_field(st, summary, spacing=spacing)
            mesh = lumen_mesh(field, summary.frame)
            charged = (charged_lumen(st, field, closure, summary, pairs,
                                     image=image)
                       if closure and field.conducts else None)
            return token, field, mesh, charged
        run_async(work, on_done=self._built,
                  on_error=lambda e: token == self._token and self._failed(e))

    def _request_reversal(self, token, st, summary) -> None:
        reading, experiment = self.box.reversal, self.box.experiment.currentData()
        self.panel.set_lumen_info(
            f"solving {st.name}'s {experiment} experiment ({reading}) to its "
            "reversal in 3-D (about a minute; pb + csc longer)…")

        def work():
            rev = reversal_lumen(st, reading, experiment, summary)
            return token, rev.lumen, lumen_mesh(rev.lumen, summary.frame), None, rev
        run_async(work, on_done=self._built,
                  on_error=lambda e: token == self._token and self._failed(e))

    def recharge(self) -> None:
        """The wall charge changed: keep the neutral field, re-solve the charge.
        Into, out of or within a reversal reading, everything is solved again."""
        if self.box.reversal is not None or self.reversal is not None:
            return self.request(self.box.show.isChecked())
        if self.field is None:
            return self.request(self.box.show.isChecked())
        self._token += 1
        token, st, summary = self._token, self.scene.structure, self.scene.summary
        field, mesh = self.field, self.mesh
        closure, pairs = self.box.closure, self.box.pairs.isChecked()
        image = self.box.image
        if not _same_grid(field, self.box.spacing):    # the image's grid
            return self.request(True)
        self.charged = None
        self.message = ""
        if closure is None or not field.conducts:
            return self._built((token, field, mesh, None))
        self.panel.set_lumen_info(f"placing {st.name}'s wall charge ({closure}"
                                  f"{' + image' if image else ''}) on the lumen "
                                  f"({self._wait(closure, image).lstrip(', ')})…")
        run_async(lambda: (token, field, mesh,
                           charged_lumen(st, field, closure, summary, pairs,
                                         image=image)),
                  on_done=self._built,
                  on_error=lambda e: token == self._token and self._failed(e))

    @staticmethod
    def _wait(closure: str, image: bool) -> str:
        if image:
            return (", and the wall charge with the image cost a minute more "
                    "if the deposit's W is cached, else about a quarter hour "
                    "on every core")
        return (", and the wall charge two or three minutes more"
                if closure == DIELECTRIC else ", and the wall charge half a "
                "minute more")

    def _built(self, result) -> None:
        token, field, mesh, charged, *rev = result
        if token != self._token:
            return
        self.field, self.mesh, self.charged = field, mesh, charged
        self.reversal = rev[0] if rev else None
        self.message = (self.reversal or charged or field).summary()
        self.recolour()
        if self.reversal is not None:
            beside = None
            if not self.reversal.candidate:
                self.own = self.reversal
            elif (self.own is not None
                  and self.own.experiment == self.reversal.experiment):
                beside = self.own
            self.panel.show_lumen_reversal(self.reversal, self.scene.summary,
                                           beside)
        else:
            self.panel.show_lumen_field(field, self.scene.summary, charged)
        self.status(self.message)

    def values(self) -> np.ndarray | None:
        """What the surface shows at each vertex, per the colour choice."""
        if self.mesh is None:
            return None
        grey = np.full(len(self.mesh.positions), np.nan)
        rk = reversal_key(self.box.colouring)
        if rk is not None:                    # grey unless drawn and carried
            kind, ion = rk
            rev = self.reversal
            if rev is None or ion not in rev.conc:
                return grey
            return self.mesh.sample(rev.conc[ion] if kind == "conc"
                                    else rev.drop[ion])
        if self.reversal is not None:         # the neutral drop on its grid
            return self.mesh.phi if self.box.colouring == "drop" else grey
        if self.box.colouring == "image":     # grey unless W was solved
            w = None if self.charged is None else self.charged.w
            return grey if w is None else self.mesh.sample(w)
        if self.box.colouring in ("wall", "energy"):
            if self.charged is None:          # the neutral pore: no wall potential
                return np.zeros(len(self.mesh.positions))
            return self.mesh.sample(self.charged.u if self.box.colouring == "wall"
                                    else self.charged.k_energy)
        if self.charged is None:
            return self.mesh.phi
        return self.mesh.sample(self.charged.mu)

    def recolour(self) -> None:
        if self.mesh is None:
            self.scene.show_lumen(None)
            return
        v = self.values()
        rk = reversal_key(self.box.colouring)
        paint = (conc_colors if rk and rk[0] == "conc" else
                 {"wall": wall_colors, "energy": wall_colors,
                  "image": image_colors}.get(self.box.colouring, ramp))
        self.mesh.colors = paint(v)
        self.scene.show_lumen(self.mesh)

    def _failed(self, error) -> None:
        self.message = f"lumen not built: {error}"
        self.scene.show_lumen(None)
        self.panel.set_lumen_info(f"<span style='color:#e06c6c'>{self.message}</span>")
        self.status(self.message)


def _same_grid(field, spacing: float | None) -> bool:
    """``field`` is cut on ``spacing`` (None: the registered grid)."""
    from ..parameters import PARAMETERS as _P
    want = _P.value("pore3d.spacing") if spacing is None else spacing
    return abs(field.volume.spacing - want) < 1e-9
