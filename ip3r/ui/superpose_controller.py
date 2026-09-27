"""The Superpose selector: a second deposit drawn on the one shown (Round 7.20).

The second deposit is loaded and fitted on a worker
(:func:`~ip3r.structure.superpose.superpose`, the transition's own fit), then
drawn as its own :class:`~ip3r.render.representations.MolecularView`
(batches ``overlay:*``) in one colour, :data:`~ip3r.render.colormaps.SUPERPOSE_COLOR`,
whatever the shown deposit is coloured by. It follows the shown deposit's
style, ligands and visible subunits (its chains are relabelled to the matched
ones). It is placed on the deposit *as deposited* and does not move with a
morph or mode frame: it is a fixed reference, which a morph toward it can be
watched approaching.

Stale results are the one hazard: a fit belongs to one pair, so a new deposit
refits (keeping the choice while it is still a candidate), and a result
arriving after a newer request is dropped.
"""

from __future__ import annotations

from ..io import loader
from ..render.colormaps import SUPERPOSE_COLOR
from ..render.representations import ColorBy, MolecularView
from ..structure.superpose import candidates, superpose
from ..structure.transition import TransitionUnavailable
from .workers import run_async

__all__ = ["SuperposeController"]


class SuperposeController:
    def __init__(self, scene, panel, style, status=None) -> None:
        self.scene, self.panel = scene, panel
        self.style = style                      # () -> the shown deposit's style kwargs
        self.status = status or (lambda msg: None)
        self.result = None                      # the drawn Superposition
        self.view: MolecularView | None = None
        self._token = 0
        #: The last outcome as text (for the smoke test and the status bar).
        self.message = ""
        panel.superpose_changed.connect(self.request)
        scene.on_clear.append(self._drop)

    def loaded(self, st, paralog: str | None) -> None:
        """A new deposit is shown: offer its candidates and refit the choice."""
        self.panel.set_superpose_candidates(candidates(st.name, paralog))
        self.request()

    def request(self) -> None:
        self._token += 1
        token, st = self._token, self.scene.structure
        other, fit = self.panel.current_superpose()
        if not other or st is None:
            self._clear()
            self.message = ""
            self.panel.set_superpose_info("")
            return
        self.panel.set_superpose_info(f"fitting {other} onto {st.name}…")
        loader.ALLOW_FETCH = True

        def work():
            try:
                return token, superpose(st, loader.load(other), fit)
            except TransitionUnavailable as exc:
                return token, str(exc)
        run_async(work, on_done=self._built,
                  on_error=lambda e: token == self._token
                  and self._refuse(f"superposition failed: {e}"))

    def _built(self, result) -> None:
        token, sup = result
        if token != self._token:
            return                  # superseded: another deposit or another choice
        if isinstance(sup, str):
            self._refuse(sup)
            return
        self._clear()
        self.result, self.message = sup, sup.summary()
        kw = self.style()
        self.view = MolecularView(
            self.scene.scene, sup.structure, name="overlay", style=kw["style"],
            color_by=ColorBy.UNIFORM, uniform_rgb=SUPERPOSE_COLOR,
            show_ligands=kw["show_ligands"], visible_chains=kw["visible_chains"])
        self.view.rebuild()
        self.panel.set_superpose_info(self.message)
        self.status(self.message)
        self.scene.viewport.update()

    def restyle(self) -> None:
        """Follow the shown deposit's style, ligands and subunits."""
        if self.view is None:
            return
        kw = self.style()
        chains_changed = kw["visible_chains"] != self.view.visible_chains
        self.view.style, self.view.show_ligands = kw["style"], kw["show_ligands"]
        self.view.visible_chains = kw["visible_chains"]
        if chains_changed:
            self.view._build_traces()
        self.view.rebuild()
        self.scene.viewport.update()

    def _refuse(self, text: str) -> None:
        self._clear()
        self.message = f"not superposed: {text.splitlines()[0]}"
        self.panel.set_superpose_info(f"<span style='color:#e06c6c'>{self.message}</span>")
        self.status(self.message)

    def _drop(self) -> None:
        """The scene's deposit was replaced: the old fit no longer applies."""
        self._token += 1
        self._clear()

    def _clear(self) -> None:
        if self.view is not None:
            self.view.clear()
        self.view, self.result = None, None
        if self.scene.scene is not None:
            self.scene.viewport.update()
