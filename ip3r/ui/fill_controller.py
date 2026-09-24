"""The Completeness selector: build an AlphaFold fill on a worker, draw it.

Owns the one decision the rest of the window should not have to make: when
the fill is stale. A fill belongs to one deposit, so a load rebuilds it (or
reports the refusal), and a result arriving after a newer request (another
deposit, another choice, a session restore) is dropped.
"""

from __future__ import annotations

from ..structure.graft import GraftRefusal, fill_structure
from .workers import run_async

__all__ = ["FillController", "fill_html"]


def fill_html(model) -> str:
    """The panel text for a built fill: the summary and every warning."""
    warn = "".join(f"<br>• {w}" for w in model.warnings())
    return f"{model.summary()}<span style='color:#d9a441'>{warn}</span>"


class FillController:
    def __init__(self, scene, panel, status=None) -> None:
        self.scene, self.panel = scene, panel
        self.status = status or (lambda msg: None)
        self.model = None
        self._token = 0                         # the latest request wins
        #: The last outcome as text (for the smoke test and the status bar).
        self.message = ""
        panel.completeness_changed.connect(self.request)

    def loaded(self, st) -> None:
        """A new deposit is shown; rebuild the fill for it if one is asked for."""
        self.model = None
        self.request(self.panel.current_completeness())

    def request(self, mode: str) -> None:
        self._token += 1
        token, st = self._token, self.scene.structure
        if mode == "none" or st is None:
            self.model, self.message = None, ""
            self.scene.show_fill(None)
            self.panel.set_fill_info("")
            return
        self.panel.set_fill_info(f"filling {st.name} from AlphaFold…")

        def work():
            try:
                return token, fill_structure(st, mode)
            except GraftRefusal as exc:
                return token, str(exc)
        run_async(work, on_done=self._built,
                  on_error=lambda e: token == self._token
                  and self._show_refusal(f"fill failed: {e}"))

    def _built(self, result) -> None:
        token, model = result
        if token != self._token:
            return                  # superseded: another deposit or another choice
        if isinstance(model, str):
            self._show_refusal(model)
            return
        self.model, self.message = model, model.summary()
        self.scene.show_fill(model)
        self.panel.set_fill_info(fill_html(model))
        self.status(model.summary())

    def _show_refusal(self, text: str) -> None:
        self.model, self.message = None, f"not filled: {text}"
        self.scene.show_fill(None)
        self.panel.set_fill_info(f"<span style='color:#e06c6c'>{self.message}</span>")
        self.status(self.message)
