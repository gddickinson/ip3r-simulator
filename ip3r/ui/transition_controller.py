"""Builds a two-state transition and plays it on the loaded structure.

The path is built from **the structure already displayed** — its own
C-alphas are frame 0 — so the drawn start is the deposit exactly and the
drawn end puts every basis C-alpha on the superposed end deposit exactly
(:func:`ip3r.structure.transition.displaced_coords`). The PIEZO1 simulator
learned this the hard way: a path built in one frame and drawn in another
landed 36 Å from its own endpoint while every interpolation test passed.

Playback shows solved frames only, never a blend of two, because a blend of
two restrained frames is not itself restrained.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..io import loader
from ..physics.transition_modes import TransitionOverlap, transition_overlap
from ..render.representations import Style
from ..structure.morph import MorphTrajectory, morph
from ..structure.transition import (Transition, atom_displacement, atom_site_index,
                                    displaced_coords, prepare_transition)

__all__ = ["TransitionController", "TransitionResult", "build_transition"]


@dataclass
class TransitionResult:
    transition: Transition
    trajectory: MorphTrajectory
    overlap: TransitionOverlap
    site_index: np.ndarray          # per atom of the start deposit
    displacement: np.ndarray        # per atom, Å; NaN off the basis


def build_transition(start, end_id: str, fit: str, method: str) -> TransitionResult:
    """Everything slow, for a worker thread. ``start`` is the displayed Structure."""
    tr = prepare_transition(start, loader.load(end_id), fit)
    traj = morph(tr.start, tr.end, method)
    idx = atom_site_index(start, tr)
    return TransitionResult(tr, traj, transition_overlap(tr, "start"), idx,
                            atom_displacement(start, tr, idx))


class TransitionController:
    def __init__(self, scene) -> None:
        self.scene = scene
        self.result: TransitionResult | None = None
        self.frame = 0

    @property
    def viewport(self):
        return self.scene.viewport

    def install(self, result: TransitionResult) -> None:
        """Accept a built transition for the structure on screen."""
        st = self.scene.structure
        if st is None or st.name != result.transition.start_id:
            raise ValueError("the loaded structure changed while the transition was built")
        self.reset()
        self.result = result
        self.scene.view.displacement = result.displacement
        self.scene.view.rebuild()
        self.viewport.update()

    def coords_at(self, i: int) -> np.ndarray:
        """Atom coordinates drawn at frame ``i`` — the single expression both
        the picture and the tests go through."""
        r = self.result
        return displaced_coords(self.scene.structure.xyz, r.trajectory.frames, i,
                                r.site_index)

    def show_frame(self, i: int, stop: bool = True) -> None:
        if self.result is None or self.scene.view is None:
            return
        if stop:
            self.scene.stop_animation()
        n = len(self.result.trajectory)
        self.frame = int(np.clip(i, 0, n - 1))
        if self.scene.view.style is Style.CARTOON and self.frame:
            self.scene.view.style = Style.TUBE       # a cartoon per frame is too slow
        xyz = self.coords_at(self.frame)
        self.scene.view.update_coords(xyz)
        self.scene.move_overlays(xyz)
        self.viewport.update()

    def play(self, on_frame=None) -> None:
        """Ping-pong through the solved frames; ``on_frame(i)`` follows along."""
        if self.result is None:
            return
        self.scene.stop_animation()
        n = len(self.result.trajectory)
        rate = 12.0                                   # frames per second of playback
        state = {"t": self.frame / rate}

        def tick(dt):
            state["t"] += dt
            k = int(state["t"] * rate) % (2 * (n - 1))
            i = k if k < n else 2 * (n - 1) - k
            if i != self.frame:
                self.show_frame(i, stop=False)
                if on_frame:
                    on_frame(i)
            return True
        self.viewport.add_animation(tick)

    def stop(self) -> None:
        self.viewport.clear_animations()

    def reset(self) -> None:
        """Forget the transition and put the deposit back as loaded."""
        self.stop()
        had = self.result is not None
        self.result, self.frame = None, 0
        view = self.scene.view
        if view is not None and had:
            view.displacement = None
            view.update_coords(self.scene.structure.xyz)
            self.scene.move_overlays(self.scene.structure.xyz)
            self.viewport.update()
