"""The AlphaFold fill drawn beside the deposit, with its seams.

The fill is its own :class:`~ip3r.render.representations.MolecularView`
(batches ``fill:*``), always coloured by pLDDT on AlphaFold's fixed bands, so
predicted residues can never be mistaken for deposited ones, whatever colouring
the deposit is in. Each seam is a bond from the deposit's flanking C-alpha to
the fill's end C-alpha: pale when it closes, red when it does not
(``graft.join_tolerance``). Both follow the deposit through a morph or a mode
frame, because :meth:`~ip3r.structure.graft.FilledModel.place` re-fits every
stretch on its own anchors.
"""

from __future__ import annotations

import numpy as np

from ..render.colormaps import SEAM_COLORS
from ..render.representations import ColorBy, MolecularView, Style

__all__ = ["FillOverlay", "SEAM_RADIUS"]

SEAM_RADIUS = 0.45


class FillOverlay:
    def __init__(self, scene, viewport) -> None:
        self.scene = scene
        self.viewport = viewport
        self.model = None
        self.view: MolecularView | None = None
        self._deposit_xyz = None

    def show(self, model, deposit_xyz, style: Style, chains) -> None:
        self.clear()
        self.model, self._deposit_xyz = model, np.asarray(deposit_xyz, np.float64)
        if model.atoms.n_atoms == 0:
            return
        # Placed on the coordinates now shown, which may be a morph or mode
        # frame rather than the deposit the fill was built on.
        atoms = model.atoms.copy_with_coords(model.place(self._deposit_xyz))
        self.view = MolecularView(self.scene, atoms, name="fill", style=style,
                                  color_by=ColorBy.PLDDT, show_ligands=False,
                                  visible_chains=chains)
        self.view.rebuild()
        self._draw_seams()

    def restyle(self, style: Style, chains) -> None:
        if self.view is None:
            return
        changed = chains != self.view.visible_chains
        self.view.style, self.view.visible_chains = style, chains
        if changed:
            self.view._build_traces()
        self.view.rebuild()
        self._draw_seams()

    def move(self, deposit_xyz) -> None:
        if self.view is None:
            return
        self._deposit_xyz = np.asarray(deposit_xyz, np.float64)
        self.view.update_coords(self.model.place(self._deposit_xyz))
        self._draw_seams()

    def _draw_seams(self) -> None:
        self.scene.remove("seams")
        dep, end, joined = self.model.seam_segments(self._deposit_xyz,
                                                    self.view.structure.xyz)
        chains = self.view.visible_chains
        if chains is not None:
            keep = np.isin([f.stretch.chain for f in self.model.fills
                            for _ in f.seam_fill], list(chains))
            dep, end, joined = dep[keep], end[keep], joined[keep]
        if not len(dep):
            return
        rgb = np.array([SEAM_COLORS[bool(j)] for j in joined], np.float32)
        self.scene.cylinders("seams").upload(
            dep.astype(np.float32), end.astype(np.float32),
            np.full(len(dep), SEAM_RADIUS, np.float32), rgb, rgb)

    def clear(self) -> None:
        if self.view is not None:
            self.view.clear()
        self.scene.remove("seams")
        self.view, self.model = None, None
