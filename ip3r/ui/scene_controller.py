"""What is drawn in the viewport, and how it changes.

Owns the :class:`~ip3r.render.representations.MolecularView` of the loaded
structure plus the extras drawn on top of it: the pore (probe spheres down
the axis), highlighted residues, and the normal-mode animation. The main
window calls these methods; nothing here builds a widget.
"""

from __future__ import annotations

import numpy as np

from ..core.annotations import constraint_at, element_of, functional_sites
from ..physics.anm import ANM, atom_displacements, tetramer_sites
from ..render import colormaps
from ..render.representations import MolecularView, Style

__all__ = ["SceneController", "matrix_to_quat"]


def matrix_to_quat(r: np.ndarray) -> np.ndarray:
    """Unit quaternion (w, x, y, z) of a rotation matrix."""
    w = np.sqrt(max(0.0, 1.0 + r[0, 0] + r[1, 1] + r[2, 2])) / 2.0
    x = np.sqrt(max(0.0, 1.0 + r[0, 0] - r[1, 1] - r[2, 2])) / 2.0
    y = np.sqrt(max(0.0, 1.0 - r[0, 0] + r[1, 1] - r[2, 2])) / 2.0
    z = np.sqrt(max(0.0, 1.0 - r[0, 0] - r[1, 1] + r[2, 2])) / 2.0
    x = np.copysign(x, r[2, 1] - r[1, 2])
    y = np.copysign(y, r[0, 2] - r[2, 0])
    z = np.copysign(z, r[1, 0] - r[0, 1])
    q = np.array([w, x, y, z])
    return q / np.linalg.norm(q)


class SceneController:
    def __init__(self, viewport) -> None:
        self.viewport = viewport
        self.view: MolecularView | None = None
        self.structure = None
        self.summary = None
        self.modes = None
        self._sites: list[str] = []
        self._base_xyz = None
        self._disp = None

    def attach(self) -> None:
        """Called once the GL context exists."""

    @property
    def scene(self):
        return self.viewport.scene

    # ------------------------------------------------------------ structure

    def set_structure(self, st, summary, style: dict) -> None:
        self.stop_animation()
        if self.view is not None:
            self.view.clear()
        self.scene.remove("pore")
        self.structure, self.summary, self.modes = st, summary, None
        paralog = summary.numbering.paralog if summary.numbering else None
        self.view = MolecularView(self.scene, st, name="model", paralog=paralog, **style)
        self.view.rebuild()
        self.viewport.set_pick_source(st.xyz)
        self.side_view()

    def restyle(self, style: dict) -> None:
        if self.view is None:
            return
        chains_changed = style.get("visible_chains") != self.view.visible_chains
        for k, v in style.items():
            setattr(self.view, k, v)
        if chains_changed:
            self.view._build_traces()
        self.view.rebuild()
        self.viewport.update()

    # --------------------------------------------------------------- camera

    def _orient(self, rot_rows: np.ndarray) -> None:
        cam = self.scene.camera
        cam.rotation = matrix_to_quat(rot_rows)
        cam.frame(self.structure.xyz)
        self.viewport.update()

    def side_view(self) -> None:
        """The four-fold axis vertical on screen, cytosolic cap up."""
        if self.summary is None or self.scene is None:
            return
        e1, e2, e3 = self.summary.frame.basis.T
        self._orient(np.array([e1, e3, np.cross(e1, e3)]))

    def top_view(self) -> None:
        """Looking down the pore from the cytosol."""
        if self.summary is None or self.scene is None:
            return
        e1, e2, e3 = self.summary.frame.basis.T
        self._orient(np.array([e1, e2, e3]))

    # --------------------------------------------------------------- extras

    def show_pore(self, on: bool) -> None:
        self.scene.remove("pore")
        if on and self.summary is not None:
            p = self.summary.profile
            keep = p.r_free > 0.3
            pts = np.column_stack([np.zeros(keep.sum()), np.zeros(keep.sum()), p.z[keep]])
            xyz = self.summary.frame.from_frame(pts)
            r = p.r_free[keep]
            col = colormaps.ramp(1.0 - np.clip((r - 1.0) / 6.0, 0, 1))
            batch = self.scene.spheres("pore")
            batch.alpha = 0.45
            batch.upload(xyz, r, col, np.zeros(len(r), np.float32))
        self.viewport.update()

    def highlight_sites(self, classes: list[str]) -> None:
        self._sites = classes
        if self.view is None:
            return
        mask = np.zeros(self.structure.n_atoms, bool)
        if self.view.paralog:
            fs = functional_sites(self.view.paralog)
            for k in classes:
                mask |= np.isin(self.structure.res_seq, fs.get(k, ()))
        self.view.highlight = mask if mask.any() else None
        self.view.rebuild()
        self.viewport.update()

    def highlight_residue(self, paralog: str, resi: int) -> str:
        if self.view is None:
            return "load a structure first"
        if self.view.paralog != paralog:
            return (f"{self.structure.name} is not in human {paralog} numbering "
                    f"(it is in {self.view.paralog or 'no human'} numbering); "
                    "residue not highlighted rather than placed on the wrong residue")
        mask = (self.structure.res_seq == resi) & ~self.structure.hetero
        if not mask.any():
            return f"{paralog} residue {resi} is not resolved in {self.structure.name}"
        self.view.highlight = mask
        self.view.highlight_color = (1.0, 0.3, 0.9)
        self.view.rebuild()
        self.viewport.update()
        return f"{paralog} {resi}: highlighted on {len(set(self.structure.chain[mask]))} subunits"

    def describe_atom(self, i: int) -> str:
        if i < 0 or self.structure is None:
            return ""
        st = self.structure
        r = int(st.res_seq[i])
        text = (f"chain {st.chain[i]} {st.res_name[i]}{r} atom {st.atom_name[i]} "
                f"(B {st.b_factor[i]:.1f})")
        par = self.view.paralog if self.view else None
        if par and not st.hetero[i]:
            jsd = constraint_at(par, np.array([r]), "deep")[0]
            text += f" — {par} element {element_of(par, r)}, deep JSD " + (
                f"{jsd:.3f}" if np.isfinite(jsd) else "not scored")
        return text

    # ---------------------------------------------------------------- modes

    def compute_modes(self):
        """Runs on a worker thread: returns (ModeSet, description)."""
        st, fr = self.structure, self.summary.frame
        coords, residues = tetramer_sites(st, fr)
        anm = ANM(coords, axis=fr.axis)
        modes = anm.label_symmetry(anm.calc_modes())
        self.modes = (modes, residues)
        return modes, f"{len(coords):,} sites ({len(residues)} per subunit)"

    def animate_mode(self, index: int, amplitude: float) -> None:
        if self.modes is None or self.view is None:
            return
        modes, residues = self.modes
        self.stop_animation()
        self._base_xyz = self.structure.xyz.astype(np.float64).copy()
        self._disp = atom_displacements(self.structure, self.summary.frame, residues,
                                        modes.mode(index, amplitude))
        if self.view.style is Style.CARTOON:
            self.view.style = Style.TUBE           # a cartoon per frame is too slow
        state = {"t": 0.0}

        def tick(dt):
            state["t"] += dt
            self.view.update_coords(self._base_xyz + np.sin(2.0 * state["t"]) * self._disp)
            return True
        self.viewport.add_animation(tick)

    def stop_animation(self) -> None:
        self.viewport.clear_animations()
        if self._base_xyz is not None and self.view is not None:
            self.view.update_coords(self._base_xyz)
            self.viewport.update()
        self._base_xyz = None
