"""What is drawn in the viewport, and how it changes.

Owns the :class:`~ip3r.render.representations.MolecularView` of the loaded
structure plus the extras drawn on top of it: the pore (probe spheres down
the axis), highlighted residues, and the normal-mode animation. The main
window calls these methods; nothing here builds a widget.
"""

from __future__ import annotations

import numpy as np

from ..core.annotations import constraint_at, element_of, functional_sites
from ..core.modules import MODULE_COLORS, MODULES_KEY, ModuleRefusal, modules
from ..physics.anm import ANM, atom_displacements, tetramer_sites
from ..render import colormaps
from ..render.representations import MolecularView, Style
from ..render.variant_spheres import variant_spheres
from .fill_overlay import FillOverlay

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
        self._variants: tuple | None = None       # (paralog, buckets, layer)
        self._variant_atoms = np.zeros(0, int)
        self._fill: FillOverlay | None = None

    @property
    def fill(self) -> FillOverlay:
        """The AlphaFold fill overlay (created once the scene exists)."""
        if self._fill is None:
            self._fill = FillOverlay(self.scene, self.viewport)
        return self._fill

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
        self.scene.remove("variants")
        self.fill.clear()
        self._variants, self._variant_atoms = None, np.zeros(0, int)
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
            if self._variants is not None:
                self.show_variants(*self._variants)
        self.view.rebuild()
        self.fill.restyle(self.view.style, self.view.visible_chains)
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
        st = self.structure
        mask = np.zeros(st.n_atoms, bool)
        rgb = np.tile(np.asarray(self.view.highlight_color, np.float32), (st.n_atoms, 1))
        if self.view.paralog:
            fs = functional_sites(self.view.paralog)
            for k in classes:
                mask |= np.isin(st.res_seq, fs.get(k, ()))
            try:        # a numbering without the modules' sites (RyR1) has none
                mods = modules(self.view.paralog) if MODULES_KEY in classes else ()
            except ModuleRefusal:
                mods = ()
            if mods:
                # Paper 6's two modules as a Cα trace, each in its own colour;
                # drawn under the site balls, which keep the default colour.
                ca = (st.atom_name == "CA") & ~st.hetero
                for m in mods:
                    sel = ca & np.isin(st.res_seq, m.residues) & ~mask
                    rgb[sel] = MODULE_COLORS[m.definition]
                    mask |= sel
        self.view.highlight = mask if mask.any() else None
        self.view.highlight_rgb = rgb
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
        self.view.highlight_rgb = None
        self.view.highlight_color = (1.0, 0.3, 0.9)
        self.view.rebuild()
        self.viewport.update()
        return f"{paralog} {resi}: highlighted on {len(set(self.structure.chain[mask]))} subunits"

    def show_variants(self, paralog: str, buckets, layer: str | None = None) -> str:
        """Variant spheres on every visible subunit; ``buckets`` empty clears."""
        self.scene.remove("variants")
        self._variants, self._variant_atoms = None, np.zeros(0, int)
        if self.view is None or not buckets:
            self.viewport.update()
            return "" if self.view is not None else "load a structure first"
        if self.view.paralog != paralog:
            self.viewport.update()
            return (f"{self.structure.name} is not in human {paralog} numbering "
                    f"(it is in {self.view.paralog or 'no human'} numbering); "
                    "no variant drawn rather than drawn on the wrong residue")
        st = self.structure
        idx, radius, rgb, labels = variant_spheres(
            st, paralog, tuple(buckets), layer, self.view._chain_ok())
        self._variants, self._variant_atoms = (paralog, tuple(buckets), layer), idx
        batch = self.scene.spheres("variants")
        # The view's coordinates, not the deposit's: a morph or mode frame may be up.
        batch.upload(self.view.structure.xyz[idx], radius, rgb,
                     np.zeros(len(idx), np.float32))
        self.viewport.update()
        n_res = len(set(st.res_seq[idx].tolist()))
        n_ch = len(set(st.chain[idx].tolist()))
        return f"{n_res} variant residues drawn on {n_ch} subunits ({len(idx)} spheres)"

    def move_overlays(self, xyz: np.ndarray) -> None:
        """Follow a morph or mode frame: the variant spheres ride their Cα, and
        each AlphaFold fill is re-fitted on its own anchors."""
        batch = self.scene.get("variants")
        if batch is not None and len(self._variant_atoms):
            batch.update_centers(np.asarray(xyz, np.float32)[self._variant_atoms])
        self.fill.move(xyz)

    def show_fill(self, model) -> None:
        """Draw an AlphaFold fill (None clears) on the coordinates now shown."""
        if model is None or self.view is None:
            self.fill.clear()
        else:
            self.fill.show(model, self.view.structure.xyz, self.view.style,
                           self.view.visible_chains)
        self.viewport.update()

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
            self.fill.restyle(Style.TUBE, self.view.visible_chains)
        state = {"t": 0.0}

        def tick(dt):
            state["t"] += dt
            xyz = self._base_xyz + np.sin(2.0 * state["t"]) * self._disp
            self.view.update_coords(xyz)
            self.move_overlays(xyz)
            return True
        self.viewport.add_animation(tick)

    def stop_animation(self) -> None:
        self.viewport.clear_animations()
        if self._base_xyz is not None and self.view is not None:
            self.view.update_coords(self._base_xyz)
            self.move_overlays(self._base_xyz)
            self.viewport.update()
        self._base_xyz = None
