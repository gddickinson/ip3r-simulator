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
from ..physics.network_checks import describe_local, local_modes
from ..physics.anm import ANM, atom_displacements, tetramer_sites
from ..render import colormaps
from ..render.representations import MolecularView, Style
from ..render.variant_spheres import variant_spheres
from ..structure.ligand import ligand_sites, neighbourhood
from .fill_overlay import FillOverlay

__all__ = ["SceneController", "matrix_to_quat", "SITE_MARGIN"]

#: The IP3 site view pulls back this much beyond S22's 15 Å pocket, for context.
SITE_MARGIN = 1.3


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
        #: (mode index, amplitude Å) while a mode animates, else None.
        self.animated_mode: tuple[int, float] | None = None
        self._variants: tuple | None = None       # (paralog, buckets, layer)
        self._variant_atoms = np.zeros(0, int)
        self._fill: FillOverlay | None = None
        #: What the camera keeps framed while the user has not moved it:
        #: "all" (the visible subunits), "site" (one IP3 site) or None.
        self.fit_target: str | None = None

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
        self.scene.remove("lumen")
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
        if chains_changed and self.fit_target is not None:
            self.refit()
        self.viewport.update()

    # --------------------------------------------------------------- camera

    def _visible_atoms(self) -> np.ndarray:
        """Atoms of the visible subunits (all of them if none is visible)."""
        ok = self.view._chain_ok() if self.view is not None else None
        if ok is None or not ok.any():
            return np.arange(self.structure.n_atoms)
        return np.flatnonzero(ok)

    def _site_atoms(self) -> tuple[np.ndarray, str] | None:
        """The pocket of the first IP3 on a visible subunit, and its subunit."""
        visible = self.view.visible_chains if self.view is not None else None
        sites = ligand_sites(self.structure)
        shown = [s for s in sites if visible is None or s.subunit in visible]
        if not shown:
            return None
        return neighbourhood(self.structure, shown[0]), shown[0].subunit

    def refit(self) -> str:
        """Frame :attr:`fit_target` in the current orientation.

        The camera's aspect is the viewport's, so a fit made before the docks
        settle is redone on every resize until the user moves the camera
        (:meth:`navigated`). Returns what was framed.
        """
        if self.structure is None or self.scene is None or self.fit_target is None:
            return ""
        xyz = self.view.structure.xyz if self.view is not None else self.structure.xyz
        everything = xyz[self._visible_atoms()]
        if self.fit_target == "site":
            found = self._site_atoms()
            if found is not None:
                atoms, subunit = found
                self.scene.camera.frame(xyz[atoms], margin=SITE_MARGIN,
                                        scene=everything, slab=True)
                self.viewport.update()
                return f"IP3 site on subunit {subunit}"
            self.fit_target = "all"
        self.scene.camera.frame(everything)
        self.viewport.update()
        return "whole structure"

    def navigated(self) -> None:
        """The user moved the camera: stop re-fitting behind their back."""
        self.fit_target = None

    def resized(self) -> None:
        if self.fit_target is not None:
            self.refit()

    def _orient(self, rot_rows: np.ndarray, target: str = "all") -> str:
        self.scene.camera.rotation = matrix_to_quat(rot_rows)
        self.fit_target = target
        return self.refit()

    def _side_rows(self) -> np.ndarray:
        e1, e2, e3 = self.summary.frame.basis.T
        return np.array([e1, e3, np.cross(e1, e3)])

    def side_view(self) -> None:
        """The four-fold axis vertical on screen, cytosolic cap up."""
        if self.summary is None or self.scene is None:
            return
        self._orient(self._side_rows())

    def top_view(self) -> None:
        """Looking down the pore from the cytosol."""
        if self.summary is None or self.scene is None:
            return
        e1, e2, e3 = self.summary.frame.basis.T
        self._orient(np.array([e1, e2, e3]))

    def fit_view(self) -> None:
        """Frame the visible subunits without turning the camera."""
        if self.summary is not None and self.scene is not None:
            self.fit_target = "all"
            self.refit()

    def site_view(self) -> str:
        """Side-on, centred on one IP3 site and its S22 pocket (15 Å).

        The site is seen from outside the tetramer: the camera looks along
        the radial direction from the axis to the ligand, cytosol up.
        Deposits without IP3 fall back to the whole structure, and say so.
        """
        if self.summary is None or self.scene is None:
            return ""
        found = self._site_atoms()
        if found is None:
            self._orient(self._side_rows())
            return f"{self.structure.name} has no IP3 on a visible subunit: whole structure framed"
        atoms, _ = found
        fr = self.summary.frame
        local = fr.to_frame(self.structure.xyz[atoms].mean(0)[None])[0]
        e1, e2, e3 = fr.basis.T
        radial = local[0] * e1 + local[1] * e2
        radial /= max(np.linalg.norm(radial), 1e-9)
        right = np.cross(e3, radial)          # screen x; e3 up; toward the viewer = radial
        return self._orient(np.array([right, e3, radial]), target="site")

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

    def show_lumen(self, mesh) -> None:
        """Draw the lumen surface (a :class:`~ip3r.render.lumen_mesh.LumenMesh`;
        None clears). It is the deposit's: hidden on a morph or mode frame."""
        self.scene.remove("lumen")
        if mesh is not None and self.scene is not None:
            from ..parameters import PARAMETERS as _P
            batch = self.scene.mesh("lumen", two_sided=True, transparent=True)
            batch.upload(mesh.positions, mesh.normals, mesh.colors, mesh.indices,
                         alpha=_P.value("display.lumen_alpha"))
            if self.view is not None:
                batch.visible = self._at_deposit(self.view.structure.xyz)
        self.viewport.update()

    def _at_deposit(self, xyz) -> bool:
        return self.structure is not None and np.array_equal(
            np.asarray(xyz, np.float32), np.asarray(self.structure.xyz, np.float32))

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
        each AlphaFold fill is re-fitted on its own anchors; the lumen, solved
        on the deposit, shows only when the deposit's coordinates are back."""
        batch = self.scene.get("variants")
        if batch is not None and len(self._variant_atoms):
            batch.update_centers(np.asarray(xyz, np.float32)[self._variant_atoms])
        lumen = self.scene.get("lumen")
        if lumen is not None:
            lumen.visible = self._at_deposit(xyz)
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
        local = local_modes(modes, residues, tetramer_sites(st, fr, stride=1)[1])
        return (modes, f"{len(coords):,} sites ({len(residues)} per subunit)",
                {m.index: describe_local([m])[0].split(") ", 1)[1] for m in local})

    def animate_mode(self, index: int, amplitude: float) -> None:
        if self.modes is None or self.view is None:
            return
        modes, residues = self.modes
        self.stop_animation()
        self.animated_mode = (int(index), float(amplitude))
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
        self.animated_mode = None
        if self._base_xyz is not None and self.view is not None:
            self.view.update_coords(self._base_xyz)
            self.move_overlays(self._base_xyz)
            self.viewport.update()
        self._base_xyz = None
