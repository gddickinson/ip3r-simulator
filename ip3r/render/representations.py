"""Turning a :class:`~ip3r.core.structure.Structure` into GPU batches.

A :class:`MolecularView` owns every batch belonging to one loaded model and
rebuilds them when the style or colouring changes, so the Qt layer only says
"cartoon, coloured by conservation" and never touches an OpenGL buffer.

Residue-keyed colourings (functional element, conservation) need to know
which human paralog's numbering the deposit is in. ``paralog`` is set from
:func:`ip3r.structure.numbering.best_numbering`; when it is ``None`` those
colourings paint grey rather than borrowing a numbering that does not fit.

Adapted from the PIEZO1 simulator's representations module (cartoon/tube
sweeps, impostor atoms and bonds), with the PIEZO-specific entity and
component machinery removed and the IP3R colourings added.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from ..core.annotations import constraint_at
from ..core.structure import Structure
from . import colormaps
from .geometry_builders import Mesh, build_cartoon, build_tube
from .scene import Scene
from .spline import assign_secondary_structure

__all__ = ["Style", "ColorBy", "MolecularView", "COLOR_LABELS", "STYLE_LABELS"]

BALL_RADIUS = 0.42
BOND_RADIUS = 0.20
STICK_RADIUS = 0.28


class Style(str, Enum):
    CARTOON = "cartoon"
    TUBE = "tube"
    BACKBONE = "backbone"
    SPHERES = "spheres"
    STICKS = "sticks"
    BALL_AND_STICK = "ball_and_stick"


class ColorBy(str, Enum):
    ELEMENT_DOMAIN = "element"        # functional element (domain map)
    CHAIN = "chain"
    CONSERVATION = "conservation"     # S17 JSD, fixed scale
    SECONDARY = "secondary"
    BFACTOR = "bfactor"
    ATOM = "atom"
    VALUE = "value"                   # an arbitrary per-atom scalar
    DISPLACEMENT = "displacement"     # between two states, fixed scale
    UNIFORM = "uniform"


STYLE_LABELS = {Style.CARTOON: "Cartoon", Style.TUBE: "Tube",
                Style.BACKBONE: "Backbone trace", Style.SPHERES: "Spheres",
                Style.STICKS: "Sticks", Style.BALL_AND_STICK: "Ball and stick"}
COLOR_LABELS = {ColorBy.ELEMENT_DOMAIN: "Functional element",
                ColorBy.CHAIN: "Subunit", ColorBy.CONSERVATION: "Conservation (JSD)",
                ColorBy.SECONDARY: "Secondary structure",
                ColorBy.BFACTOR: "B-factor / pLDDT", ColorBy.ATOM: "Atom type",
                ColorBy.VALUE: "Mode amplitude",
                ColorBy.DISPLACEMENT: "Displacement (Transition tab)",
                ColorBy.UNIFORM: "Uniform"}

SS_COLORS = np.array([[0.55, 0.58, 0.66], [0.94, 0.42, 0.42],
                      [0.98, 0.82, 0.35]], dtype=np.float32)
_RIBBONS = (Style.CARTOON, Style.TUBE, Style.BACKBONE)
_SUFFIXES = ("ribbon", "atoms", "bonds", "ligands", "highlight")


@dataclass
class ChainTrace:
    chain: str
    indices: np.ndarray
    xyz: np.ndarray
    res_seq: np.ndarray
    ss: np.ndarray


@dataclass
class MolecularView:
    scene: Scene
    structure: Structure
    name: str = "model"
    style: Style = Style.CARTOON
    color_by: ColorBy = ColorBy.ELEMENT_DOMAIN
    paralog: str | None = None
    layer: str = "deep"
    show_ligands: bool = True
    show_hydrogens: bool = False
    values: np.ndarray | None = None
    displacement: np.ndarray | None = None    # per atom, Å; NaN = not measured
    highlight: np.ndarray | None = None       # per-atom bool, drawn as balls
    highlight_color: tuple = (1.0, 0.85, 0.2)
    highlight_rgb: np.ndarray | None = None   # per-atom colours; overrides the above
    visible_chains: frozenset | None = None
    traces: list[ChainTrace] = field(default_factory=list)

    def __post_init__(self) -> None:
        self._build_traces()

    # ---------------------------------------------------------------- traces

    def _chain_ok(self) -> np.ndarray:
        st = self.structure
        if self.visible_chains is None:
            return np.ones(st.n_atoms, bool)
        return np.isin(st.chain, list(self.visible_chains))

    def _build_traces(self) -> None:
        st = self.structure
        ca = st.mask_ca() & self._chain_ok()
        self.traces = []
        for chain in st.chains:
            idx = np.flatnonzero(ca & (st.chain == chain))
            if len(idx) < 4:
                continue
            xyz = st.xyz[idx].astype(np.float64)
            self.traces.append(ChainTrace(chain, idx, xyz, st.res_seq[idx],
                                          assign_secondary_structure(xyz)))

    # --------------------------------------------------------------- colours

    def atom_colors(self) -> np.ndarray:
        st = self.structure
        cb = self.color_by
        if cb is ColorBy.ELEMENT_DOMAIN:
            return colormaps.element_colors(st, self.paralog)
        if cb is ColorBy.CHAIN:
            return colormaps.chain_colors(st)
        if cb is ColorBy.CONSERVATION:
            if self.paralog is None:
                return np.tile(colormaps.MISSING, (st.n_atoms, 1))
            vals = constraint_at(self.paralog, st.res_seq, self.layer)
            vals[st.hetero] = np.nan
            return colormaps.constraint_colors(vals)
        if cb is ColorBy.BFACTOR:
            return colormaps.bfactor_colors(st)
        if cb is ColorBy.ATOM:
            return st.element_colors()
        if cb is ColorBy.VALUE and self.values is not None:
            return colormaps.value_colors(self.values)
        if cb is ColorBy.DISPLACEMENT:
            if self.displacement is None:
                return np.tile(colormaps.MISSING, (st.n_atoms, 1))
            return colormaps.displacement_colors(self.displacement)
        if cb is ColorBy.SECONDARY:
            out = np.tile(SS_COLORS[0], (st.n_atoms, 1))
            for tr in self.traces:
                out[tr.indices] = SS_COLORS[tr.ss]
            return out
        return colormaps.uniform_color(st)

    # ------------------------------------------------------------- building

    def rebuild(self) -> None:
        self.clear()
        colors = self.atom_colors()
        if self.style in _RIBBONS:
            self._build_ribbons(colors)
        else:
            self._build_atoms(colors)
        if self.show_ligands:
            self._build_ligands()
        self._build_highlight(colors)

    def _segments(self, tr: ChainTrace, max_gap: int = 1):
        breaks = np.flatnonzero(np.diff(tr.res_seq) > max_gap) + 1
        return np.split(np.arange(len(tr.xyz)), breaks)

    def _build_ribbons(self, colors: np.ndarray) -> None:
        mesh = Mesh.empty()
        for tr in self.traces:
            for seg in self._segments(tr):
                if len(seg) < 4:
                    continue
                xyz, col = tr.xyz[seg], colors[tr.indices[seg]].astype(np.float64)
                if self.style is Style.CARTOON:
                    part = build_cartoon(xyz, tr.ss[seg], col)
                elif self.style is Style.TUBE:
                    part = build_tube(xyz, col, radius=0.85, sides=10)
                else:
                    part = build_tube(xyz, col, radius=0.35, sides=6, subdivisions=3)
                mesh = mesh.concat(part)
        if mesh.n_vertices:
            self.scene.mesh(f"{self.name}:ribbon").upload(
                mesh.positions, mesh.normals, mesh.colors, mesh.indices)

    def _atom_mask(self) -> np.ndarray:
        st = self.structure
        m = self._chain_ok() & ~st.hetero
        if not self.show_hydrogens:
            m &= st.element != "H"
        return m

    def _build_atoms(self, colors: np.ndarray) -> None:
        st = self.structure
        mask = self._atom_mask()
        if self.style is Style.SPHERES:
            radii = st.vdw_radii()
        else:
            r = STICK_RADIUS if self.style is Style.STICKS else BALL_RADIUS
            radii = np.full(st.n_atoms, r, np.float32)
        self.scene.spheres(f"{self.name}:atoms").upload(
            st.xyz[mask], radii[mask], colors[mask], np.zeros(int(mask.sum()), np.float32))
        if self.style in (Style.STICKS, Style.BALL_AND_STICK):
            self._build_bonds(colors, mask)

    def _build_bonds(self, colors: np.ndarray, mask: np.ndarray,
                     name: str = "bonds", radius: float | None = None) -> None:
        from scipy.spatial import cKDTree
        st = self.structure
        idx = np.flatnonzero(mask)
        pairs = np.asarray(sorted(cKDTree(st.xyz[idx]).query_pairs(1.95)), np.int64)
        if len(pairs) == 0:
            return
        pairs = idx[pairs]
        pairs = pairs[st.chain[pairs[:, 0]] == st.chain[pairs[:, 1]]]
        r = radius or (STICK_RADIUS if self.style is Style.STICKS else BOND_RADIUS)
        self.scene.cylinders(f"{self.name}:{name}").upload(
            st.xyz[pairs[:, 0]], st.xyz[pairs[:, 1]], np.full(len(pairs), r, np.float32),
            colors[pairs[:, 0]], colors[pairs[:, 1]])

    def _build_ligands(self) -> None:
        st = self.structure
        mask = st.hetero & self._chain_ok() & (st.element != "H") \
            & ~np.isin(st.res_name, ["HOH", "WAT"])
        if not mask.any():
            return
        colors = st.element_colors()
        # IP3 drawn full size so the four sites read at a glance.
        scale = np.where(st.res_name[mask] == "I3P", 0.9, 0.6)
        self.scene.spheres(f"{self.name}:ligands").upload(
            st.xyz[mask], st.vdw_radii()[mask] * scale, colors[mask],
            np.zeros(int(mask.sum()), np.float32))

    def _build_highlight(self, colors: np.ndarray) -> None:
        if self.highlight is None or not np.any(self.highlight):
            return
        st = self.structure
        mask = self.highlight & self._chain_ok() & (st.element != "H") & ~st.hetero
        if self.highlight_rgb is not None:
            col = np.asarray(self.highlight_rgb, np.float32)[mask]
        else:
            col = np.tile(np.asarray(self.highlight_color, np.float32), (int(mask.sum()), 1))
        self.scene.spheres(f"{self.name}:highlight").upload(
            st.xyz[mask], np.full(int(mask.sum()), 1.1, np.float32), col,
            np.ones(int(mask.sum()), np.float32))

    # ------------------------------------------------------------- animation

    def update_coords(self, xyz: np.ndarray) -> None:
        """Push new coordinates (a normal-mode frame) and redraw."""
        self.structure = self.structure.copy_with_coords(xyz)
        for tr in self.traces:
            tr.xyz = self.structure.xyz[tr.indices].astype(np.float64)
        self.rebuild()

    def clear(self) -> None:
        for s in _SUFFIXES:
            self.scene.remove(f"{self.name}:{s}")

    def set_visible(self, visible: bool) -> None:
        for s in _SUFFIXES:
            self.scene.set_visible(f"{self.name}:{s}", visible)

    def all_coords(self) -> np.ndarray:
        return self.structure.xyz[self._chain_ok()]
