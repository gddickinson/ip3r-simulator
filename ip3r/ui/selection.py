"""What is selected, what a click does with it, and distances (Round 7.18,
ported in method from PIEZO1's ``selection.py``).

A selection is a set of residues, ``(chain, residue number)`` in the
deposit's own numbering, drawn as gold spheres over the model in their own
batch, so it never disturbs the site highlights the Structure panel owns.
The sequence window and the context menu read and write the same set.

A click selects the residue under it (shift-click adds or removes). With
measuring armed, clicks collect atoms instead: every second atom closes a
distance, drawn as a thin cylinder with its length as a label. Both follow a
morph or mode frame, because the spheres and the rods ride the displayed
coordinates, and both are dropped on a new deposit: a residue number means
nothing on a different structure.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QObject, pyqtSignal

__all__ = ["SelectionController", "residue_mask", "SELECT_RGB", "MEASURE_RGB"]

SELECT_RGB = (1.0, 0.82, 0.25)
MEASURE_RGB = (0.35, 0.95, 0.95)
#: Selection spheres: this fraction of each atom's van der Waals radius.
SELECT_SCALE = 1.0


def residue_mask(st, residues) -> np.ndarray:
    """Atoms of the ``(chain, residue)`` pairs in ``residues`` (protein only)."""
    mask = np.zeros(st.n_atoms, bool)
    by_chain: dict[str, list[int]] = {}
    for ch, r in residues:
        by_chain.setdefault(ch, []).append(int(r))
    for ch, rs in by_chain.items():
        mask |= (st.chain == ch) & np.isin(st.res_seq, rs)
    return mask & ~st.hetero


class SelectionController(QObject):
    #: The selected residue set changed (the sequence window follows).
    changed = pyqtSignal()

    def __init__(self, scene, hud, status) -> None:
        super().__init__()
        self.scene_ctl = scene
        self.hud = hud
        self.status = status
        self.residues: set[tuple[str, int]] = set()
        self.measuring = False
        self._armed: list[int] = []
        #: Completed measurements: (atom i, atom j).
        self.distances: list[tuple[int, int]] = []
        self._xyz: np.ndarray | None = None
        scene.followers.append(self._follow)
        scene.on_clear.append(self.reset)

    @property
    def structure(self):
        return self.scene_ctl.structure

    # ------------------------------------------------------------ selecting

    def pick(self, index: int, extend: bool = False) -> str:
        """A click on atom ``index`` (-1 = empty space). Returns the status."""
        st = self.structure
        if st is None:
            return ""
        if index < 0:
            return "nothing under the cursor"
        if self.measuring:
            return self._measure(index)
        text = self.scene_ctl.describe_atom(index)
        if st.hetero[index]:
            return text
        key = (str(st.chain[index]), int(st.res_seq[index]))
        if extend:
            self.residues ^= {key}
        else:
            self.residues = {key}
        self._redraw()
        return text

    def select(self, residues, replace: bool = True) -> None:
        residues = {(str(c), int(r)) for c, r in residues}
        self.residues = residues if replace else self.residues | residues
        self._redraw()

    def select_everywhere(self, resi: int) -> int:
        """Residue ``resi`` on every subunit that resolves it; how many."""
        st = self.structure
        chains = sorted(set(st.chain[(st.res_seq == resi) & ~st.hetero].tolist()))
        self.select({(c, resi) for c in chains})
        return len(chains)

    def select_chain(self, chain: str) -> None:
        st = self.structure
        rs = np.unique(st.res_seq[(st.chain == chain) & ~st.hetero])
        self.select({(chain, int(r)) for r in rs})

    def clear(self) -> None:
        self.residues = set()
        self._redraw()

    def reset(self) -> None:
        """A new deposit: drop the selection and every measurement."""
        self.residues, self.distances, self._armed = set(), [], []
        self._xyz = None
        self.hud.set_readout("selection", "")
        self.hud.set_readout("distance", "")
        self.changed.emit()

    def describe(self) -> str:
        if not self.residues:
            return ""
        chains = sorted({c for c, _ in self.residues})
        nums = sorted({r for _, r in self.residues})
        span = (f"{nums[0]}" if len(nums) == 1 else
                f"{len(nums)} residues {nums[0]}–{nums[-1]}")
        return f"selected: {span} on {', '.join(chains)}"

    # ------------------------------------------------------------ measuring

    def arm(self, on: bool) -> None:
        self.measuring = on
        self._armed = []
        self.status("measuring: click two atoms (Esc or the menu stops)" if on
                    else "measuring off")

    def add_point(self, index: int) -> str:
        """Add an atom to a measurement without arming clicks (context menu)."""
        return self._measure(index)

    def _measure(self, index: int) -> str:
        self._armed.append(index)
        if len(self._armed) < 2:
            return f"measure from {self._atom_label(index)}: now click the second atom"
        i, j = self._armed[:2]
        self._armed = []
        self.distances.append((i, j))
        self._redraw()
        return self.hud.readouts.get("distance", "")

    def clear_distances(self) -> None:
        self.distances, self._armed = [], []
        self._redraw()

    def _atom_label(self, i: int) -> str:
        st = self.structure
        return f"{st.chain[i]}:{st.res_name[i]}{int(st.res_seq[i])}.{st.atom_name[i]}"

    # -------------------------------------------------------------- drawing

    def _coords(self) -> np.ndarray:
        view = self.scene_ctl.view
        return (self._xyz if self._xyz is not None
                else view.structure.xyz if view is not None else self.structure.xyz)

    def _follow(self, xyz) -> None:
        self._xyz = np.asarray(xyz, np.float64)
        if self.residues or self.distances:
            self._redraw(announce=False)

    def _redraw(self, announce: bool = True) -> None:
        scene = self.scene_ctl.scene
        st = self.structure
        if scene is None or st is None:
            return
        xyz = np.asarray(self._coords(), np.float32)
        scene.remove("selection")
        mask = residue_mask(st, self.residues)
        if mask.any():
            idx = np.flatnonzero(mask)
            batch = scene.spheres("selection")
            batch.alpha = 0.85
            batch.upload(xyz[idx], (st.vdw_radii()[idx] * SELECT_SCALE).astype(np.float32),
                         np.tile(np.float32(SELECT_RGB), (len(idx), 1)),
                         np.zeros(len(idx), np.float32))
        scene.remove("distances")
        labels, text = [], ""
        if self.distances:
            a = np.array([xyz[i] for i, _ in self.distances])
            b = np.array([xyz[j] for _, j in self.distances])
            n = len(a)
            scene.cylinders("distances").upload(
                a, b, np.full(n, 0.12, np.float32), np.tile(np.float32(MEASURE_RGB), (n, 1)))
            d = np.linalg.norm(a - b, axis=1)
            colour = tuple(int(255 * c) for c in MEASURE_RGB)
            labels = [((p + q) / 2, f"{v:.2f} Å", colour) for p, q, v in zip(a, b, d)]
            i, j = self.distances[-1]
            text = f"{self._atom_label(i)} – {self._atom_label(j)}: {d[-1]:.2f} Å"
        self.scene_ctl.viewport.set_overlay_labels(labels)
        self.hud.set_readout("distance", text)
        self.hud.set_readout("selection", self.describe())
        self.scene_ctl.viewport.update()
        if announce:
            self.changed.emit()

    # --------------------------------------------------------------- camera

    def centre_on(self, index: int, zoom: float | None = 40.0) -> None:
        """Pivot about atom ``index``; with ``zoom``, bring the camera to that
        distance (Å). The automatic fit stops: the user chose where to look."""
        scene = self.scene_ctl.scene
        if scene is None or index < 0:
            return
        cam = scene.camera
        cam.pivot = np.asarray(self._coords()[index], np.float64).copy()
        cam.pan = np.zeros(3)
        if zoom is not None:
            cam.distance = float(zoom)
        cam.slab_front = None
        self.scene_ctl.fit_target = None
        self.scene_ctl.viewport.update()

    def centre_on_selection(self) -> None:
        st = self.structure
        mask = residue_mask(st, self.residues) if st is not None else None
        if mask is None or not mask.any():
            return
        idx = np.flatnonzero(mask)
        centre = self._coords()[idx].mean(0)
        k = int(idx[np.argmin(np.linalg.norm(self._coords()[idx] - centre, axis=1))])
        spread = float(np.linalg.norm(self._coords()[idx] - centre, axis=1).max())
        self.centre_on(k, zoom=max(30.0, 3.0 * spread))
