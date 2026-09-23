"""Anisotropic network model (ANM) of the IP3R tetramer, with C4 irreps.

An ANM joins C-alpha beads with harmonic springs; the low-frequency normal
modes of the resulting Hessian are the large collective motions (Atilgan et
al. 2001). Springs are weighted ``gamma (d0/d)^2`` (Yang et al. 2009), which
softens the long-range contacts a hard cutoff over-stiffens.

**Why the symmetry label matters for this channel.** The receptor is a C4
homotetramer, so every mode transforms as one irreducible representation of
C4. The character of a mode under a 90° turn is ``+1`` (A: all four subunits
move alike), ``-1`` (B: neighbours move oppositely) or ``0`` (E: a degenerate
pair). IP3 binding at all four sites and the opening of a four-fold pore are
both A-symmetric, so only A modes can couple to them at first order — the
lowest A mode is the candidate gating coordinate, and E modes (which tilt
the cap) cannot be. The same argument the PIEZO1 simulator makes for C3,
with one more irrep.

Sites are built subunit-by-subunit over the residues **all four** subunits
resolve, in right-handed order about the axis, which is what makes the
permutation in :meth:`ANM.label_symmetry` well defined.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.symmetry import Frame, rotation_matrix, subunit_ca

__all__ = ["build_hessian", "ModeSet", "ANM", "tetramer_sites", "atom_displacements",
           "IRREP_CHARACTERS"]

#: Character of each C4 irrep under the 90° generator.
IRREP_CHARACTERS = {"A": 1.0, "B": -1.0, "E": 0.0}


def build_hessian(coords: np.ndarray, cutoff: float, gamma: float,
                  d0: float) -> sp.csr_matrix:
    """Sparse 3N x 3N ANM Hessian with inverse-square spring weights."""
    coords = np.ascontiguousarray(coords, np.float64)
    n = len(coords)
    if n < 4:
        raise ValueError("need at least four sites for an ANM")
    pairs = np.asarray(sorted(cKDTree(coords).query_pairs(cutoff)), np.int64)
    if len(pairs) == 0:
        raise ValueError(f"no contacts within {cutoff} A")
    i, j = pairs[:, 0], pairs[:, 1]
    diff = coords[j] - coords[i]
    dist = np.linalg.norm(diff, axis=1)
    unit = diff / dist[:, None]
    k = gamma * (d0 / np.maximum(dist, 1e-6)) ** 2
    blocks = -k[:, None, None] * unit[:, :, None] * unit[:, None, :]
    a, b = (x.ravel() for x in np.meshgrid(np.arange(3), np.arange(3),
                                            indexing="ij"))
    rows = np.concatenate([(3 * i[:, None] + a).ravel(), (3 * j[:, None] + a).ravel()])
    cols = np.concatenate([(3 * j[:, None] + b).ravel(), (3 * i[:, None] + b).ravel()])
    vals = np.concatenate([blocks.reshape(-1, 9).ravel(),
                           blocks.transpose(0, 2, 1).reshape(-1, 9).ravel()])
    # Diagonal super-elements: minus the sum of the off-diagonal blocks, per
    # component, which keeps the six rigid-body modes exactly at zero.
    diag = np.zeros((n, 3, 3))
    np.add.at(diag, i, -blocks)
    np.add.at(diag, j, -blocks)
    d_rows = (3 * np.arange(n)[:, None] + a).ravel()
    d_cols = (3 * np.arange(n)[:, None] + b).ravel()
    rows = np.concatenate([rows, d_rows])
    cols = np.concatenate([cols, d_cols])
    vals = np.concatenate([vals, diag.reshape(-1, 9).ravel()])
    return sp.coo_matrix((vals, (rows, cols)), shape=(3 * n, 3 * n)).tocsr()


@dataclass
class ModeSet:
    """Normal modes: ``vectors`` is ``(n_modes, n_sites, 3)``."""

    eigenvalues: np.ndarray
    vectors: np.ndarray
    symmetry: np.ndarray | None = None      # 'A', 'B', 'E' or 'mixed'
    character: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    @property
    def n_modes(self) -> int:
        return len(self.eigenvalues)

    def mode(self, index: int, amplitude: float = 1.0) -> np.ndarray:
        """One mode's displacement field scaled to a peak of ``amplitude`` Å."""
        v = self.vectors[index]
        peak = np.linalg.norm(v, axis=1).max()
        return v * (amplitude / peak) if peak > 0 else v

    def msf(self, n_modes: int | None = None) -> np.ndarray:
        """Mean-square fluctuation per site, summed over modes as 1/lambda."""
        k = self.n_modes if n_modes is None else min(n_modes, self.n_modes)
        w = 1.0 / np.maximum(self.eigenvalues[:k], 1e-12)
        return np.einsum("m,mij,mij->i", w, self.vectors[:k], self.vectors[:k])

    def overlap(self, displacement: np.ndarray) -> np.ndarray:
        """|cosine| of each mode with an observed displacement field."""
        d = np.asarray(displacement, np.float64).ravel()
        if not np.linalg.norm(d):
            return np.zeros(self.n_modes)
        flat = self.vectors.reshape(self.n_modes, -1)
        flat = flat / np.linalg.norm(flat, axis=1, keepdims=True)
        return np.abs(flat @ (d / np.linalg.norm(d)))

    def first(self, irrep: str) -> int | None:
        if self.symmetry is None:
            return None
        hits = np.flatnonzero(self.symmetry == irrep)
        return int(hits[0]) if len(hits) else None


def tetramer_sites(st: Structure, frame: Frame, stride: int | None = None
                   ) -> tuple[np.ndarray, np.ndarray]:
    """C-alpha sites over residues all four subunits resolve.

    Returns ``(coords (4m, 3), residue_numbers (m,))`` with subunit blocks in
    ``frame.chains`` order.
    """
    stride = int(_P.value("anm.stride")) if stride is None else int(stride)
    ca = subunit_ca(st)
    shared = sorted(set.intersection(*(set(ca[c]) for c in frame.chains)))
    shared = shared[::max(stride, 1)]
    coords = np.vstack([[ca[c][r] for r in shared] for c in frame.chains])
    return coords, np.asarray(shared)


@dataclass
class ANM:
    coords: np.ndarray
    n_subunits: int = 4
    axis: np.ndarray | None = None
    cutoff: float = field(default_factory=lambda: _P.value("anm.cutoff"))
    gamma: float = field(default_factory=lambda: _P.value("anm.gamma"))
    d0: float = field(default_factory=lambda: _P.value("anm.d0"))
    hessian: sp.csr_matrix | None = field(default=None, repr=False)

    def build(self) -> "ANM":
        self.hessian = build_hessian(self.coords, self.cutoff, self.gamma, self.d0)
        return self

    def n_components(self) -> int:
        """Disconnected pieces; each adds six zero modes that must be dropped."""
        tree = cKDTree(self.coords)
        adj = tree.sparse_distance_matrix(tree, self.cutoff, output_type="coo_matrix")
        return int(connected_components(adj, directed=False)[0])

    def calc_modes(self, n_modes: int | None = None) -> ModeSet:
        n_modes = int(_P.value("anm.n_modes")) if n_modes is None else n_modes
        if self.hessian is None:
            self.build()
        zero = 6 * self.n_components()
        k = n_modes + zero
        vals, vecs = spla.eigsh(self.hessian, k=k, sigma=-1e-6, which="LM")
        order = np.argsort(vals)
        vals, vecs = vals[order][zero:], vecs[:, order][:, zero:]
        return ModeSet(vals, vecs.T.reshape(len(vals), -1, 3),
                       meta={"cutoff": self.cutoff, "d0": self.d0,
                             "zero_modes_dropped": zero,
                             "n_sites": len(self.coords)})

    def label_symmetry(self, modes: ModeSet, tolerance: float | None = None) -> ModeSet:
        """Character of each mode under the 90° generator, and its irrep.

        ``S`` moves subunit k's displacement onto subunit k+1 and rotates it
        by 360°/n about the axis. For a degenerate E pair the character of
        each member is 0 only if the pair is resolved; it is always 0 for
        their sum, which is why E is identified by |chi| < tolerance.
        """
        tol = _P.value("anm.symmetry_tolerance") if tolerance is None else tolerance
        n = self.n_subunits
        per = len(self.coords) // n
        axis = self.axis if self.axis is not None else np.array([0.0, 0.0, 1.0])
        rot = rotation_matrix(axis, 2 * np.pi / n)
        chars = np.empty(modes.n_modes)
        for m in range(modes.n_modes):
            u = modes.vectors[m]
            su = np.roll(u.reshape(n, per, 3), 1, axis=0).reshape(-1, 3) @ rot.T
            chars[m] = float(np.sum(u * su) / np.sum(u * u))
        labels = np.full(modes.n_modes, "mixed", dtype="U5")
        for name, chi in IRREP_CHARACTERS.items():
            labels[np.abs(chars - chi) < tol] = name
        modes.symmetry, modes.character = labels, chars
        return modes


def atom_displacements(st: Structure, frame: Frame, residues: np.ndarray,
                       mode_vectors: np.ndarray) -> np.ndarray:
    """Spread a per-site displacement onto every atom of the structure.

    Each atom takes its subunit's displacement at the nearest sampled residue
    (the stride skips residues; a collective mode is smooth along the chain).
    Atoms of chains outside the model do not move.
    """
    out = np.zeros_like(st.xyz, dtype=np.float64)
    per = len(residues)
    for k, ch in enumerate(frame.chains):
        m = st.chain == ch
        idx = np.clip(np.searchsorted(residues, st.res_seq[m]), 0, per - 1)
        lower = np.clip(idx - 1, 0, per - 1)
        use_lower = np.abs(residues[lower] - st.res_seq[m]) < np.abs(residues[idx] - st.res_seq[m])
        idx = np.where(use_lower, lower, idx)
        out[m] = mode_vectors[k * per + idx]
    return out
