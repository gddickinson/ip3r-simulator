"""Interpolation between two experimental states of the tetramer.

A morph is an **interpolation, not a simulated trajectory**: it shows a
geometrically plausible path between two observed endpoints and says
nothing about the barrier between them or the order in which parts move.
Every :class:`MorphTrajectory` carries that sentence in ``meta["note"]``,
and the viewer prints it.

Two methods, ported from the PIEZO1 simulator's ``structure/morph.py``
(the modal method stays there; here the ANM is used to *score* the
transition, :mod:`ip3r.physics.transition_modes`, not to drive it):

``linear``
    Straight-line interpolation of each C-alpha. Exact at both endpoints,
    but a rotating domain takes chords through space, so C-alpha–C-alpha
    distances contract mid-path. The IP3R cytosolic assembly rotates as a
    body between resting and activated states, so this is not academic.

``restrained`` (default)
    Linear interpolation, then iterative restoration of every peptide
    C-alpha–C-alpha distance toward its interpolated target. Removes the
    chord artefact; the endpoints are left untouched because they are
    experimental.

The honesty metric is ``bond_error``: per frame, the worst deviation of a
peptide C-alpha–C-alpha distance from its interpolated target (Å).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["MorphTrajectory", "morph", "interpolate_linear",
           "restrained_morph", "peptide_pairs", "METHODS", "NOTE"]

METHODS = ("restrained", "linear")

NOTE = ("An interpolation between two observed states, not a simulated "
        "trajectory: it shows a plausible path, not the energy barrier or the "
        "order of events.")


@dataclass
class MorphTrajectory:
    """Frames ``(n_frames, n_sites, 3)`` between two superposed endpoints."""

    frames: np.ndarray
    method: str
    bond_error: np.ndarray
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.frames)

    def nearest(self, fraction: float) -> int:
        """Index of the solved frame nearest ``fraction`` in [0, 1].

        A blend of two frames is not itself a restrained frame, so playback
        shows solved frames only.
        """
        f = float(np.clip(fraction, 0.0, 1.0))
        return int(round(f * (len(self.frames) - 1)))

    def summary(self) -> str:
        return (f"{len(self)} frames, {self.method}, worst Cα–Cα "
                f"error {self.bond_error.max():.2f} Å")


def peptide_pairs(start: np.ndarray, end: np.ndarray, n_subunits: int = 4,
                  max_bond: float | None = None) -> np.ndarray:
    """Consecutive-site pairs within each subunit block that are real bonds.

    A pair spanning an unmodelled gap is excluded (longer than ``max_bond``
    in either endpoint): forcing a chain break to a fixed length would drag
    two genuinely distant parts of the model together.
    """
    max_bond = _P.value("morph.max_bond") if max_bond is None else max_bond
    per = len(start) // n_subunits
    idx = np.concatenate([np.arange(k * per, (k + 1) * per - 1)
                          for k in range(n_subunits)])
    pairs = np.stack([idx, idx + 1], axis=1)
    d0 = np.linalg.norm(start[pairs[:, 1]] - start[pairs[:, 0]], axis=1)
    d1 = np.linalg.norm(end[pairs[:, 1]] - end[pairs[:, 0]], axis=1)
    return pairs[(d0 < max_bond) & (d1 < max_bond)]


def _lengths(x: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    return np.linalg.norm(x[pairs[:, 1]] - x[pairs[:, 0]], axis=1)


def interpolate_linear(start: np.ndarray, end: np.ndarray, n_frames: int) -> np.ndarray:
    t = np.linspace(0.0, 1.0, n_frames)[:, None, None]
    return (1.0 - t) * start[None] + t * end[None]


def _targets(start, end, pairs, n_frames):
    d0, d1 = _lengths(start, pairs), _lengths(end, pairs)
    t = np.linspace(0.0, 1.0, n_frames)[:, None]
    return (1.0 - t) * d0[None] + t * d1[None]


def _errors(frames, pairs, targets) -> np.ndarray:
    if not len(pairs):
        return np.zeros(len(frames))
    return np.array([np.abs(_lengths(f, pairs) - tg).max()
                     for f, tg in zip(frames, targets)])


def restrained_morph(start: np.ndarray, end: np.ndarray, pairs: np.ndarray,
                     n_frames: int, iterations: int) -> np.ndarray:
    """Linear frames with peptide distances projected onto their targets.

    Symmetric SHAKE-style projection: each partner moves half the error.
    """
    frames = interpolate_linear(start, end, n_frames)
    targets = _targets(start, end, pairs, n_frames)
    i, j = pairs[:, 0], pairs[:, 1]
    for f in range(1, n_frames - 1):
        x = frames[f]
        for _ in range(iterations):
            delta = x[j] - x[i]
            d = np.maximum(np.linalg.norm(delta, axis=1), 1e-9)
            corr = delta * (0.5 * (d - targets[f]) / d)[:, None]
            np.add.at(x, i, corr)
            np.subtract.at(x, j, corr)
    return frames


def morph(start: np.ndarray, end: np.ndarray, method: str = "restrained",
          n_frames: int | None = None, n_subunits: int = 4) -> MorphTrajectory:
    """A :class:`MorphTrajectory` between two superposed site arrays."""
    start = np.ascontiguousarray(start, np.float64)
    end = np.ascontiguousarray(end, np.float64)
    if start.shape != end.shape:
        raise ValueError(f"endpoint shapes differ: {start.shape} vs {end.shape}")
    if method not in METHODS:
        raise ValueError(f"unknown morph method {method!r}; one of {METHODS}")
    n_frames = int(_P.value("morph.n_frames")) if n_frames is None else int(n_frames)
    pairs = peptide_pairs(start, end, n_subunits)
    if method == "linear":
        frames = interpolate_linear(start, end, n_frames)
    else:
        frames = restrained_morph(start, end, pairs, n_frames,
                                  int(_P.value("morph.iterations")))
    err = _errors(frames, pairs, _targets(start, end, pairs, n_frames))
    return MorphTrajectory(frames, method, err,
                           {"n_bonds": len(pairs), "note": NOTE})
