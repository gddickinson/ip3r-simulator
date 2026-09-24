"""Does the elastic network of one state point towards the other?

The question Round 2 asks of the IP3R: the ANM of the resting tetramer
(8TKG) predicts its softest collective motions; the activated deposit (8TKF)
records where the receptor actually went. The **overlap** of a mode with the
observed displacement is the |cosine| between them (Tama & Sanejouand 2001);
the **cumulative overlap** of the first k modes is the length of the
displacement's projection onto their span, sqrt(sum of squared overlaps).

Two things make the number interpretable rather than decorative:

* **Rigid body removed.** Modes are orthogonal to the six rigid-body motions,
  so any rigid component left in the displacement by the superposition only
  shrinks every cosine. It is projected out exactly, which also makes the
  result independent of the fit (pore or global).
* **A null.** A random internal direction has an expected squared overlap
  of 1/D with any one mode of a D-dimensional space. Because a C4-imposed
  reconstruction makes the displacement symmetric by construction (8TKG ->
  8TKF is 100.0 % A), the null is taken over random directions **with the
  displacement's own irrep composition**: E[cumulative^2] = sum over irreps
  of f_g k_g / D_g, where k_g counts the retained modes of irrep g and D_g is
  that irrep's internal dimension (3N/4 - 2 for A, less z-translation and
  z-rotation; 3N/4 for B; 3N/2 - 4 for E, less the x/y translations and
  rotations). The report states it beside the measured value.

The symmetry decomposition answers the C4 argument directly: only the
A-symmetric part of the displacement can be reached by A modes, so the
fraction of the transition that is A (``irrep_fraction["A"]``) bounds what
the "lowest A mode is the gating coordinate" claim can explain. For deposits
reconstructed with C4 imposed that fraction is 1 by construction and is a
property of the processing, not of the receptor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from ..structure.transition import Transition
from .anm import ANM, ModeSet, apply_generator
from .network_checks import describe_local, local_modes

__all__ = ["TransitionOverlap", "transition_overlap", "remove_rigid_body",
           "irrep_fractions", "null_cumulative"]


def remove_rigid_body(disp: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """``disp`` less its projection on the six rigid-body motions of ``coords``."""
    x = coords - coords.mean(0)
    n = len(x)
    basis = []
    for a in range(3):
        t = np.zeros((n, 3))
        t[:, a] = 1.0
        basis.append(t.ravel())
        e = np.zeros(3)
        e[a] = 1.0
        basis.append(np.cross(e, x).ravel())
    q, _ = np.linalg.qr(np.array(basis).T)
    d = disp.ravel()
    return (d - q @ (q.T @ d)).reshape(disp.shape)


def irrep_fractions(disp: np.ndarray, axis: np.ndarray, n: int = 4) -> dict[str, float]:
    """Share of |disp|^2 in each C4 isotypic component (A, B, E).

    Projectors: P_A = mean of S^k, P_B = mean of (-1)^k S^k, E the rest.
    """
    powers = [apply_generator(disp, axis, n, k) for k in range(n)]
    a = sum(powers) / n
    b = sum((-1) ** k * p for k, p in enumerate(powers)) / n
    total = float(np.sum(disp * disp))
    fa, fb = float(np.sum(a * a)) / total, float(np.sum(b * b)) / total
    return {"A": fa, "B": fb, "E": max(0.0, 1.0 - fa - fb)}


def _irrep_dims(n_sites: int, n: int = 4) -> dict[str, int]:
    per = 3 * n_sites // n
    return {"A": per - 2, "B": per, "E": 2 * per - 4}


def null_cumulative(symmetry: np.ndarray, fractions: dict, n_sites: int,
                    n: int = 4) -> np.ndarray:
    """Expected cumulative overlap of a random direction of the same irreps.

    Modes labelled 'mixed' are counted against the full internal space.
    """
    dims = _irrep_dims(n_sites, n)
    total = 3 * n_sites - 6
    step = np.array([fractions.get(g, 0.0) / dims[g] if g in dims else 1.0 / total
                     for g in symmetry])
    return np.sqrt(np.cumsum(step))


@dataclass
class TransitionOverlap:
    transition: Transition
    reference: str                  # the deposit whose network was solved
    modes: ModeSet
    residues: np.ndarray            # strided basis residues the ANM used
    stride: int
    overlap: np.ndarray             # |cos| per mode
    cumulative: np.ndarray          # sqrt(cumsum overlap^2)
    null_cumulative: np.ndarray     # random direction, same irrep make-up
    irrep_fraction: dict
    rigid_fraction: float           # share of |disp|^2 that was rigid body

    @property
    def lowest_a(self) -> int | None:
        """The lowest *collective* A mode (local fragments skipped)."""
        return self.modes.first("A")

    @property
    def best(self) -> int:
        return int(np.argmax(self.overlap))

    def report(self) -> list[str]:
        tr, k = self.transition, len(self.overlap)
        a = self.lowest_a
        f = self.irrep_fraction
        lines = [
            f"{tr.start_id} -> {tr.end_id} ({tr.paralog}): {len(tr.residues)} residues "
            f"x {tr.n_subunits} subunits; RMSD {tr.rmsd:.2f} A after the {tr.fit} fit",
            f"network of {self.reference}: {self.modes.meta['n_sites']} sites "
            f"(stride {self.stride}), {k} modes",
            f"displacement by C4 irrep: A {f['A']:.1%}, B {f['B']:.1%}, E {f['E']:.1%}; "
            f"rigid body removed {self.rigid_fraction:.1%}",
            f"best single mode #{self.best + 1} ({self.modes.symmetry[self.best]}), "
            f"overlap {self.overlap[self.best]:.3f}",
        ]
        if a is not None:
            lines.append(f"lowest collective A mode #{a + 1}: overlap {self.overlap[a]:.3f}")
        local = local_modes(self.modes, self.residues, tr.residues)
        if local:
            lines.append("local network artefacts (collectivity below threshold): "
                         + "; ".join(describe_local(local)))
        coll_a = self.modes.is_collective() & (self.modes.symmetry == "A")
        if coll_a.any():
            lines.append(f"collective A modes together ({int(coll_a.sum())}): overlap "
                         f"{np.sqrt(np.sum(self.overlap[coll_a] ** 2)):.3f}; "
                         "the split among them moves with the cutoff, the total "
                         "far less (network_checks.cutoff_scan)")
        lines.append(f"cumulative overlap of {k} modes {self.cumulative[-1]:.3f} "
                     f"(random direction of the same symmetry: "
                     f"{self.null_cumulative[-1]:.3f})")
        return lines


def transition_overlap(tr: Transition, reference: str = "start",
                       stride: int | None = None,
                       n_modes: int | None = None,
                       cutoff: float | None = None) -> TransitionOverlap:
    """Overlap of the ANM of one endpoint with the observed displacement.

    ``reference="start"`` solves the network of the start state and scores
    start -> end; ``"end"`` solves the end state and scores end -> start.
    The ANM sites are the transition's own basis (strided), so the modes and
    the displacement are indexed identically by construction. ``cutoff``
    overrides ``anm.cutoff`` (for the sensitivity scan).
    """
    if reference not in ("start", "end"):
        raise ValueError("reference is 'start' or 'end'")
    stride = int(_P.value("anm.stride")) if stride is None else int(stride)
    n = tr.n_subunits
    m = len(tr.residues)
    keep = np.arange(0, m, max(stride, 1))
    sites = np.concatenate([k * m + keep for k in range(n)])
    x0, x1 = (tr.start, tr.end) if reference == "start" else (tr.end, tr.start)
    x0, disp = x0[sites], (x1 - x0)[sites]
    internal = remove_rigid_body(disp, x0)
    rigid = 1.0 - float(np.sum(internal ** 2) / np.sum(disp ** 2))

    anm = ANM(x0, n_subunits=n, axis=tr.frame.axis,
              cutoff=_P.value("anm.cutoff") if cutoff is None else float(cutoff))
    modes = anm.label_symmetry(anm.calc_modes(n_modes))
    ov = modes.overlap(internal)
    fractions = irrep_fractions(internal, tr.frame.axis, n)
    return TransitionOverlap(
        tr, tr.start_id if reference == "start" else tr.end_id, modes,
        tr.residues[keep], max(stride, 1), ov, np.sqrt(np.cumsum(ov ** 2)),
        null_cumulative(modes.symmetry, fractions, len(x0), n), fractions, rigid)
