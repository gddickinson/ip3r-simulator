"""The four-fold axis of a tetramer, found two independent ways.

``ip3r_genes`` S0 found the 6DQN axis as the normal of the plane through the
four subunit centroids (the smallest singular direction of the centroid
cloud). This module computes it **differently** — by superposing one subunit
on its neighbour (Kabsch 1976) and reading the axis off the rotation that
does it — so that agreement between the two is a confirmation rather than a
re-run of the same arithmetic. The centroid method is also provided, because
two routes are only a check if both can be run and compared.

Coordinates in the returned :class:`Frame` put the origin on the axis and z
along it, oriented so that the cytosolic cap (the end with the larger radial
extent) is at +z — the convention S0 used.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.structure import Structure

__all__ = ["kabsch", "rotation_matrix", "rotation_axis_angle", "Frame",
           "subunit_ca", "axis_by_superposition", "axis_by_centroids",
           "c4_residual", "tetramer_frame", "order_subunits"]


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation ``R`` and translation ``t`` minimising ``|R p + t - q|``."""
    p = np.asarray(p, np.float64)
    q = np.asarray(q, np.float64)
    pc, qc = p.mean(0), q.mean(0)
    h = (p - pc).T @ (q - qc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, qc - r @ pc


def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Right-handed rotation by ``angle`` radians about ``axis``."""
    a = np.asarray(axis, np.float64)
    a = a / np.linalg.norm(a)
    k = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * (k @ k)


def rotation_axis_angle(r: np.ndarray) -> tuple[np.ndarray, float]:
    """Axis (unit) and angle (radians, 0..pi) of a rotation matrix."""
    angle = float(np.arccos(np.clip((np.trace(r) - 1.0) / 2.0, -1.0, 1.0)))
    axis = np.array([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]])
    n = np.linalg.norm(axis)
    if n < 1e-9:                                   # identity or a half-turn
        w, v = np.linalg.eigh(r + r.T)
        axis = v[:, np.argmax(w)]
    else:
        axis = axis / n
    return axis, angle


def subunit_ca(st: Structure) -> dict[str, dict[int, np.ndarray]]:
    """Per protein chain, ``residue number -> CA coordinate``."""
    out: dict[str, dict[int, np.ndarray]] = {}
    m = st.mask_ca()
    for ch, rs, xyz in zip(st.chain[m], st.res_seq[m], st.xyz[m]):
        out.setdefault(str(ch), {})[int(rs)] = xyz.astype(np.float64)
    return out


def _largest_chains(ca: dict, n: int) -> list[str]:
    return sorted(sorted(ca, key=lambda c: -len(ca[c]))[:n])


def _shared(a: dict, b: dict) -> tuple[np.ndarray, np.ndarray]:
    keys = sorted(set(a) & set(b))
    return np.array([a[k] for k in keys]), np.array([b[k] for k in keys])


def order_subunits(ca: dict, chains: list[str], axis: np.ndarray,
                   centre: np.ndarray) -> list[str]:
    """Chains in right-handed order about ``axis`` (by centroid azimuth)."""
    e1 = np.cross(axis, [1.0, 0.0, 0.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(axis, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    az = {}
    for ch in chains:
        c = np.mean(list(ca[ch].values()), axis=0) - centre
        az[ch] = np.arctan2(c @ e2, c @ e1)
    first = min(chains)
    return sorted(chains, key=lambda c: (az[c] - az[first]) % (2 * np.pi))


def axis_by_superposition(ca: dict, n: int = 4) -> tuple[np.ndarray, np.ndarray, float]:
    """Axis from the rotation mapping each subunit onto its neighbour.

    Returns ``(axis, point_on_axis, mean_angle_degrees)``. The angle is
    reported because it is itself a test: for a C4 tetramer it must be 90°.
    """
    chains = _largest_chains(ca, n)
    cents = np.array([np.mean(list(ca[c].values()), axis=0) for c in chains])
    centre = cents.mean(0)
    # Provisional axis only to order the chains; replaced below.
    _, _, vt = np.linalg.svd(cents - centre)
    order = order_subunits(ca, chains, vt[2], centre)
    axes, angles = [], []
    for k in range(n):
        p, q = _shared(ca[order[k]], ca[order[(k + 1) % n]])
        r, _ = kabsch(p, q)
        ax, ang = rotation_axis_angle(r)
        if axes and ax @ axes[0] < 0:
            ax = -ax
        axes.append(ax)
        angles.append(np.degrees(ang))
    axis = np.mean(axes, axis=0)
    axis /= np.linalg.norm(axis)
    return axis, centre, float(np.mean(angles))


def axis_by_centroids(ca: dict, n: int = 4) -> tuple[np.ndarray, np.ndarray, float]:
    """S0's method: normal of the subunit-centroid plane, plus the residual."""
    chains = _largest_chains(ca, n)
    cents = np.array([np.mean(list(ca[c].values()), axis=0) for c in chains])
    centre = cents.mean(0)
    _, sv, vt = np.linalg.svd(cents - centre)
    return vt[2] / np.linalg.norm(vt[2]), centre, float(sv[2])


@dataclass
class Frame:
    """An axis-aligned frame: z along the four-fold axis, origin on it."""

    axis: np.ndarray            # unit, in deposit coordinates, pointing +z
    centre: np.ndarray          # a point on the axis, deposit coordinates
    basis: np.ndarray           # (3, 3) columns e1, e2, e3
    chains: list                # subunits in right-handed order about +z

    def to_frame(self, xyz: np.ndarray) -> np.ndarray:
        return (np.asarray(xyz, np.float64) - self.centre) @ self.basis

    def from_frame(self, xyz: np.ndarray) -> np.ndarray:
        return np.asarray(xyz, np.float64) @ self.basis.T + self.centre


def _basis(axis: np.ndarray) -> np.ndarray:
    e3 = axis / np.linalg.norm(axis)
    seed = np.array([1.0, 0.0, 0.0])
    if abs(seed @ e3) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(seed, e3)
    e1 /= np.linalg.norm(e1)
    return np.column_stack([e1, np.cross(e3, e1), e3])


def tetramer_frame(st: Structure, n: int = 4, method: str = "superposition") -> Frame:
    """The axis frame of a tetrameric deposit, cytosolic cap at +z."""
    ca = subunit_ca(st)
    finder = axis_by_superposition if method == "superposition" else axis_by_centroids
    axis, centre, _ = finder(ca, n)
    heavy = st.xyz[(st.element != "H") & ~st.hetero].astype(np.float64)
    z = (heavy - centre) @ axis
    rad = np.linalg.norm((heavy - centre) - np.outer(z, axis), axis=1)
    lo, hi = np.percentile(z, [15, 85])
    if rad[z < lo].mean() > rad[z > hi].mean():
        axis = -axis
    chains = order_subunits(ca, _largest_chains(ca, n), axis, centre)
    return Frame(axis=axis, centre=centre, basis=_basis(axis), chains=chains)


def c4_residual(st: Structure, frame: Frame, n: int = 4) -> float:
    """Worst RMSD of subunit 1 turned by k x 90° onto subunit k+1 (Å).

    The definition S0 used: after an *ideal* rotation about the fitted axis,
    how far is each other subunit from the first. Both senses of rotation are
    tried so chain lettering cannot flip the test.
    """
    ca = subunit_ca(st)
    ref = frame.chains[0]
    worst = 0.0
    for k, ch in enumerate(frame.chains[1:], start=1):
        p, q = _shared(ca[ref], ca[ch])
        pf, qf = frame.to_frame(p), frame.to_frame(q)
        best = np.inf
        for sign in (1.0, -1.0):
            r = rotation_matrix([0, 0, 1], sign * 2 * np.pi * k / n)
            best = min(best, float(np.sqrt(((pf @ r.T - qf) ** 2).sum(1).mean())))
        worst = max(worst, best)
    return worst
