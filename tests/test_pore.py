"""The pore profile on a synthetic channel with a known constriction."""

import numpy as np

from ip3r.structure.pore import pore_profile
from ip3r.structure.symmetry import Frame


def _ring_structure():
    from conftest import make_structure
    blocks = []
    z = np.arange(-20, 21, 1.0)
    for k in range(4):
        pts = []
        for zi in z:
            radius = 3.0 if abs(zi) < 2 else 8.0      # constriction at z = 0
            t = np.pi / 2 * k
            pts.append([radius * np.cos(t), radius * np.sin(t), zi])
        blocks.append(np.array(pts))
    st = make_structure(blocks, resname="LEU", atom="CB")
    return st


def test_profile_finds_the_constriction():
    st = _ring_structure()
    frame = Frame(np.array([0, 0, 1.0]), np.zeros(3), np.eye(3), list("ABCD"))
    p = pore_profile(st, frame, -15, 15, step=0.5, slab=1.0)
    i = np.argmin(p.r_min)
    assert abs(p.z[i]) <= 2.0          # the 3 Å ring spans |z| < 2
    assert abs(p.r_min[i] - 3.0) < 1e-6
    assert abs(p.at(10.0) - 8.0) < 1e-6
    # The free radius subtracts carbon's vdW radius (1.70 Å).
    assert abs(p.r_free[i] - 1.3) < 1e-6
