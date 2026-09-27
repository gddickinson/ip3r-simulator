"""Round 7.25: the wall-model search's pieces, and 8TKF's ceilings."""

from types import SimpleNamespace

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics.wall_search import (Candidate, depths, regions,
                                      ring_positions, scales, score)


def test_score_by_hand():
    m = (0.27, 15.2)
    assert score(15.2, 0.27, m) == 0.0
    assert score(1.52, 0.027, m) == pytest.approx(2 * np.log(10))
    assert score(0.0, 0.27, m) == float("inf")


def test_grids_follow_their_parameters():
    s = scales()
    assert s[0] == pytest.approx(_P.value("wallsearch.scale_min"))
    assert s[-1] == pytest.approx(1.0)
    assert len(s) == int(_P.value("wallsearch.scale_points"))
    d = depths()
    assert d[0] == pytest.approx(_P.value("wallsearch.depth_step"))
    assert d[-1] == pytest.approx(_P.value("wallsearch.depth_max"))


def test_regions_and_labels():
    cons = {"filter": SimpleNamespace(z=-86.0), "gate": SimpleNamespace(z=-76.0)}
    summ = SimpleNamespace(constrictions=cons, span=(-106.0, -56.0))
    r = regions(summ)
    half = _P.value("lumen.constriction_half_width")
    assert r["span"] == (-106.0, -56.0)
    assert r["filter"] == (-86.0 - half, -86.0 + half)
    assert r["filter to gate"] == (-86.0 - half, -76.0 + half)
    assert Candidate("charge", 0.1).label == "charge x0.1"
    assert Candidate("charge", 0.0).label == "charge x0"
    assert Candidate("well", 0.0, "gate", 4.0).label == "Ca well 4 kT at gate"
    assert "charge x0.5" in Candidate("well + charge", 0.5, "span", 4.0).label


def test_ring_positions_are_c4_pairs_at_the_wall():
    prof = SimpleNamespace(at=lambda z: 6.0)
    p = ring_positions(SimpleNamespace(profile=prof), -80.0)
    assert p.shape == (8, 3)
    assert np.allclose(np.hypot(p[:, 0], p[:, 1]), 6.0)
    assert np.allclose(p[:, 2], -80.0)
    rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])   # 90° about z
    assert np.allclose(p[:4] @ rot.T, np.roll(p[:4], -1, axis=0))
    ang = np.degrees(np.arctan2(p[4, 1], p[4, 0]) - np.arctan2(p[0, 1], p[0, 0]))
    assert ang == pytest.approx(45.0)


@pytest.fixture(scope="module")
def pore_8tkf():
    from ip3r.io import loader
    from ip3r.physics.selectivity3d import prepare
    try:
        st = loader.load("8TKF")
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"8TKF unavailable: {e}")
    return prepare(st, spacing=_P.value("reversal3d.spacing"))


def test_8tkf_well_ceilings_fall_short_of_vais(pore_8tkf):
    """Finding: in linear response even an unbounded Ca2+-only well over
    the whole membrane span lifts P_Ca:P_K by 46x: short of 58, Vais's
    15.2 over the uncharged pore's reading at reversal (0.26). One at a
    single constriction lifts it 1.14x, because the access resistance and the
    rest of the path hold the remainder."""
    from ip3r.physics.selectivity_bound import well_ceiling
    r = regions(pore_8tkf.summary)
    span = well_ceiling(pore_8tkf, *r["span"])
    gate = well_ceiling(pore_8tkf, *r["gate"])
    assert 20.0 < span < 58.0
    assert 1.05 < gate < 3.0


def test_8tkf_deposit_wall_obeys_the_bound(pore_8tkf):
    from ip3r.physics.reversal3d import wall
    from ip3r.physics.selectivity_bound import linear_ratios
    base = linear_ratios(pore_8tkf)
    r = linear_ratios(pore_8tkf, wall(pore_8tkf).fixed)
    a, b = r.pca_pk / base.pca_pk, r.pcl_pk / base.pcl_pk
    assert a > 1.5 and b < 0.1                     # selects Ca2+, shuts Cl-
    assert b * a * a < 1.0
