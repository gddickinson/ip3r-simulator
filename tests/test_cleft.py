"""The junctional cleft and the sparks it carries.

* The solve conserves Ca2+: what each open channel releases leaves across
  the edges, to rounding. The grid is converged.
* The edge coefficient is Stern's Eq. 13, by hand.
* The couplings stand beside Stern's Fig. 9, read from the figure (per pA:
  about 185 µM at the source, 75 across the row, 73 at 30 nm along it).
  The solve sits below, most at 30 nm, so the tolerance says how far.
* The chessboard: a C channel's nearest C neighbours are the two diagonal
  ones, 42 nm away, and the V sites in between release nothing.
* The Gillespie array, uncoupled at a clamped Ca2+, reproduces the
  stationary P_open, so its sampling is right.
* Sparks in the cleft end by inactivation, far sooner than the mean-field
  cluster's, and uncoupled channels never make a multi-channel event.
"""

from dataclasses import replace

import numpy as np
import pytest

from ip3r.physics import cleft as cl
from ip3r.physics import puff_compare as pc
from ip3r.physics import sparks_cleft as sc
from ip3r.physics.ryr_gating import open_probability


@pytest.fixture(scope="module")
def geom():
    return cl.CleftGeometry.from_parameters()


def _edge_outflow(g, u):
    nx, ny, dx, dy = cl._grid(g)
    grid = u.reshape(nx, ny)
    k = cl.edge_transfer(g)
    return ((grid[0, :].sum() + grid[-1, :].sum()) * dy / (dx / 2 + 1 / k)
            + (grid[:, 0].sum() + grid[:, -1].sum()) * dx / (dy / 2 + 1 / k))


def test_release_leaves_across_the_edges(geom):
    u, src = cl._solved(geom)
    j = 7
    assert src[:, j].sum() == pytest.approx(1.0)
    assert _edge_outflow(geom, u[:, j] / cl._scale(geom)) == pytest.approx(1.0, rel=1e-9)


def test_grid_is_converged(geom):
    fine = cl.coupling_matrix(replace(geom, grid=0.5))
    assert fine == pytest.approx(cl.coupling_matrix(geom), rel=0.01)


def test_edge_coefficient_is_eq_13(geom):
    assert cl.edge_transfer(geom) == pytest.approx(
        2 * np.pi / (15.0 * np.log(1000.0 / 15.0)))


def test_couplings_beside_sterns_fig9(geom):
    g = replace(geom, current=1.0)
    u, _ = cl._solved(g)
    x0, y0 = cl.couplon(g).c_sites[7]
    pts = np.array([[x0, y0], [x0, y0 + 30], [x0 + 30, y0]])
    own, across, along = cl._sample(g, u[:, [7]], pts).ravel()
    assert own == pytest.approx(185, rel=0.15)
    assert across == pytest.approx(75, rel=0.2)
    assert along == pytest.approx(73, rel=0.3)
    assert own > across > along


def test_chessboard_neighbours(geom):
    cp = cl.couplon(geom)
    assert cp.is_c.sum() == 30 and len(cp.sites) == 60
    c = cp.c_sites
    d = np.linalg.norm(c[:, None] - c[None], axis=-1) + np.eye(len(c)) * 1e9
    assert d.min() == pytest.approx(30 * np.sqrt(2))
    g = cl.coupling_matrix(geom)
    i = 7
    assert np.argmax(g[i] - np.diag(np.diag(g))[i]) == np.argmin(d[i])


def test_uncoupled_array_samples_the_stationary_open_probability(monkeypatch):
    monkeypatch.setattr(sc, "couplings_for", lambda pp: np.zeros((30, 30)))
    pp = sc.CleftSparkParams(ca_rest=5.0)
    tr = sc.simulate_sparks_cleft(0.0, 40.0, 3, pp)
    assert tr.n_open.mean() / 30 == pytest.approx(float(open_probability(5.0)), rel=0.05)


def test_the_same_seed_gives_the_same_array():
    a = sc.simulate_sparks_cleft(0.0, 2.0, 5)
    b = sc.simulate_sparks_cleft(0.0, 2.0, 5)
    assert np.array_equal(a.peaks, b.peaks)
    assert np.array_equal(a.n_inactivated, b.n_inactivated)


def test_cleft_sparks_end_by_inactivation():
    ends = [e for seed in range(3)
            for e in pc.spark_ends(pc.simulate(pc.SPARK_CLEFT, 0.0, 20.0, seed))]
    assert len(ends) >= 20
    assert np.median([e["duration"] for e in ends]) < 0.04        # mean field: ~130 ms
    assert np.median([e["inactivated_start"] for e in ends]) <= 2
    assert np.median([e["inactivated_end"] for e in ends]) >= 10


def test_uncoupled_cleft_makes_only_blips():
    r = pc.coupling_effect(0.0, 20.0, 0, model=pc.SPARK_CLEFT)
    assert r["uncoupled"]["multi"] == 0 and r["coupled"]["large"] > 0
