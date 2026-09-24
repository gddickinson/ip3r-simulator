import numpy as np

from ip3r.physics.puff_compare import coupling_effect
from ip3r.physics.puffs import PuffParams, detect_events, simulate_cluster


def test_no_ip3_no_openings():
    assert simulate_cluster(0.0, duration=2.0).n_open.max() == 0


def test_seeded_and_reproducible():
    a = simulate_cluster(0.2, duration=2.0, seed=5)
    b = simulate_cluster(0.2, duration=2.0, seed=5)
    assert np.array_equal(a.n_open, b.n_open)


def test_coupling_clusters_openings():
    out = coupling_effect(0.2, duration=20.0)
    assert out["coupled"]["fano"] > 1.2 > out["uncoupled"]["fano"] > 0.8


def test_event_detection():
    from ip3r.physics.puffs import PuffTrace
    tr = PuffTrace(np.arange(6.0), np.array([0, 1, 3, 0, 1, 0]), np.zeros(6),
                   PuffParams(), 0.2)
    ev = detect_events(tr)
    assert [e["peak_open"] for e in ev] == [3, 1]
    # with per-bin peaks, a bin whose snapshot is 0 still counts
    tr.n_peak = np.array([0, 1, 3, 2, 1, 0])
    assert [e["peak_open"] for e in detect_events(tr)] == [3]


def test_peak_bounds_snapshot():
    tr = simulate_cluster(0.2, duration=2.0)
    assert np.all(tr.n_peak >= tr.n_open)
