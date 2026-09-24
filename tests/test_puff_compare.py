"""One ruler for both clusters: the recruitment statistics on a synthetic
trace, and the measured contrast between the two receptors."""

import numpy as np
import pytest

from ip3r.physics.puff_compare import (coupling_effect, recruitment,
                                       scan_couplings, simulate, params_for)
from ip3r.physics.puffs import PuffParams, PuffTrace


def _trace(peaks, n=20):
    pp = PuffParams()
    pp.n_channels = n
    peaks = np.asarray(peaks)
    return PuffTrace(np.arange(len(peaks)) * 1e-3, peaks, np.zeros(len(peaks)),
                     pp, 0.2, peaks)


def test_recruitment_counts_by_hand():
    # events: [1] [2] [3,12,4] [1,1] -> sizes 1, 2, 12, 1
    r = recruitment(_trace([1, 0, 2, 0, 3, 12, 4, 0, 1, 1, 0]))
    assert (r["n_events"], r["blips"], r["multi"], r["large"]) == (4, 2, 2, 1)
    assert r["large_at"] == 10 and r["large_share"] == 0.5


def test_scan_starts_uncoupled_and_is_geometric():
    cs = scan_couplings()
    assert cs[0] == 0.0 and np.allclose(cs[2:] / cs[1:-1], cs[2] / cs[1])


def test_park_drive_coupling_clusters_openings():
    out = coupling_effect(0.2, duration=20.0, model="park-drive")
    assert out["coupled"]["fano"] > 2.0 > 1.2 > out["uncoupled"]["fano"] > 0.8


def test_park_drive_is_quieter_at_rest():
    dyk = coupling_effect(0.2, duration=10.0, model="dyk")["uncoupled"]
    pd = coupling_effect(0.2, duration=10.0, model="park-drive")["uncoupled"]
    assert pd["open_fraction"] < dyk["open_fraction"] / 2


@pytest.mark.parametrize("coupling", [0.1, 0.2])
def test_only_park_drive_recruits_the_cluster(coupling):
    # Same seed, same cluster, same coupling: the DYK cluster essentially
    # never reaches half the cluster; the park/drive one repeatedly does.
    r = {m: recruitment(simulate(m, 0.2, 20.0, 0, params_for(m, coupling)))
         for m in ("dyk", "park-drive")}
    assert r["dyk"]["large"] <= 1
    assert r["park-drive"]["large"] >= 5
