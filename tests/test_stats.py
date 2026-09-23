"""Calibrate the statistics on inputs whose answer is known exactly."""

import numpy as np

from ip3r.analysis.stats import auc, mean_by_group, rank_average


def test_auc_perfect_and_reversed():
    assert auc([3, 4, 5], [0, 1, 2]) == 1.0
    assert auc([0, 1, 2], [3, 4, 5]) == 0.0


def test_auc_ties_count_one_half():
    assert auc([1, 1], [1, 1]) == 0.5
    # one pos beats one neg, ties the other: (1 + 0.5) / 2
    assert auc([2], [1, 2]) == 0.75


def test_auc_matches_brute_force():
    rng = np.random.default_rng(1)
    p, n = rng.integers(0, 5, 40), rng.integers(0, 5, 30)
    brute = np.mean([(a > b) + 0.5 * (a == b) for a in p for b in n])
    assert abs(auc(p, n) - brute) < 1e-12


def test_auc_empty_class_is_nan():
    assert np.isnan(auc([], [1.0]))


def test_rank_average_ties():
    assert list(rank_average([10, 20, 20, 30])) == [1.0, 2.5, 2.5, 4.0]


def test_mean_by_group_skips_nan():
    out = mean_by_group(["a", "a", "b"], [1.0, float("nan"), 3.0])
    assert out == {"a": (1.0, 1), "b": (3.0, 1)}
