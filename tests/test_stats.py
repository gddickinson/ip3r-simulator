"""Calibrate the statistics on inputs whose answer is known exactly."""

import numpy as np
import pytest

from ip3r.analysis.stats import (auc, mean_by_group, rank_average, sign_test,
                                 signed_rank_test)


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


def test_sign_test_exact_small_case():
    # five positives, no negatives: P = 2 * (1/2)^5
    r = sign_test([1, 2, 3, 4, 5, 0])
    assert (r["n_pos"], r["n_neg"], r["n_ties"]) == (5, 0, 1)
    assert abs(r["p"] - 2 / 32) < 1e-12
    assert sign_test([1, -1])["p"] == 1.0


def test_sign_test_large_n_does_not_overflow():
    r = sign_test([1.0] * 223 + [-1.0] * 32)
    assert 0 < r["p"] < 1e-30


def test_signed_rank_symmetric_is_null():
    r = signed_rank_test(list(range(1, 31)) + [-x for x in range(1, 31)])
    assert abs(r["z"]) < 1e-12 and abs(r["p"] - 1.0) < 1e-12


def test_signed_rank_matches_exact_null_variance():
    # the variance formula against the permutation distribution of W+,
    # with ties: every sign pattern equally likely under the null
    rng = np.random.default_rng(3)
    mags = rng.integers(1, 6, 22).astype(float)
    ranks = rank_average(mags)
    signs = rng.choice([-1.0, 1.0], size=(200_000, len(mags)))
    w = ((signs > 0) * ranks).sum(1)
    n = len(mags)
    _, c = np.unique(mags, return_counts=True)
    var = n * (n + 1) * (2 * n + 1) / 24 - (c ** 3 - c).sum() / 48
    assert abs(w.var() / var - 1) < 0.02


def test_signed_rank_too_few_is_nan():
    assert np.isnan(signed_rank_test([1, 2, -3])["p"])


def test_mann_whitney_greater_matches_scipy_asymptotic():
    from scipy.stats import mannwhitneyu
    from ip3r.analysis.stats import mann_whitney_greater
    rng = np.random.default_rng(5)
    for n1, n2 in ((12, 2700), (14, 40), (59, 300)):
        a = np.round(rng.normal(0.1, 1, n1), 1)        # rounded: plenty of ties
        b = np.round(rng.normal(0.0, 1, n2), 1)
        ref = mannwhitneyu(a, b, alternative="greater", method="asymptotic").pvalue
        assert mann_whitney_greater(a, b)["p"] == pytest.approx(ref, rel=1e-9)
    assert np.isnan(mann_whitney_greater([1.0] * 5, [0.0] * 50)["p"])


def test_spearman_matches_scipy():
    from scipy.stats import spearmanr
    from ip3r.analysis.stats import spearman
    rng = np.random.default_rng(6)
    x = rng.uniform(0, 15, 125)
    y = np.round(-0.01 * x + rng.normal(0, 0.1, 125), 2)
    ref = spearmanr(x, y)
    got = spearman(x, y)
    assert got["rho"] == pytest.approx(ref.statistic, abs=1e-12)
    assert got["p"] == pytest.approx(ref.pvalue, rel=1e-9)
    assert got["n"] == 125
    assert spearman([1, 2, 3, 4], [4, 3, 2, 1])["rho"] == -1.0
