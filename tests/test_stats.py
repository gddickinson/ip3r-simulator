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


def test_fisher_exact_known_and_scipy():
    from scipy.stats import fisher_exact as ref
    from ip3r.analysis.stats import fisher_exact
    # Fisher's lady tasting tea: 3/1 vs 1/3, two-sided p = 34/70
    assert fisher_exact([[3, 1], [1, 3]])["p"] == pytest.approx(34 / 70, rel=1e-12)
    for t in ([[140, 783], [42, 267]], [[3, 509], [137, 274]], [[0, 172], [42, 95]],
              [[10, 2], [3, 15]]):
        got, want = fisher_exact(t), ref(t)
        assert got["p"] == pytest.approx(want.pvalue, rel=1e-9)
        if np.isfinite(want.statistic) and want.statistic > 0:
            assert got["odds"] == pytest.approx(want.statistic, rel=1e-12)
    assert fisher_exact([[0, 172], [42, 95]])["odds"] == 0.0


def test_logistic_fit_binary_predictor_is_log_odds_ratio():
    from ip3r.analysis.stats import logistic_fit
    # With a 0/1 predictor the MLE slope is exactly log(ad / bc) and its
    # standard error sqrt(1/a + 1/b + 1/c + 1/d).
    a, b, c, d = 30, 10, 12, 28            # y=1|x=1, y=0|x=1, y=1|x=0, y=0|x=0
    x = [1] * (a + b) + [0] * (c + d)
    y = [1] * a + [0] * b + [1] * c + [0] * d
    f = logistic_fit(x, y)
    assert f["converged"]
    assert f["b1"] == pytest.approx(np.log(a * d / (b * c)), rel=1e-9)
    assert f["se"] == pytest.approx(np.sqrt(1 / a + 1 / b + 1 / c + 1 / d), rel=1e-9)
    assert f["b0"] == pytest.approx(np.log(c / d), rel=1e-9)


def test_logistic_fit_recovers_a_planted_slope():
    from ip3r.analysis.stats import logistic_fit
    rng = np.random.default_rng(7)
    x = rng.uniform(3, 8, 20000)
    y = rng.uniform(size=x.size) < 1 / (1 + np.exp(-(-9.0 + 2.0 * x)))
    f = logistic_fit(x, y)
    assert f["b1"] == pytest.approx(2.0, abs=4 * f["se"])
    assert f["p"] < 1e-100 or f["p"] == 0.0


def test_wilson_known_value():
    from ip3r.analysis.stats import wilson
    # Newcombe 1998, Table I: 81/263 → 0.2553–0.3662 at 95 %
    lo, hi = wilson(81, 263, 0.05)
    assert (round(lo, 4), round(hi, 4)) == (0.2553, 0.3662)
    lo, hi = wilson(0, 10, 0.05)
    assert lo == pytest.approx(0.0, abs=1e-15) and hi == pytest.approx(0.2775, abs=1e-4)


def test_mann_whitney_less_and_two_sided_match_scipy_on_zero_inflated_data():
    """FEL's β: mostly exact zeros, so nearly every value is a tie."""
    from scipy.stats import mannwhitneyu

    from ip3r.analysis.stats import mann_whitney_less, mann_whitney_two_sided
    rng = np.random.default_rng(7)
    a = np.where(rng.random(40) < 0.8, 0.0, rng.exponential(0.05, 40))
    b = np.where(rng.random(300) < 0.6, 0.0, rng.exponential(0.05, 300))
    for mine, alt in ((mann_whitney_less(a, b), "less"),
                      (mann_whitney_two_sided(a, b), "two-sided")):
        ref = mannwhitneyu(a, b, alternative=alt, method="asymptotic")
        assert mine["p"] == pytest.approx(ref.pvalue, rel=1e-9)
        assert mine["u"] == pytest.approx(ref.statistic)
    assert mann_whitney_two_sided(a, b)["cles"] == pytest.approx(
        mannwhitneyu(a, b).statistic / (len(a) * len(b)))


def test_benjamini_hochberg_by_hand():
    from ip3r.analysis.stats import benjamini_hochberg
    # m = 4 sorted: 0.01·4/1 = 0.04, 0.02·4/2 = 0.04, 0.03·4/3 = 0.04, 0.5·4/4 = 0.5;
    # input order kept, NaN takes no rank
    q = benjamini_hochberg([0.5, 0.03, float("nan"), 0.01, 0.02])
    assert np.allclose(q[[0, 1, 3, 4]], [0.5, 0.04, 0.04, 0.04])
    assert np.isnan(q[2])
    # step-up: a later small ratio pulls an earlier one down, capped at 1
    assert np.allclose(benjamini_hochberg([0.04, 0.045, 0.9, 0.9]), [0.09, 0.09, 0.9, 0.9])
    assert benjamini_hochberg([0.9, 0.95])[0] == pytest.approx(0.95)
