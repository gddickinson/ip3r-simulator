"""Small statistics the checks re-derive published numbers with.

Written here rather than imported from the publication project, so that a
re-derived number and a published one agreeing is two implementations
agreeing. Each function is calibrated in ``tests/test_stats.py`` on inputs
whose answer is known exactly.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["auc", "rank_average", "mean_by_group", "sign_test",
           "signed_rank_test"]


def rank_average(values: np.ndarray) -> np.ndarray:
    """1-based ranks with ties given their average rank."""
    v = np.asarray(values, float)
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(len(v))
    sv = v[order]
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def auc(positive, negative) -> float:
    """Area under the ROC curve = P(score_pos > score_neg) + 0.5 P(tie).

    Computed from the Mann-Whitney U statistic (Hanley & McNeil 1982), so
    ties count one half. NaN scores are dropped; an empty class gives NaN.
    """
    pos = np.asarray(positive, float)
    neg = np.asarray(negative, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    ranks = rank_average(np.concatenate([pos, neg]))
    u = ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2.0
    return float(u / (len(pos) * len(neg)))


def mean_by_group(keys, values) -> dict:
    """``key -> (mean, n)`` over finite values."""
    out: dict = {}
    for k, v in zip(keys, values):
        if v is None or not np.isfinite(v):
            continue
        s, n = out.get(k, (0.0, 0))
        out[k] = (s + v, n + 1)
    return {k: (s / n, n) for k, (s, n) in out.items()}


def sign_test(diffs) -> dict:
    """Exact two-sided sign test; zero differences are dropped and counted."""
    d = np.asarray(diffs, float)
    pos, neg = int((d > 0).sum()), int((d < 0).sum())
    n, k = pos + neg, min(pos, neg)
    if n == 0:
        return {"n_pos": 0, "n_neg": 0, "n_ties": int((d == 0).sum()), "p": float("nan")}
    # log-space binomial tail, so n in the hundreds does not overflow a float
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
            - n * math.log(2) for i in range(k + 1)]
    top = max(logs)
    tail = math.exp(top) * sum(math.exp(x - top) for x in logs)
    return {"n_pos": pos, "n_neg": neg, "n_ties": int((d == 0).sum()),
            "p": min(1.0, 2.0 * tail)}


def signed_rank_test(diffs) -> dict:
    """Wilcoxon signed-rank, two-sided, zeros dropped, normal approximation
    with the tie correction and no continuity correction.

    The approximation is the large-sample form; it is used for n >= 20 only
    (below that it returns NaN rather than a wrong p).
    """
    d = np.asarray(diffs, float)
    d = d[np.isfinite(d) & (d != 0)]
    n = len(d)
    if n < 20:
        return {"n": n, "w_plus": float("nan"), "z": float("nan"), "p": float("nan")}
    ranks = rank_average(np.abs(d))
    w_plus = float(ranks[d > 0].sum())
    mean = n * (n + 1) / 4.0
    _, counts = np.unique(np.abs(d), return_counts=True)
    var = n * (n + 1) * (2 * n + 1) / 24.0 - float((counts ** 3 - counts).sum()) / 48.0
    z = (w_plus - mean) / math.sqrt(var)
    return {"n": n, "w_plus": w_plus, "z": z,
            "p": math.erfc(abs(z) / math.sqrt(2.0))}
