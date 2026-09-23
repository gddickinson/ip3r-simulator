"""Small statistics the checks re-derive published numbers with.

Written here rather than imported from the publication project, so that a
re-derived number and a published one agreeing is two implementations
agreeing. Each function is calibrated in ``tests/test_stats.py`` on inputs
whose answer is known exactly.
"""

from __future__ import annotations

import numpy as np

__all__ = ["auc", "rank_average", "mean_by_group"]


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
