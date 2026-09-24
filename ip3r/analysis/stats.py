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
           "signed_rank_test", "mann_whitney_greater", "spearman",
           "fisher_exact", "logistic_fit", "wilson"]


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


def mann_whitney_greater(a, b) -> dict:
    """One-sided Mann-Whitney U: is ``a`` stochastically greater than ``b``?

    Normal approximation with the tie correction and a 0.5 continuity
    correction — the large-sample form, used when both groups hold at least
    ten values (NaN below that). NaN values are dropped.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if min(n1, n2) < 10:
        return {"u": float("nan"), "z": float("nan"), "p": float("nan")}
    ranks = rank_average(np.concatenate([a, b]))
    u = float(ranks[:n1].sum() - n1 * (n1 + 1) / 2.0)
    n = n1 + n2
    _, t = np.unique(np.concatenate([a, b]), return_counts=True)
    var = n1 * n2 / 12.0 * ((n + 1) - float((t ** 3 - t).sum()) / (n * (n - 1)))
    z = (u - n1 * n2 / 2.0 - 0.5) / math.sqrt(var)
    return {"u": u, "z": z, "p": 0.5 * math.erfc(z / math.sqrt(2.0))}


def spearman(x, y) -> dict:
    """Spearman's rho (Pearson on average ranks) with a two-sided p from
    Student's t on n − 2 degrees of freedom. Pairs with a NaN are dropped;
    fewer than four pairs give NaN."""
    from scipy.special import stdtr
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 4:
        return {"rho": float("nan"), "p": float("nan"), "n": n}
    rx, ry = rank_average(x) - (n + 1) / 2.0, rank_average(y) - (n + 1) / 2.0
    den = math.sqrt(float((rx ** 2).sum() * (ry ** 2).sum()))
    if den == 0:
        return {"rho": float("nan"), "p": float("nan"), "n": n}
    rho = float((rx * ry).sum() / den)
    if abs(rho) >= 1.0:
        return {"rho": rho, "p": 0.0, "n": n}
    t = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
    return {"rho": rho, "p": float(2.0 * stdtr(n - 2, -abs(t))), "n": n}


def fisher_exact(table) -> dict:
    """Two-sided Fisher exact test of a 2 × 2 table ``[[a, b], [c, d]]``.

    The p-value sums the hypergeometric probabilities of every table with the
    observed margins that is no more probable than the observed one (with a
    relative tolerance of 1e-7 for ties, as R does). ``odds`` is the sample
    odds ratio ad / bc (inf or 0 at a zero cell).
    """
    (a, b), (c, d) = [[int(v) for v in row] for row in table]
    r1, c1, n = a + b, a + c, a + b + c + d
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)

    def logp(x):
        return (math.lgamma(r1 + 1) - math.lgamma(x + 1) - math.lgamma(r1 - x + 1)
                + math.lgamma(n - r1 + 1) - math.lgamma(c1 - x + 1)
                - math.lgamma(n - r1 - c1 + x + 1)
                - math.lgamma(n + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1))

    obs = logp(a)
    p = sum(math.exp(logp(x)) for x in range(lo, hi + 1)
            if logp(x) <= obs + math.log1p(1e-7))
    odds = (a * d / (b * c)) if b * c else (math.inf if a * d else float("nan"))
    return {"odds": odds, "p": min(1.0, p)}


def logistic_fit(x, y, iterations: int = 100) -> dict:
    """Logistic regression ``logit P(y) = b0 + b1·x`` by iteratively
    reweighted least squares; Wald z and two-sided normal p for the slope.
    ``y`` is 0/1. Returns ``b0``, ``b1``, ``se``, ``z``, ``p`` and
    ``converged``."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    a = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2)
    converged = False
    for _ in range(iterations):
        p = 1.0 / (1.0 + np.exp(-(a @ beta)))
        info = a.T @ (a * (p * (1.0 - p))[:, None])
        step = np.linalg.solve(info, a.T @ (y - p))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-10:
            converged = True
            break
    p = 1.0 / (1.0 + np.exp(-(a @ beta)))
    cov = np.linalg.inv(a.T @ (a * (p * (1.0 - p))[:, None]))
    se = math.sqrt(cov[1, 1])
    z = float(beta[1] / se)
    return {"b0": float(beta[0]), "b1": float(beta[1]), "se": se, "z": z,
            "p": math.erfc(abs(z) / math.sqrt(2.0)), "converged": converged}


def wilson(k: int, n: int, alpha: float) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion k/n at level
    1 − alpha."""
    from scipy.special import ndtri
    z = float(ndtri(1.0 - alpha / 2.0))
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half
