"""Calibrate the paralog-transfer aligner against a plain reference DP."""

import numpy as np
import pytest

from ip3r.core.annotations import functional_sites
from ip3r.core.pairwise import BLOSUM62, align, paralog_transfer, transfer_map

_AA = "ARNDCQEGHILKMFPSTWYV"
_IDX = {c: i for i, c in enumerate("ARNDCQEGHILKMFPSTWYVBZX*")}


def _reference_score(a: str, b: str, go: float, ge: float) -> float:
    """Gotoh's three-state recursion cell by cell, end gaps free."""
    n, m, neg = len(a), len(b), -1e18
    M = [[neg] * (m + 1) for _ in range(n + 1)]
    X = [[neg] * (m + 1) for _ in range(n + 1)]     # a[i] against a gap
    Y = [[neg] * (m + 1) for _ in range(n + 1)]     # b[j] against a gap
    for i in range(n + 1):
        M[i][0] = 0.0
    for j in range(m + 1):
        M[0][j] = 0.0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = BLOSUM62[_IDX[a[i - 1]], _IDX[b[j - 1]]]
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + s
            X[i][j] = max(M[i - 1][j] - go, Y[i - 1][j] - go, X[i - 1][j] - ge)
            Y[i][j] = max(M[i][j - 1] - go, X[i][j - 1] - go, Y[i][j - 1] - ge)
    best = lambda i, j: max(M[i][j], X[i][j], Y[i][j])  # noqa: E731
    return max(max(best(n, j) for j in range(m + 1)), max(best(i, m) for i in range(n + 1)))


def _score_pairs(a, b, pairs, go, ge) -> float:
    """Score an alignment from its pairs: interior gaps cost, end gaps free."""
    s = sum(BLOSUM62[_IDX[a[i - 1]], _IDX[b[j - 1]]] for i, j in pairs)
    for (i0, j0), (i1, j1) in zip(pairs, pairs[1:]):
        for gap in (i1 - i0 - 1, j1 - j0 - 1):
            if gap > 0:
                s -= go + ge * (gap - 1)
    return s


@pytest.mark.parametrize("seed", range(12))
def test_score_matches_reference_dp(seed):
    rng = np.random.default_rng(seed)
    a = "".join(rng.choice(list(_AA), rng.integers(5, 30)))
    b = "".join(rng.choice(list(_AA), rng.integers(5, 30)))
    for go, ge in ((10.0, 0.5), (4.0, 1.0)):
        score, pairs = align(a, b, go, ge)
        assert score == pytest.approx(_reference_score(a, b, go, ge))
        # and the traceback is the alignment that score belongs to
        assert _score_pairs(a, b, pairs, go, ge) == pytest.approx(score)


def test_identity_and_planted_insertion():
    rng = np.random.default_rng(7)
    a = "".join(rng.choice(list(_AA), 120))
    assert transfer_map(a, a) == {i: i for i in range(1, 121)}
    b = a[:60] + "GGSGGSG" + a[60:]               # 7 residues inserted after 60
    m = transfer_map(a, b)
    assert m[60] == 60 and m[61] == 68 and m[120] == 127
    assert len(m) == 120


def test_leading_extension_is_free():
    rng = np.random.default_rng(3)
    a = "".join(rng.choice(list(_AA), 80))
    m = transfer_map(a, "".join(rng.choice(list(_AA), 40)) + a)
    assert all(m[i] == i + 40 for i in range(1, 81))


def test_itpr3_contacts_land_on_itpr1_contacts():
    """The ten IP3 contacts carried ITPR3 → ITPR1 are ITPR1's ten."""
    move = paralog_transfer("ITPR3", "ITPR1")
    c3 = functional_sites("ITPR3")["ip3_contact"]
    assert tuple(move[r] for r in c3) == functional_sites("ITPR1")["ip3_contact"]
