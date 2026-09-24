"""Pairwise protein alignment, for carrying residue numbers between paralogs.

S22 measured its ligand shells on human ITPR3 and carried them to ITPR1 and
ITPR2 through S17's ``transfer_positions`` (a MAFFT pairwise alignment). This
module is the independent route: a global Gotoh (1982) alignment with affine
gaps and BLOSUM62 (Henikoff & Henikoff 1992), end gaps free, so a longer
N- or C-terminus in one paralog costs nothing. Gap costs are registered
parameters (``align.gap_open``, ``align.gap_extend``; EMBOSS needle's).

The recursion runs one row at a time with the horizontal-gap term taken as a
running maximum, so the whole ~2,700 × 2,700 matrix is filled in under a
second without compiled code. A gap of length L costs ``open + extend·(L−1)``.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["BLOSUM62", "align", "transfer_map", "paralog_transfer"]

_AA = "ARNDCQEGHILKMFPSTWYVBZX*"
_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""
#: BLOSUM62 (NCBI's table), indexed by ``_AA``.
BLOSUM62 = np.array([[int(x) for x in r.split()] for r in _ROWS.strip().splitlines()],
                    dtype=float)

_NEG = -1e18


def _encode(seq: str) -> np.ndarray:
    lut = {c: i for i, c in enumerate(_AA)}
    return np.array([lut.get(c, lut["X"]) for c in seq.upper()], dtype=np.intp)


def align(a: str, b: str, gap_open: float | None = None,
          gap_extend: float | None = None) -> tuple[float, list[tuple[int, int]]]:
    """Global alignment of ``a`` and ``b`` with free end gaps.

    Returns ``(score, pairs)``; ``pairs`` are the aligned residue pairs as
    1-based ``(i in a, j in b)``, in order. Gapped residues are absent.
    """
    go = _P.value("align.gap_open") if gap_open is None else gap_open
    ge = _P.value("align.gap_extend") if gap_extend is None else gap_extend
    ea, eb = _encode(a), _encode(b)
    n, m = len(ea), len(eb)
    sub = BLOSUM62[ea][:, eb]                       # n × m substitution scores
    ks = ge * np.arange(m + 1)

    H = np.zeros(m + 1)                             # row 0: free leading gaps
    F = np.full(m + 1, _NEG)
    # traceback: H = E if h_horz else Hp; Hp = F if hp_vert else diagonal
    h_horz = np.zeros((n + 1, m + 1), bool)
    hp_vert = np.zeros((n + 1, m + 1), bool)
    f_open = np.zeros((n + 1, m + 1), bool)
    e_open = np.zeros((n + 1, m + 1), bool)
    last_col = np.zeros(n + 1)
    steps = ge * np.arange(m)
    for i in range(1, n + 1):
        f_open[i] = (H - go) >= (F - ge)
        F = np.maximum(H - go, F - ge)
        Hp = np.empty(m + 1)
        Hp[0] = 0.0                                 # column 0: free leading gaps
        Hp[1:] = H[:-1] + sub[i - 1]
        vert = F > Hp
        vert[0] = False
        Hp = np.where(vert, F, Hp)
        # E[j] = max_{k<j} (Hp[k] − go − ge·(j−1−k)), as a running maximum
        E = np.full(m + 1, _NEG)
        E[1:] = np.maximum.accumulate(Hp + ks)[:-1] - go - steps
        e_open[i, 1:] = Hp[:-1] - go >= E[:-1] - ge
        horz = E > Hp
        H = np.where(horz, E, Hp)
        h_horz[i], hp_vert[i] = horz, vert
        last_col[i] = H[m]

    # free trailing gaps: best cell on the last row or last column
    j_best = int(np.argmax(H))
    i_best = int(np.argmax(last_col))
    if H[j_best] >= last_col[i_best]:
        i, j, score = n, j_best, float(H[j_best])
    else:
        i, j, score = i_best, m, float(last_col[i_best])

    pairs: list[tuple[int, int]] = []
    state = "H"
    while i > 0 and j > 0:
        if state == "H":
            state = "E" if h_horz[i, j] else "Hp"
        elif state == "Hp":
            if hp_vert[i, j]:
                state = "F"
            else:
                pairs.append((i, j))
                i, j, state = i - 1, j - 1, "H"
        elif state == "F":
            state = "H" if f_open[i, j] else "F"
            i -= 1
        else:                                       # "E": opened from Hp
            state = "Hp" if e_open[i, j] else "E"
            j -= 1
    return score, pairs[::-1]


def transfer_map(src: str, dst: str) -> dict[int, int]:
    """Residue ``src`` number → residue ``dst`` number, aligned pairs only."""
    if src == dst:
        return {i: i for i in range(1, len(src) + 1)}
    return dict(align(src, dst)[1])


@lru_cache(maxsize=None)
def _cached(src: str, dst: str, go: float, ge: float) -> tuple[tuple[int, int], ...]:
    return tuple(align(src, dst, go, ge)[1])


def paralog_transfer(src: str, dst: str) -> dict[int, int]:
    """``transfer_map`` between two human paralogs' UniProt sequences,
    memoised on the sequences and the registered gap costs."""
    from .annotations import reference_sequence
    a, b = reference_sequence(src), reference_sequence(dst)
    if a == b:
        return transfer_map(a, b)
    return dict(_cached(a, b, _P.value("align.gap_open"), _P.value("align.gap_extend")))
