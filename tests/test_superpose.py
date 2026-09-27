"""Round 7.20: a second deposit drawn on the one shown.

Calibrated as the transition fit was: a deposit superposed on a rigidly moved,
chain-shuffled copy of itself must land on the copy atom for atom — side
chains and ligands included, not only the fitted C-alphas — with each
subunit relabelled to the one it sits on.
"""

from __future__ import annotations

import numpy as np
import pytest

from conftest import needs_structure
from ip3r.io.loader import load
from ip3r.structure.superpose import candidates, superpose
from ip3r.structure.symmetry import rotation_matrix
from ip3r.structure.transition import TransitionUnavailable


def _moved_copy(st, shift=1):
    """``st`` turned, shifted and with its protein chains relabelled cyclically
    (so no correspondence is 'the deposited labels' by accident)."""
    r = rotation_matrix(np.array([0.3, -0.5, 0.8]), 1.1)
    moved = st.copy_with_coords(st.xyz.astype(np.float64) @ r.T + [40.0, -25.0, 60.0],
                                name="COPY")
    ch = list(st.chains)
    relabel = dict(zip(ch, ch[shift:] + ch[:shift]))
    moved.chain = np.array([relabel[c] for c in st.chain])
    return moved, relabel


@needs_structure("8TKG")
def test_a_moved_copy_lands_on_the_deposit_atom_for_atom():
    st = load("8TKG")
    copy, relabel = _moved_copy(st)
    sup = superpose(st, copy)
    assert sup.transition.fit_rmsd < 1e-3
    # Every atom, including side chains and IP3, back where it was deposited.
    assert np.abs(sup.structure.xyz - st.xyz).max() < 1e-2
    # ...and on the subunit it sits on, whatever the copy called it.
    assert (sup.structure.chain == st.chain).all()
    assert sup.chain_map == {v: k for k, v in relabel.items()}


@needs_structure("8TKG", "8TKF")
def test_the_overlay_is_the_transitions_end():
    """The drawn C-alphas of the basis are exactly the transition's end."""
    a, b = load("8TKG"), load("8TKF")
    sup = superpose(a, b)
    tr = sup.transition
    st = sup.structure
    ca = (st.atom_name == "CA") & ~st.hetero
    pos = {int(r): i for i, r in enumerate(tr.residues)}
    m = len(tr.residues)
    for k, ch in enumerate(tr.chains_start):
        sel = np.flatnonzero(ca & (st.chain == ch) & np.isin(st.res_seq, tr.residues))
        rows = [k * m + pos[int(r)] for r in st.res_seq[sel]]
        assert np.abs(st.xyz[sel] - tr.end[rows]).max() < 1e-3
    assert st.n_atoms == b.n_atoms          # the deposit itself, every atom
    assert "8TKF on 8TKG (pore fit" in sup.summary()


@needs_structure("8TKG", "8TKF")
def test_the_fit_choice_is_carried():
    a, b = load("8TKG"), load("8TKF")
    pore, whole = superpose(a, b, "pore"), superpose(a, b, "global")
    assert whole.transition.rmsd <= pore.transition.rmsd
    assert not np.allclose(pore.structure.xyz, whole.structure.xyz)


@needs_structure("7LHF", "8TKF")
def test_a_deposit_in_no_human_numbering_is_refused():
    with pytest.raises(TransitionUnavailable, match="numbering"):
        superpose(load("8TKF"), load("7LHF"))


def test_candidates_are_the_paralogs_other_deposits():
    got = [pid for pid, _ in candidates("8TKG", "ITPR3")]
    assert "8TKF" in got and "8TKG" not in got
    assert all(pid not in got for pid, _ in candidates("9HEO", "RYR1"))
    assert candidates("7LHF", None) == []
