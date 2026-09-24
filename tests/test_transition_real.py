"""8TKG (resting) -> 8TKF (activated), measured (skipped when not downloaded).

The endpoint tests are the PIEZO1 lesson: that simulator's morph passed every
interpolation test while the drawn end landed 36 Å from the deposit it
claimed, because the path and the picture were in different frames. So the
drawn end is compared with the **deposited** 8TKF as a shape (optimal
superposition, read by a route the transition does not use), with a case
that must fail beside it.
"""

import numpy as np
import pytest

from conftest import needs_structure
from ip3r.io.loader import load
from ip3r.structure.symmetry import kabsch, subunit_ca
from ip3r.structure.transition import (TransitionUnavailable, atom_displacement,
                                       atom_site_index, displaced_coords,
                                       prepare_transition)


@pytest.fixture(scope="module")
def pair():
    return load("8TKG"), load("8TKF")


@pytest.fixture(scope="module")
def tr(pair):
    return prepare_transition(*pair)


def _shape_rmsd(a, b):
    r, t = kabsch(a, b)
    return float(np.sqrt((((a @ r.T + t) - b) ** 2).sum(1).mean()))


@needs_structure("8TKG", "8TKF")
def test_basis_and_correspondence(tr):
    assert tr.paralog == "ITPR3" and len(tr.residues) == 2194
    assert not tr.meta["correspondence_determined"]      # C4: all four fit alike
    assert np.ptp(tr.meta["fit_rmsd_by_shift"]) < 0.05
    assert tr.fit_rmsd < 4.0 < tr.rmsd
    m = tr.element_means()
    assert m["channel"] < 3.0 < 15.0 < m["RIH_N"]        # pore still, cap moves


@needs_structure("8TKG", "8TKF")
def test_drawn_path_starts_on_screen_and_ends_on_8tkf(pair, tr):
    from ip3r.structure.morph import morph
    st, end_st = pair
    traj = morph(tr.start, tr.end)
    idx = atom_site_index(st, tr)
    first = displaced_coords(st.xyz, traj.frames, 0, idx)
    last = displaced_coords(st.xyz, traj.frames, len(traj) - 1, idx)
    assert np.array_equal(first, st.xyz.astype(np.float64))
    ca_end = subunit_ca(end_st)
    deposited = np.vstack([[ca_end[c][r] for r in tr.residues] for c in tr.chains_end])
    drawn_ca = []
    for c in tr.chains_start:
        m = st.mask_ca() & (st.chain == c)
        at = dict(zip(st.res_seq[m].tolist(), np.flatnonzero(m)))
        drawn_ca.append(last[[at[r] for r in tr.residues]])
    drawn_ca = np.vstack(drawn_ca)
    assert _shape_rmsd(drawn_ca, deposited) < 1e-3
    # the case that must fail: the start drawn unmoved is 14 Å from 8TKF
    assert _shape_rmsd(tr.start, deposited) > 10.0


@needs_structure("8TKG", "8TKF")
def test_off_basis_atoms_ride_the_nearest_site(pair, tr):
    st = pair[0]
    idx = atom_site_index(st, tr)
    ip3 = st.res_name == "I3P"
    assert ip3.any()
    gap = np.linalg.norm(st.xyz[ip3] - tr.start[idx[ip3]], axis=1)
    assert gap.max() < 12.0
    disp = atom_displacement(st, tr, idx)
    assert np.isnan(disp[ip3]).all()                     # ridden, not measured
    assert np.isfinite(disp[st.mask_ca() & np.isin(st.res_seq, tr.residues)]).all()


@needs_structure("8TKG", "8TKF")
def test_lowest_collective_a_mode_points_toward_opening(tr):
    from ip3r.physics.transition_modes import transition_overlap
    ov = transition_overlap(tr)
    a = ov.lowest_a
    assert ov.irrep_fraction["A"] > 0.99                 # C4-imposed maps
    assert ov.overlap[a] > 10 * ov.null_cumulative[a]
    assert 0.3 < ov.overlap[a] < 0.55
    assert ov.cumulative[-1] > 0.55


@needs_structure("8TKG", "8TKF")
def test_overlap_independent_of_fit(pair, tr):
    from ip3r.physics.transition_modes import transition_overlap
    other = prepare_transition(*pair, fit="global")
    a, b = transition_overlap(tr), transition_overlap(other)
    assert np.allclose(a.overlap, b.overlap, atol=0.02)


@needs_structure("8TKG", "8TKF")
def test_local_fragment_is_not_called_the_lowest_a_mode(tr):
    from ip3r.physics.transition_modes import transition_overlap
    ov = transition_overlap(tr, stride=4)
    naive = ov.modes.first("A", collective=False)
    assert ov.modes.collectivity()[naive] < 0.05         # a fragment ...
    assert ov.lowest_a != naive                          # ... is skipped
    assert ov.overlap[ov.lowest_a] > 0.3


@needs_structure("7LHF", "8TKF")
def test_rat_deposit_refused():
    with pytest.raises(TransitionUnavailable):
        prepare_transition(load("7LHF"), load("8TKF"))
