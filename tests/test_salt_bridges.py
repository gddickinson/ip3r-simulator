"""Salt bridges: the pairing rule on hand-placed atoms, then 8TKF.

Each synthetic case is built so that a wrong rule gives a different answer:
many-to-one pairing would pair both aspartates with one arginine; a
centre-to-centre distance would miss the Lys at 3.9 Å N-O; ignoring the
registered cutoff would pair the Lys at 4.5 Å.
"""

import numpy as np
import pytest

from ip3r.core.structure import Structure
from ip3r.parameters import PARAMETERS as _P
from ip3r.physics.pore_charge import pore_charge
from ip3r.physics.salt_bridges import salt_bridges
from ip3r.structure.pore import PoreProfile
from ip3r.structure.symmetry import Frame
from conftest import needs_structure


def _structure(atoms):
    """``atoms``: (chain, res_seq, res_name, atom_name, xyz)."""
    n = len(atoms)
    st = Structure(
        xyz=np.array([a[4] for a in atoms], np.float32),
        element=np.array([a[3][0] for a in atoms], "U2"),
        atom_name=np.array([a[3] for a in atoms], "U6"),
        res_name=np.array([a[2] for a in atoms], "U5"),
        res_seq=np.array([a[1] for a in atoms], np.int32),
        chain=np.array([a[0] for a in atoms], "U6"),
        hetero=np.zeros(n, bool), b_factor=np.zeros(n, np.float32),
        occupancy=np.ones(n, np.float32), alt_loc=np.full(n, ".", "U2"),
        entity=np.full(n, "1", "U6"), name="SYN")
    st._build_residue_index()
    return st


def _asp(chain, seq, at):
    at = np.asarray(at, float)
    return [(chain, seq, "ASP", "OD1", at), (chain, seq, "ASP", "OD2", at + [0, 0, 2.2])]


def test_an_arginine_cancels_one_carboxylate_not_two():
    atoms = (_asp("A", 1, [0, 0, 0]) + _asp("A", 2, [6.0, 0, 0])
             + [("A", 3, "ARG", "NH1", np.array([2.5, 0, 0])),
                ("A", 3, "ARG", "NH2", np.array([3.0, 0, 0]))])
    bridges = salt_bridges(_structure(atoms))
    # Both aspartates are within 4 A of the guanidinium; only the closer pairs.
    assert [(b.acid, b.base) for b in bridges] == [(("A", 1), ("A", 3))]
    assert bridges[0].distance == pytest.approx(2.5, abs=1e-5)


def test_the_distance_is_nitrogen_to_oxygen_and_respects_the_cutoff():
    atoms = (_asp("A", 1, [0, 0, 0]) + [("B", 5, "LYS", "NZ", np.array([3.9, 0, 0]))]
             + _asp("C", 1, [30, 0, 0]) + [("D", 5, "LYS", "NZ", np.array([34.5, 0, 0]))])
    st = _structure(atoms)
    assert [b.base for b in salt_bridges(st)] == [("B", 5)]          # 4.5 A is not
    assert len(salt_bridges(st, cutoff=5.0)) == 2
    assert salt_bridges(st, cutoff=3.0) == []
    assert _P.value("pore_charge.salt_bridge_cutoff") == 4.0


def test_a_bridged_lining_group_is_dropped_and_its_partner_is_not_added():
    """An Asp lines a 5 A lumen on each subunit; on A and B an Arg of the
    next subunit, outside the lumen, bridges it. Two charges are left."""
    atoms = []
    for k, ch in enumerate("ABCD"):
        u = np.array([np.cos(np.pi / 2 * k), np.sin(np.pi / 2 * k), 0.0])
        atoms += [(ch, 10, "ASP", "OD1", u * 5.5), (ch, 10, "ASP", "OD2", u * 5.5 + [0, 0, 0.5])]
        if ch in "AB":
            nxt = "ABCD"[(k + 1) % 4]
            atoms += [(nxt, 3, "ARG", "NE", u * 8.3), (nxt, 3, "ARG", "CZ", u * 9.0),
                      (nxt, 3, "ARG", "NH1", u * 8.3 + [0, 0, 1.0])]
    st = _structure(atoms)
    z = np.arange(-20.0, 20.25, 0.5)
    frame = Frame(np.array([0, 0, 1.0]), np.zeros(3), np.eye(3), list("ABCD"))
    prof = PoreProfile(z, np.full_like(z, 5.0), np.full_like(z, 3.5))
    raw = pore_charge(st, frame, prof)
    paired = pore_charge(st, frame, prof, pair_bridges=True)
    assert raw.net_charge == -4.0 and raw.bridged == []
    assert paired.net_charge == -2.0
    assert sorted(g.chain for g in paired.groups) == ["C", "D"]
    assert {b.base for b in paired.bridged} == {("B", 3), ("C", 3)}


@needs_structure("8TKF")
def test_8tkf_filter_aspartates_are_bridged_and_pairing_does_not_rescue_the_conductance():
    from ip3r.io import loader
    from ip3r.physics.unitary import published, unitary
    u = unitary(loader.load("8TKF"))
    bridged = u.paired_charge.bridged
    assert sorted(b.acid[0] for b in bridged) == list("ABCD")
    assert {(b.acid_name, b.acid[1], b.base_name, b.base[1]) for b in bridged} == {
        ("ASP", 2478, "ARG", 2471)}
    assert all(b.acid[0] != b.base[0] for b in bridged)       # the neighbour's Arg
    assert all(b.distance < 2.7 for b in bridged)
    assert u.paired_charge.net_charge == u.charge.net_charge + 4.0
    # The finding: pairing lowers the charged reading further, and it stays
    # far below both the neutral pore and the measurements.
    assert u.paired.converged
    assert u.paired.conductance_pS < u.charged.conductance_pS < u.neutral.conductance_pS
    assert u.paired.conductance_pS < min(published().values()) / 5
