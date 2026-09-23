"""Measurements on the real deposits (skipped when not downloaded)."""

from ip3r.io.loader import load
from ip3r.structure.channel import measure_channel
from ip3r.structure.numbering import best_numbering
from conftest import needs_structure


@needs_structure("6DQN")
def test_6dqn_channel():
    s = measure_channel(load("6DQN"))
    assert abs(s.superposition_angle - 90) < 0.1
    assert s.numbering.paralog == "ITPR3"
    assert s.constrictions["gate"].residues == ["PHE2513", "ILE2517"]
    assert s.n_ip3 == 4
    # the stubbed 1434-1546 segment is reported, not hidden
    assert any(a >= 1434 and b <= 1546 for a, b in s.numbering.mismatch_segments)


@needs_structure("7LHF")
def test_rat_itpr1_fits_no_human_numbering():
    assert best_numbering(load("7LHF")) is None


@needs_structure("8TKF", "8TKG")
def test_activated_gate_is_open():
    open_ = measure_channel(load("8TKF")).constrictions["gate"].radius
    closed = measure_channel(load("8TKG")).constrictions["gate"].radius
    assert open_ > 5.0 > 3.0 > closed
