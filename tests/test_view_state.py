"""Seeded controls follow a parameter edit until typed over; restoring a
saved control reports what it could not set instead of forcing it."""

from __future__ import annotations

from ip3r.parameters import PARAMETERS
from ip3r.ui.view_state import Seeded, set_spin


class _Spin:
    """A spin box's surface: value, range, name, and rounding to decimals."""

    def __init__(self, value=0.0, lo=0.0, hi=100.0, decimals=2, name=""):
        self._v, self.lo, self.hi, self.dec, self.name = value, lo, hi, decimals, name

    def value(self):
        return self._v

    def setValue(self, v):
        self._v = round(min(max(v, self.lo), self.hi), self.dec)

    def minimum(self):
        return self.lo

    def maximum(self):
        return self.hi

    def objectName(self):
        return self.name


KEY = "puff.n_channels"


def test_a_seeded_spin_follows_and_a_typed_one_does_not():
    seeded, typed = _Spin(name="seeded"), _Spin(name="typed")
    s = Seeded()
    try:
        s.seed(seeded, lambda: PARAMETERS.value(KEY))
        s.seed(typed, lambda: PARAMETERS.value(KEY))
        typed.setValue(3.0)                               # the user's own
        PARAMETERS.set_value(KEY, PARAMETERS.default(KEY) + 5)
        assert s.follow() == ["seeded"]
        assert seeded.value() == PARAMETERS.default(KEY) + 5
        assert typed.value() == 3.0
        PARAMETERS.reset()
        s.follow()
        assert seeded.value() == PARAMETERS.default(KEY)  # back with the reset
    finally:
        PARAMETERS.reset()


def test_seeding_compares_the_value_the_spin_shows():
    """A value rounded by the widget still counts as seeded."""
    spin, s = _Spin(decimals=1), Seeded()
    s.seed(spin, lambda: 2.345)
    assert spin.value() == 2.3 and s.is_seeded(spin)


def test_the_getter_is_read_at_every_follow():
    """The cluster size depends on the receptor chosen: a getter, not a value."""
    box = {"v": 4.0}
    spin, s = _Spin(), Seeded()
    s.seed(spin, lambda: box["v"])
    box["v"] = 9.0
    s.follow()
    assert spin.value() == 9.0


def test_set_spin_clamps_and_refuses_with_a_note():
    spin, notes = _Spin(value=1.0, hi=10.0), []
    assert set_spin(spin, 50, "cluster size", notes) and spin.value() == 10.0
    assert "clamped" in notes[-1]
    assert not set_spin(spin, "big", "cluster size", notes) and spin.value() == 10.0
    assert not set_spin(spin, True, "cluster size", notes)
    assert len(notes) == 3
