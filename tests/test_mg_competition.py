"""Meissner 1997's competition law, checked against the paper's own tables.

The published rows below are test data (Meissner, Rios, Tripathy & Pasek
1997, JBC 272:1628): Table IV (0.5 M choline-Cl, no AMP) gives each
monovalent cation's constants at the activation site; Table II gives the
half-activation measured in 0.25 M of each salt. If the ions really compete
(Eq. 4), the first predicts the second.
"""

import pytest

from ip3r.parameters import PARAMETERS as P
from ip3r.physics import ryr_gating as rg
from ip3r.physics import spark_mg as sm
from ip3r.physics.mg_competition import equivalent_k_mg_a, hill_ka, ka_shift

#: Table IV, -AMP: (Ka µM, na, Ki µM, ni)
TABLE_IV = {"Na": (0.44, 3.4, 27e3, 1.8), "K": (0.43, 2.5, 42e3, 1.6),
            "Cs": (0.31, 3.3, 56e3, 1.4)}
#: Table II: Ka (µM) in 0.25 M of the salt
TABLE_II = {"Na": 1.30, "K": 0.92, "Cs": 0.73}


@pytest.mark.parametrize("ion", sorted(TABLE_II))
def test_competition_predicts_table_ii(ion):
    ka, na, ki, ni = TABLE_IV[ion]
    predicted = hill_ka(ka, na, 250e3, ki, ni)
    assert 1 / 1.6 < predicted / TABLE_II[ion] < 1.6


def test_without_competition_table_ii_is_missed():
    # The planted alternative: monovalent ions do not occupy the site.
    ka = TABLE_IV["Na"][0]
    assert TABLE_II["Na"] / ka > 2.5


def test_shift_closed_form_without_sodium():
    mg = 1000.0
    m = (mg / P.value("ryr.meissner_ki_mg")) ** P.value("ryr.meissner_ni_mg")
    expect = (1 + m) ** (1 / P.value("ryr.meissner_na_ca"))
    assert ka_shift(mg, na=0) == pytest.approx(expect)


def test_sodium_blunts_the_mg_shift():
    assert ka_shift(1000.0) < ka_shift(1000.0, na=0) / 4
    assert ka_shift(0.0) == pytest.approx(1.0)


def test_equivalent_reproduces_the_shift():
    for mg in (63.0, 250.0, 1000.0):
        k = equivalent_k_mg_a(mg)
        assert 1 + mg / k == pytest.approx(ka_shift(mg))
    with pytest.raises(ValueError):
        equivalent_k_mg_a(0.0)


def test_registered_values():
    assert equivalent_k_mg_a() == pytest.approx(769, abs=5)       # Murayama's NaCl
    assert equivalent_k_mg_a(na=0) == pytest.approx(98, abs=2)    # Meissner's medium


def test_reading_routes_to_the_equivalent():
    fit = rg.fit_to_bell()
    assert sm.k_mg_a_reading(fit, "meissner") == pytest.approx(equivalent_k_mg_a())
    assert "meissner" in sm.READINGS
